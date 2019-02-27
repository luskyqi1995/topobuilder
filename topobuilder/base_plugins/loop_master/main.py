# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import getpass
import os
import sys
from pathlib import Path
import shutil
from typing import Optional, Tuple, List, Union
from string import ascii_uppercase
from itertools import islice

# External Libraries
import pandas as pd
from SBI.structure import PDB

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
from topobuilder import plugin_source
from .core import core


def apply( cases: List[Case],
           prtid: int,
           **kwargs ) -> List[Case]:
    """
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: LOOP_MASTER ---\n')

    # Get list of PDS structures
    database = core.get_option('master', 'pds')
    pds_list = []
    database = Path(database)
    if database.is_file():
        pds_list = [line.strip() for line in open(database).readlines() if len(line.strip()) > 0]
    elif database.is_dir():
        pds_list = [str(x) for x in database.glob('*/*.pds')]
    else:
        raise ValueError('The provided MASTER database directory/list file cannot be found.')

    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('corrections', [])
        cases[i].data['metadata']['corrections'].append(case_apply(case, pds_list))
        cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case, pds_list: List[str] ) -> str:
    """
    """
    # Loop MASTER is only applied to a Case with one single connectivity and already reoriented
    if case.connectivity_count > 1:
        raise ValueError('Loop MASTER can only be applied to one connectivity.')
    if not case['configuration.reoriented'] and case.connectivity_count == 1:
        # We will need the coordinates of the secondary structures to execute this one
        # This will already cast it to absolute
        case = plugin_source.load_plugin('builder').case_apply(case, connectivity=True, overwrite=False)

    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0]
    folders.mkdir(parents=True, exist_ok=True)

    # Find steps: Each pair of secondary structure.
    it = case.connectivities_str[0].split('.')
    steps = [it[i:i + 2] for i in range(0, len(it) - 1)]

    for i, sse in enumerate(steps):
        wfolder = folders.joinpath('loop{:02d}'.format(i + 1))
        wfolder.mkdir(parents=True, exist_ok=True)

        sse1 = pd.DataFrame(case.get_sse_by_id(sse[0])['metadata']['atoms'],
                            columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'])
        sse2 = pd.DataFrame(case.get_sse_by_id(sse[1])['metadata']['atoms'],
                            columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'])
        structure = PDB(pd.concat([sse1, sse2])).renumber()
        outfile = wfolder.joinpath('loop_master.jump{:02d}.pdb'.format(i + 1))
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('-> generating structure {}\n'.format(outfile.resolve()))
        structure.write(output_file=str(outfile), format='pdb', clean=True,
                        force=TBcore.get_option('system', 'overwrite'))

    return case
