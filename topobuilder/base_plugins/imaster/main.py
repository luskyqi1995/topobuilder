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
from pathlib import Path
import shutil
from typing import Optional, Tuple, List, Union
from string import ascii_uppercase
import sys

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
from .analysis import get_steps
from .core import core


__all__ = ['apply', 'case_apply']


def apply( cases: List[Case], prtid: int, **kwargs ) -> List[Case]:
    """Use MASTER to correct secondary structure placement (smoothing).
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: IMASTER ---\n')

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
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case, pds_list: List[str] ) -> str:
    """
    """
    # Interactive MASTER smoothing is only applied to a Case with one single connectivity and already reoriented
    if case.connectivity_count > 1:
        raise ValueError('Interactive MASTER smoothing can only be applied to one connectivity.')
    if not case['configuration.reoriented'] and case.connectivity_count == 1:
        case = case.cast_absolute().apply_topologies()[0]

    # Find steps: Tops we will submit 2-layer searches
    steps = get_steps([x[-1] == 'E' for x in case.architecture_str.split('.')])

    # Work by layers
    for step in steps:
        pass


    # Make slurm file template
    # print(make_slurm_file(pds_list))
    return case

# apply(Case('test'), '/Volumes/MiniTwo/data/TopoBuilderData/database/master', 'master', 'createPDS', '.', 'serial', 700)
