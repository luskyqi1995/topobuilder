# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Optional
from pathlib import Path
import copy
import sys

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from .architects import SSEArchitect

__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           connectivity: Optional[bool] = True,
           **kwargs ) -> List[Case]:
    """Create coordinate entities from the Case data.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: BUILDER ---\n')

    for i, case in enumerate(cases):
        cases[i] = case_apply(case, connectivity, write2disc=True)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                connectivity: bool,
                write2disc: Optional[bool] = False ) -> Case:
    """
    """
    if case.connectivity_count > 1:
        raise ValueError('Only single connectivity cases can be build.')

    if connectivity:
        case = case.cast_absolute().apply_topologies()[0]
    else:
        case = case.cast_absolute()

    for i, j, sse in case:
        if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
            if TBcore.get_option('system', 'verbose'):
                sys.stdout.write('{0}.{1} already has atoms defined\n'.format(case.name, sse['id']))
            continue

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Building coordinates for {0}.{1}\n'.format(case.name, sse['id']))
        case.data['topology']['architecture'][i][j] = make_structure(sse)

    if write2disc:
        if connectivity:
            # Generate the folder tree for a single connectivity.
            folders = case.connectivities_paths[0]
            outfile = folders.joinpath('directed_sketch.pdb')
        else:
            # Generate the folder tree for the sketch.
            folders = case.main_path.joinpath('architecture')
            outfile = folders.joinpath('undirected_sketch.pdb')
        folders.mkdir(parents=True, exist_ok=True)

        structure, _ = TButil.build_pdb_object( case.ordered_structures, 2 )

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Writing structure {0}\n'.format(outfile))
        structure.write(output_file=str(outfile), format='pdb', clean=True,
                        force=TBcore.get_option('system', 'overwrite'))

    return case


def make_structure( sse: Dict ) -> Case:
    """
    """

    structure = SSEArchitect(sse, type=sse['type']).pdb
    sse['metadata'].setdefault('atoms', None)
    sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                               'Cartn_x', 'Cartn_y', 'Cartn_z']].values)

    return sse
