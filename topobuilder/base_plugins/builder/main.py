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
from .architects import SSEArchitect

__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           connectivity: Optional[bool] = True,
           overwrite: Optional[bool] = False,
           **kwargs ) -> List[Case]:
    """Create coordinate entities from the Case data.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: BUILDER ---\n')

    for i, case in enumerate(cases):
        cases[i] = case_apply(case, connectivity)
        cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                connectivity: bool,
                overwrite: bool ) -> Case:
    """
    """
    if case.connectivity_count > 1:
        raise ValueError('Only single connectivity cases can be build.')
    if connectivity:
        case = case.cast_absolute().apply_topologies()[0]
    else:
        case = case.cast_absolute()

    for i, j, sse in case:
        if 'atoms' in sse['metadata'] and not overwrite:
            if TBcore.get_option('system', 'verbose'):
                sys.stdout.write('{0}.{1} already has atoms defined\n'.format(case.name, sse['id']))
            continue

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Building coordinates for {0}.{1}\n'.format(case.name, sse['id']))
        case.data['topology']['architecture'][i][j] = make_structure(sse)

    return case


def make_structure( sse: Dict ) -> Case:
    """
    """

    structure = SSEArchitect(sse, type=sse['type']).pdb
    sse['metadata'].setdefault('atoms', None)
    sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                               'Cartn_x', 'Cartn_y', 'Cartn_z']].values)

    return sse
