# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Optional
import sys

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from .parametric import SSEArchitect

__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           connectivity: Optional[bool] = True,
           pick_aa: Optional[str] = None,
           **kwargs ) -> List[Case]:
    """Create coordinate entities from the Case data.
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i] = case_apply(case, connectivity, pick_aa, write2disc=True)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                connectivity: bool,
                pick_aa: Optional[str] = None,
                write2disc: Optional[bool] = False
                ) -> Case:
    """
    """
    # Bloc muli-connectivities.
    if case.connectivity_count > 1:
        raise ValueError('Only single connectivity cases can be build.')

    # Apply connectivity?
    if connectivity:
        case = case.cast_absolute().apply_topologies()[0]
        ofile = case.connectivities_paths[0].joinpath('directed_sketch.pdb')
    else:
        case = case.cast_absolute()
        ofile = case.main_path.joinpath('architecture').joinpath('undirected_sketch.pdb')

    for i, j, sse in case:
        if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
            if TBcore.get_option('system', 'verbose'):
                sys.stdout.write('{0}.{1} already has atoms defined\n'.format(case.name, sse['id']))
            continue

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Building coordinates for {0}.{1}\n'.format(case.name, sse['id']))
        case.data['topology']['architecture'][i][j] = make_structure(sse, pick_aa)

    if write2disc:
        ofile.parent.mkdir(parents=True, exist_ok=True)
        structure, _ = TButil.build_pdb_object( case.ordered_structures, 2 )

        TButil.plugin_filemaker('Writing structure {0}'.format(ofile))
        structure.write(output_file=str(ofile), format='pdb', clean=True,
                        force=TBcore.get_option('system', 'overwrite'))

    return case


def make_structure( sse: Dict, pick_aa: Optional[str] = None ) -> Case:
    """
    """

    structure = SSEArchitect(sse, type=sse['type'], pick_aa=pick_aa).pdb
    sse['metadata'].setdefault('atoms', None)
    sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                               'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
    return sse
