# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Tuple, Union
import sys

# External Libraries
from SBI.structure import PDB, PDBFrame
import pandas as pd

# This Library
import topobuilder.core as TBcore


__all__ = ['build_pdb_object']


def build_pdb_object( sses: List[Dict], loops: Union[List[int], int] ) -> Tuple[PDBFrame, List[int]]:
    """
    """
    if isinstance(loops, int):
        loops = [loops, ] * (len(sses) - 1)
    if len(sses) != len(loops) + 1:
        raise ValueError('Number of loops should equal number of SSE minus one.')

    pieces = []
    columns = ['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    start = 1
    for i, sse in enumerate(sses):
        start = 1 if i == 0 else int(sses[i - 1]['length']) + loops[i - 1] + start
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('PDB: Building SSE {:02d}:{} starting at {}\n'.format(i + 1, sse['id'], start))
        pieces.append(PDB(pd.DataFrame(sse['metadata']['atoms'], columns=columns)).renumber(start))

    structure = pd.concat(pieces)
    structure['id'] = list(range(1, structure.shape[0] + 1))
    return structure, [int(p.iloc[-1]['auth_seq_id']) for p in pieces]
