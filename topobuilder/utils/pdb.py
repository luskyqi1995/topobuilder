# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Tuple, Union
from pathlib import Path
import sys
import math
from string import ascii_uppercase

# External Libraries
try:
    from SBI.structure import PDB, Frame3D
    import SBI.core as SBIcr
except ImportError:
    class Frame3D():
        pass

import pandas as pd
import numpy as np
import sympy as sy

# This Library
import topobuilder.core as TBcore


__all__ = ['build_pdb_object', 'pdb_geometry_from_rules']


np.set_printoptions(precision=3)


def build_pdb_object( sses: List[Dict], loops: Union[List[int], int] ) -> Tuple[Frame3D, List[int]]:
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

    structure = pd.concat(pieces, sort=False).reset_index()
    structure['id'] = list(range(1, structure.shape[0] + 1))
    return structure, [int(p.iloc[-1]['auth_seq_id']) for p in pieces]


def pdb_geometry_from_rules( pdb_file: Union[Path, str, Frame3D], rules: List[Tuple] ) -> pd.DataFrame:
    """
    """
    if isinstance(pdb_file, (Path, str)):
        pdb_file = Path(pdb_file)
        if not pdb_file.is_file():
            raise IOError('PDB structure {} not found.'.format(pdb_file))
        pdb3d = PDB(str(pdb_file), format='pdb', clean=True, dehydrate=True, hetatms=False)['AtomTask:PROTEINBACKBONE']
    elif isinstance(pdb_file, Frame3D):
        pdb3d = pdb_file['AtomTask:PROTEINBACKBONE']
    else:
        raise ValueError('Unexpected type for pdb_file.')
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('PDB:Analyzing geometry of {}\n'.format(pdb3d.id))

    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('PDB:Available secondary structures {}\n'.format(','.join([x[0] for x in rules])))
        sys.stdout.write('PDB:With ranges {}\n'.format(','.join(['{}-{}'.format(*x[1]) for x in rules])))
        sys.stdout.write('PDB:With flip policy {}\n'.format(','.join([str(x[2]) for x in rules])))

    pieces = make_pieces(pdb3d, rules)
    pieces = make_vectors(pieces, rules)
    pieces = make_planes(pieces)
    df = make_angles_and_distances(pieces)
    df = df.assign(pdb_path=[str(pdb_file) if not isinstance(pdb_file, Frame3D) else pdb3d.id, ] * df.shape[0])
    return df


def make_pieces( pdb3d: Frame3D, rules: List[Tuple] ) -> Dict:
    """
    """
    pieces = {}
    for piece in rules:
        sse_id, ranges, _ = piece
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Individualize SSE:{} of {}\n'.format(sse_id, pdb3d.id))
        with SBIcr.on_option_value('structure', 'source', 'label'):
            segment = pdb3d['Residue:{0}-{1}'.format(ranges[0], ranges[1])]
            pieces.setdefault(sse_id, {}).setdefault('atoms', segment)
        with SBIcr.on_option_value('structure', 'source', 'auth'):
            if TBcore.get_option('system', 'verbose'):
                first, last = segment.first_compound.number, segment.last_compound.number
                sys.stdout.write('PDB:{2} - Range: {0}-{1}\n'.format(first, last, pdb3d.id))
    return pieces


def make_vectors( pieces: Dict, rules: List[Tuple] ) -> Dict:
    """
    """
    for piece in rules:
        sse_id, _, flip = piece
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Getting vector for SSE:{}{}\n'.format(sse_id, '' if not flip else ' - flipped!'))
        pieces[sse_id].setdefault('vector', list(pieces[sse_id]['atoms'].eigenvectors(40)[-1]))
        if flip:
            pieces[sse_id]['vector'] = [list(x) for x in np.flip(np.asarray(pieces[sse_id]['vector']), axis=0)]
        print('PYMOL:make_arrow("vect{0}", {1}, {2})'.format(sse_id, np.asarray(pieces[sse_id]['vector'][0]).tolist(),
                                                             np.asarray(pieces[sse_id]['vector'][-1]).tolist()))
    return pieces


def make_planes( pieces: Dict ) -> Dict:
    """
    """
    blayers = sorted(set([x[0] for x in pieces if x.endswith('E') and len(x) == 3]))
    hlayers = [x[0] for x in pieces if x.endswith('H') and len(x) == 3]
    hlayers = sorted(set([x for x in hlayers if hlayers.count(x) > 1]))

    for layer in blayers:
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Generating plane for beta layer {}\n'.format(layer))
        pieces.setdefault(layer, {'layer': [], 'floor': [], 'side': []})
        structure = pd.concat([pieces[x]['atoms'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        eign = structure.eigenvectors(30)
        eign = eigenlayers_fix(eign, [pieces[x]['vector'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        pieces[layer]['layer'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])]
        pieces[layer]['floor'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])]
        pieces[layer]['side'] = [list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])]
        for i, v in enumerate([('floor', 'red'), ('side', 'green'), ('layer', 'blue')]):
            print('PYMOL:make_arrow("b{0}", {1}, {2}, color="{3}")'.format(v[0], eign[i][0].tolist(),
                                                                           eign[i][-1].tolist(), v[1]))

    for layer in hlayers:
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Generating plane for helix layer {}\n'.format(layer))
        pieces.setdefault(layer, {'layer': [], 'floor': [], 'side': []})
        structure = pd.concat([pieces[x]['atoms'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        eign = structure.eigenvectors(30)
        eign = eigenlayers_fix(eign, [pieces[x]['vector'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        pieces[layer]['layer'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])]
        pieces[layer]['floor'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])]
        pieces[layer]['side'] = [list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])]
        for i, v in enumerate([('floor', 'red'), ('side', 'green'), ('layer', 'blue')]):
            print('PYMOL:make_arrow("h{0}", {1}, {2}, color="{3}")'.format(v[0], eign[i][0].tolist(),
                                                                           eign[i][-1].tolist(), v[1]))

    return pieces


def eigenlayers_fix( eign: np.ndarray, vectors: np.ndarray, scape: bool = False ) -> np.ndarray:
    """
    """
    eign = eign.copy()
    eign2 = sy.Line(eign[2][0], eign[2][-1])
    angles_layer = []
    for sse in vectors:
        angles_layer.append(math.degrees(sy.Line(sse[0], sse[-1]).angle_between(eign2)))

    # Fix layer direction
    if len(vectors) > 3:
        do_flip = np.isclose(angles_layer, [180, ] * len(angles_layer), atol=40)
        do_flip = sum(do_flip) >= len(do_flip) - 1
    else:
        do_flip = np.allclose(angles_layer, [180, ] * len(angles_layer), atol=40)
    # Fix vectors -> switch side and layer
    if len(vectors) > 3:
        do_switch = [not x for x in np.isclose(angles_layer, [0, ] * len(angles_layer), atol=40)]
        do_switch = sum(do_switch) >= len(do_switch) - 1
    else:
        do_switch = not np.allclose(angles_layer, [0, ] * len(angles_layer), atol=35)

    # Apply
    if do_flip:
        eign[2] = [list(x) for x in np.flip(np.asarray(eign[2]), axis=0)]
    elif do_switch and not scape:
        eign[[1, 2]] = eign[[2, 1]]
        eign = eigenlayers_fix(eign, vectors, True)

    return eign


def make_angles_and_distances( pieces: Dict ) -> pd.DataFrame:
    """
    """
    data = {'sse': [], 'layer': [],
            'angles_layer': [], 'angles_floor': [], 'angles_side': [],
            'points_layer': [], 'points_floor': [], 'points_side': [],
            'tilted_layer': [], 'tilted_floor': [], 'tilted_side': []}

    for layer in sorted(set([x[0] for x in pieces if len(x) == 1])):
        for sse in [x for x in pieces if len(x) == 3]:
            if abs(ascii_uppercase.find(layer) - ascii_uppercase.find(sse[0])) <= 1:
                data['sse'].append(sse)
                data['layer'].append(layer)
                for iplane, plane in enumerate(pieces[layer]):
                    if TBcore.get_option('system', 'debug'):
                        sys.stdout.write('PDB:{} geometry plane {} vs. sse {}\n'.format(plane, layer, sse))
                    syPlane = sy.Plane(sy.Point3D(pieces[layer][plane][0]),
                                       sy.Point3D(pieces[layer][plane][1]),
                                       sy.Point3D(pieces[layer][plane][2]))
                    syLine = sy.Line(pieces[sse]['vector'][0], pieces[sse]['vector'][-1])
                    syPoint = sy.Point3D(*pieces[sse]['vector'][1])
                    data[f'angles_{plane}'].append(math.degrees(syPlane.angle_between(syLine)))
                    data[f'points_{plane}'].append(float(syPlane.distance(syPoint)))
                    data[f'tilted_{plane}'].append(float(syPlane.distance(default_plane(iplane))))
    return pd.DataFrame(data)


def default_plane( pick: int ) -> sy.Plane:
    """
    """
    x = sy.Point3D([30, 0, 0])
    y = sy.Point3D([0, 30, 0])
    z = sy.Point3D([0, 0, 30])
    c = sy.Point3D([0, 0, 0])

    if pick == 0:
        return sy.Plane(y, c, x)
    elif pick == 1:
        return sy.Plane(x, c, z)
    elif pick == 2:
        return sy.Plane(z, c, y)
    else:
        raise ValueError('Selection must be between 0 and 2')
