# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, Tuple
import os
import sys

# External Libraries
import SBI.databases as SBIdb
import SBI.core as SBIcr
import SBI.structure as SBIstr
import pandas as pd
import numpy as np

# This Library
import topobuilder.core as TBcore
import topobuilder.utils as TButil
try:
    from .core import core
except ImportError:
    from core import core


__all__ = ['geometric_analysis', 'get_steps', 'process_master_geometries']


def get_steps( blist: List[bool] ) -> List[Tuple[int]]:
    """
    """
    def cstp(seq):
        return [seq[i] for i in range(len(seq)) if len(seq[i]) == 1 or (seq[i][0] + 1 == seq[i][1])]

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    steps = []
    plist = list(range(len(blist)))
    betas = list(np.where(np.asarray(blist))[0])
    if len(betas) == 0:
        steps.append((0, ))
        for i in range(0, len(plist) - 1):
            steps.append(tuple(plist[i: i + 2]))
    else:
        steps.append((betas[0], ))
        for i in range(1, len(betas)):
            if betas[i] > betas[i - 1] + 1:
                steps.append((betas[i], ))
        for i in range(0, len(betas) - 1):
            steps.append(tuple(betas[i: i + 2]))
        for i in range(0, len(plist) - 1):
            steps.append(tuple(plist[i: i + 2]))
        steps = f7(cstp(steps))
    return steps


# assert get_steps([True, False]) == [(0,), (0, 1)]
# assert get_steps([False, True]) == [(1,), (0, 1)]
# assert get_steps([False, False]) == [(0,), (0, 1)]
# assert get_steps([False, False, False]) == [(0,), (0, 1), (1, 2)]
# assert get_steps([False, True, False]) == [(1,), (0, 1), (1, 2)]
# assert get_steps([False, True, True, False]) == [(1,), (1, 2), (0, 1), (2, 3)]
# assert get_steps([False, True, False, True]) == [(1,), (3,), (0, 1), (1, 2), (2, 3)]
# assert get_steps([False, True, True, False, True]) == [(1,), (4,), (1, 2), (0, 1), (2, 3), (3, 4)]


def geometric_analysis( masterdf: pd.DataFrame,
                        selections: Optional[List] = None
                        ) -> dict:
    """
    """
    if core.get_option('master', 'sample') > 0:
        wdf = masterdf.sample(core.get_option('master', 'sample'))
    else:
        wdf = masterdf
    return wdf


# def geometric_stats( filename: str,
#                      selections: Optional[List] = None,
#                      directionality: Optional[str] = None,
#                      selection: Optional[str] = None
#                      ) -> pd.DataFrame:
#     """
#     """
#     df = pd.DataFrame([x.resolve() for x in Path(filename).glob('*.pdb')], columns=['pdb_path'])
#     selection = literal_eval(selection)
#     pymol = None
#
#     def execute( row ):
#         nonlocal pymol
#         # 1. Download file
#         pdb3d = SBIstr.PDB(str(row['pdb_path']), format='pdb', clean=True, dehydrate=True, hetatms=False)
#         pdb3d = pdb3d['AtomType:CA']
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.write('  Structure {0} contains {1} residues.\n'.format(row['pdb'], pdb3d.residue_count))
#         # 2. Get pieces
#         nonlocal selection
#         pieces = make_pieces(pdb3d, selection, pymol)
#         # 3. Get Vectors
#         nonlocal directionality
#         vectors = make_vectors(pieces, directionality, pymol)
#         # 4. Get Plane
#         nonlocal selections
#         planes = make_planes(pieces, selections, pymol)
#         # 5. Vectors vs. planes
#         angles = make_angles(planes, vectors)
#         # 6. Points vs. planes
#         points = make_distances(planes, vectors)
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.flush()
#         return [row['pdb_path'], vectors, planes, angles, points]
#
#     df[['pdb_path', 'vectors', 'planes', 'angles', 'points']] = df.apply(execute, axis=1, result_type='expand')
#     return df


def process_master_geometries( masterdf: pd.DataFrame,
                               structures: List[str],
                               flip: List[str]
                               ) -> pd.DataFrame:
    """Use the PDS matches provided by MASTER to obtain the geometric properties.

    Column modifications:

    =========== ========================================= =======
    column name description                               status
    =========== ========================================= =======
    pds_path    Path to the PDS file containing the match dropped
    pdb_path    Path to the PDB file containing the match new!
    sse         Identifier of the query secondary         new!
                structre.
    layer       Layer against which the data SSE is       new!
                evaluated
    vectors     List of fitting vectors to the SSE.       new!
    planes      List of :class:`tuple` with fitting       new!
                planes and edges to the SSE.
    angles      List of angles between vectors and planes new!
    points      List of distances between the center of   new!
                each vector and the planes
    =========== ========================================= =======
    """
    # Define PDB database.
    with SBIcr.on_option_value('structure', 'format', 'pdb'):
        pdbdb = SBIdb.PDBLink(core.get_option('master', 'pdb'))
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Set PDB database as: {}\n'.format(core.get_option('master', 'pdb')))

    # Prepare data: Sort by structure-chain to avoid multi-reading
    newdf = masterdf.sort_values(['pdb', 'chain'])
    pdb3d = None

    def execute( row, structures, flip ):
        # 1. Download file
        nonlocal pdbdb
        nonlocal pdb3d
        filename, pdb3d = download_pdb(pdbdb, pdb3d, row['pdb'], row['chain'])
        pdb3d = pdb3d['AtomType:CA']
        df = TButil.pdb_geometry_from_rules(pdb3d, list(zip(structures, row['match'], flip)))
        df = df.assign(pdb=[row['pdb'], ] * df.shape[0])
        df = df.assign(chain=[row['chain'], ] * df.shape[0])
        df = df.assign(rmsd=[row['rmsd'], ] * df.shape[0])
        df = df.assign(match=[row['match'], ] * df.shape[0])
        df['pdb_path'] = [filename, ] * df.shape[0]
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.flush()
        return df

    # Get the appropiate path each structure should have.
    data = []
    for _, row in newdf.iterrows():
        data.append(execute(row, structures, flip))
    return pd.concat(data)
    return newdf.apply(execute, axis=1, result_type='broadcast', args=(structures, flip))


# def geometric_properties( masterdf: pd.DataFrame,
#                           selections: Optional[List] = None,
#                           directionality: Optional[str] = None,
#                           pymol: Optional[str] = None,
#                           ) -> pd.DataFrame:
#     """Use the PDS matches provided by MASTER to obtain the geometric
#     properties.
#
#     Column modifications:
#
#     =========== ========================================= =======
#     column name description                               status
#     =========== ========================================= =======
#     pds_path    Path to the PDS file containing the match dropped
#     pdb_path    Path to the PDB file containing the match new!
#     vectors     List of fitting vectors to the SSE.       new!
#     planes      List of :class:`tuple` with fitting       new!
#                 planes and edges to the SSE.
#     angles      List of angles between vectors and planes new!
#     points      List of distances between the center of   new!
#                 each vector and the planes
#     =========== ========================================= =======
#     """
#     # Define PDB database.
#     with SBIcr.on_option_value('structure', 'format', 'pdb'):
#         pdbdb = SBIdb.PDBLink(core.get_option('master', 'pdb'))
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.write('Set PDB database as: {}\n'.format(core.get_option('master', 'pdb')))
#
#     # Prepare data: Sort by structure-chain to avoid multi-reading
#     newdf = masterdf.sort_values(['pdb', 'chain'])
#     pdb3d = None
#
#     # PyMOL
#     pymol = open(pymol, 'w') if pymol is not None else None
#
#     def execute( row ):
#         # 1. Download file
#         nonlocal pdbdb
#         nonlocal pdb3d
#         nonlocal pymol
#         if core.get_option('imaster', 'pymol'):
#             pymol.write('fetch {}\n'.format(row['pdb']))
#             pymol.write('remove solvent\n')
#         filename, pdb3d = download_pdb(pdbdb, pdb3d, row['pdb'], row['chain'])
#         pdb3d = pdb3d['AtomType:CA']
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.write('  Structure {0} contains {1} residues.\n'.format(row['pdb'], pdb3d.residue_count))
#         # 2. Get pieces
#         pieces = make_pieces(pdb3d, row['match'], pymol)
#         # 3. Get Vectors
#         nonlocal directionality
#         vectors = make_vectors(pieces, directionality, pymol)
#         # 4. Get Plane
#         nonlocal selections
#         planes = make_planes(pieces, selections, pymol)
#         # 5. Vectors vs. planes
#         angles = make_angles(planes, vectors)
#         # 6. Points vs. planes
#         points = make_distances(planes, vectors)
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.flush()
#         return [filename, vectors, planes, angles, points]
#
#     # Get the appropiate path each structure should have.
#     newdf[['pdb_path', 'vectors', 'planes', 'angles', 'points']] = newdf.apply(execute, axis=1, result_type='expand')
#     if pymol is not None:
#         pymol.close()
#     return newdf.drop(columns='pds_path')


def download_pdb( pdbdb: SBIdb.PDBLink,
                  pdb3d: SBIstr.PDBFrame,
                  pdbid: str,
                  pdbchain: str
                  ) -> Tuple[str, SBIstr.PDBFrame]:
    filename = pdbdb.store_local_path('{0}_{1}.pdb.gz'.format(pdbid, pdbchain))
    if not os.path.isfile(filename):
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('  Downloading PDB {}\n'.format(pdbid))
        pdb3d = SBIstr.PDB('fetch:{0}'.format(pdbid), format='pdb', clean=True,
                           dehydrate=True, hetatms=False)['Chain:{}'.format(pdbchain)]
        pdb3d.write(filename, format='pdb')
        return filename, pdb3d

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('  {} already exists\n'.format(filename))
    if pdb3d is None or pdb3d.id != '{0}_{1}'.format(pdbid, pdbchain):
        pdb3d = SBIstr.PDB('fetch:{0}'.format(pdbid), format='pdb', clean=True,
                           dehydrate=True, hetatms=False)['Chain:{}'.format(pdbchain)]
        return filename, pdb3d
    else:
        return filename, pdb3d

#
# def make_pieces( pdb3d: SBIstr.PDBFrame,
#                  matches: List[Tuple[int]],
#                  pymol: Optional[TextIO] = None
#                  ) -> List[SBIstr.PDBFrame]:
#     pieces = []
#     if TBcore.get_option('system', 'verbose'):
#         sys.stdout.write('  Generating {0} vectors for {1}\n'.format(len(matches), pdb3d.id))
#     for i, segment in enumerate(matches):
#         with SBIcr.on_option_value('structure', 'source', 'label'):
#             # MASTER starts match count at 0!
#             piece = pdb3d['Residue:{0}-{1}'.format(segment[0] + 1, segment[1] + 1)]
#             pieces.append(piece)
#         with SBIcr.on_option_value('structure', 'source', 'auth'):
#             if TBcore.get_option('system', 'verbose'):
#                 first, last = piece.first_compound.number, piece.last_compound.number
#                 sys.stdout.write('    {2} - Range: {0}-{1}\n'.format(first, last, pdb3d.id))
#             if core.get_option('imaster', 'pymol'):
#                 first, last = piece.first_compound.number, piece.last_compound.number
#                 pymol.write('sele sse_{0:02d}, {1} AND c. {2} AND i. {3}-{4}\n'.format(
#                     i + 1, pdb3d.id.lower().split('_')[0], piece.chain, first, last))
#     return pieces
#
#
# def make_vectors( pieces: List[SBIstr.PDBFrame],
#                   directionality: Optional[str] = None,
#                   pymol: Optional[TextIO] = None
#                   ) -> List[List[float]]:
#     vectors = []
#     if directionality is None:
#         directionality = '0' * len(pieces)
#     for i, segment in enumerate(pieces):
#         vectors.append([list(x) for x in segment.eigenvectors(40)[-1]])
#         if bool(int(directionality[i])):
#             vectors[-1] = [list(x) for x in np.flip(np.asarray(vectors[-1]), axis=0)]
#         if core.get_option('imaster', 'pymol'):
#             pymol.write('make_arrow("vect{0:02d}", {1}, {2})\n'.format(i, vectors[-1][0], vectors[-1][-1]))
#     return vectors
#
#
# def make_planes( pieces: List[SBIstr.PDBFrame],
#                  selections: Optional[List[List[int]]] = None,
#                  pymol: Optional[TextIO] = None
#                  ) -> List[List[float]]:
#     planes = []
#     if selections is None:
#         selections = [list(range(len(pieces)))]
#     for i, sele in enumerate(selections):
#         structure = pd.concat(itemgetter(*sele)(pieces))
#         eign = structure.eigenvectors(30)
#         planes.append([list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])])  # layer plane
#         planes.append([list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])])  # floor plane
#         planes.append([list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])])  # side plane
#         if core.get_option('imaster', 'pymol'):
#             cmd = 'make_plane_points("plane{0:02d}_{4}", {1}, {2}, {3}, settings={{"ALPHA": 0.3}})\n'
#             pymol.write(cmd.format(i, planes[-3][0], planes[-3][1], planes[-3][2], 'layer'))
#             pymol.write(cmd.format(i, planes[-2][0], planes[-2][1], planes[-2][2], 'floor'))
#             pymol.write(cmd.format(i, planes[-1][0], planes[-1][1], planes[-1][2], 'side'))
#             for i, v in enumerate([('ePrep', 'red'), ('eSide', 'green'), ('eMajor', 'blue')]):
#                 pymol.write('make_arrow("{0}", {1}, {2}, color="{3}")\n'.format(v[0], list(eign[i][0]), list(eign[i][-1]), v[1]))
#             pymol.flush()
#     return planes
#
#
# def make_angles( planes: List[List[float]],
#                  vectors: List[List[float]]
#                  ) -> List[List[float]]:
#     angles = []
#     for p in planes:
#         angles.append([])
#         p = sy.Plane(sy.Point3D(p[0]), sy.Point3D(p[1]), sy.Point3D(p[2]))
#         for v in vectors:
#             v = sy.Line(v[0], v[-1])
#             angles[-1].append(math.degrees(p.angle_between(v)))
#     return angles
#
#
# def make_distances( planes: List[List[float]],
#                     vectors: List[List[float]]
#                     ) -> List[List[float]]:
#     distances = []
#     for p in planes:
#         distances.append([])
#         p = sy.Plane(sy.Point3D(p[0]), sy.Point3D(p[1]), sy.Point3D(p[2]))
#         for v in vectors:
#             v = sy.Point3D(*v[1])
#             distances[-1].append(float(p.distance(v)))
#     return distances
