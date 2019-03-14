# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, TextIO, Tuple
from operator import itemgetter
from ast import literal_eval
from pathlib import Path
import os
import math
import sys

# External Libraries
import SBI.databases as SBIdb
import SBI.core as SBIcr
import SBI.structure as SBIstr
import pandas as pd
import numpy as np
import sympy as sy

# This Library
import topobuilder.core as TBcore
import topobuilder.utils as TButil
try:
    from .core import core
except ImportError:
    from core import core


__all__ = ['geometric_analysis', 'geometric_stats', 'geometric_properties',
           'parse_master_file', 'get_steps', 'process_master_geometries']


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


def geometric_stats( filename: str,
                     selections: Optional[List] = None,
                     directionality: Optional[str] = None,
                     selection: Optional[str] = None
                     ) -> pd.DataFrame:
    """
    """
    df = pd.DataFrame([x.resolve() for x in Path(filename).glob('*.pdb')], columns=['pdb_path'])
    selection = literal_eval(selection)
    pymol = None

    def execute( row ):
        nonlocal pymol
        # 1. Download file
        pdb3d = SBIstr.PDB(str(row['pdb_path']), format='pdb', clean=True, dehydrate=True, hetatms=False)
        pdb3d = pdb3d['AtomType:CA']
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('  Structure {0} contains {1} residues.\n'.format(row['pdb'], pdb3d.residue_count))
        # 2. Get pieces
        nonlocal selection
        pieces = make_pieces(pdb3d, selection, pymol)
        # 3. Get Vectors
        nonlocal directionality
        vectors = make_vectors(pieces, directionality, pymol)
        # 4. Get Plane
        nonlocal selections
        planes = make_planes(pieces, selections, pymol)
        # 5. Vectors vs. planes
        angles = make_angles(planes, vectors)
        # 6. Points vs. planes
        points = make_distances(planes, vectors)
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.flush()
        return [row['pdb_path'], vectors, planes, angles, points]

    df[['pdb_path', 'vectors', 'planes', 'angles', 'points']] = df.apply(execute, axis=1, result_type='expand')
    return df


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


def geometric_properties( masterdf: pd.DataFrame,
                          selections: Optional[List] = None,
                          directionality: Optional[str] = None,
                          pymol: Optional[str] = None,
                          ) -> pd.DataFrame:
    """Use the PDS matches provided by MASTER to obtain the geometric
    properties.

    Column modifications:

    =========== ========================================= =======
    column name description                               status
    =========== ========================================= =======
    pds_path    Path to the PDS file containing the match dropped
    pdb_path    Path to the PDB file containing the match new!
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

    # PyMOL
    pymol = open(pymol, 'w') if pymol is not None else None

    def execute( row ):
        # 1. Download file
        nonlocal pdbdb
        nonlocal pdb3d
        nonlocal pymol
        if core.get_option('imaster', 'pymol'):
            pymol.write('fetch {}\n'.format(row['pdb']))
            pymol.write('remove solvent\n')
        filename, pdb3d = download_pdb(pdbdb, pdb3d, row['pdb'], row['chain'])
        pdb3d = pdb3d['AtomType:CA']
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('  Structure {0} contains {1} residues.\n'.format(row['pdb'], pdb3d.residue_count))
        # 2. Get pieces
        pieces = make_pieces(pdb3d, row['match'], pymol)
        # 3. Get Vectors
        nonlocal directionality
        vectors = make_vectors(pieces, directionality, pymol)
        # 4. Get Plane
        nonlocal selections
        planes = make_planes(pieces, selections, pymol)
        # 5. Vectors vs. planes
        angles = make_angles(planes, vectors)
        # 6. Points vs. planes
        points = make_distances(planes, vectors)
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.flush()
        return [filename, vectors, planes, angles, points]

    # Get the appropiate path each structure should have.
    newdf[['pdb_path', 'vectors', 'planes', 'angles', 'points']] = newdf.apply(execute, axis=1, result_type='expand')
    if pymol is not None:
        pymol.close()
    return newdf.drop(columns='pds_path')


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


def make_pieces( pdb3d: SBIstr.PDBFrame,
                 matches: List[Tuple[int]],
                 pymol: Optional[TextIO] = None
                 ) -> List[SBIstr.PDBFrame]:
    pieces = []
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('  Generating {0} vectors for {1}\n'.format(len(matches), pdb3d.id))
    for i, segment in enumerate(matches):
        with SBIcr.on_option_value('structure', 'source', 'label'):
            # MASTER starts match count at 0!
            piece = pdb3d['Residue:{0}-{1}'.format(segment[0] + 1, segment[1] + 1)]
            pieces.append(piece)
        with SBIcr.on_option_value('structure', 'source', 'auth'):
            if TBcore.get_option('system', 'verbose'):
                first, last = piece.first_compound.number, piece.last_compound.number
                sys.stdout.write('    {2} - Range: {0}-{1}\n'.format(first, last, pdb3d.id))
            if core.get_option('imaster', 'pymol'):
                first, last = piece.first_compound.number, piece.last_compound.number
                pymol.write('sele sse_{0:02d}, {1} AND c. {2} AND i. {3}-{4}\n'.format(
                    i + 1, pdb3d.id.lower().split('_')[0], piece.chain, first, last))
    return pieces


def make_vectors( pieces: List[SBIstr.PDBFrame],
                  directionality: Optional[str] = None,
                  pymol: Optional[TextIO] = None
                  ) -> List[List[float]]:
    vectors = []
    if directionality is None:
        directionality = '0' * len(pieces)
    for i, segment in enumerate(pieces):
        vectors.append([list(x) for x in segment.eigenvectors(40)[-1]])
        if bool(int(directionality[i])):
            vectors[-1] = [list(x) for x in np.flip(np.asarray(vectors[-1]), axis=0)]
        if core.get_option('imaster', 'pymol'):
            pymol.write('make_arrow("vect{0:02d}", {1}, {2})\n'.format(i, vectors[-1][0], vectors[-1][-1]))
    return vectors


def make_planes( pieces: List[SBIstr.PDBFrame],
                 selections: Optional[List[List[int]]] = None,
                 pymol: Optional[TextIO] = None
                 ) -> List[List[float]]:
    planes = []
    if selections is None:
        selections = [list(range(len(pieces)))]
    for i, sele in enumerate(selections):
        structure = pd.concat(itemgetter(*sele)(pieces))
        eign = structure.eigenvectors(30)
        planes.append([list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])])  # layer plane
        planes.append([list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])])  # floor plane
        planes.append([list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])])  # side plane
        if core.get_option('imaster', 'pymol'):
            cmd = 'make_plane_points("plane{0:02d}_{4}", {1}, {2}, {3}, settings={{"ALPHA": 0.3}})\n'
            pymol.write(cmd.format(i, planes[-3][0], planes[-3][1], planes[-3][2], 'layer'))
            pymol.write(cmd.format(i, planes[-2][0], planes[-2][1], planes[-2][2], 'floor'))
            pymol.write(cmd.format(i, planes[-1][0], planes[-1][1], planes[-1][2], 'side'))
            for i, v in enumerate([('ePrep', 'red'), ('eSide', 'green'), ('eMajor', 'blue')]):
                pymol.write('make_arrow("{0}", {1}, {2}, color="{3}")\n'.format(v[0], list(eign[i][0]), list(eign[i][-1]), v[1]))
            pymol.flush()
    return planes


def make_angles( planes: List[List[float]],
                 vectors: List[List[float]]
                 ) -> List[List[float]]:
    angles = []
    for p in planes:
        angles.append([])
        p = sy.Plane(sy.Point3D(p[0]), sy.Point3D(p[1]), sy.Point3D(p[2]))
        for v in vectors:
            v = sy.Line(v[0], v[-1])
            angles[-1].append(math.degrees(p.angle_between(v)))
    return angles


def make_distances( planes: List[List[float]],
                    vectors: List[List[float]]
                    ) -> List[List[float]]:
    distances = []
    for p in planes:
        distances.append([])
        p = sy.Plane(sy.Point3D(p[0]), sy.Point3D(p[1]), sy.Point3D(p[2]))
        for v in vectors:
            v = sy.Point3D(*v[1])
            distances[-1].append(float(p.distance(v)))
    return distances


# def parse_master_file( filename: str,
#                        max_rmsd: Optional[float] = None
#                        ) -> pd.DataFrame:
#     """Load MASTER output data.
#
#     :param str filename: Output file.
#     :param float max_rmsd: Maximum RMSD value to recover.
#
#     :return: :class:`~pandas.DataFrame`
#
#     Columns of the returned :class:`~pandas.DataFrame` are:
#
#     =========== ===================================================
#     column name description
#     =========== ===================================================
#     rmsd        RMSD value between query and match
#     pds_path    Path to the PDS file containing the match
#     pdb         PDB identifier
#     chain       Chain identifier
#     match       List of ranges for the match (Rosetta count)
#     =========== ===================================================
#
#     This assumes that the PDS files basename has the standard nomenclature
#     ``<pdbid>_<chain>.pds``.
#     """
#     if TBcore.get_option('system', 'verbose'):
#         sys.stdout.write('Reading MASTER file {}\n'.format(filename))
#     df = pd.read_csv(filename,
#                      names=list(range(20)), engine='python',
#                      sep=r'\s+', header=None).dropna(axis=1, how='all')
#     df['match'] = df[df.columns[2:]].astype(str).sum(axis=1).apply(literal_eval)
#     df = df.rename(columns={0: 'rmsd', 1: 'pds_path'})
#     df = df.drop(columns=[i for i in df.columns if isinstance(i, int)])
#     df[['pdb', 'chain']] = (pd.DataFrame(list(df['pds_path'].str.replace('.pds', '')
#                             .apply(lambda x: os.path.basename(x).split('_')).values)))
#     if max_rmsd is not None:
#         df = df[(df['rmsd'] <= max_rmsd)]
#     return df[['rmsd', 'pds_path', 'pdb', 'chain', 'match']]

#
# def get_pdbs(masterdf: pd.DataFrame) -> pd.DataFrame:
#     """From the MASTER output, find the matching structures and their paths.
#
#     If the structure is not found in the database, download it and select the
#     appropiate chain.
#
#     :param masterdf: Result from the MASTER search.
#     :type masterdf: :class:`~pandas.DataFrame`
#
#     :return: :class:`~pandas.DataFrame`
#
#     Column modifications:
#     =========== ========================================= =======
#     column name description                               status
#     =========== ========================================= =======
#     pds_path    Path to the PDS file containing the match dropped
#     pdb_path    Path to the PDB file containing the match new!
#     =========== ========================================= =======
#     """
#     def get_structure(row, pdbdb):
#         filename = pdbdb.store_local_path('{0}_{1}.pdb'.format(row['pdb'], row['chain']))
#         if not os.path.isfile(filename):
#             if core.get_option('imaster', 'verbose'):
#                 print('  Downloading PDB {}'.format(row['pdb']))
#             pdb3d = SBIstr.PDB('fetch:{0}'.format(row['pdb']), format='pdb', clean=True,
#                                dehydrate=True, hetatms=False)['Chain:{}'.format(row['chain'])]
#             pdb3d.write(filename, format='pdb')
#             if core.get_option('imaster', 'verbose'):
#                 print('  Structure {0} contains {1} residues.'.format(row['pdb'], pdb3d.residue_count))
#         else:
#             if core.get_option('imaster', 'verbose'):
#                 print('  {} already exists'.format(filename))
#         return filename
#
#     # Define PDB database.
#     if core.get_option('imaster', 'verbose'):
#         print('Retrieving PDB files for analysis')
#     with SBIcr.on_option_value('structure', 'format', 'pdb'):
#         pdbdb = SBIdb.PDBLink(core.get_option('imaster', 'pdb'))
#     newdf = masterdf
#     # Get the appropiate path each structure should have.
#     newdf = newdf.assign(pdb_path=masterdf.apply(lambda row: get_structure(row, pdbdb), axis=1))
#     return newdf[['rmsd', 'pdb_path', 'pdb', 'chain', 'match']]
#
#
# def parse_master_file( filename: str,
#                        max_rmsd: Optional[float] = None
#                        ) -> pd.DataFrame:
#     """Load MASTER output data.
#
#     :param str filename: Output file.
#     :param float max_rmsd: Maximum RMSD value to recover.
#
#     :return: :class:`~pandas.DataFrame`
#
#     Columns of the returned :class:`~pandas.DataFrame` are:
#
#     =========== ===================================================
#     column name description
#     =========== ===================================================
#     rmsd        RMSD value between query and match
#     pds_path    Path to the PDS file containing the match
#     pdb         PDB identifier
#     chain       Chain identifier
#     match       List of ranges for the match (Rosetta count)
#     =========== ===================================================
#
#     This assumes that the PDS files basename has the standard nomenclature
#     ``<pdbid>_<chain>.pds``.
#     """
#     if core.get_option('imaster', 'verbose'):
#         print('Reading MASTER file {}'.format(filename))
#     df = pd.read_csv(filename,
#                      names=list(range(20)), engine='python',
#                      sep=r'\s+', header=None).dropna(axis=1, how='all')
#     df['match'] = df[df.columns[2:]].astype(str).sum(axis=1).apply(literal_eval)
#     df = df.rename(columns={0: 'rmsd', 1: 'pds_path'})
#     df = df.drop(columns=[i for i in df.columns if isinstance(i, int)])
#     df[['pdb', 'chain']] = (pd.DataFrame(list(df['pds_path'].str.replace('.pds', '')
#                             .apply(lambda x: os.path.basename(x).split('_')).values)))
#     if max_rmsd is not None:
#         df = df[(df['rmsd'] <= max_rmsd)]
#     return df[['rmsd', 'pds_path', 'pdb', 'chain', 'match']]
#
#
# def fit_vector_to_data(df: pd.DataFrame,
#                        backbone: bool = False
#                        ) -> pd.DataFrame:
#     """Add fitting lines to secondary structure segments.
#
#     :param df: Data derived from the MASTER search. Requires the
#         ``pdb_path`` column.
#     :type df: :class:`~pandas.DataFrame`
#     :param bool backbone: If :data:`True`, use full backbone atoms;
#         otherwise use only C-alpha.
#
#     :return: :class:`~pandas.DataFrame`
#
#     :raises:
#         :ValueError: If ``pdb_path`` column is not present.
#
#     Column modifications:
#     =========== ========================================= =======
#     column name description                               status
#     =========== ========================================= =======
#     vectors     List of fitting vectors to the SSE.       new!
#     =========== ========================================= =======
#     """
#     if 'pdb_path' not in df.columns:
#         raise ValueError('The pdb_path column is necessary. Get it with src.io.get_pdbs')
#
#     def row_vector(row, backbone):
#         pdb3d = PDB(row['pdb_path'])
#         pieces = []
#         backbone = 'AtomTask:PROTEINBACKBONE' if backbone else 'AtomType:CA'
#         for segment in row['match']:
#             with SBIcr.on_option_value('structure', 'source', 'label'):
#                 piece = pdb3d['Residue:{0}-{1}'.format(segment[0], segment[1])][backbone]
#                 pieces.append(fit_vector(piece.coordinates))
#         if core.get_option('imaster', 'verbose'):
#             print('  TOTAL: {} vectors generated'.format(len(pieces)))
#         return pieces
#     if core.get_option('imaster', 'verbose'):
#         print('Obtaining vectors for SSE')
#     return df.assign(vectors=df.apply(lambda row: row_vector(row, backbone), axis=1))
#
#
# def fit_plane_to_data(df: pd.DataFrame,
#                       backbone: bool = False,
#                       order: int = 1,
#                       selections: Optional[List] = None,
#                       ) -> pd.DataFrame:
#     """Add fitting plane to secondary structure segments.
#
#     :param df: Data derived from the MASTER search. Requires the
#         ``pdb_path`` column.
#     :type df: :class:`~pandas.DataFrame`
#     :param bool backbone: If :data:`True`, use full backbone atoms;
#         otherwise use only C-alpha.
#     :param int order: 1 for linear fitting, 2 for quadratic.
#     :param selections: If provided, only use the given SSE for fitting.
#         **Count the structures starting at 0!**.
#     :type selections: :class:`List`[:class:`List`]
#
#     :return: :class:`~pandas.DataFrame`
#
#     :raises:
#         :ValueError: If ``pdb_path`` column is not present.
#
#     Column modifications:
#     =========== ========================================== =======
#     column name description                               status
#     =========== ========================================== =======
#     planes      List of :class:`tuple` with fitting planes new!
#                 and edges to the SSE.
#     =========== ========================================== =======
#     """
#     if 'pdb_path' not in df.columns:
#         raise ValueError('The pdb_path column is necessary.')
#
#     def row_plane(row, backbone, order, selections):
#         pdb3d = PDB(row['pdb_path'])
#         pieces = []
#         backbone = 'AtomTask:PROTEINBACKBONE' if backbone else 'AtomType:CA'
#         for segment in row['match']:
#             with SBIcr.on_option_value('structure', 'source', 'label'):
#                 piece = pdb3d['Residue:{0}-{1}'.format(segment[0], segment[1])][backbone]
#                 pieces.append(piece)
#         if selections is None:
#             selections = [list(range(len(pieces)))]
#         if core.get_option('imaster', 'verbose'):
#             print('  planes generated from selections {}'.format(selections))
#         planes = []
#         # pc=pd.concat(pieces)
#         # for e in pc.eigenvectors(30):
#         #     print([list(x) for x in e])
#         for sele in selections:
#             planes.append(fit_plane(pd.concat(itemgetter(*sele)(pieces)).coordinates, order))
#         if core.get_option('imaster', 'verbose'):
#             print('  {} planes generated'.format(len(planes)))
#         return planes
#     if core.get_option('imaster', 'verbose'):
#         print('Obtaining planes for selected SSE')
#     return df.assign(planes=df.apply(lambda row: row_plane(row, backbone, order, selections), axis=1))
#
#
# def fit_vector(coordinates: np.ndarray) -> sy.geometry.Line:
#     """
#     Fits a vector over a set of 3D points.
#
#     Adapted from https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
#     """
#     # Calculate the mean of the points, i.e. the 'center' of the cloud
#     datamean = coordinates.mean(axis=0)
#
#     # Do an SVD on the mean-centered data.
#     uu, dd, vv = np.linalg.svd(coordinates - datamean)
#
#     # Now vv[0] contains the first principal component, i.e. the direction
#     # vector of the 'best fit' line in the least squares sense.
#
#     # Now generate some points along this best fit line, for plotting.
#     linepts = vv[0] * np.mgrid[-20:20:2j][:, np.newaxis]
#
#     # shift by the mean to get the line in the right place
#     linepts += datamean
#     return sy.geometry.Line(*linepts)
#
#
# def fit_plane(coordinates: np.ndarray,
#               order: int = 1
#               ) -> sy.geometry.Plane:
#     """Fits a plane in either 1st or 2nd order over a set of 3D points.
#
#     Adapted from https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
#     """
#     # Compute Meshgrid
#     x = coordinates[:, 0]
#     y = coordinates[:, 1]
#     z = coordinates[:, 2]
#
#     # Find edges to trim meshgrid
#     x_max, x_min, y_max, y_min = max(x), min(x), max(y), min(y)
#
#     # Create regular grid covering the domain of the data with some addition
#     X, Y = np.meshgrid(np.arange(x_min - 10, x_max + 10, 1.0),
#                        np.arange(y_min - 10, y_max + 10, 1.0))
#     XX = X.flatten()
#     YY = Y.flatten()
#
#     # For 2nd order plane
#     xy = np.column_stack((x, y))
#
#     # Chose order of approximation plane fitting
#     # 1: linear, 2: quadratic
#     order = 1
#     if order == 1:
#         # Best-fit linear plane
#         A = np.c_[x, y, np.ones(x.shape[0])]
#
#         # Find coefficients
#         C, _, _, _ = sc.linalg.lstsq(A, z)
#
#         # Evaluate it on grid
#         Z = C[0] * X + C[1] * Y + C[2]
#
#     elif order == 2:
#         # best-fit quadratic curve
#         A = np.c_[np.ones(x.shape[0]), xy, np.prod(xy, axis=1), xy**2]
#         C, _, _, _ = sc.linalg.lstsq(A, z)
#
#         # Evaluate it on a grid
#         Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX * YY, XX ** 2, YY ** 2], C).reshape(X.shape)
#
#     # For plane creation
#     edge0 = np.array([X[0][0], Y[0][0], Z[0][0]]).flatten()
#     edge1 = np.array([X[0][-1:], Y[0][-1:], Z[0][-1:]]).flatten()
#     edge2 = np.array([X[-1][-1:], Y[-1][-1:], Z[-1][-1:]]).flatten()
#     edges = np.array([edge0, edge1, edge2])
#
#     return sy.Plane(sy.Point3D(edges[0]), sy.Point3D(edges[1]), sy.Point3D(edges[2])), edges
#
#
# def vectors2planes( df: pd.DataFrame ) -> pd.DataFrame:
#     """Adds angle between each :class:`~scipy.Line` and each :class:`~scipy.Plane`.
#
#     :param df: Data derived from the MASTER search. Requires the
#         ``vectors`` and the ``planes`` columns.
#     :type df: :class:`~pandas.DataFrame`
#
#     :return: :class:`~pandas.DataFrame`
#
#     :raises:
#         :ValueError: If ``vectors`` and ``planes`` columns are not present.
#
#     Column modifications:
#
#     =========== ========================================== =======
#     column name description                                status
#     =========== ========================================== =======
#     angles      List of angles between vectors and planes. new!
#     =========== ========================================== =======
#     """
#     if 'vectors' not in df.columns or 'planes' not in df.columns:
#         raise ValueError('The vectors and plane columns are necessary.')
#
#     def vector2plane( row ):
#         angles = []
#         for p in row['planes']:
#             angles.append([])
#             for v in row['vectors']:
#                 angles[-1].append(math.degrees(p[0].angle_between(v)))
#         return angles
#
#     if core.get_option('imaster', 'verbose'):
#         print('Calculating vector vs. plane angles.')
#     return df.assign(angles=df.apply(vector2plane, axis=1))
#
#
# def pymol_commands( df: pd.DataFrame ) -> pd.DataFrame:
#     """Adds the PyMOL commands to draw the vectors and plane.
#
#     :param df: Data derived from the MASTER search. Requires the
#         ``vectors`` and the ``planes`` columns.
#     :type df: :class:`~pandas.DataFrame`
#
#     :return: :class:`~pandas.DataFrame`
#
#     :raises:
#         :ValueError: If ``vectors`` and ``planes`` columns are not present.
#
#     Column modifications:
#
#     =========== ======================= =======
#     column name description             status
#     =========== ======================= =======
#     pymol       List of PyMOL commands. new!
#     =========== ======================= =======
#     """
#     if 'vectors' not in df.columns or 'planes' not in df.columns:
#         raise ValueError('The vectors and plane columns are necessary.')
#
#     def pymolcom( row ):
#         comm = []
#         dots = [list(x) for x in row['planes'][0][1]]
#         # dots = []
#
#         for i, v in enumerate(row['vectors']):
#             a = list(v.p1.evalf())
#             b = list(v.p2.evalf())
#             # if i == 0 or i == len(row['vectors']) - 1:
#             #     dots.extend([list(plane.projection(a).evalf()), list(plane.projection(b).evalf())])
#             comm.append('make_arrow vector{0}, {1}, {2}'.format(i + 1, a, b))
#         comm.append('make_plane_points plane, {0}, {1}, {2}'.format(*dots))
#         return comm
#
#     return df.assign(pymol=df.apply(pymolcom, axis=1))
