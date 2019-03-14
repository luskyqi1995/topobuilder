# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys
import os
from typing import Optional
from ast import literal_eval

# External Libraries
import pandas as pd
import numpy as np

# This Library
import topobuilder.core as TBcore

__all__ = ['parse_master_file']


def parse_master_file( filename: str,
                       max_rmsd: Optional[float] = None,
                       piece_count: Optional[int] = 18,
                       shift_0: Optional[bool] = False
                       ) -> pd.DataFrame:
    """Load MASTER output data.

    :param str filename: Output file.
    :param float max_rmsd: Maximum RMSD value to recover.
    :param int piece_count: If known, specify the number of structural pieces
        in the MASTER search. Otherwise, it is assumed.
    :param bool shift_0: MASTER matches start counting in 0. If the flag
        is :data:`.True`, change it to start with 1.

    :return: :class:`~pandas.DataFrame`

    Columns of the returned :class:`~pandas.DataFrame` are:

    =========== ===================================================
    column name description
    =========== ===================================================
    rmsd        RMSD value between query and match
    pds_path    Path to the PDS file containing the match
    pdb         PDB identifier
    chain       Chain identifier
    match       List of ranges for the match (Rosetta count)
    =========== ===================================================

    This assumes that the PDS files basename has the standard nomenclature
    ``<pdbid>_<chain>.pds``.
    """
    def shift(x):
        return list(np.asarray(x) + 1)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Reading MASTER file {}\n'.format(filename))
    df = pd.read_csv(filename,
                     names=list(range(piece_count + 2)), engine='python',
                     sep=r'\s+', header=None).dropna(axis=1, how='all')
    df['match'] = df[df.columns[2:]].astype(str).sum(axis=1).apply(literal_eval)
    if shift_0:
        df['match'] = df['match'].apply(shift)
    df = df.rename(columns={0: 'rmsd', 1: 'pds_path'})
    df = df.drop(columns=[i for i in df.columns if isinstance(i, int)])
    df[['pdb', 'chain']] = (pd.DataFrame(list(df['pds_path'].str.replace('.pds', '')
                            .apply(lambda x: os.path.basename(x).split('_')).values)))
    if max_rmsd is not None:
        df = df[(df['rmsd'] <= max_rmsd)]
    return df[['rmsd', 'pds_path', 'pdb', 'chain', 'match']]
