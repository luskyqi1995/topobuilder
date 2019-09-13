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
import shlex
from pathlib import Path
from typing import Optional, Tuple, List, Union
from ast import literal_eval
from tempfile import NamedTemporaryFile

# External Libraries
import pandas as pd
import numpy as np

# This Library
import topobuilder.core as TBcore
from .plugins import plugin_filemaker

__all__ = ['createPDS', 'master_best_each', 'pds_database']

pds_file = None
pds_list = None


def get_master_exes() -> Tuple[str, str]:
    """Provide the path to the ``MASTER`` executables.

    .. note::
        Depends on the ``master.master`` configuration option
        Depends on the ``master.create`` configuration option

    :return: The ``master`` and ``createPDS`` executable paths. In this order.
    """
    exes = []
    for name, option_name in [('master', 'master'), ('createPDS', 'create')]:
        # Get the executable and check that it works
        exe = TBcore.get_option('master', option_name)
        if not os.access(exe, os.X_OK):
            raise IOError('Unable to find a proper executable for {} at {}'.format(name, exe))
        exes.append(exe)
    return exes


def pds_database( force: Optional[bool] = False ) -> Tuple[Path, List]:
    """Provide the list of available PDS as a file and a list.

    :param force: When :data:`.True`, recheck the PDS database even if one is
        already assigned.

    .. note::
        Depends onf the ``master.pds`` configuration option
    """
    global pds_file
    global pds_list

    if not force:
        if pds_file is not None and pds_list is not None:
            return pds_file, pds_list

    pds_file = TBcore.get_option('master', 'pds')
    pds_list = []
    pds_file = Path(pds_file)
    if pds_file.is_file():
        pds_list = [line.strip() for line in open(pds_file).readlines() if len(line.strip()) > 0]
        return pds_file, pds_list
    elif pds_file.is_dir():
        pds_list = [str(x.resolve()) for x in pds_file.glob('*/*.pds')]
        f = NamedTemporaryFile(mode='w', delete=False)
        plugin_filemaker('Temporary file for PDS database: {}'.format(f.name))
        [f.write(x + '\n') for x in pds_list]
        f.close()
        pds_file = Path(f.name)
        return pds_file, pds_list
    else:
        raise ValueError('The provided MASTER database directory/list file cannot be found.')


def createPDS( infile: Union[Path, str], outfile: Optional[str] = None ) -> List[str]:
    """Make the createPDS command call.

    .. note::
        Depends onf the ``master.create`` configuration option
    """
    _, createPDS = get_master_exes()
    createbash = '{0} --type query --pdb {1} --pds {2}'
    infile = Path(infile)
    if not infile.is_file():
        raise IOError('Unable to find structure file {}'.format(infile))
    outfile = outfile if outfile is not None else infile.with_suffix('.pds')
    return shlex.split(createbash.format(createPDS, infile, outfile))


def master_best_each( infile: Union[Path, str],
                      outdir: Union[Path, str],
                      rmsd: Optional[float] = 5.0
                      ) -> List[List[str]]:
    """Create one MASTER call for each PDS file with the --topN 1 flag.

    .. note::
        Depends onf the ``master.master`` configuration option
        Depends onf the ``master.pds`` configuration option
    """
    master, _ = get_master_exes()
    _, pds_list = pds_database()

    infile = Path(infile)
    if not infile.is_file():
        raise IOError('Unable to find PDS file {}'.format(infile))
    outdir = Path(outdir)
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)
    createbash = '{0} --query {1} --target {2} --rmsdCut {3}  --topN 1 --matchOut {4}'

    cmds = []
    for pds in pds_list:
        tid = str(Path(Path(pds).name).with_suffix(''))
        outfile = outdir.joinpath('{}.master'.format(tid))
        cmds.append(shlex.split(createbash.format(master, infile, pds, rmsd, outfile)))
    return cmds
