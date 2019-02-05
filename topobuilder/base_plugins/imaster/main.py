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

# External Libraries

# This Library
from topobuilder.case import Case
from .analysis import get_steps
from .slurm import make_slurm_file
from . import core


__all__ = ['requirements', 'apply']


def requirements( **kwargs ) -> dict:
    """Extra information that needs to be provided in order to run the ``imaster plugin``.

    :param str database: Path to the root directory of the MASTER database or to a file listing the PDS.
        If it does not exist a ``pdb database`` for each master match will be created.
    :param str master: Full path of the MASTER executable.
    :param str createPDS: Full path of the createPDS executable.

    :param str slurm_logs: Path to the directory where the SLURM log files should be placed.
        If not provided defaults to ``/scratch/$USER/logs``.
    :param str slurm_partition: Name of the partition to use in SLURM. If not provided defaults to ``serial``.
    :param int slurm_array: Number of splits for cluster submission. If not provided defaults to ``700``.

    :param steps: Steps in which the construction/correction loops will progress. Each position contains 2 entries:
        a :ref:`topology` definition and a tuple of the structures to use for the plane (all if :data:`None` are defined).
    :type steps: :class:`List`[:class:`Tuple`[:class:`str`, Optional[:class:`Tuple`]]]

    :return: :class:`dict`
    """
    # DATABASE
    if 'database' not in kwargs:
        raise KeyError('A path pointing to the root folder of the MASTER database or to a list file is required.')
    if not Path(kwargs['database']).is_dir() and not Path(kwargs['database']).is_file():
        raise ValueError('The provided MASTER database/list file cannot be found.')

    # MASTER AND CREATEPDS
    if 'master' not in kwargs:
        kwargs['master'] = shutil.which('master')
        if shutil.which('master') is None:
            raise KeyError('Full path of the executable MASTER command needs to be provided.')
    if not Path(kwargs['master']).is_file():
        raise ValueError('Provided MASTER path cannot be found.')
    if 'createPDS' not in kwargs:
        kwargs['createPDS'] = shutil.which('createPDS')
        if shutil.which('createPDS') is None:
            raise KeyError('Full path of the executable createPDS command needs to be provided.')
    if not Path(kwargs['createPDS']).is_file():
        raise ValueError('Provided createPDS path cannot be found.')

    # SLURM
    if 'slurm_logs' not in kwargs:
        try:
            kwargs['slurm_logs'] = str(Path(os.environ['SCRATCH']).joinpath('logs'))
        except KeyError:
            kwargs['slurm_logs'] = str(Path('/scratch').joinpath(getpass.getuser()).joinpath('logs'))
    if 'slurm_partition' not in kwargs:
        kwargs['slurm_partition'] = 'serial'
    if 'slurm_array' not in kwargs:
        kwargs['slurm_array'] = 700

    # STEPS
    if 'steps' not in kwargs:
        raise KeyError('Steps for the interactive master calculation need to be specified.')

    return kwargs


def apply( cases: List[Case], prtid: int, **kwargs ) -> List[Case]:
    """
    """
    # Get list of PDS structures
    database = core.get_option('imaster', 'pds')
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
        cases[i].set_protocol_done(prtid)
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
    print(steps)


    # Make slurm file template
    # print(make_slurm_file(pds_list))
    return case

# apply(Case('test'), '/Volumes/MiniTwo/data/TopoBuilderData/database/master', 'master', 'createPDS', '.', 'serial', 700)
