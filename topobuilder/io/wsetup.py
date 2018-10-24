# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from collections import OrderedDict
from pathlib import Path
import shutil

# External Libraries


# This Library
from topobuilder.io import write_case

__all__ = ['setup_build']


def setup_build( case: OrderedDict, overwrite: bool = False ) -> dict:
    """Generate the folder structure needed to proceed with a building experiment.

    :param OrderedDict case: data of the :class:`.CaseSchema`.
    :param bool overwrite: If :data:`True`, delete previous content of the target
        directory and start again from scratch.

    :return: :class:`dict` of :class:`Path` to the different places of interest.
    """
    name = case['configuration']['name']

    paths = {'wdir': Path(name)}

    if overwrite:
        shutil.rmtree(paths['wdir'], True)

    paths.setdefault('architecture', paths['wdir'].joinpath('architecture'))
    paths['architecture'].mkdir(parents=True, exist_ok=overwrite)
    paths.setdefault('arch_sketch', paths['architecture'].joinpath('sketch.pdb'))

    paths.setdefault('casefile', paths['wdir'].joinpath(name))
    write_case(case, str(paths['casefile']), format='yaml')

    return paths
