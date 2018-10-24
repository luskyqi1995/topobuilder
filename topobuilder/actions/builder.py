# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: build
"""
# Standard Libraries
import argparse
import os
import shutil
from collections import OrderedDict

# External Libraries


# This Library
from topobuilder.io import read_case, write_case
from topobuilder.io import setup_build
from topobuilder.coordinates import GeneralArchitect

__all__ = ['build']


def build( case: str, overwrite: bool = False ):
    """
    """
    # 1. Load case and make sure it is absolute.
    data = read_case(case, True)

    # 2. Create output working directory tree
    paths = setup_build(data, overwrite)

    # 3. Generate Sketch
    arch = GeneralArchitect(data, paths)
    arch.build_sketch()

