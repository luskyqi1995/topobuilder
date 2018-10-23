# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import argparse
import os
import shutil
from collections import OrderedDict

# External Libraries


# This Library
from topobuilder.io import read_case, write_case

__all__ = ['build']


def build( case: str, overwrite: bool = False ):
    """
    """
    # 1. Load case
    data = read_case(case, True)

    # 2. Create output working directory
    wdir = case['description']['name']
    if os.path.isdir(wdir):
        if not overwrite:
            raise IOError('Unable to overwrite output directory {}'.format(wdir))
        else:
            shutil.rmtree(wdir)
    os.mkdir(wdir)
    write_case(data, os.path.join(wdir, wdir), format)

    # 2.


def cli_build():
    """
    """
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_case_template.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-case', dest='case', action='store',
                        help='Case file.', required=True)
    parser.add_argument('-overwrite', dest='overwrite', action='store_true',
                        default=False, help='Ignore previously existing data. '
                        'Will delete previously existing work folder.')
    options = parser.parse_args()

    build(**vars(options))
