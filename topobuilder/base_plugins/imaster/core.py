# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import shutil

# External Libraries
from libconfig import *

# This Library


with ifndef():
    register_option('imaster', 'pymol', False, 'bool', 'Print PyMOL commands.')
    # register_option('imaster', 'sample', 200, 'int', 'Number of structures to sample per partition. 0 or lower means all.')
    register_option('imaster', 'master', shutil.which('master'), 'path_in', 'MASTER executable.')
    register_option('imaster', 'createPDS', shutil.which('createPDS'), 'path_in', 'createPDS executable.')
    register_option('imaster', 'pds', None, 'path_in', 'Local PDS database.')
    register_option('imaster', 'pdb', None, 'path_in', 'Local PDB database.')
    register_option('imaster', 'slurm.partition', 'serial', 'string', 'Name of the available SLURM partition.')
    register_option('imaster', 'slurm.array', 700, 'int', 'Into how may nodes is the search splitted.')
    register_option('imaster', 'slurm.logs', os.getcwd(), 'path_in', 'Path on were to dump the log files.')

    config_file = get_local_config_file('.topobuilder.cfg')
    # Either make or read from the file.
    if config_file is not None:
        set_options_from_YAML( config_file )

for name in user_forbidden:
    del globals()[name]
