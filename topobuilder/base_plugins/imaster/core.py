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
from libconfig import Config

# This Library

core = Config()
with core.ifndef():
    core.register_option('imaster', 'pymol', False, 'bool', 'Print PyMOL commands.')
    # core.register_option('imaster', 'sample', 200, 'int', 'Number of structures to sample per partition. 0 or lower means all.')
    core.register_option('master', 'master', shutil.which('master'), 'path_in', 'MASTER executable.')
    core.register_option('master', 'create', shutil.which('createPDS'), 'path_in', 'createPDS executable.')
    core.register_option('master', 'pds', None, 'path_in', 'Local PDS database.')
    core.register_option('master', 'pdb', None, 'path_in', 'Local PDB database.')
    core.register_option('slurm', 'use', True, 'bool', 'Use SLURM cluster submission system.')
    core.register_option('slurm', 'partition', 'serial', 'string', 'Name of the available SLURM partition.')
    core.register_option('slurm', 'array', 700, 'int', 'Into how may nodes is the search splitted.')
    core.register_option('slurm', 'logs', os.getcwd(), 'path_in', 'Path on were to dump the log files.')

    # There are different levels of configuration files that can be picked.
    # If any configuration file is set up, the priority goes as follows:
    #   1) Local config file (in the actual executable directory)
    #   2) Root of the current working repository (if any)
    #   3) User's home path
    config_file = core.get_local_config_file('.topobuilder.cfg')
    if config_file is not None:
        core.set_options_from_YAML( config_file )

    core.lock_configuration()
