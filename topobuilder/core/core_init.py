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
    # Register IO control options
    register_option('topobuilder', 'verbose', False, 'bool', 'Makes topobuilder chatty.')

    config_file = get_local_config_file('.topobuilder.cfg')
    # Either make or read from the file.
    if not os.path.isfile(config_file):
        write_options_to_YAML( config_file )
    else:
        set_options_from_YAML( config_file )

for name in user_forbidden:
    del globals()[name]
