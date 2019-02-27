# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries

# External Libraries
from libconfig import Config
core = Config()

# This Library

with core.ifndef():
    # Register IO control options
    core.register_option('system', 'verbose', False, 'bool', 'Makes topobuilder chatty.')
    core.register_option('system', 'debug', False, 'bool', 'Makes topobuilder VERY chatty.')
    core.register_option('system', 'overwrite', False, 'bool', 'Overwrite existing files.')

    # There are different levels of configuration files that can be picked.
    # If any configuration file is set up, the priority goes as follows:
    #   1) Local config file (in the actual executable directory)
    #   2) Root of the current working repository (if any)
    #   3) User's home path
    config_file = core.get_local_config_file('.topobuilder.cfg')
    if config_file is not None:
        core.set_options_from_YAML( config_file )

    core.lock_configuration()
