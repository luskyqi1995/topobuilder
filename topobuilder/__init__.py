# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. obj:: CaseSchema
"""
# Standard Libraries
import os
from pathlib import Path

# External Libraries
from pluginbase import PluginBase

# This Library
from .core import core
from . import case
from . import utils

# Plugin System
# Load all present plugins
tp_plugin_dirs = [str(Path(__file__).parent.joinpath('base_plugins')),
                  str(Path.home().joinpath('.topobuilder'))]
try:
    tp_plugin_dirs.append(os.environ['TOPOBULDERPLUGINS'])
except KeyError:
    pass
plugin_base = PluginBase(package='topobuilder.plugins')
plugin_source = plugin_base.make_plugin_source(searchpath=tp_plugin_dirs)

# This Library
from . import interface

