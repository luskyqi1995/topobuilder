# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. obj:: CaseSchema
"""

from . import io
from . import case

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
