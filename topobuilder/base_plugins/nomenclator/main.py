# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union
import sys
import copy

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore

__all__ = ['apply']

_PLT_TYPES_ = ['sketchXZ', 'sketchXY']


def apply( cases: List[Case],
           prtid: int,
           subnames: Union[str, List[str]],
           **kwargs ) -> List[Case]:
    """Add subnames to the case.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: NOMENCLATOR ---\n')

    if not isinstance(subnames, list):
        subnames = [subnames, ]

    for i, case in enumerate(cases):
        sn = copy.deepcopy(subnames)
        sn.insert(0, case.name)
        cases[i].data['configuration']['name'] = '_'.join(sn)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases
