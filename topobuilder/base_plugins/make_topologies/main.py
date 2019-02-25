# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List
import sys

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore

__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           **kwargs ) -> List[Case]:
    """
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: MAKE_TOPOLOGIES ---\n')

    new_cases = []
    for case in cases:
        new_cases.extend(case.apply_topologies())

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('From initial {0} architectures, {1} topologies were generated\n'.format(len(cases), len(new_cases)))

    for i, _ in enumerate(new_cases):
        new_cases[i].set_protocol_done(prtid)
    return new_cases
