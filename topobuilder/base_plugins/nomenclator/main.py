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
import topobuilder.utils as TButil

__all__ = ['apply', 'apply_case']


def apply( cases: List[Case],
           prtid: int,
           subnames: Union[str, List[str]],
           **kwargs ) -> List[Case]:
    """Add subnames to a list of Case.
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i] = apply_case(case, subnames)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


def apply_case( case: Case,
                subnames: Union[str, List[str]]
                ) -> Case:
    """Add subnames to a Case.

    :param case: Target :class:`.Case`.
    :param subnames: List of or subname to append.

    :return: :class:`.Case` with ``configuration.name`` modified
    """
    c = Case(case)

    # Unify subnames behaviour for 1 to N
    if not isinstance(subnames, list):
        subnames = [subnames, ]

    # Check name was not already added.
    sn = copy.deepcopy(subnames)
    if case.name.endswith('_'.join(sn)):
        TButil.plugin_warning('Seems the subnames {} already existed.'.format('_'.join(sn)))
        TButil.plugin_warning('Will NOT re-append.')
        return c

    # Add new names
    sn.insert(0, case.name)
    c.data['configuration']['name'] = '_'.join(sn)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Renamed case {} to {}\n'.format(case.name, c.name))
    return c
