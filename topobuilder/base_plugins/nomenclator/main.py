# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict
import sys
import copy

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil


__all__ = ['metadata', 'apply', 'case_apply']


def metadata() -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: Funtion on the expected input connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    """
    def isngl( count ):
        return True

    return {'name': 'nomenclator',
            'Itags': [],
            'Otags': [],
            'Isngl': isngl,
            'Osngl': False}


def apply( cases: List[Case],
           prtid: int,
           subnames: Union[str, List[str]],
           **kwargs ) -> List[Case]:
    """Add subnames to a list of Case, defining the folder configuration.

    :param cases: List of :class:`.Case` to execute.

    :return: :func:`list` of :class:`.Case`
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i] = case_apply(case, subnames)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


@TButil.plugin_conditions(metadata())
def case_apply( case: Case,
                subnames: Union[str, List[str]]
                ) -> Case:
    """Add subnames to a Case, defining the folder configuration.

    :param case: Target :class:`.Case`.
    :param subnames: List of or subname to append.

    :return: :class:`.Case` with ``configuration.name`` modified
    """
    kase = Case(case)

    # Unify subnames behaviour for 1 to N
    if not isinstance(subnames, list):
        subnames = [subnames, ]

    # Reserved keywords
    for i, sn in enumerate(subnames):
        if sn == 'architecture':
            subnames[i] = kase.architecture_str.replace('.', '')

    # Check name was not already added.
    sn = copy.deepcopy(subnames)
    if kase.name.endswith('_'.join(sn)):
        TButil.plugin_warning('Seems the subnames {} already existed.'.format('_'.join(sn)))
        TButil.plugin_warning('Will NOT re-append.')
        return kase

    # Add new names
    sn.insert(0, kase.name)
    kase.data['configuration']['name'] = '_'.join(sn)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Renamed case {} to {}\n'.format(case.name, kase.name))

    return kase
