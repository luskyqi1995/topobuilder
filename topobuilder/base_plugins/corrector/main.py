# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict, Optional
from pathlib import Path
import copy
import sys

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
           corrections: Optional[Union[str, Dict, Path, List]] = None,
           **kwargs ) -> List[Case]:
    """Apply corrections to the Case.
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i] = case_apply(case, corrections)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


@TButil.plugin_conditions(metadata())
def case_apply( case: Case,
                corrections: Optional[Union[str, Dict, Path, List]] = None
                ) -> Case:
    """
    """
    kase = Case(case)

    # Make sure we have a list of corrections.
    if corrections is None and case['metadata.corrections'] is None:
        return kase
    if corrections is not None:
        if not isinstance(corrections, list):
            corrections = [corrections, ]
        for i, c in enumerate(corrections):
            if isinstance(c, str):
                corrections[i] = Path(c)
    else:
        corrections = []

    crr = copy.deepcopy(corrections)
    krr = case['metadata.corrections']
    if krr is not None:
        crr.extend(krr)

    for c in crr:
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Applying correction: {0}\n'.format(c))
        kase = kase.apply_corrections(c)
    return kase
