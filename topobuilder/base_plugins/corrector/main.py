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

__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           corrections: Optional[Union[str, Dict, Path, List]] = None,
           **kwargs ) -> List[Case]:
    """Apply corrections to the Case.
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i] = apply_case(case, corrections)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def apply_case( case: Case,
                corrections: Optional[Union[str, Dict, Path, List]] = None
                ) -> Case:
    """
    """
    kase = Case(case)

    # Make sure we have a list of corrections.
    if corrections is None:
        return kase
    if not isinstance(corrections, list):
        corrections = [corrections, ]
    for i, c in enumerate(corrections):
        if isinstance(c, str):
            corrections[i] = Path(c)

    crr = copy.deepcopy(corrections)
    krr = case['metadata.corrections']
    if krr is not None:
        crr.extend(krr)

    for c in crr:
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Applying correction: {0}\n'.format(c))
        kase = kase.apply_corrections(c)
    return kase
