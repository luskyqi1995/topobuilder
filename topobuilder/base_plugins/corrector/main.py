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

__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           corrections: Optional[Union[str, Dict, Path, List]] = None,
           **kwargs ) -> List[Case]:
    """Apply corrections to the Case.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: CORRECTOR ---\n')

    if corrections is None:
        corrections = []
    if not isinstance(corrections, list):
        corrections = [corrections, ]
    for i, c in enumerate(corrections):
        if isinstance(c, str):
            corrections[i] = Path(c)

    new_cases = []
    for case in cases:
        crr = copy.deepcopy(corrections)
        krr = case['metadata.corrections']
        if krr is not None:
            crr.extend(krr)
        kase = Case(case)
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Applying {0} corrections to case {1}\n'.format(len(crr), kase.name))
        for c in crr:
            if TBcore.get_option('system', 'verbose'):
                sys.stdout.write('Applying correction: {0}\n'.format(c))
            kase = kase.apply_corrections(c)
        new_cases.append(kase)
        new_cases[-1].set_protocol_done(prtid)
    return new_cases
