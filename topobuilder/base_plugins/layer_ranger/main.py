# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict
import copy
import sys

# External Libraries
import numpy as np

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore

__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           ranger: Dict,
           **kwargs ) -> List[Case]:
    """
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: LAYER_RANGER ---\n')

    expected = []
    for layer in sorted(ranger):
        a, b = ranger[layer]
        expected.append(b - a + 1)
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Applying ranges {0}-{1} to layer {2}\n'.format(a, b, layer))
    expected = np.prod(expected)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Starting from {} cases\n'.format(len(cases)))

    for case in cases:
        if case.connectivity_count != 0:
            raise ValueError('Cannot apply layer_ranger to topologies with defined connectivity')

    ncases = len(cases)
    new_cases = copy.deepcopy(cases)
    for ilayer, layer in enumerate(sorted(ranger.keys())):
        lcases = []
        ncases = len(new_cases)
        for i in range(ranger[layer][0], ranger[layer][1] + 1):
            for j in range(len(new_cases)):
                lcases.append(new_cases[j].set_type_for_layer(layer, i))
        new_cases = new_cases[ncases:]
        new_cases.extend(lcases)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Generated a total of {} cases\n'.format(len(new_cases)))
    return new_cases
