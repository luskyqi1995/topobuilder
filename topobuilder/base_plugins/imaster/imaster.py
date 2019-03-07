# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import argparse
from pathlib import Path

# External Libraries

# This Library
from analysis import parse_master_file, geometric_properties, geometric_stats
from core import core


def options():
    """
    """
    # Parse Arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-in', dest='input', action='store', required=True)
    parser.add_argument('-out', dest='out', action='store', required=True)
    parser.add_argument('-mode', dest='mode', action='store', default='master', choices=['master', 'stats'])
    parser.add_argument('-selection', dest='selection', action='store', default=None)
    parser.add_argument('-directionality', dest='directionality', action='store', required=True)
    parser.add_argument('-pymol', dest='pymol', action='store_true', default=False)
    parser.add_argument('-planepick', dest='planepick', action='store', nargs='+', default=[])

    options = parser.parse_args()

    if options.mode == 'master' and not Path(options.input).is_file():
        raise IOError('Unable to find MASTER file {}.'.format(options.input))
    if options.mode == 'stats' and not Path(options.input).is_dir():
        raise IOError('Unable to find PDB dir {}.'.format(options.input))

    if options.mode == 'stats' and options.selection is None:
        raise AttributeError('In stats mode, a selection must be provided.')

    if core.get_option('imaster', 'pymol') or options.pymol:
        core.set_option('imaster', 'pymol', True)
        options.pymol = options.out + '.pymol'
    else:
        core.set_option('imaster', 'pymol', False)
        options.pymol = None

    if len(options.planepick) == 0:
        options.planepick = None
    else:
        options.planepick = [int(x) for x in options.planepick]
        options.planepick = [options.planepick, ]

    return options


def main( options ):
    """
    """
    if options.mode == 'master':
        # Load MASTER search data.
        masterdf = parse_master_file(options.input)
        # Geometric properties retrieval
        masterdf = geometric_properties(masterdf, options.planepick, options.directionality, options.pymol)
    if options.mode == 'stats':
        masterdf = geometric_stats(options.input, options.planepick, options.directionality, options.selection)
    # Output data
    masterdf.to_csv(options.out + '.csv', index=False)


if __name__ == '__main__':
    main(options())
