# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys
import argparse
from pathlib import Path
from itertools import cycle

# External Libraries
import pandas as pd

# This Library
import topobuilder.utils as TButil
import topobuilder.core as TBcore
from topobuilder.case import Case


def options():
    """
    """
    # Parse Arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-case', dest='case', action='store', required=True)
    parser.add_argument('-indir', dest='indir', action='store', required=True)
    parser.add_argument('-out', dest='out', action='store', required=True)

    return parser.parse_args()


def main( options ):
    """
    """
    data = []
    case = Case(Path(options.case))
    # Get connectivities
    sse = case.connectivities_str[0].split('.')
    # Get flips
    flip = cycle([case['configuration.flip_first'], not case['configuration.flip_first']])
    flip = [next(flip) for _ in range(len(sse))]
    # Get ranges
    loops = case['metadata.loop_lengths']
    lsses = [x['length'] for x in case.ordered_structures]
    ranges = []
    start = 1
    for i, s in enumerate(lsses):
        ranges.append([start, start + s - 1])
        start += s
        if i < len(loops):
            start += (loops[i])

    rules = list(zip(sse, ranges, flip))
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Applying rules:\n')
        sys.stdout.write(','.join(sse) + '\n')
        sys.stdout.write(','.join(['{}-{}'.format(x, y) for x, y in ranges]) + '\n')
        sys.stdout.write(','.join([str(x) for x in flip]) + '\n')
        sys.stdout.write('\n')

    for pdbf in Path(options.indir).glob('*pdb'):
        data.append(TButil.pdb_geometry_from_rules(pdbf, rules))
        sys.stdout.flush()

    df = pd.concat(data)
    df.to_csv(options.out, index=False)


if __name__ == '__main__':
    main(options())
