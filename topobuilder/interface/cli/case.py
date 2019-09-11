# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: cli_case_template
.. func:: cli_absolute_case
"""
# Standard Libraries
import os
import argparse
from pathlib import Path

# External Libraries

# This Library
from topobuilder.case import Case, case_template
from topobuilder import plugin_source
import topobuilder.utils as TButil


def cli_case_template():
    """Generate a :class:`.Case`.
    """
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_case_template.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-name', dest='name', action='store',
                        help='Job Name.', required=True)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-architecture', dest='architecture', action='store',
                       help='Architecture string definition.', default=None)
    group.add_argument('-topology', dest='topology', action='store',
                       help='Topology string definition.', default=None)

    parser.add_argument('-format', dest='format', action='store', choices=['json', 'yaml'],
                        help='Format for the case file.',
                        default='yaml')

    options = parser.parse_args()

    _, outfile = case_template(**vars(options))
    TButil.plugin_filemaker('New case file created at: {}'.format(os.path.abspath(outfile)))


def cli_absolute_case():
    """Transform a relative :class:`.Case` to absolute.
    """
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_case_template.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-case', dest='case', action='store',
                        help='Relative case file.', required=True)
    parser.add_argument('-corrections', dest='corrections', action='store', nargs='+',
                        help='Correction file for the case.', default=None)
    parser.add_argument('-caseout', dest='caseout', action='store',
                        help='Output filename (optional).', default=None)
    options = parser.parse_args()

    # Process naming system
    prefix = options.case.split('.')
    format = 'yaml' if prefix[-1] == 'yml' else 'json'
    if options.caseout is None:
        prefix[-1] = 'absolute'
        options.caseout = '.'.join(prefix)

    # Read, transform and write
    case = plugin_source.load_plugin('corrector').apply([Case(Path(options.case)), ], -1, options.corrections)[0]
    case = case.cast_absolute()
    outfile = case.write(options.caseout, format)
    TButil.plugin_filemaker('New case file created at: {}'.format(outfile.resolve()))
