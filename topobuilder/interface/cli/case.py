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
import argparse

# External Libraries

# This Library
from topobuilder.io import case_template, read_case, write_case


def cli_case_template():
    """Generate a :class:`.CaseSchema`.
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

    case_template(**vars(options))


def cli_absolute_case():
    """Transform a relative :class:`.CaseSchema` to absolute.
    """
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_case_template.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-case', dest='case', action='store',
                        help='Relative case file.', required=True)
    options = parser.parse_args()

    # Process naming system
    prefix = options.case.split('.')
    format = 'yaml' if prefix[-1] == 'yml' else 'json'
    prefix[-1] = 'absolute'
    prefix = '.'.join(prefix)

    # Read, transform and write
    data = read_case(options.case, True)
    write_case(data, prefix, format)
