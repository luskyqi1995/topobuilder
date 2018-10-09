# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. obj:: CaseSchema
.. obj:: CaseError
.. func:: case_template
"""
# Standard Libraries
import getpass
import re
import argparse
import json
import string
from collections import OrderedDict

# External Libraries
import yaml
from yaml import CLoader as Loader, CDumper as Dumper
from yaml.representer import SafeRepresenter
from marshmallow import Schema, fields, pprint
from marshmallow.validate import Regexp

# This Library
from ..__init__ import __version__ as __version__

__all__ = ['CaseSchema', 'CaseError', 'case_template']

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

# Default Values
_DEFAULT_HELIX_LENGTH_        = 13
_DEFAULT_BETA_LENGTH_         = 7
_DEFAULT_HELIX_DISTANCE_      = 10
_DEFAULT_HELIX_BETA_DISTANCE  = 11
_DEFAULT_BETA_PAIR_DISTANCE_  = 5
_DEFAULT_BETA_STACK_DISTANCE_ = 8
_DEFAULT_LOOP_DISTANCE_       = 18.97

_ACCEPTED_SSE_TYPES_          = '^[HE]$|^S[2-9]\d*$'
_ACCEPTED_SSE_PATTERN_        = re.compile(_ACCEPTED_SSE_TYPES_)
_ACCEPTED_SSE_ERROR_          = "Structure type should meet " \
                                "the pattern: '{}'".format(_ACCEPTED_SSE_TYPES_)

_ACCEPTED_SSE_ID_             = '^\w\d+$'
_ACCEPTED_SSE_ID_PATTERN_     = re.compile(_ACCEPTED_SSE_ID_)
_ACCEPTED_SSE_ID_ERROR_       = "Secondary structure id should meet " \
                                "the pattern: '{}'".format(_ACCEPTED_SSE_ID_)


# Global Configuration
class LengthsSchema( Schema ):
    class Meta:
        ordered = True

    alpha = fields.Integer(default=_DEFAULT_HELIX_LENGTH_,
                           metadata='Number of amino acids in unspecified alpha helix.')
    beta = fields.Integer(default=_DEFAULT_BETA_LENGTH_,
                          metadata='Number of amino acids in unspecified beta strand.')


class DistanceSchema( Schema ):
    class Meta:
        ordered = True

    aa = fields.Number(default=_DEFAULT_HELIX_DISTANCE_)
    ab = fields.Number(default=_DEFAULT_HELIX_BETA_DISTANCE)
    bb_pair = fields.Number(default=_DEFAULT_BETA_PAIR_DISTANCE_)
    bb_stack = fields.Number(default=_DEFAULT_BETA_STACK_DISTANCE_)
    max_loop = fields.Number(default=_DEFAULT_LOOP_DISTANCE_)


class DefaultSchema( Schema ):
    class Meta:
        ordered = True

    length = fields.Nested(LengthsSchema(), default=LengthsSchema().dump({}))
    distance = fields.Nested(DistanceSchema(), default=DistanceSchema().dump({}))


class ConfigurationSchema( Schema ):
    class Meta:
        ordered = True

    name = fields.String(required=True, default='<name>',
                         error_messages={'required': 'A case identifier is required'},
                         metadata='Case identifier.')
    user = fields.String(default=getpass.getuser(),
                         metadata='User identifier.')
    defaults = fields.Nested(DefaultSchema(), default=DefaultSchema().dump({}),
                             metadata='Default parameters.')
    relative = fields.Boolean(default=True,
                              metadata='Relative vs. absolute coordinates.')
    comments = fields.List(fields.String, default=[__version__],
                           metadata='Relative vs. absolute coordinates.')


# Topology
class CoordinateSchema( Schema ):
    class Meta:
        ordered = True

    x = fields.Number(metadata='Value for the X axis.')
    y = fields.Number(metadata='Value for the Y axis.')
    z = fields.Number(metadata='Value for the Z axis.')


class StructureSchema( Schema ):
    class Meta:
        ordered = True

    id = fields.String(validate=Regexp(_ACCEPTED_SSE_ID_, error=_ACCEPTED_SSE_ID_ERROR_),
                       metadata='Secondary structure identifier.')
    type = fields.String(required=True, default='<type>',
                         validate=Regexp(_ACCEPTED_SSE_PATTERN_, error=_ACCEPTED_SSE_ERROR_),
                         metadata='Type of secondary structure.')
    length = fields.Number(metadata='Amino acid length of the secondary structure.')
    coordinates = fields.Nested(CoordinateSchema())
    tilt = fields.Nested(CoordinateSchema())
    reference = fields.String()


class TopologySchema( Schema ):
    class Meta:
        ordered = True

    architecture = fields.List(fields.List(fields.Nested(StructureSchema())),
                               metadata='Positions and definitions for the secondary structures.')
    connectivity = fields.List(fields.List(fields.String()),
                               metadata='Sequence order of the secondary structures.')


# Case
class CaseSchema( Schema ):
    class Meta:
        ordered = True

    configuration = fields.Nested(ConfigurationSchema(), required=True,
                                  default=ConfigurationSchema().dump({}),
                                  error_messages={'required': 'Configuration data is required'},
                                  metadata='TopoBuilder case definition.')
    topology = fields.Nested(TopologySchema(), required=True, default=TopologySchema().dump({}),
                             error_messages={'required': 'A topological definition is required'},
                             metadata='Topology Definition.')


# Error
class CaseError( Exception ):
    """Errors referring to :class:`.CaseSchema` processing"""


#
# CLI command to create default case files (starting point)
#
def describe_architecture( architecture: str ) -> OrderedDict:
    """
    """
    expression = re.compile('^(\d+)([EH])$')
    architecture = architecture.upper()
    result = {'architecture': []}

    for layer in architecture.split('.'):
        m = re.match(expression, layer)
        if not m:
            raise CaseError('Architecture format not recognized.')
        result['architecture'].append([])
        for i in range(int(m.group(1))):
            name = '{0}{1}{2}'.format(string.ascii_uppercase[len(result['architecture']) - 1],
                                      i + 1, m.group(2))
            result['architecture'][-1].append({'type': m.group(2),
                                               'id': name})

    schema = TopologySchema()
    return schema.dump(result)


def describe_topology( topology: str) -> OrderedDict:
    """
    """
    expression = re.compile('^([A-Z])(\d+)([EH])$')
    topology = topology.upper()
    tp = {}
    result = {'architecture': [], 'connectivity': [[]]}
    for sse in topology.split('.'):
        m = re.match(expression, sse)
        if not m:
            raise CaseError('Topology format not recognized.')
        result['connectivity'][0].append(sse)
        tp.setdefault(string.ascii_uppercase.find(m.group(1)) + 1,
                      {}).setdefault(int(m.group(2)), (m.group(3), sse))

    if list(sorted(tp.keys())) != list(range(min(tp.keys()), max(tp.keys()) + 1)):
        raise CaseError('Topology format skips layers.')
    for k1 in sorted(tp.keys()):
        result['architecture'].append([])
        if list(sorted(tp[k1].keys())) != list(range(min(tp[k1].keys()), max(tp[k1].keys()) + 1)):
            raise CaseError('Topology format skips positions in layer {}.'.format(k1))
        for k2 in sorted(tp[k1].keys()):
            result['architecture'][-1].append({'type': tp[k1][k2][0],
                                               'id': tp[k1][k2][1]})

    schema = TopologySchema()
    return schema.dump(result)


def case_template():
    """Generate an empty :class:`.CaseSchema`.

    :param str topology:
    :param str connectivity:

    """

    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=case_template.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-name', dest='name', action='store',
                        help='Job Name.', required=True)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-architecture', dest='architecture', action='store',
                       help='Initial expected architecture.', default=None)
    group.add_argument('-topology', dest='topology', action='store',
                       help='Expected connected topology.',
                       default=None)

    parser.add_argument('-absolute', dest='relative', action='store_false',
                        help='Position is absolute instead of relative (default).',
                        default=True)
    parser.add_argument('-format', dest='format', action='store', choices=['json', 'yaml'],
                        help='Format for the case file.',
                        default='yaml')

    options = parser.parse_args()

    # Create the case
    case = {'configuration': {'name': options.name}}
    schema = CaseSchema()

    # Architecture-defined case
    if options.architecture:
        case.setdefault('topology', describe_architecture(options.architecture))

    # Topology-defined case (architecture + connectivity)
    if options.topology:
        case.setdefault('topology', describe_topology(options.topology))

    result = schema.dump(case)

    # Output
    if options.format == 'yaml':
        # This is required for YAML to properly print the Schema as an OrderedDict
        # Adapted from https://gist.github.com/oglops/c70fb69eef42d40bed06 to py3
        def dict_representer(dumper, data):
            return dumper.represent_dict(data.items())

        def dict_constructor(loader, node):
            return OrderedDict(loader.construct_pairs(node))

        Dumper.add_representer(OrderedDict, dict_representer)
        Loader.add_constructor(_mapping_tag, dict_constructor)
        Dumper.add_representer(str, SafeRepresenter.represent_str)

        with open('{}.yml'.format(options.name), 'w') as fd:
            yaml.dump(result, fd, Dumper=Dumper, default_flow_style=False)
    else:
        with open('{}.json'.format(options.name), 'w') as fd:
            json.dump(result, fd, indent=2)
    pprint(result)
    print(type(result))
