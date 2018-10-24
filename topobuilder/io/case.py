# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: CaseSchema
.. class:: CaseError
.. func:: read_case
.. func:: write_case
.. func:: case_template
"""
# Standard Libraries
import getpass
import re
import os
import argparse
import json
import string
import copy
from typing import Optional
from collections import OrderedDict

# External Libraries
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
from yaml.representer import SafeRepresenter

from marshmallow import Schema, ValidationError, fields, pprint
from marshmallow import (pre_load, validates_schema, post_load)
from marshmallow.validate import Regexp

# This Library
from ..__init__ import __version__ as __version__

__all__ = ['CaseSchema', 'CaseError', 'read_case', 'write_case', 'case_template']

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

_ACCEPTED_SSE_ID_             = '^\w\d+[HE]$'
_ACCEPTED_SSE_ID_PATTERN_     = re.compile(_ACCEPTED_SSE_ID_)
_ACCEPTED_SSE_ID_ERROR_       = "Secondary structure id should meet " \
                                "the pattern: '{}'".format(_ACCEPTED_SSE_ID_)


# Global Configuration
class LengthsSchema( Schema ):
    class Meta:
        ordered = True

    H = fields.Integer(default=_DEFAULT_HELIX_LENGTH_,
                       metadata='Number of amino acids in unspecified alpha helix.')
    E = fields.Integer(default=_DEFAULT_BETA_LENGTH_,
                       metadata='Number of amino acids in unspecified beta strand.')


class DistanceSchema( Schema ):
    class Meta:
        ordered = True

    aa = fields.Number(default=_DEFAULT_HELIX_DISTANCE_)
    ab = fields.Number(default=_DEFAULT_HELIX_BETA_DISTANCE)
    bb_pair = fields.Number(default=_DEFAULT_BETA_PAIR_DISTANCE_)
    bb_stack = fields.Number(default=_DEFAULT_BETA_STACK_DISTANCE_)
    max_loop = fields.Number(default=_DEFAULT_LOOP_DISTANCE_)

    def get_x_distance( self, data: dict, type1: str, type2: str ) -> float:
        """Provide the x distance between 2 secondary structures depending on their type.
        """
        if type1 is None:
            return 0
        if type1 == 'H' and type2 == 'H':
            return data['aa']
        if type1 == 'H' and type2 == 'E':
            return data['ab']
        if type1 == 'E' and type2 == 'H':
            return data['ab']
        if type1 == 'E' and type2 == 'E':
            return data['bb_pair']

    def get_z_distance( self, data: dict, type1: str, type2: str ) -> float:
        """Provide the z distance between 2 secondary structures depending on their type.
        """
        if type1 is None:
            return 0
        if type1 == 'H' and type2 == 'H':
            return data['aa']
        if type1 == 'H' and type2 == 'E':
            return data['ab']
        if type1 == 'E' and type2 == 'H':
            return data['ab']
        if type1 == 'E' and type2 == 'E':
            return data['bb_stack']


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

    def check_completeness( self, data: dict, info: dict ):
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``,
        :class:`.CoordinateSchema` should have explicitely defined ``x``, ``y`` and ``z``.
        """
        if 'x' not in data:
            raise CaseError('X value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'x')
        if 'y' not in data:
            raise CaseError('Y value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'y')
        if 'z' not in data:
            raise CaseError('Z value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'z')

    def fill_missing( self, data: dict, value: float ) -> dict:
        """Fill non-specified coordinates with the provided value
        """
        for i in ['x', 'y', 'z']:
            data.setdefault(i, value)
        return data

    def append_values( self, data: dict, value: dict ) -> dict:
        """Fill non-specified coordinates with the provided value
        """
        for i in ['x', 'y', 'z']:
            data[i] += value[i]
        return data


class StructureSchema( Schema ):
    class Meta:
        ordered = True

    id = fields.String(validate=Regexp(_ACCEPTED_SSE_ID_, error=_ACCEPTED_SSE_ID_ERROR_),
                       metadata='Secondary structure identifier.')
    type = fields.String(required=True, default='<type>',
                         validate=Regexp(_ACCEPTED_SSE_PATTERN_, error=_ACCEPTED_SSE_ERROR_),
                         metadata='Type of secondary structure.')
    length = fields.Integer(metadata='Amino acid length of the secondary structure.')
    coordinates = fields.Nested(CoordinateSchema())
    tilt = fields.Nested(CoordinateSchema())
    reference = fields.String()

    def get_position( self, data: dict ) -> dict:
        """Shortcut to the x, y, z coordinates of the :class:`.StructureSchema`.
        """
        return copy.deepcopy({'x': data['coordinates']['x'], 'y': data['coordinates']['y'], 'z': data['coordinates']['z']})

    def check_completeness( self, data: dict ) -> OrderedDict:
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``,
        :class:`.StructureSchema` should have explicitely defined ``coordinates``, ``angles`` and
        ``length``.
        """
        schema = CoordinateSchema()

        if 'length' not in data:
            raise CaseError('The length of the secondary structures must be provided in absolute mode.', 'length')
        if 'coordinates' not in data:
            raise CaseError('Coordinates of the secondary structures must be provided in absolute mode.', 'coordinates')
        else:
            schema.check_completeness(data['coordinates'], 'coordinates')
        if 'tilt' not in data:
            raise CaseError('Tilt of the secondary structures must be provided in absolute mode.', 'tilt')
        else:
            schema.check_completeness(data['tilt'], 'tilt')

    def cast_absolute( self, data: dict, position: dict, defaults: dict ) -> dict:
        """Transform a ``relative`` :class:`.StructureSchema` into an ``absolute`` one.
        """
        cschema = CoordinateSchema()

        # length
        if 'length' not in data:
            data.setdefault('length', defaults['length'][data['type']])

        # coordinates
        if 'coordinates' not in data:
            data.setdefault('coordinates', {'x': 0, 'y': 0, 'z': 0})
        data['coordinates'] = cschema.fill_missing(data['coordinates'], 0)
        data['coordinates'] = cschema.append_values(data['coordinates'], position)

        # tilt
        if 'tilt' not in data:
            data.setdefault('tilt', {'x': 0, 'y': 0, 'z': 0})
        else:
            data['tilt'] = cschema.fill_missing(data['tilt'], 0)

        return data


class TopologySchema( Schema ):
    class Meta:
        ordered = True

    architecture = fields.List(fields.List(fields.Nested(StructureSchema())),
                               default=[[]], required=True,
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

    @validates_schema
    def validates_absolute( self, data: dict ):
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``, all
        secondary structures from the defined ``architecture`` should have explicitely defined
        ``coordinates``, ``angles`` and ``length``.
        """
        if data['configuration']['relative'] == False:
            schema = StructureSchema()
            for layer in data['topology']['architecture']:
                for sse in layer:
                    schema.check_completeness(sse)

    def cast_absolute( self, data: dict ) -> dict:
        """Transform a ``relative`` :class:`.CaseSchema` into an ``absolute`` one.
        """
        if not data['configuration']['relative']:
            return data

        data['configuration']['relative'] = False
        sschema = StructureSchema()
        dschema = DistanceSchema()

        position = {'x': 0, 'y': 0, 'z': 0}
        defaults = data['configuration']['defaults']
        for i, layer in enumerate(data['topology']['architecture']):
            position['x'] = 0
            back = None if i == 0 else data['topology']['architecture'][i - 1][0]['type']
            here = layer[0]['type']
            position['z'] += dschema.get_z_distance(defaults['distance'], back, here)
            for j, sse in enumerate(layer):
                left = None if j == 0 else layer[j - 1]['type']
                here = sse['type']
                position['x'] += dschema.get_x_distance(defaults['distance'], left, here)
                data['topology']['architecture'][i][j] = sschema.cast_absolute(sse, position, defaults)
                position = sschema.get_position(data['topology']['architecture'][i][j])

        return data


# Error
class CaseError( ValidationError ):
    """Errors referring to :class:`.CaseSchema` processing"""


# Functions
def read_case( filename: str, make_absolute: bool = False ) -> OrderedDict:
    """Read a :class:`.CaseSchema` file into a :class:`OrderedDict`

    :param str filename: Name of the :class:`.CaseSchema` file.
    :param bool make_absolute: If :data:`True`, coordinates and shifts are
        defined as absolute positions.

    :return: :class:`OrderedDict`

    :raises:
        :IOError: If ``filename`` is not found.
    """
    result = None
    schema = CaseSchema()

    if not os.path.isfile(filename):
        raise IOError('Unable to find case file {}'.format(filename))

    try:
        result = json.loads("".join([x.strip() for x in open(filename).readlines()]))
    except json.JSONDecodeError as je:
        result = yaml.load(open(filename))

    result = schema.dump(result)
    if make_absolute:
        result = schema.cast_absolute(result)
    return schema.load(result)


def write_case( data: OrderedDict, prefix: Optional[str] = None, format: str = 'yaml' ):
    """Write :class:`.CaseSchema` into a file.

    :param OrderedDict data: :class:`.CaseSchema` content.
    :param str prefix: If provided, used as prefix for the filename. Otherwise the
        :class:`.CaseSchema` ``configuration.name`` is used.
    :param str format: Output format.

    :raises:
        :ValueError: If format is not ``yaml`` or ``json``.

    """
    if format not in ['yaml', 'json']:
        raise ValueError('Available formats are yaml or json.')

    prefix = data['configuration']['name'] if prefix is None else prefix

    if format == 'yaml':
        # This is required for YAML to properly print the Schema as an OrderedDict
        # Adapted from https://gist.github.com/oglops/c70fb69eef42d40bed06 to py3
        _mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

        def dict_representer(dumper, data):
            return dumper.represent_dict(data.items())

        def dict_constructor(loader, node):
            return OrderedDict(loader.construct_pairs(node))

        Dumper.add_representer(OrderedDict, dict_representer)
        Loader.add_constructor(_mapping_tag, dict_constructor)
        Dumper.add_representer(str, SafeRepresenter.represent_str)

        with open('{}.yml'.format(prefix), 'w') as fd:
            yaml.dump(data, fd, Dumper=Dumper, default_flow_style=False)
    else:
        with open('{}.json'.format(prefix), 'w') as fd:
            json.dump(data, fd, indent=2)


def describe_architecture( architecture: str ) -> OrderedDict:
    """Transform an architecture definition into a :class:`.TopologySchema`.

    According to [CATH](http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath), an
    architecture defines *structures that are classified according to their overall
    shape as determined by the orientations of the secondary structures in 3D space
    but ignores the connectivity between them*.

    In the **TopoSuite**, we have adopted and adapted this definition to a layer-based
    [FORM]() definition of secondary structure placements.

    The format of an architecture string definition is as follows (lower or upper case
    are both accepted)::

        2h.4e.2h

    This defines a 3-layer topology (each layer separated by points) formed by a first
    layer with 2 helices, a mid-layer with 4 beta strands and a final layer with 2 helices.

    Normally, the length of each secondary structure will be determined by the appropriate
    default setting. If length want to be defined, this can be done by providing length data
    to each layer using `:`::

        2h:13:10.4e:5:5:5:5.2h:7:13

    Notice that, if secondary structure residue length is provided, **it has to be provided
    for all secondary structures**. Unexpected behaviour might arise from doing otherwise.

    .. note::
        All architecture string definitions can be processed by the **TopoSuite**, but not
        all constructs generated by the **TopoSuite** can be minimized into an architecture
        string definition.

    :param str architecture: Architecture string definition.

    :return: :class:`OrderedDict` definition of the :class:`.TopologySchema`.

    :raises:
        :CaseError: If the string cannot be properly parsed.

    .. seealso::
        :func:`.describe_topology`
    """
    expression = re.compile('^(\d+)([EH])$')
    architecture = architecture.upper()
    result = {'architecture': []}

    for layer in architecture.split('.'):
        layer = layer.split(':')
        m = re.match(expression, layer[0])
        if not m:
            raise CaseError('Architecture format not recognized.')
        result['architecture'].append([])
        for i in range(int(m.group(1))):
            name = '{0}{1}{2}'.format(string.ascii_uppercase[len(result['architecture']) - 1],
                                      i + 1, m.group(2))
            result['architecture'][-1].append({'type': m.group(2), 'id': name})
            if len(layer) > 1:
                try:
                    result['architecture'][-1][-1].setdefault('length', int(layer[i + 1]))
                except Exception as e:
                    print(e)

    schema = TopologySchema()
    return schema.dump(result)


def describe_topology( topology: str ) -> OrderedDict:
    """Transform a topology definition into a :class:`.TopologySchema`.

    According to [CATH](http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath), a
    topology defines *structures are grouped into fold groups at this level depending
    on both the overall shape and connectivity of the secondary structures*.

    In the **TopoSuite**, a topology string defines the connectivity of the secondary
    structures and, by using the [FORM]() systematic naming system (a ``row==letters``,
    ``columns==numbers`` grid system), it also defines their position. As such, a definition
    such as this::

        B2E.C1H.B1E.A1H.B3E.A2H.B4E.C2H

    would represent a 3-layered topology (layers A, B and C) with two structures in the first
    and third layers and 4 in the middle one.

    By default, residue length of each secondary structure is defined by the default configuration,
    but it can be set up **on an individual basis** by providing the length after the secondary
    structure type::

        B2E5.C1H12.B1E4.A1H16.B3E5.A2H13.B4E5.C2H13

    :param str topology: Topology string definition.

    :return: :class:`OrderedDict` definition of the :class:`.TopologySchema`.

    :raises:
        :CaseError: If the string cannot be properly parsed.

    .. seealso::
        :func:`.describe_architecture`
    """
    expression = re.compile('^([A-Z])(\d+)([EH])(\d*)$')
    topology = topology.upper()
    tp = {}
    result = {'architecture': [], 'connectivity': [[]]}
    for sse in topology.split('.'):
        m = re.match(expression, sse)
        if not m:
            raise CaseError('Topology format not recognized.')
        sse_id = '{0}{1}{2}'.format(m.group(1), m.group(2), m.group(3))
        result['connectivity'][0].append(sse_id)
        tp.setdefault(string.ascii_uppercase.find(m.group(1)) + 1,
                      {}).setdefault(int(m.group(2)), (m.group(3), sse_id, m.group(4)))

    if list(sorted(tp.keys())) != list(range(min(tp.keys()), max(tp.keys()) + 1)):
        raise CaseError('Topology format skips layers.')
    for k1 in sorted(tp.keys()):
        result['architecture'].append([])
        if list(sorted(tp[k1].keys())) != list(range(min(tp[k1].keys()), max(tp[k1].keys()) + 1)):
            raise CaseError('Topology format skips positions in layer {}.'.format(k1))
        for k2 in sorted(tp[k1].keys()):
            scaffold = {'type': tp[k1][k2][0], 'id': tp[k1][k2][1]}
            if tp[k1][k2][2] != '':
                scaffold.setdefault('length', tp[k1][k2][2])
            result['architecture'][-1].append(scaffold)

    schema = TopologySchema()
    return schema.dump(result)


def case_template( name: str, architecture: Optional[str] = None, topology: Optional[str] = None,
                   format: str = 'yaml' ) -> OrderedDict:
    """Generate a :class:`.CaseSchema`.

    :param str name: Identifier of the case.
    :param str architecture: Definition of unconnected, unordered secondary structure.
    :param str topology: Definition of connected,ordered secondary structure. If provided, it will
        overwrite ``architecture``.
    :param str format: Format of the output file (``yaml`` or ``json``).

    :return: :class:`OrderedDict` definition of the :class:`.CaseSchema`.
    """

    # Create the case
    case = {'configuration': {'name': name}}
    schema = CaseSchema()

    # Architecture-defined case
    if architecture and not topology:
        case.setdefault('topology', describe_architecture(architecture))

    # Topology-defined case (architecture + connectivity)
    if topology and not architecture:
        case.setdefault('topology', describe_topology(topology))

    result = schema.dump(case)
    result = schema.load(result)

    # Output
    write_case(result, name, format)

    return result
