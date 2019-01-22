# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import getpass
import re
import os
import argparse
import json
import string
import copy
from typing import Optional, Tuple, Dict, List
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

__all__ = ['CaseSchema', 'CaseError']

# Default Values
_DEFAULT_HELIX_LENGTH_        = 13
_DEFAULT_BETA_LENGTH_         = 7
_DEFAULT_HELIX_DISTANCE_      = 10
_DEFAULT_HELIX_BETA_DISTANCE  = 11
_DEFAULT_BETA_PAIR_DISTANCE_  = 4.85
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
    reoriented = fields.Boolean(default=False,
                                metadata='Has connectivity directions been applied?')
    protocols = fields.List(fields.String, default=[],
                            metadata='Pipeline of protocols to run.')
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
            if i in value:
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
    layer_tilt = fields.Nested(CoordinateSchema())
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

        # layer_tilt
        if 'layer_tilt' not in data:
            data.setdefault('layer_tilt', {'x': 0, 'y': 0, 'z': 0})
        else:
            data['layer_tilt'] = cschema.fill_missing(data['layer_tilt'], 0)

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


# Error
class CaseError( ValidationError ):
    """Errors referring to :class:`.CaseSchema` processing"""


def case_corrections( case: OrderedDict, corrections: Dict) -> OrderedDict:
    """ Add positional corrections to a :class:`.CaseSchema`.

    .. note::
        This function requires a :class:`.CaseSchema` for which a ``topology.architecture``
        has already been defined.

    :param case: Current :class:`.CaseSchema` to modify.
    :type case: :class:`OrderedDict`
    :param dict corrections: Corrections to apply to the SSE, defined by SSE id.

    :return: :class:`OrderedDict` - Corrected :class:`.CaseSchema`.

    :raises:
        :CaseError: If ``topology.architecture`` has not been defined.
    """
    if 'topology' not in case or 'architecture' not in case['topology']:
        raise CaseError('Unable to apply corrections to non-existing fields.')

    schema = CaseSchema()
    case = copy.deepcopy(case)

    for j, layer in enumerate(case['topology']['architecture']):
        for i, sse in enumerate(layer):
            if sse['id'] in corrections:
                for c in corrections[sse['id']]:
                    if not isinstance(corrections[sse['id']][c], (dict, OrderedDict)):
                        case['topology']['architecture'][j][i].setdefault(c, None)
                        case['topology']['architecture'][j][i][c] = corrections[sse['id']][c]
                    else:
                        case['topology']['architecture'][j][i].setdefault(c, {})
                        for k in corrections[sse['id']][c]:
                            case['topology']['architecture'][j][i][c].setdefault(k, None)
                            case['topology']['architecture'][j][i][c][k] = corrections[sse['id']][c][k]
    return schema.load(schema.dump(case))



def case_template( name: str, architecture: Optional[str] = None, topology: Optional[str] = None,
                   corrections: Optional[Dict] = None, format: str = 'yaml', make_absolute: bool = False
                   ) -> Tuple[OrderedDict, str]:
    """Generate a :class:`.CaseSchema`.

    :param str name: Identifier of the case.
    :param str architecture: Definition of unconnected, unordered secondary structure.
    :param str topology: Definition of connected,ordered secondary structure. If provided, it will
        overwrite ``architecture``.
    :param dict corrections: Corrections to apply to the default case, identified by the SSE id.
    :param str format: Format of the output file (``yaml`` or ``json``).
    :param bool make_absolute: If :data:`True`, coordinates and shifts are
        defined as absolute positions.

    :return: :class:`OrderedDict` definition of the :class:`.CaseSchema` and
        :class:`str` generated filename.
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

    if corrections is not None:
        result = case_corrections(result, corrections)

    if make_absolute:
        result = schema.cast_absolute(result)

    # Output
    outfile = write_case(result, name, format)

    return result, outfile
