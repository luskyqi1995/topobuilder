# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import re
import json
import string
from copy import deepcopy
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Union, TypeVar
from collections import OrderedDict

# External Libraries
import yaml
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
from yaml.representer import SafeRepresenter

# This Library
from .schema import (CaseSchema, CaseError, TopologySchema, CoordinateSchema,
                     StructureSchema, DistanceSchema)

__all__ = ['Case']

C = TypeVar('C', bound='Case')


class Case( object ):
    """
    """
    def __init__( self, init: Optional[Union[str, dict, Path, C]] = None ):
        self.data = OrderedDict()
        self.schema = CaseSchema()

        if isinstance(init, str):
            self.data = OrderedDict({'configuration': {'name': init}})
        elif isinstance(init, Case):
            self.data = init.check().data
        elif isinstance(init, (dict, OrderedDict)):
            self.data = OrderedDict(deepcopy(init))
            self.data = self.schema.load(self.data)
        elif isinstance(init, Path):
            if not init.is_file():
                raise IOError('Unable to find case file {}'.format(init.resolve()))
            if init.suffix == '.gz':
                raise IOError('Unable to manage gzipped file case {}'.format(init.resolve()))
            try:
                self.data = json.loads("".join([x.strip() for x in open(init).readlines()]))
            except json.JSONDecodeError:
                self.data = yaml.load(open(init))

        self.check()

    @property
    def name( self ) -> str:
        """Returns the :class:`.Case` identifier.

        :return: class:`str`
        """
        return self['configuration.name']

    @property
    def shape( self ) -> Tuple[int]:
        """Returns a :class:`tuple` with the number of secondary structures in each layer.

        :return: :class:`tuple`
        """
        if 'architecture' not in self:
            return tuple()
        return tuple([len(layer) for layer in self['topology.architecture']])

    @property
    def shape_len( self ) -> Tuple[Tuple[int]]:
        """Returns the length of each secondary structure in a :class:`tuple` with
        the same shape as the ``architecture``.

        :return: :class:`tuple`
        """
        if 'architecture' not in self:
            return tuple()

        # make a copy absolute first
        c = Case(self.data).cast_absolute()
        out = []
        for layer in c['topology.architecture']:
            out.append([])
            for sse in layer:
                out[-1].append(sse['length'])
            out[-1] = tuple(out[-1])
        return tuple(out)

    @property
    def center_shape( self ) -> Dict:
        """Returns the limit positions for each layer and in total, taking
        only into account the center of each sse.

        :return: :class:`dict`
        """
        if 'architecture' not in self:
            return {}

        asciiU = string.ascii_uppercase

        # make a copy absolute first
        c = Case(self.data).cast_absolute()
        result = {}
        for il, layer in enumerate(c['topology.architecture']):
            result.setdefault(asciiU[il], {'top': 0, 'bottom': 0, 'left': 0, 'right': 0})
            for sse in layer:
                if sse['coordinates']['y'] > result[asciiU[il]]['top']:
                    result[asciiU[il]]['top'] = sse['coordinates']['y']
                if sse['coordinates']['y'] < result[asciiU[il]]['bottom']:
                    result[asciiU[il]]['bottom'] = sse['coordinates']['y']
                if sse['coordinates']['x'] < result[asciiU[il]]['left']:
                    result[asciiU[il]]['left'] = sse['coordinates']['x']
                if sse['coordinates']['x'] > result[asciiU[il]]['right']:
                    result[asciiU[il]]['right'] = sse['coordinates']['x']
            result[asciiU[il]]['width'] = result[asciiU[il]]['right'] - result[asciiU[il]]['left']
            result[asciiU[il]]['hight'] = result[asciiU[il]]['top'] - result[asciiU[il]]['bottom']
        return result

    @property
    def connectivity_count( self ) -> int:
        """Returns the number of connectivities defined in the ``topology``.

        :return: :class:`int`
        """
        if 'connectivity' not in self:
            return 0
        else:
            return len(self['topology.connectivity'])

    @property
    def architecture_str( self ) -> str:
        """Returst a string representation of the architecture.

        :return: :class:`str`
        """
        if 'architecture' not in self:
            return ''
        return architecture_cast(self['topology.architecture'])

    @property
    def connectivities_str( self ) -> Tuple[str]:
        """Returns a list of string representations of topologies.

        :return: :class:`tuple` of :class:`str`
        """
        result = []
        if self.connectivity_count == 0:
            return result

        for i, _ in enumerate(self['topology.connectivity']):
            result.append(topology_cast(self.data, i))
        return tuple(result)

    @property
    def directionality_profile( self ) -> str:
        """Returns a binary-type signature representing the directionality
        of all secondary structures without taking connectivity into account.

        .. note::
            Proper values for this function requires the :class:`Case` to be
            **reoriented**.

        :return: :class:`str`
        """
        result = ''
        if 'architecture' not in self:
            return result

        c = Case(self).cast_absolute()
        for layer in c['topology.architecture']:
            for sse in layer:
                result += '1' if sse['tilt']['x'] > 90 and sse['tilt']['x'] < 270 else '0'
        return result

    @property
    def is_absolute( self ) -> bool:
        cr = self['configuration.relative']
        if cr is not None:
            return not cr
        return False

    @property
    def is_relative( self ) -> bool:
        return not self.is_absolute

    @property
    def is_reoriented( self ) -> bool:
        cr = self['configuration.reoriented']
        if cr is not None:
            return cr
        return False

    def get_type_for_layer( self, layer: Union[int, str] ) -> str:
        layerint = layer_int(layer)

        if layerint > len(self['topology.architecture']):
            raise IndexError('Requested layeer is bigger than any available.')

        if len(self['topology.architecture'][layerint]) == 0:
            return 'X'

        return self['topology.architecture'][layerint][0]['type']

    def set_type_for_layer( self, layer: Union[int, str], sse_count: int) -> C:
        sschema = StructureSchema()
        sse_type = self.get_type_for_layer(layer)
        layerint = layer_int(layer)
        layerstr = layer_str(layer)

        ks = Case(self)
        sse = []
        for i in range(sse_count):
            sse.append(sschema.dump({'id': '{0}{1}{2}'.format(layerstr, i + 1, sse_type),
                                     'type': sse_type}))
        ks.data['topology']['architecture'][layerint] = sse
        return ks

    def check( self ) -> C:
        """Evaluate the :class:`.Case` content thourhg the :class:`.CaseSchema`.
        """
        self.data = self.schema.load(self.schema.dump(self.data))
        return self

    def add_architecture( self, architecture: Optional[str] = None ) -> C:
        """Adds an architecture definition to the :class:`.Case`.

        According to `CATH <http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath>`_, an
        architecture defines *structures that are classified according to their overall
        shape as determined by the orientations of the secondary structures in 3D space
        but ignores the connectivity between them*.

        In the **TopoSuite**, we have adopted and adapted this definition to a layer-based
        `FORM <https://doi.org/10.1016/j.str.2009.07.012>`_ definition of secondary structure
        placements.

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

        :return: :class:`.Case` - updated case.

        :raises:
            :CaseOverwriteError: If an architecture is already defined in this :class:`.Case`
            :CaseError: If the string cannot be properly parsed.

        .. seealso::
            :func:`.add_topology`
        """
        if architecture is None:
            return Case(self.data).check()

        if 'architecture' in self:
            raise CaseOverwriteError('An arquitecture is already defined.')

        c = Case(self.data)
        c.data['topology']['architecture'] = architecture_cast(architecture)['architecture']
        return c.check()

    def add_topology( self, topology: Optional[str] = None ) -> C:
        """Adds a topology definition to the :class:`.Case`.

        According to `CATH <http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath>`_, a
        topology defines *structures are grouped into fold groups at this level depending
        on both the overall shape and connectivity of the secondary structures*.

        In the **TopoSuite**, a topology string defines the connectivity of the secondary
        structures and, by using the `FORM <https://doi.org/10.1016/j.str.2009.07.012>`_ systematic
        naming system (a ``row==letters``, ``columns==numbers`` grid system), it also defines their
        position. As such, a definition such as this::

            B2E.C1H.B1E.A1H.B3E.A2H.B4E.C2H

        would represent a 3-layered topology (layers A, B and C) with two structures in the first
        and third layers and 4 in the middle one.

        By default, residue length of each secondary structure is defined by the default configuration,
        but it can be set up **on an individual basis** by providing the length after the secondary
        structure type::

            B2E5.C1H12.B1E4.A1H16.B3E5.A2H13.B4E5.C2H13

        :param str topology: Topology string definition.

        :return: :class:`.Case` - updated case.

        :raises:
            :CaseError: If the string cannot be properly parsed.
            :CaseOverwriteError: If an architecture is already defined in this :class:`.Case` and the
                provided topology does not match

        .. seealso::
            :func:`.describe_architecture`
        """
        if topology is None:
            return Case(self.data).check()

        t = Case('temp')
        t.data['topology'] = topology_cast(topology)

        if 'architecture' in self:
            if self.architecture_str != t.architecture_str or self.shape_len != t.shape_len:
                raise CaseOverwriteError('Provided topology does not match existing architecture.')

        c = Case(self.data)
        if 'connectivity' not in c:
            c.data['topology'] = topology_cast(topology)
        else:
            conn = topology_cast(topology)['connectivity'][0]
            if ".".join(conn) not in set(self.connectivities_str):
                c.data['topology']['connectivity'].append(conn)
        return c.check()

    def cast_absolute( self ) -> C:
        """Transform a ``relative`` :class:`.CaseSchema` into an ``absolute`` one.
        """
        c = Case(self.data)
        if self.is_absolute:
            return c

        c.data['configuration']['relative'] = False
        sschema = StructureSchema()
        dschema = DistanceSchema()

        position = {'x': 0, 'y': 0, 'z': 0}
        defaults = c['configuration.defaults']
        for i, layer in enumerate(c['topology.architecture']):
            position['x'] = 0
            back = None if i == 0 else c['topology.architecture'][i - 1][0]['type']
            here = layer[0]['type']
            position['z'] += dschema.get_z_distance(defaults['distance'], back, here)
            for j, sse in enumerate(layer):
                left = None if j == 0 else layer[j - 1]['type']
                here = sse['type']
                # X shift is inherited in the following structure.
                position['x'] += dschema.get_x_distance(defaults['distance'], left, here)
                position['y'] = 0  # Reset Y coordinate
                c.data['topology']['architecture'][i][j] = sschema.cast_absolute(sse, position, defaults)
                position = sschema.get_position(c['topology.architecture'][i][j])
        return c.check()

    def apply_topologies( self ) -> List[C]:
        """Generates a :class:`List` of :class:`.Case` in which the different available connectivities
        have been applied.

        :return: :class:`List` of :class:`.Case`

        :raises:
            :CaseIncompleteError: If there is not enough data to apply topologies.
        """
        results = []

        if 'connectivity' not in self or 'architecture' not in self:
            raise CaseIncompleteError('Unable to apply corrections to non-existing fields.')

        for i, cn in enumerate(self['topology.connectivity']):
            corrections = {}
            c = Case(self.data)
            c['topology']['connectivity'] = [self['topology.connectivity'][i]]
            for turn in c['topology.connectivity'][0][1::2]:
                corrections.setdefault(turn, {'tilt': {'y': 180, 'x': 180}})
            results.append(c.apply_corrections(corrections))
            results[-1].data['configuration']['reoriented'] = True
        return results

    def apply_corrections( self, corrections: Optional[Union[Dict, str, Path]] = None ) -> C:
        """
        """
        if corrections is None:
            return Case(self.data).check()

        if isinstance(corrections, str):
            corrections = Path(corrections)

        if isinstance(corrections, Path):
            if not corrections.is_file():
                raise IOError('Unable to find corrections file {}'.format(corrections.resolve()))
            if corrections.suffix == '.gz':
                raise IOError('Unable to manage gzipped file case {}'.format(corrections.resolve()))
            try:
                crr = json.loads("".join([x.strip() for x in open(corrections).readlines()]))
            except json.JSONDecodeError:
                crr = yaml.load(open(corrections))
            return Case(self.data).apply_corrections(crr)

        # APPLY LAYER CORRECTIONS
        corrections = layer_corrections(corrections, self)

        # APPLY SSE CORRECTIONS
        return sse_corrections(corrections, self)

    def set_protocol_done( self, protocol_id: int ) -> C:
        """Label a protocol as executed.

        :param int protocol_id: Identifier of the protocol according to its position.
        """
        if protocol_id == -1:
            return
        if self['configuration.protocols'] is None or len(self['configuration.protocols']) < protocol_id:
            raise IndexError('Trying to access an unspecified protocol.')

        c = Case(self.data)
        c.data['configuration']['protocols'][protocol_id].setdefault('status', True)
        c.data['configuration']['protocols'][protocol_id]['status'] = True
        return c

    def assign_protocols( self, protocols: List[Dict] ) -> C:
        """Overwrite current protocols in the :class:`.Case` with the provided ones.

        :param protocols: New protocols for the :class:`.Case`
        """
        c = Case(self.data)
        if c['configuration.protocols'] is None:
            c.data['configuration'].setdefault('protocols', [])
        c.data['configuration']['protocols'] = protocols
        return c

    def write( self, prefix: Optional[Union[str, Path]] = None, format: str = 'yaml' ) -> Path:
        """Write :class:`.Case` into a file.

        :param str prefix: If not specified, output file is created in the **current working directory**,
            using :class:`.Case` ``configuration.name`` as prefix. A :class:`str` will modify the name of
            the file. A directory :class:`Path` will change dir and use the default prefix, while a
            non-directory :class:`Path` will set up output directory and file prefix.
        :type prefix: Union[str, Path]
        :param str format: Output format. Options are [``yaml``, ``json``]

        :return: :class:`Path` - filename

        :raises:
            :ValueError: If format is not ``yaml`` or ``json``.

        """
        if format not in ['yaml', 'json']:
            raise ValueError('Available formats are yaml or json.')

        if prefix is None:
            prefix = self['configuration.name']
        if isinstance(prefix, Path):
            if prefix.is_dir():
                prefix = prefix.joinpath(self['configuration.name'])

        if format == 'yaml':
            YAML_Dumper()
            with open('{}.yml'.format(str(prefix)), 'w') as fd:
                yaml.dump(self.data, fd, Dumper=Dumper, default_flow_style=False)
            return Path('{}.yml'.format(str(prefix)))
        else:
            with open('{}.json'.format(str(prefix)), 'w') as fd:
                json.dump(self.data, fd, indent=2)
            return Path('{}.json'.format(str(prefix)))

    def __contains__( self, item ):
        if isinstance(item, str):
            if item == 'architecture':
                ta = self['topology.architecture']
                return False if ta is None else ta != [[]]
            if item == 'connectivity':
                ta = self['topology.connectivity']
                return False if ta is None else ta != [[]]

        raise NotImplementedError()

    def __getitem__( self, key ):
        r = self.data
        for k in key.split('.'):
            r = r.get(k, None)
            if r is None:
                break
        return r

    def __eq__( self, value ):
        if not isinstance(value, Case):
            raise NotImplementedError()
        return self.data == value.data


def layer_corrections( corrections: dict, case: Case) -> dict:
    """
    """
    asciiU = string.ascii_uppercase
    sizes = case.center_shape
    maxwidth = max(sizes[l]['width'] for l in sizes)
    for j, layer in enumerate(case['topology.architecture']):
        layer_id = asciiU[j]
        if layer_id in corrections:
            c = corrections[layer_id]
            if 'xalign' in c:
                diffw = maxwidth - sizes[layer_id]['width']
                if diffw != 0 and c['xalign'].lower() != 'left':
                    sse_id = layer[0]['id']
                    corrections.setdefault(sse_id, {}).setdefault('coordinates', {}).setdefault('x', 0)
                    if c['xalign'].lower() == 'right':
                        corrections[sse_id]['coordinates']['x'] += diffw
                    elif c['xalign'].lower() == 'center':
                        corrections[sse_id]['coordinates']['x'] += diffw / 2
    return corrections


def sse_corrections( corrections: dict, case: Case) -> Case:
    """
    """
    cs = CoordinateSchema()
    case = deepcopy(case.data)
    for j, layer in enumerate(case['topology']['architecture']):
        for i, sse in enumerate(layer):
            if sse['id'] in corrections:
                for c in corrections[sse['id']]:
                    if not isinstance(corrections[sse['id']][c], (dict, OrderedDict)):
                        case['topology']['architecture'][j][i].setdefault(c, None)
                        case['topology']['architecture'][j][i][c] = corrections[sse['id']][c]
                    else:  # has to be in ['coordinates', 'tilt', 'layer_tilt']
                        case['topology']['architecture'][j][i].setdefault(c, {})
                        ks = cs.fill_missing(case['topology']['architecture'][j][i][c], 0)
                        case['topology']['architecture'][j][i][c] = cs.append_values(ks, corrections[sse['id']][c])

                        # for k in corrections[sse['id']][c]:
                        #     case['topology']['architecture'][j][i][c].setdefault(k, None)
                        #     case['topology']['architecture'][j][i][c][k] = corrections[sse['id']][c][k]
    return Case(case).check()


def layer_cast( layer: Union[int, str] ) -> Union[int, str]:
    if isinstance(layer, int):
        return string.ascii_uppercase[layer]
    elif isinstance(layer, str):
        return string.ascii_uppercase.find(layer.upper())
    else:
        raise ValueError('Layer is defined by integer or string.')


def layer_int( layer: Union[str, int] ) -> int:
    if isinstance(layer, int):
        return layer
    else:
        return layer_cast(layer)


def layer_str( layer: Union[str, int] ) -> str:
    if isinstance(layer, str):
        return layer.upper()
    else:
        return layer_cast(layer)


def architecture_cast( architecture: Union[str, List, dict] ) -> Union[dict, str]:
    """
    """
    if isinstance(architecture, str):
        expression = re.compile(r'^(\d+)([EH])$')
        architecture = architecture.upper()
        asciiU = string.ascii_uppercase
        result = {'architecture': []}

        for layer in architecture.split('.'):
            layer = layer.split(':')
            m = re.match(expression, layer[0])
            if not m:
                raise CaseError('Architecture format not recognized.')
            result['architecture'].append([])
            for i in range(int(m.group(1))):
                name = '{0}{1}{2}'.format(asciiU[len(result['architecture']) - 1],
                                          i + 1, m.group(2))
                result['architecture'][-1].append({'type': m.group(2), 'id': name})
                if len(layer) > 1:
                    try:
                        result['architecture'][-1][-1].setdefault('length', int(layer[i + 1]))
                    except IndexError:
                        print('Lengths were not provided for all defined secondary structures.')
                    except ValueError:
                        print('Length values MUST BE integers.')
                    except Exception as e:
                        print(e)

        schema = TopologySchema()
        return schema.dump(result)

    if isinstance(architecture, list):
        architecture = {'architecture': architecture}

    result = []
    for layer in architecture['architecture']:
        result.append('{0}{1}'.format(len(layer), layer[0]['type']))
    result = '.'.join(result)
    return result


def topology_cast( topology: Union[str, dict], count: Optional[int] = 0 ) -> Union[dict, str]:
    """
    """
    if isinstance(topology, str):
        expression = re.compile(r'^([A-Z])(\d+)([EH])(\d*)$')
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

    return ".".join(topology['topology']['connectivity'][count])


def YAML_Dumper():
    # This is required for YAML to properly print the Schema as an OrderedDict
    # Adapted from https://gist.github.com/oglops/c70fb69eef42d40bed06 to py3
    def dict_representer(dumper, data):
        return dumper.represent_dict(data.items())

    Dumper.add_representer(OrderedDict, dict_representer)
    Dumper.add_representer(str, SafeRepresenter.represent_str)


# Errors
class CaseOverwriteError( CaseError ):
    """Error raised when trying to override :class:`.Case` data in an unexpected way
    """


class CaseIncompleteError( CaseError ):
    """Errors raise when trying to execute a process without enough information
    """
