# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: Architect
"""
# Standard Libraries
from collections import OrderedDict
import os

# External Libraries


# This Library
from topobuilder.io import CaseSchema
from .parametric import ParametricStructure
from .virtual.VirtualMaker import VirtualMaker
from ..form.Form import Form

__all__ = ['GeneralArchitect']


class GeneralArchitect( object ):
    """
    """
    def __init__( self, case: OrderedDict, paths: dict ):
        schema = CaseSchema()
        self.case = schema.cast_absolute(case)
        self.path = paths

    def build_sketch( self ):
        """
        """
        sselist = []
        for layer in self.case['topology']['architecture']:
            for ss in layer:
                vs = SSEArchitect(ss, type=ss['type'])
                coordinates = [ss['coordinates']['x'], ss['coordinates']['y'], ss['coordinates']['z']]
                vs = VirtualMaker(ss['length'], coordinates, type=ss['type'])
                vs.tilt_y_degrees(ss['tilt']['y'])
                vs.tilt_degrees(ss['tilt']['x'], 0, ss['tilt']['z'])
                sselist.append(vs)
        print(sselist)
        shapeForm = Form("shapesketch", sselist, None)
        shapeForm.prepare_coords()
        with open(self.path['arch_sketch'], "w") as fd:
            fd.write(shapeForm.to_pdb())


class SSEArchitect( object ):
    """Decides the correct type of secondary structure to build.
    """
    def __new__( cls, *args, **kwargs ):
        sse_type = kwargs.pop('type', None)
        if sse_type is None:
            raise AttributeError('A secondary structure type must be provided.')
        sse_type = sse_type.upper()

        if sse_type == 'H':
            return AlphaHelixArchitect(*args, **kwargs)
        elif sse_type == 'G':
            return Helix310Architect(*args, **kwargs)
        elif sse_type == 'I':
            return HelixPiArchitect(*args, **kwargs)
        elif sse_type == 'E':
            return FlatBetaArchitect(*args, **kwargs)
        else:
            raise ValueError('Unrecognized secondary structure type.')


class AlphaHelixArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 1.5
    _ROTATION = 100
    _MONO = {'N': [1.321, 0.841, -0.711],
             'CA': [2.300, 0.000, 0.000],
             'C': [1.576, -1.029, 0.870],
             'O': [1.911, -2.248, 0.871]}


class Helix310Architect( ParametricStructure ):
    """
    """
    _PERIODE = 2.0


class HelixPiArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 1.1


class FlatBetaArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 3.2
