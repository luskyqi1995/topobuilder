# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: Architect
"""
# Standard Libraries
import string
from pathlib import Path
from typing import Union, Dict

# External Libraries

# This Library
from topobuilder.case import Case
from .parametric import ParametricStructure


__all__ = ['GeneralArchitect', 'SSEArchitect']


class GeneralArchitect( object ):
    """
    """
    def __init__( self,
                  case: Union[str, Dict, Case, Path],
                  paths: Dict ):
        self.case = Case(case).cast_absolute()
        self.path = paths

    def build_sketch( self ):
        """
        """
        shapeForm = self.build_structure()
        with open(self.path['arch_sketch'], "w") as fd:
            fd.write(shapeForm.to_pdb())
        return self.path['arch_sketch']

    def build_connectivity( self, subdir ):
        """
        """
        shapeForm = self.build_structure(connectivity=True)
        outd = self.path['connectivity'].joinpath(subdir)
        outd.mkdir(parents=True, exist_ok=True)
        outf = outd.joinpath('sketch.pdb')
        with open(outf, "w") as fd:
            fd.write(shapeForm.to_pdb())
        return outf

    def build_connectivities( self ):
        """
        """
        cases = self.case.apply_topologies()
        outfls = []
        for c in cases:
            subdir = '.'.join(c['topology']['connectivity'][0])
            newarch = GeneralArchitect(c, self.path)
            outfls.append(newarch.build_connectivity(subdir))
        return outfls

    # def build_structure( self, connectivity=False ):
    #     """
    #     """
    #     sselist = []
    #     for ilayer, layer in enumerate(self.case.data['topology']['architecture']):
    #         print('building layer {}'.format(ilayer + 1))
    #         for iss, ss in enumerate(layer):
    #             print('  building SSE {}'.format(iss + 1))
    #             # sselist.append(SSEArchitect(ss, type=ss['type']).pdb)
    #             # sselist[-1].write('test{}.pdb'.format(iss), format='pdb')
    #             coordinates = [ss['coordinates']['x'], ss['coordinates']['y'], ss['coordinates']['z']]
    #             vs = VirtualMaker(ss['length'], coordinates, type=ss['type'])
    #             vs.tilt_y_degrees(ss['tilt']['y'])
    #             vs.tilt_degrees(ss['tilt']['x'], 0, ss['tilt']['z'])
    #             sselist.append(vs)
    #
    #     if connectivity:
    #         mintp = '.'.join([_[:2] for _ in self.case.data['topology']['connectivity'][0]])
    #         crval = list(filter((0).__ne__, [mintp.count(cr) for cr in string.ascii_uppercase]))
    #         crval.insert(0, 0)
    #
    #         new_order = []
    #         for pid in [_[:2] for _ in mintp.split('.')]:
    #             layerv = crval[string.ascii_uppercase.find(pid[0])]
    #             new_order.append(layerv + int(pid[1]) - 1)
    #         sselist = [sselist[i] for i in new_order]
    #
    #     shapeForm = Form("shapesketch", sselist, None)
    #     shapeForm.prepare_coords()
    #     return shapeForm


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
    # CHOP780201 alpha-helix propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.57; P -> 0.57
    _AA_STAT = [("A", 1.42), ("L", 1.21), ("R", 0.98), ("K", 1.16), ("N", 0.67),
                ("M", 1.45), ("D", 1.01), ("F", 1.13), ("C", 0.70), ("P", 0.00),
                ("Q", 1.11), ("S", 0.77), ("E", 1.51), ("T", 0.83), ("G", 0.00),
                ("W", 1.08), ("H", 1.00), ("Y", 0.69), ("I", 1.08), ("V", 1.06)]


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
    _ROTATION = -180
    _MONO = {'N': [-0.440, -1.200, 0.330],
             'CA': [0.000, 0.000, 1.210],
             'C': [-0.550, 1.200, 0.330],
             'O': [-2.090, 1.300, 0.220]}
    # CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.75; P -> 0.55
    _AA_STAT = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
                ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 1.19), ("P", 0.00),
                ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.00),
                ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]
