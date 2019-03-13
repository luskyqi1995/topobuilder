# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: ParametricStructure
"""
# Standard Libraries
import sys
from random import random
from bisect import bisect

# External Libraries
import numpy as np
import pandas as pd
from SBI.structure import PDB, Frame3D
import SBI.structure.geometry as SBIgeo
from SBI.data import alphabet

# This Library
import topobuilder.core as TBcore


__all__ = ['ParametricStructure']


class ParametricStructure( object ):

    _MONO = None
    _PERIODE = None
    _ROTATION = None
    _AA_STAT = None

    def __init__( self, indata ):
        """
        """
        if isinstance(indata, Frame3D):
            self.pdb = indata
            self.desc = None
            self.reverse()
        elif isinstance(indata, dict):
            self.pdb = []
            self.desc = indata
            self.build()

    def build( self ):
        """
        """
        def weighted_choice(choices):
            values, weights = zip(*choices)
            total = 0
            cum_weights = []
            for w in weights:
                total += w
                cum_weights.append(total)
            x = random() * total
            i = bisect(cum_weights, x)
            return values[i]

        if self._MONO is None or self._PERIODE is None:
            raise NotImplementedError()

        # 1. Locate center point for each residue we need to build
        vector_module = float(self._PERIODE * (self.desc['length'] - 1))
        upper_bound = np.copy(np.array([0., 0., 0.], dtype='float64')) + np.array([0, vector_module / 2, 0])
        points = [np.copy(upper_bound) - np.array([0, self._PERIODE * x, 0]) for x in range(self.desc['length'])]

        # 2. Build. For each point, we build one periode at [0, 0, 0]. Then, we rotate and then shift.
        self.pdb = []
        _MONO = pd.DataFrame(self._MONO).T
        for i, p in enumerate(points):
            coords = SBIgeo.rotate_degrees(_MONO.values, y=self._ROTATION * i)
            coords = SBIgeo.translate(coords, p)
            self.pdb.append(coords)
        self.pdb = np.vstack(self.pdb)
        # We want undirected structures to start always looking up
        self.pdb = SBIgeo.rotate_degrees(self.pdb, x=180)

        # Apply the case-defined placements for each structure
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('tilt: ' + str(self.desc['tilt']) + '\n')
            sys.stdout.write('move: ' + str(self.desc['coordinates']) + '\n')

        self.pdb = SBIgeo.rotate_degrees(self.pdb, x=self.desc['tilt']['x'],
                                         y=self.desc['tilt']['y'],
                                         z=self.desc['tilt']['z'])
        self.pdb = SBIgeo.translate(self.pdb, [self.desc['coordinates']['x'],
                                               self.desc['coordinates']['y'],
                                               self.desc['coordinates']['z']])

        # Prepare other data to create a coordinate entity
        resis = np.repeat(list(range(1, i + 2)), _MONO.shape[0])
        atoms = np.asarray([_MONO.index.values, ] * (i + 1)).flatten()

        # Prepare sequence
        sequence = []
        for _ in range(self.desc['length']):
            sequence.append(alphabet.aminoacids1to3(weighted_choice(self._AA_STAT)))
        sequence = np.repeat(np.asarray(sequence), _MONO.shape[0])

        self.pdb = PDB(pd.DataFrame(self.pdb, columns=["Cartn_x", "Cartn_y", "Cartn_z"])
                         .assign(auth_comp_id=sequence)
                         .assign(auth_atom_id=atoms).assign(auth_seq_id=resis)
                         .assign(id=list(range(1, self.pdb.shape[0] + 1))))

    def reverse( self ):
        """
        """
        raise NotImplementedError()
