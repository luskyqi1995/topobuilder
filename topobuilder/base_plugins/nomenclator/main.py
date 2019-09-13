# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict
import copy

# External Libraries

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError
import topobuilder.utils as TButil


__all__ = ['nomenclator']


class nomenclator( Node ):
    """Alters the ``configuration.name`` of a :class:`.Case` by adding subnames.

    This affects on the creation of the subfolders where the rest of the :class:`.Pipeline`
    will be executed.
    """

    RESERVED_KEYWORDS = ['architecture', 'connectivity', 'images', 'summary']

    def __init__( self, tag: int, subnames: Union[List, str] ):
        super(nomenclator, self).__init__(tag)

        self.subnames = subnames
        # Unify subnames behaviour for 1 to N
        if not isinstance(self.subnames, list):
            self.subnames = [self.subnames, ]

        # Words used as main folders should not be used for subnaming
        if len(set([x.lower() for x in self.subnames]).intersection(set(self.RESERVED_KEYWORDS))) > 0:
            raise NodeOptionsError(f'Keywords {self.RESERVED_KEYWORDS} cannot be used for subnaming.')

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in ('configuration.name', ):
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)

        # Check name was not already added.
        sn = copy.deepcopy(self.subnames)
        if kase.name.endswith('_'.join(sn)):
            TButil.plugin_warning('Seems the subnames {} already existed.'.format('_'.join(sn)))
            TButil.plugin_warning('Will NOT re-append.')
            return kase

        # Add new names
        oname = kase.name
        sn.insert(0, oname)
        kase.data['configuration']['name'] = '_'.join(sn)

        self.log.debug(f'Renamed case {oname} to {kase.name}')
        return kase.data
