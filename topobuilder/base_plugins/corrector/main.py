# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict, Optional
from pathlib import Path
import copy
# import sys

# External Libraries

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeDataError
# import topobuilder.core as TBcore
# import topobuilder.utils as TButil


__all__ = ['corrector']


class corrector( Node ):
    """Applies corrections to the placements of the secondary structures in a :term:`FORM`.

    This affects on the creation of the subfolders where the rest of the :class:`.Pipeline`
    will be executed.

    .. note::
        On **execution**, the plugin will not append new subnames when those already exist. For example,
        if ``configuration.name`` is ``1QYS_experiment1_naive`` and ``subnames=['experiment1', 'naive']``,
        the final ``configuration.name`` will still be ``1QYS_experiment1_naive`` and not
        ``1QYS_experiment1_naive_experiment1_naive``. This is to avoid folder recursion generation when
        re-running a previous :class:`.Pipeline`.

    .. caution::
        There are some keywords that cannot be used as a subname due to them generating their own
        **first level** subfolders. These keywords are ``architecture``, ``connectivity``, ``images``
        and ``summary``. Trying to add one of those terms as subname will generate a :class:`.NodeDataError`
        on **check** time.

    .. admonition:: To Developers

        When developing a new plugin, if it is expected to create new **first level** subfolders, they should
        be listed in the class attribute :attr:`.nomenclator.RESERVED_KEYWORDS`. See more on how to
        :ref:`develop your own plugins <make_plugin>`.

    :param subnames: Subnames that will be added to the :class:`.Case` initial name.

    :raises:
        :NodeOptionsError: On **initialization**. If a reserved key is provided as a subname.
        :NodeDataError: On **check**. If the required fields to be executed are not there.

    """
    REQUIRED_FIELDS = ('architecture', )
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  corrections: Optional[Union[str, Dict, Path, List[Union[str, Path]]]] = None ):
        super(corrector, self).__init__(tag)

        # Make sure we have a list of corrections.
        self.corrections = corrections
        if self.corrections is not None:
            if not isinstance(self.corrections, list):
                self.corrections = [self.corrections, ]
            for i, c in enumerate(self.corrections):
                if isinstance(c, str):
                    self.corrections[i] = Path(c)
        else:
            corrections = []

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)

        crr = copy.deepcopy(self.corrections)

        # See if there are extra corrections attached to the case itself
        krr = kase['metadata.corrections']
        krr = [] if krr is None else krr
        crr.extend(krr)

        # Apply each set of corrections
        for c in crr:
            self.log.info('Applying correction: {0}\n'.format(c))
            kase = kase.apply_corrections(c)

        self.log.debug(f'Applied a total of {len(crr)} corrections.')
        self.log.debug(f'{len(krr)} from within the Case definition.')
        self.log.debug(f'{len(self.corrections)} from protocol-provided data.')
        return kase.data
