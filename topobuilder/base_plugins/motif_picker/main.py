# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict, Optional
from path import Path
import re
import copy

# External Libraries

# This Library
from topobuilder.case import Case
from topobuilder.case.schema import _ACCEPTED_SSE_ID_
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError


__all__ = ['motif_picker']


class motif_picker( Node ):
    """Recovers a motif of interest from a protein structure to map upon a :term:`FORM`.

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

    :param source: File name to the PDB file of interest containing the motif.

    :raises:
        :NodeOptionsError: On **initialization**. If the provided structure source file does not exist.
        :NodeOptionsError: On **initialization**. If the number of requested segments does not match the
            number of secondary structures to attach them.
        :NodeOptionsError: On **initialization**. If the number of requested segments does not match the
            shape definition.
        :NodeOptionsError: On **initialization**. If any secondary structure identifier in ``attach`` does
            not match *TopoBuilder* naming system.
        :NodeOptionsError: On **initialization**. If the ``shape`` definition format is unrecognized.
        :NodeOptionsError: On **initialization**. If amount of segments specified in ``extend_n`` and ``extend_c`` differ.
        :NodeOptionsError: On **initialization**. If amount of segments specified in ``extend_n`` and ``extend_c`` differ
            from the number of segments according to ``selection``.
        :NodeDataError: On **check**. If the required fields to be executed are not there.

    """
    SHAPE_PATTERN = r'^[HE]([lx][HE])*$'
    REQUIRED_FIELDS = ('architecture', )
    RETURNED_FIELDS = ('metadata.motif_picker')
    VERSION = 'v1.0'

    def __init__( self, tag: int, source: Union[Path, str], selection: str,
                  attach: List[str], shape: Optional[str] = None, binder: Optional[str] = None,
                  extend_n: Optional[List[int]] = None, extend_c: Optional[List[int]] = None,
                  identifier: Optional[str] = None ):
        super(motif_picker, self).__init__(tag)

        self.identifier = identifier if identifier is not None else self.tag
        # Obtaining structure input file
        self.source = source.resolve() if isinstance(source, Path) else Path(source).resolve()
        if not source.is_file():
            raise NodeOptionsError(f'Provided structure file {self.source} cannot be found.')

        # Check that the number of selection, attach and shape make sense
        self.selection = selection.split(',')
        self.attach = attach
        self.shape = shape.split()

        # Each segments belongs to a SSE
        if len(self.selection) != len(self.attach):
            err = 'The number of selection segments should match the number of structures to attach. There are '
            err += f'{len(self.selection)} segments requested to be assigned to {len(self.attach)} structures.'
            raise NodeOptionsError(err)
        # The shape defines propely the selection
        if len(self.shape) != (len(self.selection) * 2) - 1:
            raise NodeOptionsError(f'Shape definition {self.shape} does not match with the provided number of segments.')
        # The names of the SSE to attach the motif can exist.
        err = []
        for sse in self.attach:
            if not re.match(_ACCEPTED_SSE_ID_, sse):
                err.append(f'Unrecognized SSE identifier in {sse}.')
        if len(err) > 0:
            raise NodeOptionsError('\n'.join(err))
        # Check that the shape definition is sound
        if not re.match(self.SHAPE_PATTERN, self.shape):
            raise NodeOptionsError(f'The shape definition {self.shape} does not match the expected format.')

        # Capture binder selection
        self.binder = binder

        # Manage extension requests
        self.extend_n = [0, ] * len(self.selection) if extend_n is None else extend_n
        self.extend_c = [0, ] * len(self.selection) if extend_c is None else extend_c

        if len(self.extend_n) != len(self.extend_c):
            err = f'Number of segments expected by extend_n ({len(self.extend_n)}) does not match the number '
            err += f'of segments expected by extend_c ({len(self.extend_c)}).'
            raise NodeOptionsError(err)
        if len(self.extend_n) != len(self.segments):
            err = f'Number of extensions ({len(self.extend_n)}) does not match number of '
            err += f'requested segments ({len(self.segments)}).'
            raise NodeOptionsError(err)

        # Chatty
        self.log.debug(f'Picking a {len(self.segments)}-segment motif from {self.source} of shape {self.shape}.')
        self.log.debug(f'Motif segments are to be assigned to SSE(s) {",".join(self.attach)}.')
        if self.binder is not None:
            self.log.debug(f'A binder on selection {self.binder} has been included.')
        self.log.debug(f'Segments will be extended on their C-terminal by {",".join([str(x) for x in self.extend_c])}.')
        self.log.debug(f'Segments will be extended on their N-terminal by {",".join([str(x) for x in self.extend_n])}.')

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds
        kase.data.setdefault('metadata', {}).setdefault('motif_picker', [])
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)
        result = {'id': f'motif_{self.identifier}'}

        # Create a working folder
        folders = kase.undirected_path.joinpath(result['id'])
        folders.mkdir(parents=True, exist_ok=True)
        result['data_dir'] = str(folders)

        # Check name was not already added.
        sn = copy.deepcopy(self.subnames)
        if kase.name.endswith('_'.join(sn)):
            self.log.notice(f'Seems the subnames {"_".join(sn)} already existed.')
            self.log.notice('Will NOT re-append.')
            return kase

        # Add new names
        oname = kase.name
        sn.insert(0, oname)
        kase.data['configuration']['name'] = '_'.join(sn)

        self.log.debug(f'Renamed case {oname} to {kase.name}')

        # Attach data and return
        kase.data.setdefault('metadata', {}).setdefault('motif_picker', []).append(result)
        return kase.data
