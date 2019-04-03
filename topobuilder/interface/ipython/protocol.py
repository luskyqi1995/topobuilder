# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Union, List, Dict
import json

# External Libraries
import yaml

# This Library
from topobuilder import plugin_source
from topobuilder.case import Case
from topobuilder.utils import IpyExit
import topobuilder

__all__ = ['InteractiveProtocol']


class InteractiveProtocol( object ):
    def __init__( self, case: Union[Case, List[Case]], protocol: Union[Path, Dict] ):
        self.case = case
        if isinstance(protocol, Path):
            try:
                protocol = json.loads("".join([x.strip() for x in open(protocol).readlines()]))
            except json.JSONDecodeError:
                protocol = yaml.load(open(protocol))
        self.protocols = protocol
        self.current = -1

    def __next__( self ):
        if self.current + 1 >= len(self.protocols):
            raise StopIteration('Finished interactive protocols.')
        self.current += 1
        print(self.protocols[self.current])
