# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, Union, Dict
from pathlib import Path, PosixPath
import json
import sys
from types import GeneratorType

# External Libraries

# This Library
import topobuilder.core as TBcore

___all__ = ['checkpoint_in', 'checkpoint_out']


class GeneralEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, PosixPath):
            return str(o)
        if isinstance(o, GeneratorType):
            return list(o)
        return json.JSONEncoder.default(self, o)


def checkpoint_in( filename: Union[Path, str] ) -> Optional[Dict]:
    """
    """
    if TBcore.get_option('system', 'forced'):
        return None

    filename = Path(filename)
    if filename.is_file():
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('CHECKPOINT: Reloading from {}\n'.format(filename))
        with Path(filename).open() as fd:
            data = json.load(fd)
        return data

    return None


def checkpoint_out( filename: Union[Path, str], data: Dict ):
    """
    """
    filename = Path(filename)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('CHECKPOINT: Creating at {}\n'.format(filename))
    with filename.open('w') as fd:
        json.dump(data, fd, cls=GeneralEncoder)
