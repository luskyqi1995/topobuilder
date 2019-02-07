# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: protocol
"""
# Standard Libraries
from typing import Union, Dict, Optional
from pathlib import Path
import json

# External Libraries
import yaml

# # This Library
from topobuilder import plugin_source
from topobuilder.case import Case

__all__ = ['protocol']


def protocol( case: Union[str, Path, Dict, Case],
              protocol: Optional[Union[str, Path]] = None,
              overwrite: Optional[bool] = False ):
    """
    """
    c = Case(case)

    protocols = c['configuration.protocols']
    if protocols is not None:
        if len(protocols) == 1 and not bool(protocols[0]):
            protocols = None
    if protocols is None and protocol is None:
        raise AttributeError('There are no protocols to run')
    if protocol is not None and protocols is not None:
        raise AttributeError('Protocols are provided both through file and in the Case. '
                             'Pick one.')
    if protocol is not None:
        protocol = str(Path(protocol).resolve())
        try:
            protocols = json.loads("".join([x.strip() for x in open(protocol).readlines()]))
        except json.JSONDecodeError:
            protocols = yaml.load(open(protocol))

    # Check requested plugins (avoid time is something is wrong)
    for i, ptcl in enumerate(protocols):
        if 'name' not in ptcl:
            raise ValueError('All protocols require a "name" field')
        if ptcl['name'] not in plugin_source.list_plugins():
            raise ValueError('Requested protocol {} cannot be found.'.format(ptcl['name']))
        protocols[i].setdefault('status', False)

    cases = [c.assign_protocols(protocols), ]
    for i, ptcl in enumerate(protocols):
        if not ptcl['status']:
            cases = plugin_source.load_plugin(ptcl['name']).apply(cases, prtid=i, **ptcl)

    # for c in cases:
    #     c.write()
