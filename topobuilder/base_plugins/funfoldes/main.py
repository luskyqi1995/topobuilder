# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
from pathlib import Path
from typing import Optional, List, Dict
import sys
# import textwrap

# External Libraries

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from . import utils


__all__ = ['metadata', 'apply', 'case_apply']


def metadata() -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: Funtion on the expected input connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    """
    def isngl( count ):
        return count == 1

    return {'name': 'funfoldes',
            'Itags': ['loop_lengths', 'fragments'],
            'Otags': ['funfoldes'],
            'Isngl': isngl,
            'Osngl': True}


def apply( cases: List[Case],
           prtid: int,
           folding_nstruct: Optional[int] = 2000,
           design_nstruct: Optional[int] = 10,
           natbias: Optional[float] = 2.5,
           layer_design: Optional[bool] = True,
           #folding_script: Optional[str] = '',
           design_script = None,
           **kwargs ) -> List[Case]:
    """Execute the FunFolDes Rosetta protocol.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: FUNFOLDES ---\n')

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('funfoldes', '')
        cases[i] = case_apply(case, folding_nstruct, design_nstruct, natbias, design_script)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


@TButil.plugin_conditions(metadata())
def case_apply( case: Case,
                design_script,
                folding_nstruct: Optional[int] = 2000,
                design_nstruct: Optional[int] = 10,
                natbias: Optional[float] = 2.5,
                layer_design: Optional[bool] = True,
                #folding_script: Optional[str] = '',
                ) -> Case:
    """Execute the FunFolDes Rosetta protocol.
    """
    # Securing data and preparing metadata output structure
    case = Case(case)
    data = {'script': {'folding': '', 'design': ''},
            'cmd':
                {'folding': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol'],
                 'design': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol']},
            'silent_files': {'folding': [], 'design': []},
            'minisilent': {'folding': '', 'design': ''}
            }

    # Generate the folder tree for a single connectivity.
    wpaths = utils.folder_structure(case)

    # Check if checkpoint exists, retrieve and skip
    reload = TButil.checkpoint_in(wpaths['checkpoint'])
    if reload is not None:
        case.data['metadata']['funfoldes'] = reload
        return case

    # We need to check that the rosetta_scripts executable is available
    if not data['cmd']['folding'][0].is_file() or not os.access(str(data['cmd']['folding'][0]), os.X_OK):
        raise IOError('Cannot find executable {}'.format(data['cmd']['folding'][0]))

    # Build the structure
    utils.build_template_sketch( case, wpaths['pdb'] )

    # Make the Folding and Design RScripts
    data = utils.make_scripts(case, wpaths, design_script, data, natbias, layer_design)

    # Finish command
    data = utils.commands(case, folding_nstruct, design_nstruct, data, wpaths)

    # Execute
    data = utils.execute(data, wpaths)

    # Update metadata
    data = utils.update_data( data, wpaths )

    # Checkpoint save
    TButil.checkpoint_out(wpaths['checkpoint'], data)
    case.data['metadata']['funfoldes'] = data

    return case
