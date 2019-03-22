# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import math
from pathlib import Path
import gzip
from typing import Optional, List, Union
import sys
from itertools import chain, zip_longest

# External Libraries
from rstoolbox.io import open_rosetta_file

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder.utils import build_pdb_object
from topobuilder import plugin_source, PluginOrderError


__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           nstruct: int = 2000,
           script: Optional[Union[Path, str]] = None,
           **kwargs ) -> List[Case]:
    """Execute the FunFolDes Rosetta protocol.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: FUNFOLDES ---\n')

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('funfoldes', '')
        cases[i] = case_apply(case, script)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


def case_apply( case: Case,
                script: Optional[Union[Path, str]] = None,
                nstruct: int = 2000
                ) -> Case:
    """
    """
    case = Case(case)
    data = {'script': '', 'cmd': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol'],
            'silent_files': [], 'minisilent': ''}

    # FunFolDes is only applied to a Case with one single connectivity and already reoriented
    if case.connectivity_count > 1:
        raise ValueError('FunFolDes can only be applied to one connectivity.')

    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0].joinpath('funfoldes')
    folders.mkdir(parents=True, exist_ok=True)
    outdir = folders.joinpath('outputs')
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_file = folders.joinpath('template_sketch.pdb')
    ffd_file = folders.joinpath('funfoldes.xml')
    checkpoint = folders.joinpath('checkpoint.json')

    # Check if checkpoint exists, retrieve and skip
    reload = TButil.checkpoint_in(checkpoint)
    if reload is not None:
        case.data['metadata']['funfoldes'] = reload
        return case

    # We need to check that the rosetta_scripts executable is available
    if not data['cmd'][0].is_file() or not os.access(str(data['cmd'][0]), os.X_OK):
        raise IOError('Cannot find executable {}'.format(data['cmd'][0]))

    # 1. Build the structure
    # 1.1 Make sure the loop_length are specified.
    loop_lengths = case['metadata.loop_lengths']
    if TBcore.get_option('system', 'verbose'):
        looplist = ', '.join([str(x) for x in loop_lengths])
        sys.stdout.write('Gathering a total of {} loops of length: {}\n'.format(len(loop_lengths), looplist))
        sys.stdout.write('To apply over {} secondary structures: {}\n'.format(len(case), case.connectivities_str[0]))

    if loop_lengths is None or len(loop_lengths) != len(case) - 1:
        raise ValueError('Specified length of the loops is necessary to execute FunFolDes.')

    # Get the structure.
    case = plugin_source.load_plugin('builder').case_apply(case, connectivity=True, pick_aa='V')
    sse_list = case.ordered_structures
    pdb, edges = build_pdb_object(sse_list, loop_lengths)
    pdb.write(str(pdb_file), format='pdb', clean=True,
              force=TBcore.get_option('system', 'overwrite'))

    # Get the fragments.
    frags = case['metadata.fragments']
    if frags is None:
        raise PluginOrderError('Fragments to proceed to FunFolDes cannot be found. Check fragment_maker.')

    # Make the RScript
    dssp = ''.join([x for x in chain(*zip_longest([''.join([x['type'], ] * x['length']) for x in sse_list],
                                                  [''.join(['L', ] * x) for x in loop_lengths])) if x is not None])
    ffd = TButil.funfoldes_script(('SSE', dssp),
                                  ('looper', zip(edges[:-1], loop_lengths)),
                                  (str(frags['files'][0]), str(frags['files'][1])))
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing the RosettaScript file: {}\n'.format(ffd_file))
    with ffd_file.open('w') as fd:
        fd.write(ffd)
    data['script'] = ffd_file
    data['cmd'].append(data['script'])

    # Finish command
    out_prefix = case.name if not TBcore.get_option('slurm', 'use') else '_'.join([case.name, '${SLURM_ARRAY_TASK_ID}'])
    nstruct = nstruct if not TBcore.get_option('slurm', 'use') else math.ceil(nstruct / TBcore.get_option('slurm', 'array'))
    commons = ['-overwrite', '-in:ignore_unrecognized_res', '-in:ignore_waters', '-out:file:silent_struct_type',
               'binary', '-out:mute', 'protocols.abinitio', 'protocols.moves', 'core.optimization']

    data['cmd'].extend(['-in:file:s', str(pdb_file), '-out:prefix', out_prefix,
                        '-out:file:silent', str(outdir.joinpath(out_prefix)) + '.silent'])
    data['cmd'].extend(['-nstruct', str(nstruct)])
    data['cmd'].extend(commons)

    if TBcore.get_option('slurm', 'use'):
        slurm_file = folders.joinpath('submit_funfoldes.sh')
        cmd = ['srun', ]
        cmd.extend(data['cmd'])
        with slurm_file.open('w') as fd:
            fd.write(TButil.slurm_header() + '\n' )
            fd.write(' '.join([str(x) for x in cmd]) + '\n')

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Submiting jobs to SLURM... this might take a while\n')
        TButil.submit_slurm(slurm_file)

    data['silent_files'] = list(outdir.glob('*silent'))

    data['minisilent'] = folders.joinpath('output.minisilent.gz')
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Generating minisilent file at {}\n'.format(data['minisilent']))
    fd = gzip.open( data['minisilent'], "wb" )
    for line, _, _, _ in open_rosetta_file( [str(x) for x in data['silent_files']], True, check_symmetry=False ):
        fd.write(line.encode('utf-8'))

    # Checkpoint save
    TButil.checkpoint_out(checkpoint, data)
    case.data['metadata']['funfoldes'] = data

    return case
