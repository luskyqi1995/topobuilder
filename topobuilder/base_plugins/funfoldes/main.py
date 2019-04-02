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
from typing import Optional, List, Union, Dict, Tuple
import sys
import textwrap
from subprocess import run, DEVNULL

# External Libraries
from rstoolbox.io import open_rosetta_file

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder.utils import build_pdb_object
from topobuilder import plugin_source, PluginOrderError


__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           nstruct: Optional[int] = 2000,
           natbias: Optional[float] = 2.5,
           layer_design: Optional[bool] = True,
           **kwargs ) -> List[Case]:
    """Execute the FunFolDes Rosetta protocol.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: FUNFOLDES ---\n')

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('funfoldes', '')
        cases[i] = case_apply(case, nstruct, natbias)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


def case_apply( case: Case,
                nstruct: Optional[int] = 2000,
                natbias: Optional[float] = 2.5,
                layer_design: Optional[bool] = True
                ) -> Case:
    """
    """
    # Securing data and preparing metadata output structure
    case = Case(case)
    data = {'script': {'folding': '', 'design': ''},
            'cmd':
                {'folding': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol'],
                 'design': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol']},
            'silent_files': {'folding': '', 'design': ''},
            'minisilent': {'folding': '', 'design': ''}
            }

    # FunFolDes is only applied to a Case with:
    #   (1) one single connectivity
    #   (2) provided loop lengths
    #   (3) provided fragments to fold
    if case.connectivity_count > 1:
        raise ValueError('FunFolDes can only be applied to one connectivity.')
    if case['metadata.loop_lengths'] is None:
        raise PluginOrderError(textwrap.dedent("""\
            metadata.loop_lengths needs to be specified.
            Check plugins capable of providing this information."""))
    if case['metadata.fragments'] is None:
        raise PluginOrderError(textwrap.dedent("""\
            Fragments to proceed to FunFolDes cannot be found at metadata.fragments.
            Check fragment_maker."""))

    # Generate the folder tree for a single connectivity.
    wpaths = folder_structure(case)

    # Check if checkpoint exists, retrieve and skip
    reload = TButil.checkpoint_in(wpaths['checkpoint'])
    if reload is not None:
        case.data['metadata']['funfoldes'] = reload
        return case

    # We need to check that the rosetta_scripts executable is available
    if not data['cmd']['folding'][0].is_file() or not os.access(str(data['cmd']['folding'][0]), os.X_OK):
        raise IOError('Cannot find executable {}'.format(data['cmd']['folding'][0]))

    # Build the structure
    build_template_sketch( case, wpaths['pdb'] )

    # Make the Folding and Design RScripts
    data = make_scripts(case, wpaths, data, natbias, layer_design)

    # Finish command
    data = commands(case, nstruct, data, wpaths)

    # Execute
    execute(data, wpaths)

    data['silent_files']['folding'] = list(wpaths['outdir'].glob('*_funfol.silent'))
    data['silent_files']['design'] = list(wpaths['outdir'].glob('*_des.silent'))

    data['minisilent']['folding'] = wpaths['main'].joinpath('output_funfol.minisilent.gz')
    data['minisilent']['design'] = wpaths['main'].joinpath('output_des.minisilent.gz')
    for k in data['minisilent']:
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Generating minisilent file at {}\n'.format(data['minisilent'][k]))
        fd = gzip.open( data['minisilent'][k], "wb" )
        for line, _, _, _ in open_rosetta_file([str(x) for x in data['silent_files'][k]], True, check_symmetry=False ):
            fd.write(line.encode('utf-8'))

    # Checkpoint save
    TButil.checkpoint_out(wpaths['checkpoint'], data)
    case.data['metadata']['funfoldes'] = data

    return case


def folder_structure( case: Case ) -> Dict:
    """
    """
    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0].joinpath('funfoldes')
    folders.mkdir(parents=True, exist_ok=True)
    outdir = folders.joinpath('outputs')
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_file = folders.joinpath('template_sketch.pdb')
    ffd_fold_file = folders.joinpath('funfoldes_fold.xml')
    ffd_design_file = folders.joinpath('funfoldes_design.xml')
    checkpoint = folders.joinpath('checkpoint.json')

    return {'main': folders, 'outdir': outdir, 'pdb': pdb_file,
            'foldRS': ffd_fold_file, 'designRS': ffd_design_file,
            'checkpoint': checkpoint}


def build_template_sketch( case: Case, pdb_file: Union[Path, str] ):
    """
    """
    # 1.1 Make sure the loop_length are specified.
    loop_lengths = case['metadata.loop_lengths']
    if TBcore.get_option('system', 'verbose'):
        looplist = ', '.join([str(x) for x in loop_lengths])
        sys.stdout.write('Gathering a total of {} loops of length: {}\n'.format(len(loop_lengths), looplist))
        sys.stdout.write('To apply over {} secondary structures: {}\n'.format(len(case), case.connectivities_str[0]))

    # Get the structure.
    case = plugin_source.load_plugin('builder').case_apply(case, connectivity=True, pick_aa='V')
    sse_list = case.ordered_structures
    pdb, _ = build_pdb_object(sse_list, loop_lengths)
    pdb.write(str(pdb_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))


def make_scripts( case: Case,
                  wpaths: Dict,
                  data: Dict,
                  natbias: float = 2.5,
                  layer_design: bool = True
                  ) -> Tuple[str, str]:
    """
    """
    fld = TButil.rosettascript(TButil.funfoldes(case))
    dsg = TButil.rosettascript(TButil.constraint_design(case, natbias, layer_design))

    if TBcore.get_option('system', 'jupyter'):
        ifold = os.getenv('TB_FUNFOLDES_FOLD_FILE', None)
        idsgn = os.getenv('TB_FUNFOLDES_DSGN_FILE', None)

        if ifold is None:
            print('-' * 80)
            print('\n\nExpected FOLDING script will be:\n')
            print(fld)
            print(textwrap.dedent("""\n\n If a different script wants to be provided, do so by
    assigning a RosettaScript to os.environ['TB_FUNFOLDES_FOLD_FILE'] or assign '' to
    it if you are ok with the default script."""))
        elif ifold:
            ifold = Path(ifold)
            if not ifold.is_file():
                raise IOError('Unknown file {}'.format(ifold))
            fld = ''.join(list(ifold.open().readlines()))

        if idsgn is None:
            print('-' * 80)
            print('\n\nExpected DESIGN script will be:\n')
            print(dsg)
            print(textwrap.dedent("""\n\n If a different script wants to be provided, do so by
    assigning a RosettaScript to os.environ['TB_FUNFOLDES_DSGN_FILE'] or assign '' to
    it if you are ok with the default script."""))
        elif idsgn:
            idsgn = Path(idsgn)
            if not idsgn.is_file():
                raise IOError('Unknown file {}'.format(idsgn))
            dsg = ''.join(list(idsgn.open().readlines()))

        if ifold is None or idsgn is None:
            TButil.exit()

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing the folding RosettaScript file: {}\n'.format(wpaths['foldRS']))
    with wpaths['foldRS'].open('w') as fd:
        fd.write(fld)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing the design RosettaScript file: {}\n'.format(wpaths['designRS']))
    with wpaths['designRS'].open('w') as fd:
        fd.write(dsg)

    data['script']['folding'] = wpaths['foldRS']
    data['script']['design'] = wpaths['designRS']
    data['cmd']['folding'].append(wpaths['foldRS'])
    data['cmd']['design'].append(wpaths['designRS'])
    return data


def commands( case: Case, nstruct: int, data: Dict, wpaths: Dict ) -> Dict:
    """
    """
    out_prefix = (case.name if not TBcore.get_option('slurm', 'use') else '_'.join([case.name, '${SLURM_ARRAY_TASK_ID}'])) + '_'
    nstruct = nstruct if not TBcore.get_option('slurm', 'use') else math.ceil(nstruct / TBcore.get_option('slurm', 'array'))

    commons = ['-overwrite', '-in:ignore_unrecognized_res', '-in:ignore_waters', '-out:file:silent_struct_type',
               'binary', '-out:mute', 'protocols.abinitio', 'protocols.moves', 'core.optimization']

    flded = str(wpaths['outdir'].joinpath(out_prefix + 'funfol')) + '.silent'
    dsgnd = str(wpaths['outdir'].joinpath(out_prefix + 'des')) + '.silent'
    prefix1 = out_prefix + 'funfol_'
    prefix2 = out_prefix + 'des_'
    data['cmd']['folding'].extend(['-in:file:s', str(wpaths['pdb']), '-out:prefix', prefix1, '-out:file:silent', flded])
    data['cmd']['design'].extend(['-in:file:silent', flded, '-out:prefix', prefix2, '-out:file:silent', dsgnd])
    data['cmd']['folding'].extend(['-nstruct', str(nstruct)])
    data['cmd']['design'].extend(['-nstruct', str(10)])
    data['cmd']['folding'].extend(commons)
    data['cmd']['design'].extend(commons)
    return data


def execute(data: Dict, wpaths: Dict):
    """
    """
    if TBcore.get_option('slurm', 'use'):
        slurm_file = wpaths['main'].joinpath('submit_funfoldes.sh')
        TButil.plugin_filemaker('Submission file at {}'.format(slurm_file))
        with slurm_file.open('w') as fd:
            fd.write(TButil.slurm_header() + '\n' )
            for k in ['folding', 'design']:
                cmd = ['srun', ]
                cmd.extend(data['cmd'][k])
                fd.write(' '.join([str(x) for x in cmd]) + '\n')

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Submiting jobs to SLURM... this might take a while\n')
        TButil.submit_slurm(slurm_file)
    else:
        for k in ['folding', 'design']:
            TButil.plugin_bash(data['cmd'][k])
            run([str(x) for x in data['cmd'][k]], stdout=DEVNULL)
