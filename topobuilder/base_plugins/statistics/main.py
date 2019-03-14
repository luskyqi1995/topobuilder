# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import sys
from pathlib import Path
from typing import List
from subprocess import run

# External Libraries
import pandas as pd

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder import PluginOrderError


__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           source: str,
           analysis: str,
           **kwargs ) -> List[Case]:
    """Perform a define statistic analysis on a given source.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: STATISTICS ---\n')

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i] = case_apply(case, source, analysis)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                source: str,
                analysis: str ) -> str:
    """
    """
    case = Case(case)

    # Generate the folder tree for a single connectivity.
    wfolder = case.connectivities_paths[0].joinpath('statistic')
    wfolder.mkdir(parents=True, exist_ok=True)
    # Generate internal folder
    thisfolder = wfolder.joinpath('_pdb_files')
    thisfolder.mkdir(parents=True, exist_ok=True)

    # Commands
    commands = []

    # Get data by source
    if source == 'funfoldes':
        commands.extend(funfoldes2pdb(case, thisfolder))

    # Load analysis commands
    if analysis == 'geometry':
        commands.append(geometry(case, wfolder, thisfolder))

    # Execute
    execute(commands, wfolder)

    # Postprocess
    postprocess(analysis, wfolder)

    if analysis == 'geometry':
        case['metadata'].setdefault('statistic',
                                    {}).setdefault('geometry',
                                                   str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))

    return case


def funfoldes2pdb( case: Case, wfolder: Path ) -> List:
    """
    """
    silent_files = case['metadata.funfoldes.silent_files']
    if silent_files is None:
        raise PluginOrderError('There is no output data from the funfoldes plugin.')

    extract_pdb = Path(str(Path(TBcore.get_option('rosetta', 'scripts')).resolve()).replace('rosetta_scripts.', 'extract_pdbs.'))
    if not extract_pdb.is_file() or not os.access(str(extract_pdb), os.X_OK):
        raise IOError('Cannot find executable {}'.format(extract_pdb))

    if not TBcore.get_option('slurm', 'use'):
        cmd = [extract_pdb, '-in:file:silent']
        cmd.extend(silent_files)
        cmd.extend(['-out:prefix', str(wfolder)])
    else:
        indir = str(wfolder.joinpath('${SLURM_ARRAY_TASK_ID}'))
        cmd = ['srun', extract_pdb, '-in:file:silent']
        cmd.append(os.path.commonprefix([str(x) for x in silent_files]) + '${SLURM_ARRAY_TASK_ID}.silent')
        cmd.extend(['-out:prefix', str(indir) + '/'])
    return [['mkdir', '-p', indir] if TBcore.get_option('slurm', 'use') else '', cmd]


def geometry( case: Case, wfolder: Path, thisfolder: Path ) -> List:
    """
    """
    cfile = case.write(wfolder.joinpath('current_case'))
    cmd = ['python', str(Path(__file__).parent.joinpath('geometry.py')), '-case', str(cfile), '-indir']
    if not TBcore.get_option('slurm', 'use'):
        cmd.append(str(thisfolder))
        cmd.extend(['-out', wfolder.joinpath('geometry.csv')])
    else:
        cmd.append(str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))
        cmd.extend(['-out', str(wfolder.joinpath('_geometry.${SLURM_ARRAY_TASK_ID}.csv'))])
    return cmd


def execute( cmd: List, wfolder: Path ):
    """
    """
    if not TBcore.get_option('slurm', 'use'):
        run(cmd)
    else:
        slurm_file = wfolder.joinpath('submit_analytics.sh')
        with slurm_file.open('w') as fd:
            fd.write(TButil.slurm_header() + '\n')
            fd.write(TButil.slurm_pyenv() + '\n')
            for c in cmd:
                fd.write(' '.join([str(x) for x in c]) + '\n')
        TButil.submit_slurm(slurm_file)


def postprocess( analysis: str, wfolder: Path ):
    """
    """
    if analysis == 'geometry':
        df = pd.concat([pd.read_csv(x) for x in Path(wfolder).glob('_geometry.*.csv')])
        df.to_csv(wfolder.joinpath('geometry.csv'), index=False)
