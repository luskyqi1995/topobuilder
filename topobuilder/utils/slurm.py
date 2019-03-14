# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys
import os
import time
import math
import textwrap
import tempfile
from pathlib import Path
from typing import Union, Optional
import subprocess

# External Libraries

# This Library
import topobuilder.core as TBcore


__all__ = ['slurm_header', 'slurm_pyenv', 'submit_slurm', 'submit_nowait_slurm']


def submit_slurm( slurm_file: Union[Path, str],
                  condition_file: Optional[Union[Path, str]] = None ):
    """
    """
    slurm_control_file = (Path(tempfile.mkdtemp('slurm_control'))
                          .joinpath('slurm_control.{}.sh'.format(os.getpid())))
    condition_file = control_slurm_file(slurm_control_file, condition_file)

    main_id = submit_nowait_slurm(slurm_file)
    submit_nowait_slurm(slurm_control_file, 'afterany', main_id)

    wait_for(condition_file)
    os.unlink(str(condition_file))


def submit_nowait_slurm( slurm_file: Union[Path, str],
                         dependency_mode: Optional[str] = None,
                         dependency_id: Optional[int] = None
                         ) -> int:
    """
    """
    command = ['sbatch']
    if dependency_mode is not None and dependency_id is not None:
        command.append('--depend={}:{}'.format(dependency_mode, dependency_id))
    command.append('--parsable')
    command.append(str(slurm_file))
    p = subprocess.run(command, stdout=subprocess.PIPE)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write(' '.join(command) + '\n')
    return int(str(p.stdout.decode("utf-8")).strip())


def wait_for( condition_file: Optional[Union[Path, str]] ):
    """
    """
    waiting_time = 0
    while not Path(condition_file).is_file():
        time.sleep(30)
        waiting_time += 1

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Total waiting time: {}h:{}m\n'.format(math.floor(waiting_time / 60),
                                                                waiting_time % 60))


def control_slurm_file( slurm_file: Union[Path, str],
                        condition_file: Optional[Union[Path, str]] = None
                        ) -> Path:
    """SLURM submission file that will touch a control file.

    The submission file will be dependant on another running job, thus indicating when
    the job has finished.

    :param slurm_file: Name of the expected SLURM submission file.
    :param condition_file: Name of the file to create in order to assert completion of
        a previously dependant job.

    :return: Name of the ``condition_file``.
    """
    with TBcore.on_option_value('slurm', 'array', 0, 'slurm', 'memory', 2046):
        header = slurm_header()
    condition_file = condition_file if condition_file is not None else 'touch_control.{}'.format(os.getpid())
    condition_file = Path().cwd().joinpath(condition_file)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('A SLURM control file will be generated at {}\n'.format(condition_file))

    with open(slurm_file, 'w') as fd:
        fd.write(header + '\n\n')
        fd.write('echo \'finished\' > {}\n'.format(condition_file.resolve()))
    return condition_file


def slurm_header() -> str:
    """
    """
    config = {
        'slurm_array':     TBcore.get_option('slurm', 'array'),
        'partition':       TBcore.get_option('slurm', 'partition'),
        'logpath':         TBcore.get_option('slurm', 'logs'),
        'memory':          TBcore.get_option('slurm', 'memory'),
        'nodes':           TBcore.get_option('slurm', 'nodes'),
        'time':            TBcore.get_option('slurm', 'time'),
        'ntasks-per-node': TBcore.get_option('slurm', 'ntasks-per-node'),
        'cpus-per-task':   TBcore.get_option('slurm', 'cpus-per-task'),
        'sublog':          '.%a'
    }

    config['slurm_array'] = '' if config['slurm_array'] == 0 else "#SBATCH --array=1-{}\n".format(config['slurm_array'])
    config['sublog'] = '' if config['slurm_array'] == 0 else config['sublog']

    if Path(config['logpath']).is_dir():
        config['logpath'] = str(Path(config['logpath']).resolve().joinpath('_output'))
    return textwrap.dedent("""\
        #!/bin/bash
        #SBATCH --nodes {nodes}
        #SBATCH --partition={partition}
        #SBATCH --ntasks-per-node {ntasks-per-node}
        #SBATCH --cpus-per-task {cpus-per-task}
        #SBATCH --mem {memory}
        #SBATCH --time {time}
        {slurm_array}#SBATCH --output={logpath}.%A.{sublog}.out
        #SBATCH --error={logpath}.%A.{sublog}.err
    """).format(**config)


def slurm_pyenv() -> str:
    """
    """
    pypath = Path(sys.executable).parent.joinpath('activate')
    if not os.access(pypath, os.X_OK):
        return ''
    else:
        return '\n'.join(['source {}'.format(pypath), "export PYTHONPATH=''"])

