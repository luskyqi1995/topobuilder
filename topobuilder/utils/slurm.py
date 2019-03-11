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


__all__ = ['submit_slurm', 'submit_nowait_slurm']


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
        time.sleep(60)
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
    header = slurm_header(0, 2046)
    condition_file = condition_file if condition_file is not None else 'touch_control.{}'.format(os.getpid())
    condition_file = Path().cwd().joinpath(condition_file)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('A SLURM control file will be generated at {}\n'.format(condition_file))

    with open(slurm_file, 'w') as fd:
        fd.write(header + '\n\n')
        fd.write('echo \'finished\' > {}\n'.format(condition_file.resolve()))
    return condition_file


def slurm_header( slurm_array: int,
                  memory: Optional[int] = 4096
                  ) -> str:
    """
    """
    partition = TBcore.get_option('slurm', 'partition')
    logpath = TBcore.get_option('slurm', 'logs')

    slurm_array = "" if slurm_array == 0 else "#SBATCH --array=1-{2}\n".format(slurm_array)

    if Path(logpath).is_dir():
        logpath = Path(logpath).resolve().joinpath('output')
    return textwrap.dedent("""\
        #!/bin/bash
        #SBATCH --nodes 1
        #SBATCH --partition={0}
        #SBATCH --ntasks-per-node 1
        #SBATCH --cpus-per-task 1
        #SBATCH --mem {3}
        #SBATCH --time 10:00:00
        {2}#SBATCH --output={1}.%A.out
        #SBATCH --error={1}.%A.err
    """).format(partition, logpath, slurm_array, memory)
