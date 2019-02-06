# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Union, Dict, List, Optional
from pathlib import Path
import textwrap
import math
import copy

# External Libraries

# This Library
from topobuilder.case import Case
from .core import core


def make_slurm_file( pds_list: List ) -> str:
    """
    """
    createPDS = str(Path(core.get_option('imaster', 'createPDS')).absolute())
    master = str(Path(core.get_option('imaster', 'master')).absolute())
    slurm_array = core.get_option('imaster', 'slurm.array')

    outstr = [slurm_header(slurm_array)]
    outstr.append('if [ ! -f {0}.pds ]; then')
    outstr.append('  if (( ${SLURM_ARRAY_TASK_ID} == 0 )); then')
    outstr.append('    {0} --type query --pdb {{0}}.pdb --pds {{0}}.pds'.format(createPDS))
    outstr.append('  else')
    outstr.append('    sleep 30s')
    outstr.append('  fi')
    outstr.append('fi')
    outstr.append('mkdir -p {0}_master')
    sp = int(math.ceil(len(pds_list) / slurm_array))
    mcom = [master, '--query {0}.pds', '--target {0}', '--rmsdCut 5',
            '--matchOut {{0}}_master/{{0}}.{0}.master', '--topN 1']
    for i, c in enumerate(pds_list[::sp]):
        outstr.append('if (( ${{SLURM_ARRAY_TASK_ID}} == {} )); then'.format(i))
        outfiles = []
        for cc in pds_list[sp * i: (sp * i) + sp]:
            pds = cc.split('/')[-1].split('.')[0]
            outstr.append(copy.deepcopy(mcom))
            outstr[-1][2] = outstr[-1][2].format(cc)
            outstr[-1][4] = outstr[-1][4].format(pds)
            outfiles.append(outstr[-1][4].split()[-1])
            outstr[-1] = '  ' + ' '.join(outstr[-1])
        outstr.append('  cat {0} > {{0}}_master/step{1:03d}'.format(' '.join(outfiles), i))
        outstr.append('  rm {0}'.format(' '.join(outfiles)))

        outstr.append('fi')

    return '\n'.join(outstr)


def slurm_header( slurm_array: int ) -> str:

    partition = core.get_option('imaster', 'slurm.partition')
    logpath = core.get_option('imaster', 'slurm.logs')

    if Path(logpath).is_dir():
        logpath = Path(logpath).resolve().joinpath('output')
    return textwrap.dedent("""\
        #!/bin/bash
        #SBATCH --nodes 1
        #SBATCH --partition={0}
        #SBATCH --ntasks-per-node 1
        #SBATCH --cpus-per-task 1
        #SBATCH --mem 4096
        #SBATCH --time 10:00:00
        #SBATCH --array=0-{2}
        #SBATCH --output={1}.%A.out
        #SBATCH --error={1}.%A.err
    """).format(partition, logpath, slurm_array)
