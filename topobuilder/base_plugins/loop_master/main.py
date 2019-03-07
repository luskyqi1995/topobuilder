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
from typing import List, Tuple, Dict
import math
from tempfile import NamedTemporaryFile
import shlex
from subprocess import run
import gzip
import itertools


# External Libraries
import pandas as pd
from pandas.compat import StringIO
from SBI.structure import PDB, PDBFrame, ChainFrame
import SBI.structure.geometry as SBIgeo
import seaborn as sns
import matplotlib.pyplot as plt
import rstoolbox

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
from topobuilder import plugin_source
from .core import core


__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           loop_range: int = 3,
           top_loops: int = 20,
           harpins_2: bool = True,
           **kwargs ) -> List[Case]:
    """Use MASTER to cover the transitions between secondary structures.

    And something else.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: LOOP_MASTER ---\n')

    # Get list of PDS structures
    tempdb = False
    database = core.get_option('master', 'pds')
    database = Path(database)
    if database.is_file():
        pass
    elif database.is_dir():
        f = NamedTemporaryFile(mode='w', delete=False)
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('Temporary file for PDS database: {}\n'.format(f.name))
        [f.write(str(x.resolve()) + '\n') for x in database.glob('*/*.pds')]
        f.close()
        database = Path(f.name)
        tempdb = True
    else:
        raise ValueError('The provided MASTER database directory/list file cannot be found.')

    # Get ABEGOS
    abegodata = get_abegos()

    # Get FragFiles
    fragfiles = get_fragfiles()

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('loop_fragments', [])
        cases[i] = case_apply(case, database, loop_range, top_loops, abegodata, harpins_2, fragfiles)
        cases[i].set_protocol_done(prtid)

    if tempdb:
        os.unlink(f.name)

    return cases


def case_apply( case: Case,
                pds_list: Path,
                loop_range: int,
                top_loops: int,
                abego: pd.DataFrame,
                harpins_2: bool,
                fragfiles: pd.DataFrame ) -> str:
    """
    """
    # Loop MASTER is only applied to a Case with one single connectivity and already reoriented
    if case.connectivity_count > 1:
        raise ValueError('Loop MASTER can only be applied to one connectivity.')
    if not case['configuration.reoriented'] and case.connectivity_count == 1:
        # We will need the coordinates of the secondary structures to execute this one
        # This will already cast it to absolute
        case = plugin_source.load_plugin('builder').case_apply(case, connectivity=True, overwrite=False)

    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0]
    folders.mkdir(parents=True, exist_ok=True)

    # Find steps: Each pair of secondary structure.
    it = case.connectivities_str[0].split('.')
    steps = [it[i:i + 2] for i in range(0, len(it) - 1)]
    loop_step = case.cast_absolute()['configuration.defaults.distance.loop_step']
    lengths = case.connectivity_len[0]
    print(lengths)
    start = 1

    for i, sse in enumerate(steps):
        # 1. Make folders
        wfolder = folders.joinpath('loop{:02d}'.format(i + 1))
        wfolder.mkdir(parents=True, exist_ok=True)
        outfile = wfolder.joinpath('loop_master.jump{:02d}.pdb'.format(i + 1))
        masfile = outfile.with_suffix('.master')

        # 2. Check hairpin
        sse1 = case.get_sse_by_id(sse[0])
        sse1_name = sse1['id']
        sse2 = case.get_sse_by_id(sse[1])
        sse2_name = sse2['id']
        is_hairpin = check_hairpin(sse1_name, sse2_name)

        if not masfile.is_file():
            # 3. Generate structures
            sse1, sse2 = make_structure(sse1, sse2, outfile)

            # 4. calculate expected loop length by loop_step
            Mdis, mdis = get_loop_length(sse1, sse2, loop_step, loop_range)

            # 5. Run master
            execute_master(outfile, pds_list, mdis, Mdis)

            # 6. Minimize master data (pick top_loopsx3 lines to read and minimize the files)
            minimize_master_file(masfile, top_loops, 3)

        # 7. Retrieve master data
        dfloop = process_master_data(masfile, sse1_name, sse2_name, abego, fragfiles, top_loops, is_hairpin and harpins_2)
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write(dfloop.to_string())

        # 8. Make Fragments


        # 8. Make files
        file3 = outfile.with_suffix('.3mers')
        file9 = outfile.with_suffix('.9mers')

        # 9. Attach files
        case.data['metadata']['loop_fragments'].append({'length': dfloop.iloc[0].values[0],
                                                        'abego': list(dfloop['abego'].values),
                                                        'fragfiles': [str(file3.resolve()), str(file9.resolve())]})

    return case


def make_fragment_files( dfloop: pd.DataFrame, fragfiles: pd.DataFrame ) -> List[Dict]:
    """
    """
    pass


def get_fragfiles():
    """
    """
    fragpath = Path(core.get_option('master', 'fragments'))
    if not fragpath.is_dir():
        raise IOError('MASTER fragments folder cannot be found.')
    return pd.DataFrame([(x.name[:4], x.name[5:6], x, y) for x, y in zip(sorted(fragpath.glob('*/*3mers.gz')),
                                                                         sorted(fragpath.glob('*/*9mers.gz')))],
                        columns=['pdb', 'chain', '3mers', '9mers'])


def get_abegos():
    """
    """
    abegos = core.get_option('loop_master', 'abego')
    abegos = Path(abegos)
    if not abegos.is_file():
        raise IOError('The ABEGO fasta file has to be provided')
    else:
        doopen = gzip.open if abegos.suffix == '.gz' else open
        abegodata = []
        with doopen(abegos, 'rt') as fd:
            for line1, line2 in itertools.zip_longest(*[fd] * 2):
                line2 = line2 if len(line2.strip()) != 0 else 'NON\n'
                line1 = line1.strip().lstrip('>').split('_')
                abegodata.append('{},{},{}'.format(line1[0], line1[1], line2))
        abegodata = pd.read_csv(StringIO(''.join(abegodata)), names=['pdb', 'chain', 'abego'], header=None)
        abegodata = abegodata[abegodata['abego'] != 'NON']

    return abegodata


def make_structure(sse1: dict, sse2: dict, outfile: Path) -> Tuple[PDBFrame, PDBFrame]:
    """
    """
    sse1 = PDB(pd.DataFrame(sse1['metadata']['atoms'],
                            columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                     'Cartn_x', 'Cartn_y', 'Cartn_z'])).renumber(1)

    sse2 = PDB(pd.DataFrame(sse2['metadata']['atoms'],
                            columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                     'Cartn_x', 'Cartn_y', 'Cartn_z'])).renumber(sse1.iloc[-1]['auth_seq_id'] + 5)
    structure = pd.concat([sse1, sse2])
    structure['id'] = list(range(1, structure.shape[0] + 1))

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> generating structure {}\n'.format(outfile.resolve()))
    structure.write(output_file=str(outfile), format='pdb', clean=True,
                    force=TBcore.get_option('system', 'overwrite'))

    return sse1, sse2


def get_loop_length(sse1: PDB, sse2: PDB, loop_step: int, loop_range: int) -> Tuple[int, int]:
    """
    """
    res1 = ChainFrame(PDB(sse1)).last_compound
    res2 = ChainFrame(PDB(sse2)).first_compound
    distance = SBIgeo.point_distance(res1[res1['label_atom_id'] == 'N'].coordinates,
                                     res2[res2['label_atom_id'] == 'N'].coordinates)
    distance = math.ceil(distance / loop_step)
    distance = [x for x in range(distance - loop_range - 1, distance + loop_range + 1) if x > 0]
    return max(distance), min(distance)


def execute_master(outfile: Path, pds_list: Path, mdis: int, Mdis: int):
    """
    """
    createPDS = core.get_option('master', 'create')
    createbash = '{0} −−type query −−pdb {1} --pds {2}'
    master = core.get_option('master', 'master')
    masterbash = '{0} −−query {1} −−targetList {2} −−rmsdCut 5 −−matchOut {3} --gapLen {4}-{5}'

    createcmd = shlex.split(createbash.format(createPDS, outfile, outfile.with_suffix('.pds')))
    mastercmd = shlex.split(masterbash.format(master, outfile.with_suffix('.pds'),
                                              pds_list, outfile.with_suffix('.master'), mdis, Mdis))
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> Execute: {}\n'.format(' '.join(createcmd)))
    run(createcmd)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> Execute: {}\n'.format(' '.join(mastercmd)))
    run(mastercmd)


def minimize_master_file( masfile: Path, top_loops: int, multiplier: int ):
    """
    """
    try:
        with open(masfile) as fd:
            head = [next(fd) for x in range(top_loops * multiplier)]
        with open(masfile, 'w') as fd:
            fd.write(''.join(head))
    except StopIteration:
        pass


def check_hairpin( name1: str, name2: str) -> bool:
    """
    """
    if name1[0] != name2[0]:
        return False
    if name1[-1] != 'E':
        return False
    if int(name1[1]) == int(name2[1]) + 1:
        return True
    if int(name1[1]) == int(name2[1]) - 1:
        return True
    return False


def process_master_data( masfile: Path,
                         name1: str,
                         name2: str,
                         abego: pd.DataFrame,
                         fragfiles: pd.DataFrame,
                         top_loops: int,
                         hairpin: bool ) -> pd.DataFrame:
    """
    """
    def cutter(row):
        match = row['match']
        # MASTER starts match count at 0!
        loop = row['abego'][match[0][1] + 1: match[1][0]]
        return row['abego'][match[0][0]: match[1][1] + 1], loop, len(loop)

    dfloop = plugin_source.load_plugin('imaster').parse_master_file(masfile)
    dfloop = dfloop.merge(abego, on=['pdb', 'chain']).merge(fragfiles, on=['pdb', 'chain']).dropna()
    dfloop[['abego', 'loop', 'loop_length']] = dfloop.apply(cutter, axis=1, result_type='expand')
    dfloop = dfloop.iloc[:top_loops]
    dfloop['length_count'] = dfloop.loop_length.map(dfloop.loop_length.value_counts())
    finaldf = dfloop.drop_duplicates(['loop'])

    pick = 0
    if hairpin and 2 in finaldf['loop_length']:
        pick = 2
    else:
        pick = finaldf[finaldf['length_count'] == finaldf['length_count'].max()]['loop_length'].min()
    finaldf = finaldf[finaldf['loop_length'] == pick]

    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    sns.countplot(x='loop_length', data=dfloop, ax=ax, palette="PRGn")
    cnt = dfloop[dfloop['loop_length'] == pick]['length_count'].values[0]
    pick = [int(x.get_text()) for x in ax.get_xticklabels()].index(pick)
    ax.plot(pick, cnt, marker=11, color='black')
    ax.set_title('loop {} <-> {}'.format(name1, name2))
    imagename = masfile.with_suffix(TBcore.get_option('system', 'image'))
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> loop image summary at {}\n'.format(imagename))
    plt.savefig(imagename, dpi=300)

    return finaldf
