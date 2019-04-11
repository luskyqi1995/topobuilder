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
from ast import literal_eval
from subprocess import run
import gzip
import itertools


# External Libraries
import pandas as pd
from pandas.compat import StringIO
from SBI.structure import PDB, PDBFrame, ChainFrame
import SBI.structure.geometry as SBIgeo
from rstoolbox.io import parse_rosetta_fragments, write_rosetta_fragments

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder import plugin_source
from .core import core


__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           loop_range: int = 3,
           top_loops: int = 20,
           harpins_2: bool = True,
           rmsd_cut: float = 5.0,
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
        cases[i].data.setdefault('metadata', {}).setdefault('loop_lengths', [])
        cases[i] = case_apply(case, database, loop_range, top_loops, rmsd_cut, abegodata, harpins_2, fragfiles)
        cases[i] = cases[i].set_protocol_done(prtid)

    if tempdb:
        os.unlink(f.name)

    return cases


def case_apply( case: Case,
                pds_list: Path,
                loop_range: int,
                top_loops: int,
                rmsd_cut: float,
                abego: pd.DataFrame,
                harpins_2: bool,
                fragfiles: pd.DataFrame ) -> str:
    """
    """
    # Loop MASTER is only applied to a Case with one single connectivity and already reoriented
    if case.connectivity_count > 1:
        raise ValueError('Loop MASTER can only be applied to one connectivity.')

    # We will need the coordinates of the secondary structures to execute this one
    # This will already cast it to absolute
    with TBcore.on_option_value('system', 'overwrite', False):
        case = plugin_source.load_plugin('builder').case_apply(case, connectivity=True)

    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0].joinpath('loop_master')
    folders.mkdir(parents=True, exist_ok=True)

    # Find steps: Each pair of secondary structure.
    it = case.connectivities_str[0].split('.')
    steps = [it[i:i + 2] for i in range(0, len(it) - 1)]
    loop_step = case.cast_absolute()['configuration.defaults.distance.loop_step']
    lengths = case.connectivity_len[0]
    start = 1

    for i, sse in enumerate(steps):
        # 1. Make folders and files
        wfolder = folders.joinpath('loop{:02d}'.format(i + 1))
        wfolder.mkdir(parents=True, exist_ok=True)
        outfile = wfolder.joinpath('loop_master.jump{:02d}.pdb'.format(i + 1))
        masfile = outfile.with_suffix('.master')
        checkpoint = wfolder.joinpath('checkpoint.json')

        # 2. Check if checkpoint exists, retrieve and skip
        reload = TButil.checkpoint_in(checkpoint)
        if reload is not None:
            case.data['metadata']['loop_fragments'].append(reload)
            case.data['metadata']['loop_lengths'].append(int(reload['edges']['loop']))
            start += (int(reload['edges']['sse1']) + int(reload['edges']['loop']))
            continue

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
            execute_master_fixedgap(outfile, pds_list, mdis, Mdis, rmsd_cut)

            # 6. Minimize master data (pick top_loopsx3 lines to read and minimize the files)
            match_count = minimize_master_file(masfile, top_loops, 3)

        # 7. Retrieve master data
        dfloop = process_master_data(masfile, sse1_name, sse2_name, abego, fragfiles, top_loops, is_hairpin and harpins_2)
        sse1l, loopl, sse2l = lengths[i], int(dfloop['loop_length'].values[0]), lengths[i + 1]
        total_len = sse1l + loopl + sse2l
        end_edge = total_len + start - 1
        edges = {'ini': int(start), 'end': int(end_edge), 'sse1': int(sse1l), 'loop': int(loopl), 'sse2': int(sse2l)}
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('\nINI: {}; END: {}; SSE1: {}; LOOP: {}; SSE2: {}\n\n'.format(start, end_edge, sse1l, loopl, sse2l))
            sys.stdout.write(dfloop.to_string() + '\n')

        # 8. Make Fragments
        loop_data = make_fragment_files(dfloop, edges, masfile)
        loop_data['match_count'] += match_count
        case.data['metadata']['loop_fragments'].append(loop_data)
        case.data['metadata']['loop_lengths'].append(int(loopl))

        start += (sse1l + loopl)

        # Checkpoint save
        TButil.checkpoint_out(checkpoint, loop_data)

    return case


def make_fragment_files( dfloop: pd.DataFrame, edges: Dict, masfile: Path ) -> Dict:
    """
    """
    data = {'loop_length': int(dfloop.iloc[0]['loop_length']), 'abego': list(dfloop['loop'].values),
            'edges': edges, 'fragfiles': [], 'match_count': 0}

    dfs3 = []
    dfs9 = []
    sample = math.ceil(200 / dfloop.shape[0])
    for i, row in dfloop.iterrows():
        # Remember: MASTER match starts with 0!
        dfs3.append((parse_rosetta_fragments(str(row['3mers']), source='{}_{}'.format(row['pdb'], row['chain']))
                     .slice_region(row['match'][0][0] + 1, row['match'][1][1] + 1).sample_top_neighbors(sample)
                     .renumber(edges['ini']).top_limit(edges['end'])))
        dfs9.append((parse_rosetta_fragments(str(row['9mers']), source='{}_{}'.format(row['pdb'], row['chain']))
                     .slice_region(row['match'][0][0] + 1, row['match'][1][1] + 1).sample_top_neighbors(sample)
                     .renumber(edges['ini']).top_limit(edges['end'])))

    # Merge Fragments
    dfs3all = dfs3[0]
    dfs9all = dfs9[0]
    for i in range(1, len(dfs3)):
        dfs3all = dfs3all.add_fragments(dfs3[i], ini=edges['ini'], how='append')
        dfs9all = dfs9all.add_fragments(dfs9[i], ini=edges['ini'], how='append')
    dfs3all = dfs3all.sample_top_neighbors(200)
    dfs9all = dfs9all.sample_top_neighbors(200)

    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('Writing 3mers fragfile\n')
    data['fragfiles'].append(write_rosetta_fragments(dfs3all, prefix=str(masfile.with_suffix('')), strict=True))
    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('3mers fragfile: {}\n'.format(data['fragfiles'][-1]))
        sys.stdout.write('Writing 9mers fragfile\n')
    data['fragfiles'].append(write_rosetta_fragments(dfs9all, prefix=str(masfile.with_suffix('')), strict=True))
    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('9mers fragfile: {}\n'.format(data['fragfiles'][-1]))

    dfs3all.drop(columns=['pdb', 'frame', 'neighbors', 'neighbor',
                          'aa', 'sse', 'phi', 'psi', 'omega']).to_csv(data['fragfiles'][0] + '.csv', index=False)
    dfs9all.drop(columns=['pdb', 'frame', 'neighbors', 'neighbor',
                          'aa', 'sse', 'phi', 'psi', 'omega']).to_csv(data['fragfiles'][1] + '.csv', index=False)
    imageprefix = masfile.with_suffix('.fragprofile')
    TButil.plot_fragment_templates(dfs3all, dfs9all, imageprefix)

    return data


def get_fragfiles():
    """
    """
    fragpath = Path(core.get_option('master', 'fragments'))
    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('Listing available fragment files at: {}\n'.format(fragpath.name))
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

    if TBcore.get_option('system', 'debug'):
        sys.stdout.write('Loading ABEGO data from: {}\n'.format(abegos.name))
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


def execute_master_fixedgap(outfile: Path, pds_list: Path, mdis: int, Mdis: int, rmsd_cut: float):
    """
    """
    createPDS = core.get_option('master', 'create')
    createbash = '{0} --type query --pdb {1} --pds {2}'
    master = core.get_option('master', 'master')
    masterbash = '{0} --query {1} --targetList {2} --rmsdCut {6} --matchOut {3} --gapLen {4}-{5}'

    createcmd = shlex.split(createbash.format(createPDS, outfile, outfile.with_suffix('.pds')))
    mastercmd = shlex.split(masterbash.format(master, outfile.with_suffix('.pds'),
                                              pds_list, outfile.with_suffix('.master'), mdis, Mdis, rmsd_cut))
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> Execute: {}\n'.format(' '.join(createcmd)))
    with open(os.devnull, 'w') as devnull:
        run(createcmd, stdout=devnull)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('-> Execute: {}\n'.format(' '.join(mastercmd)))
    with open(os.devnull, 'w') as devnull:
        run(mastercmd, stdout=devnull)


def minimize_master_file( masfile: Path, top_loops: int, multiplier: int ) -> int:
    """
    """
    try:
        with open(masfile) as fd:
            num_lines = sum(1 for line in fd if line.rstrip())
        with open(masfile) as fd:
            head = [next(fd) for x in range(top_loops * multiplier)]
        with open(masfile, 'w') as fd:
            fd.write(''.join(head))
    except StopIteration:
        pass
    return num_lines


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

    if masfile.with_suffix('.csv').is_file():
        df = pd.read_csv(masfile.with_suffix('.csv'))
        df['match'] = df['match'].apply(literal_eval)
        return df

    dfloop = TButil.parse_master_file(masfile)
    dfloop = dfloop.merge(abego, on=['pdb', 'chain']).merge(fragfiles, on=['pdb', 'chain']).dropna()
    dfloop[['abego', 'loop', 'loop_length']] = dfloop.apply(cutter, axis=1, result_type='expand')
    dfloop = dfloop.iloc[:top_loops]
    dfloop['length_count'] = dfloop.loop_length.map(dfloop.loop_length.value_counts())
    dfloop.drop(columns=['pds_path']).to_csv(masfile.with_suffix('.all.csv'), index=False)
    finaldf = dfloop.sort_values('rmsd').drop_duplicates(['loop'])

    pick = 0
    if hairpin and 2 in finaldf['loop_length'].values:
        pick = 2
    else:
        pick = finaldf[finaldf['length_count'] == finaldf['length_count'].max()]['loop_length'].min()
    finaldf = finaldf[finaldf['loop_length'] == pick]

    TButil.plot_loop_length_distribution(dfloop, pick, masfile.with_suffix(''), 'loop {} <-> {}'.format(name1, name2))

    df = finaldf.drop(columns=['pds_path'])
    df.to_csv(masfile.with_suffix('.csv'), index=False)
    return df
