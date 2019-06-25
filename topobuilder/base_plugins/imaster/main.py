# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, Dict, Set
from string import ascii_uppercase
from operator import itemgetter
from subprocess import run, DEVNULL
from pathlib import Path
from itertools import cycle
import shlex
import math
import os
import sys

# External Libraries
import numpy as np
import pandas as pd
import scipy as sc
import networkx as nx

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder import plugin_source
from .analysis import get_steps


__all__ = ['apply', 'case_apply']


def apply( cases: List[Case],
           prtid: int,
           rmsd: Optional[float] = 5.0,
           bin: Optional[str] = 'mid',
           **kwargs,
           ) -> List[Case]:
    """Use MASTER to correct secondary structure placement (smoothing).
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('corrections', [])
        cases[i] = case_apply(case, rmsd, bin)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                rmsd: Optional[float] = 5.0,
                bin: Optional[str] = 'mid',
                step: Optional[int] = None,
                corrections: Optional[Dict] = dict() ) -> Case:
    """
    """
    kase = Case(case)
    corrections = corrections if TBcore.get_option('system', 'jupyter') else {}
    kase.data.setdefault('metadata', {}).setdefault('imaster', {})
    kase.data.setdefault('metadata', {}).setdefault('corrections', [])
    kase.data.setdefault('metadata', {}).setdefault('bin', bin)

    # Interactive MASTER smoothing is only applied to a Case with one single connectivity and already reoriented
    if kase.connectivity_count > 1:
        raise ValueError('Interactive MASTER smoothing can only be applied to one connectivity.')
    if not kase['configuration.reoriented'] and kase.connectivity_count == 1:
        kase = kase.cast_absolute().apply_topologies()[0]

    # Generate the folder tree for a single connectivity.
    wfolder = kase.connectivities_paths[0].joinpath('imaster')
    wfolder.mkdir(parents=True, exist_ok=True)
    current_case_file = kase.cast_absolute().write(wfolder.joinpath('current'))

    # Find steps: Tops we will submit 2-layer searches
    steps = get_steps([x[-1] == 'E' for x in kase.architecture_str.split('.')])
    steps = [steps[step], ] if step is not None and TBcore.get_option('system', 'jupyter') else steps

    # Work by layers
    done_l = set()
    for i, step in enumerate(steps):
        # Step working directory
        stepfolder = wfolder.joinpath('step{:02d}'.format(i + 1))
        stepfolder.mkdir(parents=True, exist_ok=True)
        query = stepfolder.joinpath('imaster.query{:02d}.pdb'.format(i + 1))
        checkpoint = stepfolder.joinpath('checkpoint.json')

        reload = TButil.checkpoint_in(checkpoint)
        if reload is not None:
            kase.data['metadata']['imaster'].setdefault('step{:02d}'.format(i + 1), reload)
            kase.data['metadata']['corrections'].append(reload['corrections'])
            corrections.update(reload['corrections'])
            done_l.update(reload['layers'])
            # CKase = CKase.apply_corrections(corrections)
            continue

        # Apply corrections from previous steps and rebuild
        CKase = Case(kase).apply_corrections(corrections)
        with TBcore.on_option_value('system', 'overwrite', True):
            CKase = plugin_source.load_plugin('builder').case_apply(CKase, connectivity=True)

        # Generate structure query and get layer displacements
        layers = set(itemgetter(*step)(ascii_uppercase))
        sses = [sse for sse in CKase.ordered_structures if sse['id'][0] in layers]
        structure, _ = TButil.build_pdb_object(sses, 3)
        TButil.plugin_filemaker('Writing structure {0}'.format(query))
        structure.write(output_file=str(query), format='pdb', clean=True, force=True)

        flip = cycle([CKase['configuration.flip_first'], not CKase['configuration.flip_first']])
        counts = np.asarray([sse['length'] for sse in CKase.ordered_structures])
        cends = np.cumsum(counts)
        cstrs = cends - cends + 1

        rules = list(zip([sse['id'] for sse in CKase.ordered_structures],
                         list(zip(cstrs, cends)),
                         list(next(flip) for _ in range(len(CKase.ordered_structures)))))
        print(rules)

        # MASTER search
        createpds = TButil.createPDS(query)
        TButil.plugin_bash(createpds)
        run(createpds, stdout=DEVNULL)
        masters = TButil.master_best_each(query.with_suffix('.pds'), stepfolder.joinpath('_master'), rmsd)
        data = submit_searches(masters, stepfolder, current_case_file, '.'.join([x['id'] for x in sses]))
        data = calc_corrections(data, kase, set(data['layers']), done_l, bin)

        kase.data['metadata']['imaster'].setdefault('step{:02d}'.format(i + 1), data)
        TButil.checkpoint_out(checkpoint, data)
        kase.data['metadata']['corrections'].append(data['corrections'])
        done_l.update(data['layers'])
        corrections.update(data['corrections'])

    return kase


def submit_searches( cmd: List[str], wdir: Path, current_case_file: Path, current_sse: str ) -> Dict:
    """
    """
    unimaster = wdir.joinpath('match.master')
    imaster = Path(__file__).parent.joinpath('imaster.py')
    unidata = wdir.joinpath('geometry.csv')
    if unimaster.is_file() and unidata.is_file():
        return {'matches': unimaster, 'stats': unidata, 'corrections': None,
                'layers': list(set([x[0] for x in current_sse.split('.')]))}
    if not TBcore.get_option('slurm', 'use'):
        no_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata.with_suffix(''))
    else:
        with_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata)

    return {'matches': unimaster, 'stats': unidata, 'corrections': None,
            'layers': list(set([x[0] for x in current_sse.split('.')]))}


def calc_corrections( data: Dict, case: Case, qlayers: Set, dlayers: Set, bin: Optional[str] = 'mid' ) -> Dict:
    """
    """
    tocorrect = qlayers.difference(dlayers)
    toreference = qlayers.difference(tocorrect)
    if len(tocorrect) > 1:
        raise ValueError('Layers are corrected one by one.')
    tocorrect = list(tocorrect)[0]

    if len(toreference) == 1:
        toreference = list(toreference)[0]
    elif len(toreference) == 0:
        toreference = None
    else:
        for x in toreference:
            found = False
            if abs(ascii_uppercase.find(x) - ascii_uppercase.find(tocorrect)) == 1:
                toreference = x
                found = True
                break
            if not found:
                toreference = None

    # Load Data, bin, show and addapt if no matches for the given bin.
    bins = ["close", "mid", "far", "extreme"]
    df = pd.read_csv(data['stats'])
    df = df.assign(bin=pd.cut(df['rmsd'], bins=[0, 2, 2.5, 3, 5], labels=bins))
    _, _, isBin = TButil.plot_match_bin(df, Path(data['stats']).parent.joinpath('match_count'),
                                        len(TButil.pds_database()[1]), ['pdb', 'chain'])
    while not isBin[bin]:
        binInt = bins.index(bin)
        if binInt >= len(bins):
            if TBcore.get_option('system', 'verbose'):
                sys.stdout.write('No matches found to use for correction.')
            data['corrections'] = {}
            return data
        bin = bins[binInt]
        data['bin'] = bin

    if toreference is None:
        if case.get_type_for_layer(tocorrect) == 'E':
            data['corrections'] = first_layer_beta_correction(df, bin, Path(data['stats']).parent)
    elif case.get_type_for_layer(toreference) == 'E':
        if case.get_type_for_layer(tocorrect) == 'H':
            data['corrections'] = alpha_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Found corrections {}\n'.format(data['corrections']) )
    return data


def alpha_on_beta_correction(df: pd.DataFrame, bin: str, wdir: Path, qlayer: str, rlayer: str, case: Case ) -> Dict:
    """
    """
    # Report data
    stats = make_mode_stats(df, wdir).reset_index()
    clms = ['measure', 'layer', 'sse', bin]
    stats = stats[(stats['layer'] == rlayer)][clms]
    for layer in sorted(df.layer.unique()):
        ofile = 'geometric_distributions_layer{}'.format(layer)
        TButil.plot_geometric_distributions(df[df['layer'] == layer], Path(wdir).joinpath(ofile))

    data = {}
    for sse in [x for x in stats.sse.unique() if x.startswith(qlayer)]:
        ddf = stats[(stats['sse'] == sse)]
        data.setdefault(sse, {}).setdefault('tilt', {'x': ddf[ddf['measure'] == 'angles_layer'][bin].values[0],
                                                     'z': ddf[ddf['measure'] == 'angles_side'][bin].values[0]})
        pc = ddf[ddf['measure'] == 'points_layer'][bin].values[0] - case['configuration.defaults.distance.ab']
        if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
            pc = pc * -1

        data.setdefault(sse, {}).setdefault('coordinates', {'z': pc})
    return data


def first_layer_beta_correction( df: pd.DataFrame, bin: str, wdir: Path ) -> Dict:
    """
    """
    # Report data
    make_mode_stats(df, wdir)
    for layer in sorted(df.layer.unique()):
        ofile = 'geometric_distributions_layer{}'.format(layer)
        TButil.plot_geometric_distributions(df[df['layer'] == layer], Path(wdir).joinpath(ofile))

    # Make network
    def reshape(df):
        def reshape(g):
            g = g.drop(columns=['pdb', 'chain']).T
            g = g.rename(columns=g.loc['sse'])
            g = g.reindex(g.index.drop('sse'))
            g = g.assign(bin=g.loc['bin'].values[0])
            g = g.assign(A0E=0)
            g = g.drop(index='bin')
            return g

        df = df.copy()[['pdb', 'chain', 'angles_layer', 'sse', 'bin']]
        bins = list(range(-100, 105, 5))
        labels = list(np.arange(-97.5, 100, 5))
        df = df.assign(anglebin=pd.cut(df['angles_layer'], bins=bins, labels=labels))
        df = df.drop(columns=['angles_layer']).groupby(['pdb', 'chain']).apply(reshape)
        df.index = list(range(df.shape[0]))
        return df.drop(columns=['bin'])
    ddf = reshape(df[df['bin'] == bin])
    sses = sorted(list(ddf.columns))
    netwk = []
    for i in range(0, len(sses) - 1):
        idx = [sses[i], sses[i + 1]]
        tmp = pd.DataFrame(ddf.groupby(idx).size()).rename(columns={0: 'count'}).reset_index()
        netwk.append(tmp.reset_index())
        netwk[-1][idx[0]] = netwk[-1][idx[0]].apply(lambda v: "_".join([idx[0], str(v)]))
        netwk[-1][idx[1]] = netwk[-1][idx[1]].apply(lambda v: "_".join([idx[1], str(v)]))
        netwk[-1] = nx.from_pandas_edgelist(netwk[-1], idx[0], idx[1], ['count'], create_using=nx.DiGraph)

    posk = {}
    try:
        network = nx.compose_all(netwk)
        for n in network.nodes:
            dd = n.split('_')
            posk.setdefault(n, (int(dd[0][1]), float(dd[1])))
    except Exception:
        network = nx.Graph()

    # Image summary
    TButil.plot_angle_network(network, posk, sses, Path(wdir).joinpath('network_{}'.format(bin)))

    # Return corrections
    data = {}
    for x in nx.dag_longest_path(network, 'count', default_weight=0):
        x = x.split('_')
        if x[0] == 'A0E':
            continue
        else:
            data.setdefault(x[0], {}).setdefault('tilt', {}).setdefault('x', float(x[1]))
    return data


def make_mode_stats( df: pd.DataFrame, wdir: Path ) -> pd.DataFrame:
    """
    """
    data = {'close': [], 'mid': [], 'far': [], 'extreme': []}
    topi, midi, boti = [], [], []

    for angle in [x for x in df.columns if x.startswith('angles_') or x.startswith('points_')]:
        for part in ['close', 'mid', 'far', 'extreme']:
            dfp = df[df['bin'] == part]
            for lay in df.layer.unique():
                for sse in df.sse.unique():
                    dfs = dfp[(dfp['sse'] == sse) & (dfp['layer'] == lay)]
                    try:
                        kde = sc.stats.gaussian_kde(dfs[angle])
                        x = np.linspace(dfs[angle].min(), dfs[angle].max(), 200)
                        kde = kde(x)
                        mode = x[np.argsort(kde)[-1]]
                    except ValueError:
                        mode = 0
                    data[part].append(mode)
                    if part == 'close':
                        topi.append(angle)
                        midi.append(lay)
                        boti.append(sse)
    stats = pd.DataFrame(data, index=pd.MultiIndex.from_tuples(list(zip(*[topi, midi, boti])),
                                                               names=['measure', 'layer', 'sse']))
    stats.reset_index().to_csv(wdir.joinpath('mode_stats.csv'), index=False)
    TButil.plugin_filemaker('Mode stats stored at {}'.format(wdir.joinpath('mode_stats.csv')))
    return stats


def with_slurm( cmd: List[str],
                current_case_file: Path,
                current_sse: str,
                unimaster: Path,
                imaster: Path,
                unidata: Path ):
    """
    """
    # Make bashfile
    bashcont = []
    createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
    parts = math.ceil(len(cmd) / TBcore.get_option('slurm', 'array'))

    wwd = unimaster.parent.parent
    cwd = Path().cwd()
    os.chdir(str(wwd))

    for i, com in enumerate(cmd):
        cmd[i][2] = str(Path(com[2]).relative_to(wwd))
        cmd[i][-1] = str(Path(com[-1]).relative_to(wwd))

    for ii, cp in enumerate(cmd):
        cmd[ii][-1] = cp[-1] + '_${SLURM_ARRAY_TASK_ID}'
    for j, i in enumerate(range(0, len(cmd), parts)):
        sumfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.master').relative_to(wwd)
        datfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.geo').relative_to(wwd)
        bashcont.append('if (( ${{SLURM_ARRAY_TASK_ID}} == {} )); then'.format(j + 1))
        bashcont.extend([' '.join(x) for x in cmd[i:i + parts]])
        bashcont.append('cat {0} > {1}'.format(Path(cmd[-1][-1]).parent.joinpath('*_${SLURM_ARRAY_TASK_ID}'), sumfile))
        bashcont.append(createbash.format(imaster, current_case_file.relative_to(wwd),
                                          sumfile, current_sse, datfile))
        bashcont.append('fi')
    with unimaster.parent.joinpath('submit.sh').relative_to(wwd).open('w') as fd:
        fd.write(TButil.slurm_header())
        fd.write(TButil.slurm_pyenv())
        fd.write('\n'.join(bashcont))

    TButil.submit_slurm(unimaster.parent.joinpath('submit.sh').relative_to(wwd))
    TButil.plugin_filemaker('Creating geometric coordinate file {}'.format(unidata))
    allCSV = [str(x) for x in unimaster.parent.relative_to(wwd).glob('_*.geo.csv')]
    pd.concat([pd.read_csv(x) for x in allCSV]).to_csv(unidata.relative_to(wwd), index=False)
    TButil.plugin_filemaker('Creating MASTER search file {}'.format(unimaster))
    with unimaster.relative_to(wwd).open('w') as fd:
        for x in unimaster.parent.glob('_*.master'):
            with x.relative_to(wwd).open() as fi:
                fd.write(''.join(fi.readlines()))
    os.chdir(str(cwd))


def no_slurm( cmd: List[str],
              current_case_file: Path,
              current_sse: str,
              unimaster: Path,
              imaster: Path,
              unidata: Path ):
    """
    """
    # Search on MASTER
    result = []
    for com in cmd:
        TButil.plugin_bash(com)
        run(com, stdout=DEVNULL)
        outf = Path(com[-1])
        if outf.is_file():
            result.append(str(outf))
    result.insert(0, 'cat')
    TButil.plugin_filemaker('Unify matches at {0}'.format(unimaster))
    with unimaster.open('w') as fd:
        run(result, stdout=fd)
    result[0] = 'rm'
    run(result[:-2], stdout=DEVNULL)

    # Analyze
    createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
    cmd = shlex.split(createbash.format(imaster, current_case_file, unimaster, current_sse, unidata))
    TButil.plugin_bash(cmd)
    run(cmd, stdout=DEVNULL)
