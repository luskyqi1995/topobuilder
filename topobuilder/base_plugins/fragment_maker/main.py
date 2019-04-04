# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Optional, Tuple, List, Union
import sys

# External Libraries
import pandas as pd

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from rstoolbox.io import parse_rosetta_fragments, write_rosetta_fragments


__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           protocol: str,
           **kwargs ) -> List[Case]:
    """Generates final fragment files for the full structure.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: FRAGMENT_MAKER ---\n')

    # Execute for each case
    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('fragments', {})
        cases[i] = case_apply(case, protocol)
        cases[i] = cases[i].set_protocol_done(prtid)

    return cases


def case_apply( case: Case,
                protocol: str,
                script: Optional[Union[Path, str]] = None,
                ) -> Case:
    """
    """
    case = Case(case)
    data = {'protocol': protocol, 'files': []}

    # Fragments can only be made for a full, reoriented Case.
    if case.connectivity_count > 1:
        raise ValueError('FunFolDes can only be applied to one connectivity.')

    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0].joinpath('fragment_maker')
    folders.mkdir(parents=True, exist_ok=True)
    checkpoint = folders.joinpath('checkpoint.json')

    # Check if checkpoint exists, retrieve and skip
    reload = TButil.checkpoint_in(checkpoint)
    if reload is not None and reload['protocol'] == protocol:
        case.data['metadata']['fragments'] = reload
        return case

    # Switch depending on the fragment_protocol
    if protocol == 'loop_master':
        data['files'] = loop_master_protocol(case, folders)

    # Store data
    case.data['metadata']['fragments'] = data

    # Checkpoint save
    TButil.checkpoint_out(checkpoint, data)

    return case


def loop_master_protocol( case: Case, folders: Path ) -> Tuple[str, str]:
    """
    """
    lf = case['metadata.loop_fragments']
    if lf is None:
        raise TButil.PluginOrderError('Data that should be loaded through loop_master is not found.')

    for i, loop in enumerate(lf):
        if i == 0:
            ff3 = parse_rosetta_fragments(loop['fragfiles'][0])
            ff9 = parse_rosetta_fragments(loop['fragfiles'][1])
            df3 = [pd.read_csv(str(loop['fragfiles'][0]) + '.csv'), ]
            df9 = [pd.read_csv(str(loop['fragfiles'][1]) + '.csv'), ]
        else:
            df3.append(pd.read_csv(str(loop['fragfiles'][0]) + '.csv'))
            df9.append(pd.read_csv(str(loop['fragfiles'][1]) + '.csv'))
            ff3 = ff3.add_fragments(parse_rosetta_fragments(loop['fragfiles'][0]), ini=int(loop['edges']['ini']), how='append')
            ff9 = ff9.add_fragments(parse_rosetta_fragments(loop['fragfiles'][1]), ini=int(loop['edges']['ini']), how='append')

    TButil.plot_fragment_templates(pd.concat(df3), pd.concat(df9), folders.joinpath('template_fragment_profile'))

    small_file = write_rosetta_fragments(ff3.top_limit(lf[-1]['edges']['end']), prefix=folders.joinpath('small'), strict=True)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing small fragment file: {}\n'.format(small_file))
    large_file = write_rosetta_fragments(ff9.top_limit(lf[-1]['edges']['end']), prefix=folders.joinpath('large'), strict=True)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing large fragment files: {}\n'.format(large_file))

    return small_file, large_file
