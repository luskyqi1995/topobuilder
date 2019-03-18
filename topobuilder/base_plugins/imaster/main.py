# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, Dict
from string import ascii_uppercase
from operator import itemgetter
from subprocess import run

# External Libraries

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
           **kwargs,
           ) -> List[Case]:
    """Use MASTER to correct secondary structure placement (smoothing).
    """
    TButil.plugin_title(__file__, len(cases))

    for i, case in enumerate(cases):
        cases[i].data.setdefault('metadata', {}).setdefault('corrections', [])
        cases[i] = case_apply(case)
        cases[i] = cases[i].set_protocol_done(prtid)
    return cases


def case_apply( case: Case,
                rmsd: Optional[float] = 5.0,
                step: Optional[int] = None,
                corrections: Optional[Dict] = None ) -> Case:
    """
    """
    kase = Case(case)
    corrections = corrections if TBcore.get_option('system', 'jupyter') else None
    kase.data.setdefault('metadata', {}).setdefault('imaster', {}).setdefault('corrections', {})
    kase.data.setdefault('metadata', {}).setdefault('corrections', [])

    # Interactive MASTER smoothing is only applied to a Case with one single connectivity and already reoriented
    if kase.connectivity_count > 1:
        raise ValueError('Interactive MASTER smoothing can only be applied to one connectivity.')
    if not kase['configuration.reoriented'] and kase.connectivity_count == 1:
        kase = kase.cast_absolute().apply_topologies()[0]

    # Generate the folder tree for a single connectivity.
    wfolder = kase.connectivities_paths[0].joinpath('imaster')
    wfolder.mkdir(parents=True, exist_ok=True)

    # Find steps: Tops we will submit 2-layer searches
    steps = get_steps([x[-1] == 'E' for x in kase.architecture_str.split('.')])
    steps = [steps[step], ] if step is not None and TBcore.get_option('system', 'jupyter') else steps

    # Work by layers
    CKase = Case(kase)
    for i, step in enumerate(steps):
        # Step working directory
        stepfolder = wfolder.joinpath('step{:02d}'.format(i + 1))
        stepfolder.mkdir(parents=True, exist_ok=True)
        query = stepfolder.joinpath('imaster.query{:02d}.pdb'.format(i + 1))

        # Apply corrections from previous steps and rebuild
        CKase = CKase.apply_corrections(corrections)
        with TBcore.on_option_value('system', 'overwrite', True):
            CKase = plugin_source.load_plugin('builder').case_apply(CKase, connectivity=True)

        # Generate structure query
        layers = set(itemgetter(*step)(ascii_uppercase))
        sses = [sse for sse in CKase.ordered_structures if sse['id'][0] in layers]
        structure, _ = TButil.build_pdb_object(sses, 3)
        TButil.plugin_filemaker('Writing structure {0}'.format(query))
        structure.write(output_file=str(query), format='pdb', clean=True, force=True)

        # MASTER search
        run(TButil.createPDS(query))
        print(TButil.master_best_each(query.with_suffix('.pds'), stepfolder.joinpath('_master'), rmsd))

    kase.data['metadata']['corrections'].append(corrections)
    return kase
