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
from typing import Union, Optional, Tuple, List
import subprocess

# External Libraries

# This Library
import topobuilder.core as TBcore


__all__ = ['funfoldes_script']


def peptidestub_mover( mover_id: str,
                       rules: List[Tuple[int, int]]
                       ) -> Tuple[str, str]:
    """
    """
    body = textwrap.dedent("""\
        <PeptideStubMover name="{}" reset="false">
        {}
        </PeptideStubMover>""")
    insert = '<Insert resname="GLY" repeat="1" jump="false" anchor_rsd="{}" anchor_atom="C" connecting_atom="N" />'
    allins = []
    for x, y in rules:
        for i in range(y):
            allins.append(insert.format(x + i))
    return mover_id, body.format(mover_id, "\n".join(allins))


def secondarystructure_selector( sele_id: str, sse: str ) -> str:
    """
    """
    return sele_id, textwrap.dedent("""\
        <SecondaryStructure name="{}" overlap="0" minE="1"  minH="1"
            ss="HE" include_terminal_loops="false" use_dssp="false"
            pose_secstruct="{}" />""").format(sele_id, sse)


def funfoldes_script( sse_select: Tuple[str, str],
                      loop_mover: Tuple[str, str],
                      frag_files: Tuple[str, str]
                      ) -> str:
    """
    """
    sse_select = secondarystructure_selector(*sse_select)
    loop_mover = peptidestub_mover(*loop_mover)
    return textwrap.dedent("""\
    <ROSETTASCRIPTS>
        <SCOREFXNS>
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
            <Index name="piece" resnums="28-30" />
            {sse_selector}
        </RESIDUE_SELECTORS>
        <FILTERS>
            <RmsdFromResidueSelectorFilter name="rmsd" reference_selector="{sse_selector_id}"
                reference_name="sketchPose" query_selector="{sse_selector_id}" confidence="0." />
        </FILTERS>
        <MOVERS>
            {peptidestub_mover}
            <SavePoseMover name="inSketch" reference_name="sketchPose" restore_pose="0" />
            <StructFragmentMover name="makeFrags" prefix="frags"
                small_frag_file="{small_frag_file}" large_frag_file="{large_frag_file}"
            />
            <AddConstraints name="foldingCST" >
                <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="{sse_selector_id}" >
                    <Inner sd="1.2" weight="1." ca_only="1"
                        use_harmonic="true" unweighted="false" min_seq_sep="4" />
                    <Outer sd="2" weight="2." ca_only="1"
                        use_harmonic="true" unweighted="false"  max_distance="40" />
                </SegmentedAtomPairConstraintGenerator>
                <AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1" />
            </AddConstraints>
            <NubInitioMover name="FFL" fragments_id="frags" template_motif_selector="piece"
                        rmsd_threshold="10" residue_type="V" >
                <Nub reference_name="sketchPose" residue_selector="piece" >
                    <Segment order="1" n_term_flex="2" c_term_flex="1"/>
                </Nub>
            </NubInitioMover>
            <WriteSSEMover name="structure" dssp="true" />
        </MOVERS>
        <PROTOCOLS>
            <Add mover="{peptidestub_mover_id}" />
            <Add mover="inSketch" />
            <Add mover="makeFrags" />
            <Add mover="foldingCST" />
            <Add mover="FFL" />
            <Add mover="structure" />
            <Add filter="rmsd" />
        </PROTOCOLS>
        <OUTPUT />
    </ROSETTASCRIPTS>""").format(sse_selector_id=sse_select[0], sse_selector=sse_select[1],
                                 peptidestub_mover_id=loop_mover[0], peptidestub_mover=loop_mover[1],
                                 small_frag_file=frag_files[0], large_frag_file=frag_files[1])
