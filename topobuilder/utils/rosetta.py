# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import textwrap
import math
from typing import Dict, Union

# External Libraries
from jinja2 import Template
from bs4 import BeautifulSoup

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore


__all__ = ['rosettascript', 'funfoldes', 'constraint_design']


class ScriptPieces( dict ):
    def __add__( self, other ):
        data = ScriptPieces()
        for k in ['scorefxns', 'residueselectors', 'packerpalette', 'taskoperations', 'movemapfactory',
                  'simplemetrics', 'filters', 'movers', 'protocols', 'output']:
            if k in self or k in other:
                data.setdefault(k, [x for x in self.get(k, ['', ]) + other.get(k, ['', ]) if len(x) > 0])
        return data


def constraint_minimization(  case: Case, natbias: float ) -> ScriptPieces:
    """
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstmin" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [SELECTOR_SecondaryStructure('sse_cstmin', case), ]

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstmin" reference_selector="sse_cstmin"
            reference_name="eminPose_cstmin" query_selector="sse_cstmin" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstmin" reference_name="eminPose_cstmin" restore_pose="0" />
    <AddConstraints name="cst_cstmin" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstmin" residue_selector="sse_cstmin" >
            <Outer sd="2.0" weight="1." ca_only="1"
             use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <AutomaticSheetConstraintGenerator name="cst_sheet_cstmin" sd="2.0" distance="6.1" />
    </AddConstraints>
    <MinMover name="fast_cstmin" scorefxn="sfxn_cstmin" chi="1" bb="1" />
    """), MOVER_SetSecStructEnergies( 'ssse_cstmin', 'sfxn_cstmin', natbias, case )]

    protocols = [textwrap.dedent("""\
    <Add mover="spose_cstmin" />
    <Add mover="cst_cstmin" />
    <Add mover="ssse_cstmin" />
    <Add mover="fast_cstmin" />
    <Add filter="rmsd_cstmin" />
    """), ]

    bf = PROTOCOL_BasicFilters(case, '_cstmin')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def constraint_design( case: Case, natbias: float, layer_design: bool = True ) -> ScriptPieces:
    """
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstdes" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [SELECTOR_SecondaryStructure('sse_cstdes', case), ]

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstdes" reference_selector="sse_cstdes"
            reference_name="eminPose_cstdes" query_selector="sse_cstdes" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstdes" reference_name="eminPose_cstdes" restore_pose="0" />
    <AddConstraints name="cst_cstdes" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstdes" residue_selector="sse_cstdes" >
            <Outer sd="2.0" weight="0.5" ca_only="1" use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <AutomaticSheetConstraintGenerator name="cst_sheet_cstdes" sd="2.0" distance="6.1" />
    </AddConstraints>"""), MOVER_SetSecStructEnergies( 'ssse_cstdes', 'sfxn_cstdes', natbias, case ),
              textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes" relaxscript="MonomerDesign2019" taskoperations="layer_design"/>
    """) if layer_design else textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes" relaxscript="MonomerDesign2019"/>""")]

    protocols = [textwrap.dedent("""\
    <Add mover="spose_cstdes" />
    <Add mover="cst_cstdes" />
    <Add mover="ssse_cstdes" />
    <Add mover="design_cstdes" />
    <Add filter="rmsd_cstdes" />
    """), ]

    ld = PROTOCOL_LayerDesign(case) if layer_design else ScriptPieces()
    bf = PROTOCOL_BasicFilters(case, '_cstdes')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + ld + bf


def funfoldes( case: Case ) -> str:
    """
    """
    mid = math.floor(len(case.secondary_structure) / 2)
    residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}-{}" />""").format(mid - 1, mid + 1),
                        SELECTOR_SecondaryStructure('sse_ffd', case)]

    filters = [textwrap.dedent("""\
        <RmsdFromResidueSelectorFilter name="rmsd_ffd" reference_selector="sse_ffd"
            reference_name="sketchPose_ffd" query_selector="sse_ffd" confidence="0." />""")]

    movers = [MOVER_PeptideStubMover('add_loops_ffd', case), textwrap.dedent("""\
        <SavePoseMover name="save_ffd" reference_name="sketchPose_ffd" restore_pose="0" />
        <StructFragmentMover name="makeFrags_ffd" prefix="frags" small_frag_file="{}" large_frag_file="{}" />
        <AddConstraints name="foldingCST_ffd" >
            <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="sse_ffd" >
                <Inner sd="1.2" weight="1." ca_only="1"
                    use_harmonic="true" unweighted="false" min_seq_sep="4" />
                <Outer sd="2" weight="2." ca_only="1"
                    use_harmonic="true" unweighted="false"  max_distance="40" />
            </SegmentedAtomPairConstraintGenerator>
            <AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1" />
        </AddConstraints>
        <NubInitioMover name="FFL_ffd" fragments_id="frags" template_motif_selector="piece_ffd" rmsd_threshold="10" >
        """).format(*case['metadata.fragments.files']), textwrap.dedent("""\
            <Nub reference_name="sketchPose_ffd" residue_selector="piece_ffd" >
                <Segment order="1" n_term_flex="2" c_term_flex="1" editable="1,2,3"/></Nub>
            """), textwrap.dedent("""</NubInitioMover>""")]

    protocols = [textwrap.dedent("""\
        <Add mover="add_loops_ffd" />
        <Add mover="save_ffd" />
        <Add mover="makeFrags_ffd" />
        <Add mover="foldingCST_ffd" />
        <Add mover="FFL_ffd" />
        <Add filter="rmsd_ffd" />""")]

    with TBcore.on_option_value('psipred', 'script', None):
        bf = PROTOCOL_BasicFilters(case, '_ffd')
    return ScriptPieces({'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def rosettascript( data: Dict ) -> str:
    """
    """
    if 'output' in data and isinstance(data['output'], list):
        data['output'] = data['output'][0]
    content = BeautifulSoup(Template(textwrap.dedent("""\
    <ROSETTASCRIPTS>
        {% if scorefxns %}<SCOREFXNS>{% for item in scorefxns %}{{ item }}{% endfor %}</SCOREFXNS>{% endif %}
        {% if residueselectors %}
            <RESIDUE_SELECTORS>{% for item in residueselectors %}{{ item }}{% endfor %}</RESIDUE_SELECTORS>
        {% endif %}
        {% if packerpalettes %}<PACKER_PALETTES>{% for item in packerpalettes %}{{item}}{% endfor %}</PACKER_PALETTES>{% endif %}
        {% if taskoperations %}<TASKOPERATIONS>{% for item in taskoperations %}{{ item }}{% endfor %}</TASKOPERATIONS>{% endif %}
        {% if movemapfactory %}
            <MOVE_MAP_FACTORIES>{% for item in movemapfactory %}{{ item }}{% endfor %}</MOVE_MAP_FACTORIES>
        {% endif %}
        {% if simplemetrics %}<SIMPLE_METRICS>{% for item in simplemetrics %}{{ item }}{% endfor %}</SIMPLE_METRICS>{% endif %}
        {% if filters %}<FILTERS>{% for item in filters %}{{ item }}{% endfor %}</FILTERS>{% endif %}
        {% if movers %}<MOVERS>{% for item in movers %}{{ item }}{% endfor %}</MOVERS>{% endif %}
        {% if protocols %}<PROTOCOLS>{% for item in protocols %}{{ item }}{% endfor %}</PROTOCOLS>{% endif %}
        {% if output %}<OUTPUT scorefxn="{{ output }}"/>{% endif %}
    </ROSETTASCRIPTS>
    """)).render(data), 'xml').contents
    return '\n'.join(x.prettify() for x in content)


def SELECTOR_SecondaryStructure( name: str,
                                 sse: Union[str, Case],
                                 sse_type: str = 'HE',
                                 terminal_loops: bool = False
                                 ) -> str:
    """
    """
    if isinstance(sse, Case):
        sse = sse.secondary_structure

    return textwrap.dedent("""\
        <SecondaryStructure name="{}" overlap="0" minE="1"  minH="1" ss="{}" include_terminal_loops="{}"
        use_dssp="0" pose_secstruct="{}" />""").format(name, sse_type, int(terminal_loops), sse)


def MOVER_SetSecStructEnergies( name: str, score: str, natbias: float, case: Case ) -> str:
    """
    """
    data = dict(zip(['ss_pair', 'hh_pair', 'hss_triplets'], case.sse_pairing))
    data['sse'] = case.secondary_structure
    data['score'] = score
    data['name'] = name
    data['natbias'] = natbias
    return Template(textwrap.dedent("""\
        <SetSecStructEnergies name="{{name}}" scorefxn="{{score}}"
            secstruct="{{sse}}" use_dssp="0"
            {% if hh_pair|length > 0 %}hh_pair="{{hh_pair}}"{% endif %}
            {% if ss_pair|length > 0 %}ss_pair="{{ss_pair}}"{% endif %}
            {% if hss_triplets|length > 0 %}hss_triplets="{{hss_triplets}}"{% endif %}
            ss_from_blueprint="0"
            {% if ss_pair|length > 0 %}natbias_ss="{{natbias}}"{% endif %}
            {% if hh_pair|length > 0 %}natbias_hh="{{natbias}}"{% endif %}
            {% if hss_triplets|length > 0 %}natbias_hs="{{natbias}}"{% endif %}
        />""")).render(data)


def MOVER_PeptideStubMover( name: str, case: Case ) -> str:
    """
    """
    return Template(textwrap.dedent("""\
        <PeptideStubMover name="{{name}}" reset="false">
        {% for item in insert %}
        <Insert resname="GLY" repeat="1" jump="false" anchor_rsd="{{item}}" anchor_atom="C" connecting_atom="N" />
        {% endfor %}
        </PeptideStubMover>""")).render({'insert': [i for i, ltr in enumerate(case.secondary_structure) if ltr == 'L'],
                                         'name': name})


def PROTOCOL_LayerDesign( case: Case ) -> ScriptPieces:
    """
    """
    residueselectors = [textwrap.dedent("""\
        <Layer name="surface" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>"""),
                        SELECTOR_SecondaryStructure('sheet', case, 'E'),
                        SELECTOR_SecondaryStructure('entire_helix', case, 'H'),
                        SELECTOR_SecondaryStructure('entire_loop', case, 'L', True), textwrap.dedent("""\
        <And name="helix_cap" selectors="entire_loop">
            <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/></And>
        <And name="helix_start" selectors="entire_helix">
            <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/></And>
        <And name="helix" selectors="entire_helix"><Not selector="helix_start"/></And>
        <And name="loop" selectors="entire_loop"><Not selector="helix_cap"/></And>
        """)]

    taskoperations = [textwrap.dedent("""\
        <DesignRestrictions name="layer_design">
            <Action selector_logic="surface AND helix_start" aas="DEHKPQR"/>
            <Action selector_logic="surface AND helix" aas="EHKQR"/>
            <Action selector_logic="surface AND sheet" aas="EHKNQRST"/>
            <Action selector_logic="surface AND loop" aas="DEGHKNPQRST"/>
            <Action selector_logic="boundary AND helix_start" aas="ADEHIKLMNPQRSTVWY"/>
            <Action selector_logic="boundary AND helix" aas="ADEHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND sheet" aas="DEFHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND loop" aas="ADEFGHIKLMNPQRSTVWY"/>
            <Action selector_logic="core AND helix_start" aas="AFILMPVWY"/>
            <Action selector_logic="core AND helix" aas="AFILMVWY"/>
            <Action selector_logic="core AND sheet" aas="FILMVWY"/>
            <Action selector_logic="core AND loop" aas="AFGILMPVWY"/>
            <Action selector_logic="helix_cap" aas="DNST"/>
        </DesignRestrictions>"""), ]

    return ScriptPieces({'residueselectors': residueselectors, 'taskoperations': taskoperations})


def PROTOCOL_BasicFilters( case: Case, suffix: str = '' ) -> ScriptPieces:
    """
    """
    sse = case.secondary_structure

    residueselectors = textwrap.dedent("""\
        <Layer name="surface{suffix}" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary{suffix}" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core{suffix}" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>
        """).format(suffix=suffix)

    filters = textwrap.dedent("""\
    <PackStat name="pack{suffix}" confidence="0." />
    <CavityVolume name="cav_vol{suffix}" confidence="0." />
    <SecondaryStructure name="sse_match{suffix}" ss="{sse1}" compute_pose_secstruct_by_dssp="true" confidence="0." />
    """).format(sse1=sse, suffix=suffix)

    movers = [textwrap.dedent("""\
    <LabelPoseFromResidueSelectorMover name="labelcore{suffix}" property="CORE" residue_selector="core{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelboundary{suffix}" property="BOUNDARY" residue_selector="boundary{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelsurface{suffix}" property="SURFACE" residue_selector="surface{suffix}" />
    <DisplayPoseLabelsMover name="labeldump{suffix}" use_dssp="1" write="1" />
    """).format(suffix=suffix), ]
    if TBcore.get_option('psipred', 'script', in_path_none=True) is not None:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" cmd="{psipred}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix, psipred=TBcore.get_option('psipred', 'script')))
    else:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix))

    protocols = textwrap.dedent("""\
    <Add mover="labelcore{suffix}"/>
    <Add mover="labelboundary{suffix}"/>
    <Add mover="labelsurface{suffix}"/>
    <Add mover="labeldump{suffix}"/>
    <Add mover="sse_report{suffix}"/>
    <Add filter="pack{suffix}" />
    <Add filter="cav_vol{suffix}" />
    <Add filter="sse_match{suffix}" />
    """).format(suffix=suffix)

    return ScriptPieces({'residueselectors': [residueselectors, ], 'filters': [filters, ],
                         'movers': movers, 'protocols': [protocols, ]})
