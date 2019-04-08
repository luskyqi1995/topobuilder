# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict
import itertools
import sys

# External Libraries
import networkx as nx

# This Library
from topobuilder.case import Case
from topobuilder.case.schema import CoordinateSchema
import topobuilder.core as TBcore
import topobuilder.utils as TButil

__all__ = ['metadata', 'apply', 'case_apply']


def metadata() -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: Funtion on the expected input connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    """
    def isngl( count ):
        return True

    return {'name': 'nomenclator',
            'Itags': [],
            'Otags': ['equivalent_connectivities'],
            'Isngl': isngl,
            'Osngl': True}


def apply( cases: List[Case],
           prtid: int,
           representatives: bool = False,
           **kwargs ) -> List[Case]:
    """Generate all possible connectivities in the Case.
    """
    TButil.plugin_title(__file__, len(cases))

    new_cases = []
    for i, case in enumerate(cases):
        new_cases.extend(case_apply(case, representatives))

    for i, case in enumerate(new_cases):
        new_cases[i] = case.set_protocol_done(prtid)

    TButil.plugin_case_count(len(new_cases))
    return new_cases


def case_apply( case: Case, representatives: bool = False ) -> List[Case]:
    """Generate all possible connectivities in the Case.
    """
    case = Case(case)

    new_cases = []
    # If connectivities are pre-specified, only make those.
    if case.connectivity_count > 0:
        new_cases.extend(eval_representatives(case.apply_topologies(), representatives))
    else:
        new_cases.extend(eval_representatives(explore_connectivities(case), representatives))
    return new_cases


def explore_connectivities( case: Case ) -> List[Case]:
    """
    """
    case = Case(case)
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Exploring connectivities for {}\n'.format(case.architecture_str))

    # Create Graph
    G = make_graph(case)

    # Run connectivities
    topologies = []
    for n1, n2 in itertools.combinations(G.nodes(), 2):
        for p in search_paths(G, n1, n2):
            topologies.append(p)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Explored a total of {} unique connectivities\n'.format(len(topologies)))

    return case.add_secured_topologies(topologies).apply_topologies()


def eval_representatives( cases: List[Case], representatives: bool ) -> List[Case]:
    """
    """
    # Calculate representatives
    reps = {}
    for c in cases:
        reps.setdefault(c.directionality_profile, []).append(c.connectivities_str[0])
    reps = [list(sorted(set(x))) for x in reps.values()]

    for c in cases:
        c.data.setdefault('metadata', {}).setdefault('equivalent_connectivities', {})
        cn = c.connectivities_str[0]
        gr = [group for group in reps if cn in group][0]
        c.data['metadata']['equivalent_connectivities']['others'] = gr
        c.data['metadata']['equivalent_connectivities']['representative'] = (gr[0] == cn)

    if representatives:
        cases = [c for c in cases if c['metadata.equivalent_connectivities.representative']]
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Selecting {} representatives\n'.format(len(cases)))
    return cases


def make_graph( case: Case ) -> nx.Graph:
    """
    """
    # We will need distances between SSE
    case = case.cast_absolute()
    schema = CoordinateSchema()
    maxl = case['configuration.defaults.distance.max_loop']

    a = case['topology.architecture']
    G = nx.Graph()

    # Inner layer links
    for layer in a:
        for pair in itertools.combinations(layer, 2):
            # Mix types in layers might need special considerations
            if pair[0]['type'] == pair[1]['type']:
                # Betas can connect to +2 (for greek keys)
                if pair[0]['type'] in ['E', ]:
                    if abs(int(pair[0]['id'][1:-1]) - int(pair[1]['id'][1:-1])) > 2:
                        continue
                # Alphas can connect to +1
                if pair[0]['type'] in ['H', 'G']:
                    if abs(int(pair[0]['id'][1:-1]) - int(pair[1]['id'][1:-1])) > 1:
                        continue
            if schema.distance(pair[0]['coordinates'], pair[1]['coordinates']) <= maxl:
                G.add_edge(pair[0]['id'], pair[1]['id'])

    # Layer to Layer links (consecutive layers only)
    for ly1, ly2 in zip(a[:-1], a[1:]):
        for sse1, sse2 in itertools.product(*[ly1, ly2]):
            if schema.distance(sse1['coordinates'], sse2['coordinates']) <= maxl:
                G.add_edge(sse1['id'], sse2['id'])
    return G


def search_paths(G: nx.Graph, n1: str, n2: str) -> List:
    conns = []
    for path in nx.all_simple_paths(G, n1, n2):
        if len(path) == nx.number_of_nodes(G):
            conns.append(path)
            conns.append(list(reversed(path)))
    return conns
