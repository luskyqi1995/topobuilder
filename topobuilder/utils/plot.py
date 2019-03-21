# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Union, Tuple, Optional, List, Dict
import math

# External Libraries
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from rstoolbox.components import FragmentFrame
import pandas as pd
import networkx as nx

# This Library
import topobuilder.core as TBcore
from .plugins import plugin_imagemaker


__all__ = ['plot_fragment_templates', 'plot_loop_length_distribution', 'plot_match_bin',
           'plot_geometric_distributions', 'plot_angle_network']


def plot_fragment_templates( dfsmall: Union[FragmentFrame, pd.DataFrame],
                             dflarge: Union[FragmentFrame, pd.DataFrame],
                             prefix: Union[Path, str],
                             write: bool = True
                             ) -> Tuple[plt.Figure, Path]:
    """
    """
    fig = plt.figure(figsize=(20, 8))
    ax0 = plt.subplot2grid((2, 1), (0, 0), fig=fig)
    pp = dfsmall.groupby(['position', 'source'])['position'].count() / dfsmall['size'].values[0]
    pp = pp.unstack().fillna(0)
    color = ListedColormap(sns.diverging_palette(220, 20, n=len(pp.columns)).as_hex())
    pp.plot(kind='bar', stacked=True, ax=ax0, colormap=color)
    ax0.set_xlim(-0.5, dfsmall.position.max() - dfsmall.position.min() + 1)
    ax0.set_ylabel('3mers coverage')
    ax0.get_xaxis().set_visible(False)
    ax0.get_legend().remove()

    ax1 = plt.subplot2grid((2, 1), (1, 0), sharex=ax0, fig=fig)
    pp = dflarge.groupby(['position', 'source'])['position'].count() / dflarge['size'].values[0]
    pp = pp.unstack().fillna(0)
    pp.plot(kind='bar', stacked=True, ax=ax1, colormap=color)
    ax1.set_xlim(-0.5, dflarge.position.max() - dflarge.position.min() + 1)
    ax1.set_ylabel('9mers coverage')
    ax1.set_xticks(range(0, dflarge.position.max() - dflarge.position.min() + 1, 5))
    ax1.set_xticklabels(range(dflarge.position.min(), dflarge.position.max() + 1, 5), rotation=0)
    ax1.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=len(pp.columns))
    plt.subplots_adjust(wspace=0, hspace=0.025)

    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
    if write:
        plugin_imagemaker('fragment templates image summary at {}'.format(imagename))
        plt.savefig(imagename, dpi=300)
    return fig, imagename


def plot_loop_length_distribution( dfloop: pd.DataFrame,
                                   pick: int,
                                   prefix: Union[Path, str],
                                   title: str,
                                   write: bool = True
                                   ) -> Tuple[plt.Figure, Path]:
    """
    """
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    sns.countplot(x='loop_length', data=dfloop, ax=ax, palette="PRGn")
    cnt = dfloop[dfloop['loop_length'] == pick]['length_count'].values[0]
    pick = [int(x.get_text()) for x in ax.get_xticklabels()].index(pick)
    ax.plot(pick, cnt, marker=11, color='black')
    ax.set_title(title)

    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
    if write:
        plugin_imagemaker('loop image summary at {}'.format(imagename))
        plt.savefig(imagename, dpi=300)
    return fig, imagename


def plot_match_bin( master_match: pd.DataFrame,
                    prefix: Union[Path, str],
                    expected: Union[int, float],
                    groupby: Optional[List] = None,
                    write: bool = True
                    ) -> Tuple[plt.Figure, Path]:
    """
    """
    fig = plt.figure(figsize=[15, 5])
    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)

    master_match = master_match.sort_values('rmsd')
    df = master_match.groupby(groupby).head(1) if groupby is not None else master_match

    sns.countplot(x="bin", data=df, ax=ax, order=['close', 'mid', 'far', 'extreme'], palette="Set3")
    ax.set_xticklabels(['close\n[0, 2)', 'mid\n[2, 2.5)', 'far\n[2.5, 3)', 'extreme\n[3, 5)'])
    ax.set_ylim(top=expected)
    for p in ax.patches:
        pp = int(0 if math.isnan(p.get_height()) else p.get_height())
        ypos = pp + 1000 if pp + 1000 < expected - 500 else pp - 1000
        ax.annotate('{:d}'.format(pp), (p.get_x() + 0.37, ypos))

    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
    if write:
        plugin_imagemaker('MASTER match image summary at {}'.format(imagename))
        plt.savefig(imagename, dpi=300)
    return fig, imagename


def plot_geometric_distributions( df: pd.DataFrame,
                                  prefix: Union[Path, str],
                                  write: bool = True
                                  ) -> Tuple[plt.Figure, Path]:
    """
    """
    ordering = sorted(df.sse.unique())
    fig = plt.figure(figsize=[15, 15])
    grid = (3, 2)

    for i, l in enumerate(['layer', 'floor', 'side']):
        ax = plt.subplot2grid(grid, (i, 0), fig=fig)
        sns.violinplot(x='sse', y='angles_{}'.format(l), hue='bin', data=df, palette="Set3", order=ordering, ax=ax, cut=1)
        ax.legend().remove()
        ax.set_ylabel('angle')
        ax.set_title('angles_{}'.format(l))
        ax = plt.subplot2grid(grid, (i, 1), fig=fig)
        sns.violinplot(x='sse', y='points_{}'.format(l), hue='bin', data=df, palette="Set3", order=ordering, ax=ax, cut=0)
        ax.set_ylabel('distance')
        ax.set_title('points_{}'.format(l))
        if i != 0:
            ax.legend().remove()
        else:
            ax.legend(ncol=len(df.bin.unique()))
    plt.tight_layout()

    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
    if write:
        plugin_imagemaker('Geometric distributions image summary at {}'.format(imagename))
        plt.savefig(imagename, dpi=300)
    return fig, imagename


def plot_angle_network( network: nx.DiGraph,
                        node_positions: Dict,
                        sse_list: List[str],
                        prefix: Union[Path, str],
                        write: bool = True
                        ) -> Tuple[plt.Figure, Path]:
    """
    """
    fig = plt.figure(figsize=[25, 25])
    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)

    nx.draw(network, pos=node_positions, ax=ax)
    ax.set_axis_on()
    ax.set_xticks(range(len(sse_list)))
    ax.set_xticklabels(sse_list)
    ax.set_xlabel('SSE')
    ax.set_ylabel('zeta angle')

    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
    if write:
        plugin_imagemaker('Layer angle network image summary at {}'.format(imagename))
        plt.savefig(imagename, dpi=300)
    return fig, imagename
