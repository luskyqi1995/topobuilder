# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys
from pathlib import Path
from typing import Union, Tuple

# External Libraries
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from rstoolbox.components import FragmentFrame
import pandas as pd

# This Library
import topobuilder.core as TBcore


__all__ = ['plot_fragment_templates', 'plot_loop_length_distribution']


def plot_fragment_templates( dfsmall: FragmentFrame,
                             dflarge: FragmentFrame,
                             prefix: Union[Path, str],
                             write: bool = True
                             ) -> Tuple[plt.Figure, Path]:
    """
    """
    imagename = Path(str(prefix) + TBcore.get_option('system', 'image'))
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
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('PLOT: fragment templates image summary at {}\n'.format(imagename))
    if write:
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
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('PLOT: loop image summary at {}\n'.format(imagename))
    if write:
        plt.savefig(imagename, dpi=300)
    return fig, imagename
