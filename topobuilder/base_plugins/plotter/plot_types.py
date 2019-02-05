# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Tuple
import sys

# External Libraries
import matplotlib.pyplot as plt

# This Library
import topobuilder.core as TBcore
from topobuilder.case import Case, plot_case_sketch

__all__ = ['sketch2d']


def sketch2d( cases: List[Case], **kwargs ) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    """
    ncases = len(cases)
    columns = kwargs.pop('columns', 2 if ncases > 1 else 1)
    grid = (int(ncases / columns) + ncases % columns, columns)
    if TBcore.get_option('topobuilder', 'verbose'):
        sys.stdout.write('Generating an image grid of: {0}x{1}\n'.format(grid[0], grid[1]))

    fsize = (kwargs.pop('width', 7.5 * grid[1]),
             kwargs.pop('hight', 7.5 * grid[0]))
    fig = plt.figure(figsize=fsize)
    axs = []
    ylim, xlim = [0, 0], [0, 0]
    for i, case in enumerate(cases):
        position = (int(i / grid[1]), i % grid[1])
        title = '{0}_{1:03d}'.format(case['configuration.name'], i + 1)
        if TBcore.get_option('topobuilder', 'verbose'):
            sys.stdout.write('Showing {0}-{3} in position: {1}x{2}\n'.format(title, position[0], position[1],
                                                                             case.architecture_str))

        ax = plt.subplot2grid(grid, position, fig=fig)
        axs.append(ax)
        plot_case_sketch(case, ax,
                         kwargs.pop('connections', False),
                         kwargs.pop('beta_fill', 'red'),
                         kwargs.pop('beta_edge', 'black'),
                         kwargs.pop('alpha_fill', 'blue'),
                         kwargs.pop('alpha_edge', 'black'),
                         kwargs.pop('connection_edge', 'blak'))
        ax.set_title(title)
        cy = ax.get_ylim()
        cx = ax.get_xlim()
        ylim = [ylim[0] if cy[0] < ylim[0] else cy[0],
                ylim[1] if cy[1] > ylim[1] else cy[1]]
        xlim = [xlim[0] if cx[0] > xlim[0] else cx[0],
                xlim[1] if cx[1] < xlim[1] else cx[1]]

    for ax in axs:
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlim(xlim[0], xlim[1])

    return fig, axs
