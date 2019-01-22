# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
from pathlib import Path
from typing import Optional, Tuple, Dict

# External Libraries
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
from matplotlib.patches import PathPatch
from matplotlib.colors import ColorConverter, to_hex

# This Library
from .case import Case

__all__ = ['case_template', 'plot_case_sketch']

def case_template( name: str,
                   architecture: Optional[str] = None,
                   topology: Optional[str] = None,
                   corrections: Optional[Dict] = None,
                   format: str = 'yaml',
                   make_absolute: bool = False
                   ) -> Tuple[Case, Path]:
    """Generate a :class:`.Case`.

    :param str name: Identifier of the case.
    :param str architecture: Definition of unconnected, unordered secondary structure.
    :param str topology: Definition of connected,ordered secondary structure. If provided, it will
        overwrite ``architecture``.
    :param dict corrections: Corrections to apply to the default case, identified by the SSE id.
    :param str format: Format of the output file (``yaml`` or ``json``).
    :param bool make_absolute: If :data:`True`, coordinates and shifts are
        defined as absolute positions.

    :return: :class:`.Case` and :class:`Path` generated filename.
    """

    # Create the case
    case = Case(name)

    # Architecture-defined case
    case = case.add_architecture(architecture)

    # Topology-defined case (architecture + connectivity)
    # If match, it is appended to the pre-defined architecture
    case = case.add_topology(topology)

    # Apply corrections if any
    case = case.apply_corrections(corrections)

    if make_absolute:
        case = case.cast_absolute()

    # Output
    outfile = case.write(name, format)

    return case, outfile


def plot_case_sketch( case: Case,
                      ax: Optional[plt.Axes]=None,
                      beta_fill: Optional[str]='red',
                      beta_edge: Optional[str]='black',
                      alpha_fill: Optional[str]='blue',
                      alpha_edge: Optional[str]='black'
                      ) -> Tuple[plt.Figure, plt.Axes]:
    """
    """
    def make_triangle(y, x, rot_deg, fcolor, ecolor, scale):
        unit_triangle = Path.unit_regular_polygon(3)
        path = Path(unit_triangle.vertices * scale, unit_triangle.codes)
        trans = Affine2D().translate(x, y).rotate_deg_around(x, y, rot_deg)
        t_path = path.transformed(trans)
        patch = PathPatch(t_path, facecolor=fcolor, edgecolor=ecolor, zorder=2)
        return patch

    if ax is None:
        fig = plt.figure()
        ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    ax.set_aspect('equal', adjustable='box')

    shp = case.center_shape
    xmax = max([shp[l]['right'] for l in shp]) + 3
    xmin = min([shp[l]['left'] for l in shp]) - 3
    ymax = 0
    ymin = 0

    for layer in case.cast_absolute()['topology.architecture']:
        for sse in layer:
            ymax = sse['coordinates']['z'] if sse['coordinates']['z'] > ymax else ymax
            ymin = sse['coordinates']['z'] if sse['coordinates']['z'] < ymin else ymin
            rotation = 180 if sse['tilt']['x'] > 90 and sse['tilt']['x'] < 270 else 0
            if sse['type'] == 'H':
                c = plt.Circle((sse['coordinates']['x'], sse['coordinates']['z']), radius=3,
                               facecolor=alpha_fill, edgecolor=alpha_edge, zorder=2)
                ax.add_artist(c)
                p = make_triangle(sse['coordinates']['z'], sse['coordinates']['x'], rotation,
                                  lighten_color(alpha_fill, 0.5), alpha_edge, 2)
                ax.add_artist(p)
            if sse['type'] == 'E':
                p = make_triangle(sse['coordinates']['z'], sse['coordinates']['x'], rotation,
                                  beta_fill, beta_edge, 2)
                ax.add_artist(p)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax + 4, ymin - 4)
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.grid(zorder=0)
    return ax


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    https://stackoverflow.com/a/49601444/2806632

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except Exception:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
