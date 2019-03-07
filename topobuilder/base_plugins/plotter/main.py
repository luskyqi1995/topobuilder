# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Optional
from pathlib import Path
import os
import sys

# External Libraries
import matplotlib.pyplot as plt

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
from . import plot_types as pts

__all__ = ['apply']

_PLT_TYPES_ = ['sketchXZ', 'sketchXY']


def apply( cases: List[Case],
           prtid: int,
           outfile: Optional[Union[str, Path]] = None,
           plot_types: Optional[List[str]] = None,
           **kwargs ) -> List[Case]:
    """Generate visual representations of the Case.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: PLOTTER ---\n')

    if outfile is None:
        outfile = cases[0].main_path.joinpath('images').resolve()
        outfile.mkdir(parents=True, exist_ok=True)
    if isinstance(outfile, str):
        outfile = Path(outfile).resolve()

    outformat = TBcore.get_option('system', 'image')

    plot_types = [_PLT_TYPES_[0], ] if plot_types is None else plot_types
    if len(set(plot_types).difference(_PLT_TYPES_)) > 0:
        raise ValueError('Requested unknown plot format. '
                         'Available are: {}'.format(','.join(_PLT_TYPES_)))

    for ptype in plot_types:
        if outfile.is_dir():
            thisoutfile = outfile.joinpath(".".join([str(os.getppid()), '{:02d}'.format(prtid), ptype + outformat]))
        else:
            thisoutfile = Path(str(outfile) + '.' + ptype + outformat)
        thisoutfile.parent.mkdir(parents=True, exist_ok=True)
        if not TBcore.get_option('system', 'overwrite') and thisoutfile.is_file():
            sys.stderr.write('Unable to overwrite file {}: Already exists\n'.format(thisoutfile))
            continue

        fig, ax = getattr(pts, ptype)(cases, **kwargs.pop(ptype, {}))
        plt.tight_layout()
        plt.savefig(str(thisoutfile), dpi=300)
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Creating new image at: {}\n'.format(str(thisoutfile)))

    for i, _ in enumerate(cases):
        cases[i].set_protocol_done(prtid)
    return cases
