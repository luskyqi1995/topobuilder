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
from .plot_types import sketch2d

__all__ = ['apply']


def apply( cases: List[Case],
           prtid: int,
           outfile: Optional[Union[str, Path]] = None,
           outformat: Optional[str] = 'png',
           plot_types: Optional[List[str]] = None,
           **kwargs ) -> List[Case]:
    """
    """
    if outfile is None:
        outfile = Path().cwd().resolve()
    if isinstance(outfile, str):
        outfile = Path(outfile).resolve()

    if outformat not in ['png', 'svg']:
        raise ValueError('Output formats are limited to "png" and "svg"')

    plot_types = ['sketch2d'] if plot_types is None else plot_types
    plot_types = list(set([x.lower() for x in plot_types]))
    plot_avail = ['sketch2d']
    if len(set(plot_types).difference(plot_avail)) > 0:
        raise ValueError('Requested unknown plot format. '
                         'Available are: {}'.format(','.join(plot_avail)))

    for ptype in plot_types:
        if outfile.is_dir():
            outfile = outfile.joinpath(".".join([str(os.getppid()), str(prtid), outformat]))
        else:
            outfile = outfile.stem + '.' + outformat
        if ptype == 'sketch2d':
            fig, ax = sketch2d(cases, **kwargs.pop('sketch2d', {}))
        plt.tight_layout()
        plt.savefig(str(outfile), dpi=300)
        if TBcore.get_option('topobuilder', 'verbose'):
            sys.stdout.write('Creating new image at: {}\n'.format(str(outfile)))

    for i, _ in enumerate(cases):
        cases[i].set_protocol_done(prtid)
    return cases
