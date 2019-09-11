# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Optional, Dict
from pathlib import Path
import os
import sys

# External Libraries
import matplotlib.pyplot as plt

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from . import plot_types as pts

__all__ = ['metadata', 'apply']

_PLT_TYPES_ = ['sketchXZ', 'sketchXY']


def metadata() -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: When :data:`True`, input requires single connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    """
    def isngl( count ):
        return True

    return {'name': 'plotter',
            'Itags': [],
            'Otags': [],
            'Isngl': isngl,
            'Osngl': False}


def apply( cases: List[Case],
           prtid: int,
           outfile: Optional[Union[str, Path]] = None,
           prefix: Optional[str] = None,
           plot_types: Optional[List[str]] = None,
           **kwargs ) -> List[Case]:
    """Generate visual representations of the Case.
    """
    TButil.plugin_title(__file__, len(cases))

    # File management
    if outfile is None:
        outfile = cases[0].main_path.joinpath('images').resolve()
        outfile.mkdir(parents=True, exist_ok=True)
    if isinstance(outfile, str):
        outfile = Path(outfile).resolve()

    outformat = TBcore.get_option('system', 'image')

    # Checking available plot_types
    plot_types = [_PLT_TYPES_[0], ] if plot_types is None else plot_types
    if len(set(plot_types).difference(_PLT_TYPES_)) > 0:
        raise ValueError('Requested unknown plot format. '
                         'Available are: {}'.format(','.join(_PLT_TYPES_)))

    for ptype in plot_types:
        if outfile.is_dir():
            prefix = prefix if prefix is not None else ".".join([str(os.getppid()), '{:02d}'.format(prtid)])
            thisoutfile = outfile.joinpath(".".join([prefix, ptype + outformat]))
        else:
            thisoutfile = Path(str(outfile) + '.' + ptype + outformat)
        thisoutfile.parent.mkdir(parents=True, exist_ok=True)
        if not TBcore.get_option('system', 'overwrite') and thisoutfile.is_file():
            sys.stderr.write('Unable to overwrite file {}: Already exists\n'.format(thisoutfile))
            continue

        fig, ax = getattr(pts, ptype)(cases, **kwargs.pop(ptype, {}))
        plt.tight_layout()
        plt.savefig(str(thisoutfile), dpi=300)

        TButil.plugin_imagemaker('Creating new image at: {}'.format(str(thisoutfile)))

    for i, case in enumerate(cases):
        cases[i] = case.set_protocol_done(prtid)
    return cases
