# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import sys

# External Libraries
import colorama as cl

# This Library
import topobuilder.core as TBcore

cl.init(autoreset=True)


def plugin_title( plugin_path: str, cases: int ):
    """Print on-screen the plugin's name.

    :param str plugin_path: ``__file__`` of the plugin's main file.
    :param int cases: Number of cases to which the plugin will be applied.
    """
    if TBcore.get_option('system', 'verbose'):
        name = os.path.basename(os.path.dirname(plugin_path)).upper()
        name = ' '.join(['  TB PLUGIN:', name, ' '])
        bord = ''.join(['_', ] * len(name))

        sys.stdout.write('\n')
        sys.stdout.write(cl.Style.BRIGHT + bord + '\n')
        sys.stdout.write(cl.Style.BRIGHT + name + '\n')
        sys.stdout.write(cl.Style.BRIGHT + bord + '\n')
        sys.stdout.write(cl.Style.DIM + '* batch applied to {:03d} cases\n\n'.format(cases))


def plugin_warning( text: str ):
    """Format a warning and manage follow up behaviour.

    :param str text: Warning text.
    """
    sys.stdout.write(cl.Fore.GREEN + text + '\n')
    if TBcore.get_option('system', 'strict'):
        sys.exit(-1)
