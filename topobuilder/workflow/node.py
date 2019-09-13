# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import abc
import sys
from typing import List, Dict

# External Libraries
from logbook import Logger, StreamHandler

# This Library
from .log import logger_group, logger_level, LoggedError

StreamHandler(sys.stdout).push_application()

__all__ = ['Node', 'NodeDataError', 'NodeOptionsError', 'NodeMissingError']


class Node( abc.ABC ):
    """Defines each individual step of the :class:`.Pipeline`.

    This is the base class from which all plugins must derive.
    """

    def __init__( self, tag: int ):
        super(Node, self).__init__()
        self.tag = tag
        self.log = Logger(f'Node => {self.tag:02d} - {type(self).__name__}')
        logger_group.add_logger(self.log)
        logger_level(logger_group)
        self.log.info(f'Starting new work node.')

    def check( self, dummy: List[Dict] ) -> List[Dict]:
        """Evaluates the feasability of executing the :class:`.Node`.

        :param dummy: Shape of the data provided by the previous :class:`.Node`.

        :return: Shape of the data after being modified by the current :class:`.Node`.

        :raises:
            :NodeDataError: If the input data shape does not match the expected.
        """
        self.log.info('Checking node viability according to input data.')
        back = []
        self.log.debug(f'Checking a total of {len(dummy):04d} cases.')
        for dt in dummy:
            back.append(self.single_check(dt))
        return back

    @abc.abstractmethod
    def single_check( self, dummy: Dict ) -> Dict:
        """Evaluates the feasability of executing the :class:`.Node`.

        :param dummy: Shape of the data provided by the previous :class:`.Node`.

        :return: Shape of the data after being modified by the current :class:`.Node`.

        :raises:
            :NodeDataError: If the input data shape does not match the expected.
        """
        raise NotImplementedError()

    def execute( self, data: List[Dict] ) -> List[Dict]:
        """Process all the data.

        :param data: Data to execute.

        :return: Modified data.
        """
        self.log.info('Executing node.')
        self.log.debug(f'Executing a total of {len(data):04d} cases.')
        back = []
        for dt in data:
            back.append(self.single_execute(dt))
        self.log.debug(f'Generated a total of {len(back):04d} cases.')
        return back

    @abc.abstractmethod
    def single_execute( self, data: Dict ) -> Dict:
        """Individually process each data entry.

        :param data: Data to execute.

        :return: Modified data.
        """
        raise NotImplementedError()


class NodeDataError(LoggedError):
    """Raises when the :class:`.Node` is not being provided the rigth data shape.
    """


class NodeOptionsError(LoggedError):
    """Raises when the :class:`.Node` is not provided the necessary options to be executed.
    """


class NodeMissingError(LoggedError):
    """Raises when a requested :class:`.Node` cannot be found or is not the rigth type.
    """
