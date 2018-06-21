# coding: utf-8

import abc

import cached_property
import six
from Bio.Seq import Seq

from ..record import CircularRecord
from ._structured import StructuredRecord


class AbstractModule(StructuredRecord):
    """An abstract molecular cloning module.
    """

    _level = None

    def overhang_start(self):
        return Seq(self._match.group(1))

    def overhang_end(self):
        return Seq(self._match.group(3))

    def target_sequence(self):
        start, end = self._match.span(2)
        if isinstance(self.record, CircularRecord):
            return (self.record << start)[:end-start]
        else:
            return self.record[start:end]


class Product(AbstractModule):
    """Level -1 module.
    """

    _level = -1


class Entry(AbstractModule):
    """Level 0 module.
    """

    _level = 0


class Cassette(AbstractModule):
    """Level 1 module.
    """

    _level = 1


class Multigene(AbstractModule):
    """Level 2 module.
    """

    _level = 2
