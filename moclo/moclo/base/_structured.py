# coding: utf-8

import abc
import typing

import cached_property
import six

from ..regex import DNARegex


@six.add_metaclass(abc.ABCMeta)
class _StructuredRecord(object):

    _structure = NotImplemented
    _regex = None

    def __init__(self, record):
        self.record = record
        self.seq = record.seq

    @cached_property.cached_property
    def _match(self):
        if self._regex is None:
            self._regex = DNARegex(self._structure)
        match = self._regex.search(self.record)
        if match is None:
            msg = "'{}' does not match the expected structure !"
            raise ValueError(msg.format(self.record.name or self.record.id))
        return match

    def is_valid(self):
        try:
            return self._match is not None
        except ValueError:
            return False
