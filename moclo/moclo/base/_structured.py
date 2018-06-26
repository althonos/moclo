# coding: utf-8

import abc
import typing

import cached_property
import six

from .. import errors
from ..regex import DNARegex


@six.add_metaclass(abc.ABCMeta)
class StructuredRecord(object):

    _structure = NotImplemented
    _regexes = {}

    def __init__(self, record):
        self.record = record
        self.seq = record.seq

    @cached_property.cached_property
    def _match(self):
        regex = self._regexes.get(type(self))
        if regex is None:
            topology = self.record.annotations.get('topology', 'circular').lower()
            regex = self._regexes[type(self)] = DNARegex(self._structure, linear=topology != 'circular')
        match = regex.search(self.record)
        if match is None:
            details = "does not match '{}' structure".format(type(self).__name__)
            raise errors.InvalidSequence(self.record, details=details)
        return match

    def is_valid(self):
        try:
            return self._match is not None
        except errors.InvalidSequence:
            return False
