# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import abc
import typing

import six
from property_cached import cached_property

from .. import errors
from ..regex import DNARegex

if typing.TYPE_CHECKING:
    from typing import Text, Type  # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401
    from ..regex import SeqMatch  # noqa: F401


@six.add_metaclass(abc.ABCMeta)
class StructuredRecord(object):
    """A DNA record with a specific structure.
    """

    _regex = None

    def __init__(self, record):
        # type: (SeqRecord) -> None
        self.record = record
        self.seq = record.seq

    @classmethod
    @abc.abstractmethod
    def structure(cls):
        # type: () -> Text
        """Get the expected record structure as a DNA pattern in regex syntax.
        """
        return NotImplemented

    @classmethod
    def _get_regex(cls):
        if cls._regex is None:
            cls._regex = DNARegex(cls.structure())
        return cls._regex

    @cached_property
    def _match(self):
        # type: () -> SeqMatch[SeqRecord]
        topology = self.record.annotations.get("topology", "circular").lower()
        match = self._get_regex().search(self.record, linear=topology != "circular")
        if match is None:
            details = "does not match '{}' structure".format(type(self).__name__)
            raise errors.InvalidSequence(self.record, details=details)
        return match

    def is_valid(self):
        # type: () -> bool
        """Check if the wrapped record follows the required class structure.

        Returns:
            `bool`: `True` if the record is valid, `False` otherwise.

        """
        try:
            return self._match is not None
        except errors.InvalidSequence:
            return False
