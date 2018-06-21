# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import re
import typing

import Bio.Seq
import Bio.SeqRecord
import six

from .record import CircularRecord

if typing.TYPE_CHECKING:
    from typing import Mapping, Text


class DNARegex(object):
    """A DNA pattern matcher.
    """

    # dict: single-letter code to regex syntax for nucleotides
    _lettermap = {
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'K': '[GT]',
        'M': '[AC]',
        'N': '[ACGT]',
        'R': '[AG]',
        'S': '[CG]',
        'V': '[ACG]',
        'W': '[AT]',
        'Y': '[CT]'
    } # type: Mapping[Text, Text]

    @classmethod
    def _transcribe(cls, pattern):
        # type: (Text) -> Text
        target = ['(?i)']
        for letter in pattern:
            target.append(cls._lettermap.get(letter, letter))
        return ''.join(target)

    def __init__(self, pattern, linear=True):
        # type: (Text) -> None
        self.pattern = pattern
        self.regex = re.compile(self._transcribe(pattern))
        self.linear = True

    def search(self, string, pos=0, endpos=six.MAXSIZE):
        if isinstance(string, CircularRecord):
            return self.regex.search(str(string.seq) * 2)
        elif isinstance(string, Bio.SeqRecord.SeqRecord):
            return self.regex.search(str(string.seq) * (2 - self.linear))
        elif isinstance(string, Bio.Seq.Seq):
            return self.regex.search(str(string) * (2 - self.linear))
        else:
            return self.regex.search(string * (2 - self.linear))
