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
    from typing import Mapping, Match, Optional, Sequence, Text, Tuple  # noqa: F401


# typing.TypeVar: a generic sequence
_S = typing.TypeVar("_S", bound=typing.Sequence)


class SeqMatch(typing.Generic[_S]):
    def __init__(self, match, rec, shift=0):
        # type: (Match, _S, int) -> None  # type: ignore
        self.match = match
        self.shift = shift
        self.rec = rec  # type: _S

    def end(self):
        # type: () -> int
        return self.match.end()

    def start(self):
        # type: () -> int
        return self.match.start()

    def span(self, index=0):
        # type: (int) -> Tuple[int, int]
        return self.match.span(index)

    def group(self, index=0):
        # type: (int) -> _S
        span = self.match.span(index)
        if span[1] >= span[0] >= len(self.rec):
            return self.rec[
                span[0] % len(self.rec) : span[1] % len(self.rec)
            ]  # type: ignore
        elif span[1] >= len(self.rec) > span[0]:
            return (
                self.rec[: span[1] % len(self.rec)] + self.rec[span[0] :]
            )  # type: ignore
        else:
            return self.rec[span[0] : span[1]]  # type: ignore


class DNARegex(object):
    """A DNA pattern matcher.
    """

    # dict: single-letter code to regex syntax for nucleotides
    _lettermap = {
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "K": "[GT]",
        "M": "[AC]",
        "N": "[ACGT]",
        "R": "[AG]",
        "S": "[CG]",
        "V": "[ACG]",
        "W": "[AT]",
        "Y": "[CT]",
    }  # type: Mapping[Text, Text]

    @classmethod
    def _transcribe(cls, pattern):
        # type: (Text) -> Text
        target = ["(?i)"]
        for letter in pattern:
            target.append(cls._lettermap.get(letter, letter))
        return "".join(target)

    def __init__(self, pattern):
        # type: (Text) -> None
        self.pattern = pattern
        self.regex = re.compile(self._transcribe(pattern))

    def search(self, string, pos=0, endpos=six.MAXSIZE, linear=True):
        # type: (_S, int, int) -> Optional[SeqMatch[_S]]

        if not isinstance(string, (Bio.Seq.Seq, Bio.SeqRecord.SeqRecord)):
            t = type(string).__name__
            raise TypeError("can only match Seq or SeqRecord, not '{}'".format(t))

        if isinstance(string, Bio.SeqRecord.SeqRecord):
            data = str(string.seq)
        else:
            data = str(string)

        if not linear or isinstance(string, CircularRecord):
            data *= 2

        for i in range(pos, min(len(string), endpos)):
            match = self.regex.match(data, i, i + len(string))
            if match is not None:
                return SeqMatch(match, string)

        return None
