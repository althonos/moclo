# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import copy
import functools
import typing

import six
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

if typing.TYPE_CHECKING:
    from typing import Dict, List, Text, Union  # noqa: F401
    from Bio.Seq import Seq  # noqa: F401


def _ambiguous(func):
    @functools.wraps(func)
    def newfunc(self, *args, **kwargs):
        raise TypeError("ambiguous operation: '{}'".format(func.__name__))

    return newfunc


class CircularRecord(SeqRecord):
    """A derived ``SeqRecord`` that contains a circular DNA sequence.

    It handles the `in` operator as expected, and removes the implementation
    of the `+` operator since circular DNA sequence do not have an end to
    append more nucleotides to. In addition, it overloads the `>>` and `<<`
    operators to allow rotating the sequence and its annotations, effectively
    changing the *0* position.

    See Also:
        ``Bio.SeqRecord.SeqRecord`` documentation on the `Biopython wiki
        <https://biopython.org/wiki/SeqRecord>`_.
    """

    # Patch constructor to allow constructing a CircularRecord from
    # a SeqRecord argument (other arguments are ignored).

    def __init__(
        self,
        seq,  # type: Union[Seq, SeqRecord]
        id="<unknown id>",  # type: Text
        name="<unknown name>",  # type: Text
        description="<unknown description>",  # type: Text
        dbxrefs=None,  # type: List[Text]
        features=None,  # type: List[SeqFeature]
        annotations=None,  # type: Dict[Text, Any]
        letter_annotations=None,
    ):
        # type: (...) -> None
        """Create a new `CircularRecord` instance.

        If given a `SeqRecord` as the first argument, it will simply copy all
        attributes of the record. This allows using `Bio.SeqIO.read` to open
        records, then loading them into a `CircularRecord`.
        """
        if isinstance(seq, SeqRecord):
            self.__init__(  # noqa: T484
                seq.seq,
                seq.id,
                seq.name,
                seq.description,
                copy.deepcopy(seq.dbxrefs),
                copy.deepcopy(seq.features),
                copy.deepcopy(seq.annotations),
                copy.deepcopy(seq.letter_annotations),
            )
        else:
            if annotations is not None:
                topology = annotations.get("topology", "circular")
                if topology.lower() != "circular":
                    raise ValueError("record does not describe a circular sequence")
            super(CircularRecord, self).__init__(
                seq,
                id,
                name,
                description,
                dbxrefs,
                features,
                annotations,
                letter_annotations,
            )

    # Patch methods that are ambiguous with a non-linear sequence.

    @_ambiguous
    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Since adding an arbitrary sequence to a plasmid is ambiguous (there is
        no sequence end), trying to add a sequence to a ``CircularRecord``
        will raise a `TypeError`.
        """
        return NotImplemented

    @_ambiguous
    def __radd__(self, other):  # noqa: D105
        """Add another sequence or string to this sequence (from the left).

        Since adding an arbitrary sequence to a plasmid is ambiguous (there is
        no sequence end), trying to add a sequence to a ``CircularRecord``
        will raise a `TypeError`.
        """
        return NotImplemented

    # Patch other methods to work as intended

    def __contains__(self, char):  # noqa: D105
        """Implement the `in` keyword, searches the sequence.
        """
        return len(char) <= len(self) and char in str(self.seq) * 2

    def __getitem__(self, index):  # noqa: D105
        """Return a sub-sequence or an individual letter.

        The sub-sequence is always returned as a ``SeqRecord``, since it is
        probably not circular anymore.
        """
        rec = super(CircularRecord, self).__getitem__(index)
        if isinstance(index, slice):
            annotations = copy.deepcopy(rec.annotations)
            if "topology" in annotations:
                annotations["topology"] = "linear"
            return SeqRecord(
                rec.seq,
                rec.id,
                rec.name,
                rec.description,
                copy.deepcopy(rec.dbxrefs),
                copy.deepcopy(rec.features),
                annotations,
                copy.deepcopy(rec.letter_annotations),
            )
        else:
            return rec

    def reverse_complement(
        self,
        id=False,
        name=False,
        description=False,
        features=True,
        annotations=False,
        letter_annotations=True,
        dbxrefs=False,
    ):
        """Return a new ``CircularRecord`` with reverse complement sequence.
        """
        return type(self)(
            super(CircularRecord, self).reverse_complement(
                id=id,
                name=name,
                description=description,
                features=features,
                annotations=annotations,
                letter_annotations=letter_annotations,
                dbxrefs=dbxrefs,
            )
        )

    # Additional methods

    def __lshift__(self, index):
        """Rotate the sequence counter-clockwise, preserving annotations.
        """
        return self >> (-index % len(self.seq))

    def __rshift__(self, index):
        """Rotate the sequence clockwise, preserving annotations.
        """
        index %= len(self.seq)  # avoid unnecessary cycles

        if index == 0:
            return self
        elif index < 0:
            return self << -index

        newseq = self.seq[-index:] + self.seq[:-index]
        newfeats = []
        newletan = {
            k: v[index:] + v[:index] for k, v in six.iteritems(self.letter_annotations)
        }

        for feature in self.features:
            loc = feature.location
            if loc is None:
                newloc = None
            elif feature.type == "source" and loc.start == 0 and loc.end == len(self):
                newloc = loc
            else:
                _newloc = []
                for part in (loc + index).parts:
                    if part.end >= len(newseq) and part.start >= len(newseq):
                        r = part.start // len(newseq)  # remainder is used to
                        _newloc.append(
                            FeatureLocation(  # make sure that part.end
                                start=part.start
                                - r * len(newseq),  # is always after part.start
                                end=part.end
                                - r * len(newseq),  # even on additional end
                                strand=part.strand,  # overlap
                                ref=part.ref,
                                ref_db=part.ref_db,
                            )
                        )
                    else:
                        _newloc.append(part)
                newloc = _newloc[0] if len(_newloc) == 1 else CompoundLocation(_newloc)
            newfeats.append(
                SeqFeature(
                    location=newloc,
                    type=feature.type,
                    id=feature.id,
                    qualifiers=feature.qualifiers,
                )
            )

        return type(self)(
            seq=newseq,
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs,
            features=newfeats,
            annotations=self.annotations,
            letter_annotations=newletan,
        )
