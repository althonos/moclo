# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import copy
import functools

import six
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

def _ambiguous(func):
    @functools.wraps(func)
    def newfunc(self, *args, **kwargs):
        raise TypeError("ambiguous operation: '{}'".format(func.__name__))
    return newfunc


class CircularRecord(SeqRecord):

    # Patch constructor to allow constructing a CircularRecord from
    # a SeqRecord argument (other arguments are ignored).

    def __init__(self,
                 seq,
                 id='<unknown id>',
                 name='<unknown name>',
                 description='<unknown description>',
                 dbxrefs=None,
                 features=None,
                 annotations=None,
                 letter_annotations=None):

        if isinstance(seq, SeqRecord):
            self.__init__(
                seq.seq,
                seq.id,
                seq.name,
                seq.description,
                copy.deepcopy(seq.dbxrefs),
                copy.deepcopy(seq.features),
                copy.deepcopy(seq.annotations),
                copy.deepcopy(seq.letter_annotations))
        else:
            if annotations is not None:
                topology = annotations.get('topology', 'circular')
                if topology.lower() != 'circular':
                    raise ValueError("record does not describe a circular sequence")
            super(CircularRecord, self).__init__(
                seq,
                id,
                name,
                description,
                dbxrefs,
                features,
                annotations,
                letter_annotations)


    # Patch methods that are ambiguous with a non-linear sequence.

    @_ambiguous
    def __add__(self, other):
        pass

    @_ambiguous
    def __radd__(self, other):
        pass

    # Patch other methods to work as intended

    def __contains__(self, char):
        return char in str(self.seq)*2

    # def __getitem__(self, index):  # TODO
    #     return super(CircularRecord, self).__getitem__(index)


    # Additional methods

    def __lshift__(self, index):
        """Rotate the sequence to the left, preserving annotations.
        """
        return self >> (-index % len(self.seq))

    def __rshift__(self, index):
        """Rotate the sequence to the right, preserving annotations.
        """
        # FIXME: possible issue with CompoundLocation when operator
        #        is not join

        index %= len(self.seq)

        if index == 0:
            return self
        elif index < 0:
            return self << -index

        newseq = self.seq[-index:] + self.seq[:-index]
        newfeats = []
        newletan = {
            k: v[index:] + v[:index]
                for k, v in six.iteritems(self.letter_annotations)
        }

        for feature in self.features:
            if feature.location is None:
                newloc = None
            else:
                _newloc = []
                for part in (feature.location + index).parts:
                    if part.end > len(newseq) and part.start > len(newseq):
                        _newloc.append(FeatureLocation(
                            start=part.start % len(newseq),
                            end=part.end % len(newseq),
                            strand=part.strand,
                            ref=part.ref,
                            ref_db=part.ref_db))
                    elif part.start < len(newseq) and part.end > len(newseq):
                        _newloc.append(FeatureLocation(
                            start=part.start,
                            end=len(newseq),
                            strand=part.strand,
                            ref=part.ref,
                            ref_db=part.ref_db))
                        _newloc.append(FeatureLocation(
                            start=0,
                            end=part.end % len(newseq),
                            strand=part.strand,
                            ref=part.ref,
                            ref_db=part.ref_db))
                    else:
                        _newloc.append(part)
                newloc = _newloc[0] if len(_newloc) == 1 else CompoundLocation(_newloc)
            newfeats.append(SeqFeature(
                location=newloc,
                type=feature.type,
                id=feature.id,
                qualifiers=feature.qualifiers))

        return CircularRecord(
            seq=newseq,
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs,
            features=newfeats,
            annotations=self.annotations,
            letter_annotations=newletan)
