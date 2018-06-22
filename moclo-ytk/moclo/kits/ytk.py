# coding: utf-8

import abc

import six

from ..base import modules, vectors


### VECTORS

class YTKEntryVector(vectors.EntryVector):
    """A MoClo Yeast ToolKit entry vector.

    Any plasmid with two BsmBI restriction sites can be used to create a YTK
    entry, although the toolkit-provided entry vector (*pYTK001*) is probably
    the most appropriate plasmid to use.

    Warning:
        To the contrary of the usual MoClo entry vectors described in the
        *Weber et al.* paper, the YTK entry vectors do not provide another
        BsaI restriction site enclosing the placeholder sequence. As such,
        YTK Level -1 modules must embed the a BsaI binding site.

    See Also:
        `YTKProduct`.
    """
    _structure = (
        "(NNNN)" # Vector overhang (start)
        "(NN"
        "GAGACG" # BsmBI
        "N*?"    # Placeholder sequence
        "CGTCTC" # BsmBI
        "NN)"
        "(NNNN)" # Vector overhang (end)
    )


class YTKCassetteVector(vectors.CassetteVector):
    """A MoClo Yeast ToolKit cassette vector.
    """
    _structure = (
        "(NNNN)" # Restriction site (start)
        "(N"
        "GAGACC" # BsaI
        "N*?"    # Placeholder sequence
        "GGTCTC" # BsaI
        "N)"
        "(NNNN)" # Restriction site (end)
    )



### MODULES

class YTKProduct(modules.Product):  # FIXME ?
    _structure = (
        "CGTCTC"  # BsmBI
        "N"
        "(NNGG)"  # Product overhangs (start)
         "(TCTC"  # BsaI (first 2 nucleotides in the overhang)
        "N"
        "NNNN"   # Type specific overhang (start)
        "N*?"    # Template
        "NNNN"   # Type specific overhang (end)
        "N"
         "GAGA)" # BsaI (last 2 nucleotides in the overhang)
        "(CCNN)" # Entry overhangs (end)
        "N"
        "GAGACG" # BsmBI
    )


class YTKEntry(modules.Entry):
    _structure = (
        "GGTCTC"  # BsaI
        "N"
        "(NNNN)"  # Type specific overhang (start)
        "(N*?)"   # Target sequence
        "(NNNN)"  # Type specific overhang (end)
        "N"
        "GAGACC"  # BsaI
    )


class YTKCassette(modules.Cassette):
    _structure = (
        "CGTCTC"  # BsmBI
        "N"
        "(NNNN)"  # Cassette overhang (start)
        "(N*?)"   # Target sequence
        "(NNNN)"  # Type specific overhang (end)
        "N"
        "GAGACG"  # BsmBI
    )



### PARTS

_ent = 'GGTCTCN({start})(N*?)({end})NGAGACC'
_vec = '({end})(NGAGACCN*?GGTCTCN)({start})'

@six.add_metaclass(abc.ABCMeta)
class YTKPart(object):
    """A Yeast MoClo ToolKit standard part.

    A part is a plasmid with standardized flanking overhangs sequences
    that allow for type recognition.

    Note:
        The base MoClo framework defines entries and vectors for each level
        of assembly, but the Yeast ToolKit authors instead identify their
        sequences as "parts" of different, strongly-defined types. In order
        to follow the base MoClo frameworks, parts of Type 8, 8a, and 678 are
        handled as cassette vectors, while other parts are handled as entries.

    References:
        1. Lee, M. E., DeLoache, W. C., Cervantes, B., & Dueber, J. E. (2015).
           *A Highly Characterized Yeast Toolkit for Modular, Multipart
           Assembly.* ACS Synthetic Biology, 4(9), 975–986.
           https://doi.org/10.1021/sb500366v
    """

    @property
    @abc.abstractmethod
    def _structure(self):
        return NotImplemented


class YTKPart1(YTKPart, YTKEntry):
    _structure = _ent.format(start='CCCT', end='AACG')


class YTKPart2(YTKPart, YTKEntry):
    _structure = _ent.format(start='AACG', end='TATG')


class YTKPart3(YTKPart, YTKEntry):
    _structure = _ent.format(start='TATG', end='ATCC')


class YTKPart3a(YTKPart, YTKEntry):
    _structure = _ent.format(start='TATG', end='TTCT')


class YTKPart3b(YTKPart, YTKEntry):
    _structure = _ent.format(start='TTCT', end='ATCC')


class YTKPart4(YTKPart, YTKEntry):
    _structure = _ent.format(start='ATCC', end='GCTG')


class YTKPart4a(YTKPart, YTKEntry):
    _structure = _ent.format(start='ATCC', end='TGGC')


class YTKPart4b(YTKPart, YTKEntry):
    _structure = _ent.format(start='TGGC', end='GCTG')


class YTKPart234(YTKPart, YTKEntry):
    _structure = _ent.format(start='AACG', end='GCTG')


class YTKPart234r(YTKPart, YTKEntry):
    _structure = _vec.format(start='GCTG', end='AACG')

class YTKPart5(YTKPart, YTKEntry):
    _structure = _ent.format(start='GCTG', end='TACA')


class YTKPart6(YTKPart, YTKEntry):
    _structure = _ent.format(start='TACA', end='GAGT')


class YTKPart7(YTKPart, YTKEntry):
    _structure = _ent.format(start='GAGT', end='CCGA')


class YTKPart8(YTKPart, YTKCassetteVector):
    _structure = _vec.format(start='CCGA', end='CCCT')


class YTKPart8a(YTKPart, YTKCassetteVector):
    _structure = _vec.format(start='CCGA', end='CAAT')


class YTKPart8b(YTKPart, YTKEntry):
    _structure = _ent.format(start='CAAT', end='CCCT')


class YTKPart678(YTKPart, YTKCassetteVector):
    _structure = _vec.format(start='TACA', end='CCCT')


del _vec
del _ent
