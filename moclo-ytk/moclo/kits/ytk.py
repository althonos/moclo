# coding: utf-8
"""An implementation of the Yeast ToolKit for the Python MoClo library.

References:
    1. `Lee, M. E., DeLoache, W. C., Cervantes, B., Dueber, J. E. (2015).
       A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly.
       ACS Synthetic Biology, 4(9), 975–986. <https://doi.org/10.1021/sb500366v>`_

    2. `Weber, E., Engler, C., Gruetzner, R., Werner, S., Marillonnet, S. (2011).
       A Modular Cloning System for Standardized Assembly of Multigene Constructs.
       PLOS ONE, 6(2), e16765. <https://doi.org/10.1371/journal.pone.0016765>`_

"""

import abc

import six

from ..base import modules, vectors


### VECTORS

class YTKEntryVector(vectors.EntryVector):
    """A MoClo Yeast ToolKit entry vector.

    Any plasmid with two BsmBI restriction sites can be used to create a YTK
    entry, although the toolkit-provided entry vector (*pYTK001*) is probably
    the most appropriate plasmid to use.

    Caution:
        To the contrary of the usual MoClo entry vectors described in the
        *Weber et al.* paper, the YTK entry vectors do not provide another
        BsaI restriction site enclosing the placeholder sequence. As such,
        YTK Level -1 modules must embed the a BsaI binding site.

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

    The YTK provides a canonical integration plasmid, preassembled from
    several other parts, that can be used as a cassette vector for an
    assembly of parts 2, 3 and 4. Part 8, 8a and 678 are also considered
    to be cassette vectors.

    References:
        *Lee et al.*, Figure 2.

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
    """A MoClo Yeast ToolKit cassette vector.

    As the YTK entry vector does not contain the required *BsaI* restriction
    site, the site must be contained in the product sequence.

    Caution:
        The standard construction describe in the *Lee et al.* paper directly
        inserts the beginning of the *BsaI* recognition site inside of the two
        BsmBI overhangs at both ends of the product. Other valid constructions
        that do not proceed like so won't be considered a valid product,
        although they contain the required *BsaI* site.

    References:
        *Lee et al.*, Supplementary Figure S19.

    """

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
    """A MoClo Yeast ToolKit entry.

    YTK entries are stored and shared as plasmids flanked by *BsaI* binding
    sites at both ends of the target sequence.

    Caution:
        Although the *BsaI* binding sites is not located within the target
        sequence for almost all the standard toolkit parts, special **234r**
        parts have these sites reversed, because these parts are used to
        assemble cassette vectors and require the final construct to contain
        a *BsaI* site to allow assembly with other parts. Those parts will
        **not match the default `YTKEntry`, and must be used as `YTKPart234r`
        instances** for the assembly logic to work as expected.

    """

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
    """A MoClo Yeast ToolKit cassette.
    """

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
    that allows immediate type recognition.

    Note:
        The base MoClo framework defines entries and vectors for each level
        of assembly, but the Yeast ToolKit authors instead identify their
        sequences as "parts" of different, strongly-defined types. In order
        to follow the base MoClo frameworks, parts of Type 8, 8a, and 678 are
        handled as cassette vectors, while other parts are handled as entries.

    """

    @property
    @abc.abstractmethod
    def _structure(self):
        return NotImplemented


class YTKPart1(YTKPart, YTKEntry):
    """A YTK type 1 part (**Assembly Connector**).

    Typically, parts of this type contains non-coding and non-regulatory
    sequences that are used to direct assembly of multigene plasmids, such
    as ligation sites for other Type IIS endonucleases (e.g. *BsmBI*).

    Hint:
        | Start overhang: **CCCT**
        | End overhang:   **AACG**

    """
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
