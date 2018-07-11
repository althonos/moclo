# coding: utf-8
"""An implementation of the EcoFlex ToolKit for the Python MoClo library.

References:
    1. `Moore, S. J., Lai, H.-E., Kelwick, R. J. R., Chee, S. M., Bell, D. J.,
       Polizzi, K. M., Freemont, P. S. (2016).
       EcoFlex: A Multifunctional MoClo Kit for E. coli Synthetic Biology.
       ACS Synthetic Biology, 5(10), 1059–1069.
       <https://doi.org/10.1021/acssynbio.6b00031>`_
    2. `Weber, E., Engler, C., Gruetzner, R., Werner, S., Marillonnet, S. (2011).
       A Modular Cloning System for Standardized Assembly of Multigene Constructs.
       PLOS ONE, 6(2), e16765.
       <https://doi.org/10.1371/journal.pone.0016765>`_

"""

import abc

import six

from ..base import modules, vectors

__author__ = 'Martin Larralde'
__author_email__ = 'martin.larralde@ens-paris-saclay.fr'
__version__ = (
    __import__('pkg_resources')
        .resource_string(__name__, 'ecoflex.version')
        .strip()
        .decode('ascii')
)


### VECTORS

class EcoFlexEntryVector(vectors.EntryVector):
    """An EcoFlex MoClo entry vector.
    """

    _structure = (
        '(NNNN)'  # Vector overhang (start)
        '(N'
        'GAGACG'  # BsmBI
        'N*?'     # Placeholder sequence
        'CGTCTC'  # BsmBI
        'N)'
        '(NNNN)'  # Vector overhang (end)
    )


class EcoFlexCassetteVector(vectors.CassetteVector):
    """An EcoFlex MoClo cassette vector.
    """

    _structure = (
        'CGTCTC'  # BsmBI
        'N'
        'NNNN'    # Cassette overhang (start)
        '(NNNN)'  # Vector overhang (start)
        '(N'
        'GAGACC'  # BsaI
        'N*?'     # Placeholder sequence
        'GGTCTC'  # BsaI
        'N)'
        '(NNNN)'  # Vector overhang (end)
        'NNNN'    # Cassette overhang (end)
        'N'
        'GAGACG'  # BsmBI
    )


class EcoFlexDeviceVector(vectors.DeviceVector):
    """An EcoFlex MoClo device vector.
    """

    _structure = (
        'GGTCTC'  # BsaI
        'NNNN'    # Device overhang (start)
        '(NNNN)'  # Vector overhang (start)
        '(N'
        'GAGACG'  # BsmBI
        'N*'      # Placeholder sequence
        'CGTCTC'  # BsmBI
        'N)'
        '(NNNN)'  # Vector overhang (end)
        'NNNN'    # Device overhang (end)
        'GAGACC'  # BsaI
    )


### MODULES

class EcoFlexProduct(modules.Product):
    """An EcoFlex MoClo product.
    """

    _structure = (
        'GGTCTC'  # BsaI
        'N'
        '(NNNN)'  # Type specific overhang (start)
        '(N*?)'   # Target sequence
        '(NNNN)'  # Type specific overhang (end)
        'N'
        'GAGACC'  # BsaI
    )


class EcoFlexEntry(modules.Entry):
    """An EcoFlex MoClo entry.

    EcoFlex entries are stored and shared as plasmids flanked by *BsaI*
    binding sites at both ends of the target sequence.
    """

    _structure = (
        'GGTCTC'  # BsaI
        'N'
        '(NNNN)'  # Type specific overhang (start)
        '(N*?)'   # Target sequence
        '(NNNN)'  # Type specific overhang (end)
        'N'
        'GAGACC'  # BsaI
    )


class EcoFlexCassette(modules.Cassette):
    """An EcoFlex MoClo cassette.
    """

    _structure = (
        'CGTCTC'  # BsmBI
        'N'
        '(NNNN)'  # Cassette overhang (start)
        '(N*?)'   # Target sequence
        '(NNNN)'  # Type specific overhang (end)
        'N'
        'GAGACG'  # BsmBI
    )


class EcoFlexDevice(modules.Device):
    """An EcoFlex MoClo device.
    """

    _structure = (
        'GGTCTC'  # BsaI
        'N'
        '(NNNN)'  # Device overhang (start)
        '(N*?)'   # Target sequence
        '(NNNN)'  # Device overhang (end)
        'N'
        'GAGACC'  # BsaI
    )


### PARTS

_ent = 'GGTCTCN({start})(N*?)({end})NGAGACC'
_vec = '({end})(NGAGACCN*?GGTCTCN)({start})'


@six.add_metaclass(abc.ABCMeta)
class EcoFlexPart(object):
    """An EcoFlex MoClo standard part.
    """

    @property
    @abc.abstractmethod
    def _structure(self):
        return NotImplemented


class EcoFlexPromoter(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo promoter.

    .. image: promoter.svg
       :align: center

    """

    _structure = _ent.format(start='CTAT', end='GTAC')


class EcoFlexRBS(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo ribosome binding site.

    .. image:: rbs.svg
       :align: center

    Parts of this type contain a ribosome binding site (RBS). The last
    adenosine serves as the beginning of the start codon of the following CDS.
    """

    _structure = _ent.format(start='GTAC', end='CATA')


class EcoFlexTagLinker(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo tag linker.

    .. image:: linker.svg
       :align: center

    Parts of this type also contain a RBS, but they allow adding a N-terminal
    tag sequence before the CDS.
    """

    _structure = _ent.format(start='GTAC', end='TAAA')


class EcoFlexTag(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo N-terminal tag.

    .. image:: tag.svg
       :align: center

    Parts of this type typically contain tags that are added to the N-terminus
    of the translated protein, such as a *hexa histidine* or a *Strep(II)* tag.
    """

    _structure = _ent.format(start='TAAA', end='CATA')


class EcoFlexCodingSequence(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo coding sequence.

    .. image:: cds.svg
       :align: center

    Parts of this type contain a coding sequence (CDS), with the start codon
    beginning on the upstream overhang.

    Caution:
        Although the start codon is located on the upstream overhang, a STOP
        codon is expected to be found within this part target sequence before
        the downstream overhang.
    """

    _structure = _ent.format(start='CATA', end='TCGA')


class EcoFlexTerminator(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex MoClo terminator.

    .. image:: terminator.svg
       :align: center

    """

    _structure = _ent.format(start='TCGA', end='TGTT')


class EcoFlexPromoterRBS(EcoFlexPart, EcoFlexEntry):
    """An EcoFlex Moclo promoter followed by an RBS.

    .. image:: promoter-rbs.svg
       :align: center

    These parts contain a promoter (official parts use the T7 consensus)
    followed by a ribosome binding site, and possibly a proteic tag.
    """

    _structure = _ent.format(start='CTAT', end='CATA')
