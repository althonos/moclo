# coding: utf-8
"""An implementation of the CIDAR ToolKit for the Python MoClo library.

References:
    1. `Iverson, S. V., Haddock, T. L., Beal, J., & Densmore, D. M. (2016).
       CIDAR MoClo: Improved MoClo Assembly Standard and New E. coli Part Library
       Enable Rapid Combinatorial Design for Synthetic and Traditional Biology.
       ACS Synthetic Biology, 5(1), 99–103.
       <https://doi.org/10.1021/acssynbio.5b00124>`_
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
        .resource_string(__name__, 'cidar.version')
        .strip()
        .decode('ascii')
)

### VECTORS

class CIDAREntryVector(vectors.EntryVector):
    """A CIDAR MoClo entry vector.
    """

    _structure = (
        'GGTCTC'  # BsaI
        'N'
        '(NNNN)'  # Downstream overhang
        '(NN'
        'GTCTTC'  # BbsI
        'N*'      # Placeholder sequence
        'GAAGAC'  # BbsI
        'NN)'
        '(NNNN)'  # Upstream overhang
        'N'
        'GAGACC'  # BsaI
    )


class CIDARCassetteVector(vectors.CassetteVector):
    """A CIDAR Moclo cassette vector.

    References:
        *Iverson et al.*, Figure 1.

    """

    _structure = (
        'GAAGAC'  # BbsI
        'NN'
        '(NNNN)'  # Downstream overhang
        '(N'
        'GAGACC'  # BsaI
        'N*'      # Placeholder sequence
        'CGTCTC'  # BsaI
        'N)'
        '(NNNN)'  # Upstream overhang
        'NN'
        'GTCTTC'  # BbsI
    )


class CIDARDeviceVector(vectors.DeviceVector):
    """A CIDAR Moclo device vector.

    References:
        *Iverson et al.*, Figure 1.

    """

    _structure = (
        'GGTCTC'  # BsaI
        'N'
        '(NNNN)'  # Downstream overhang
        '(NN'
        'GTCTTC'  # BbsI
        'N*'      # Placeholder sequence
        'GAAGAC'  # BbsI
        'NN)'
        '(NNNN)'  # Upstream overhang
        'N'
        'GAGACC'  # BsaI
    )


### MODULES

class CIDARProduct(modules.Product):
    """A CIDAR MoClo product.
    """
    # TODO #


class CIDAREntry(modules.Entry):
    """A CIDAR MoClo entry.
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


class CIDARCassette(modules.Cassette):
    """A CIDAR MoClo cassette.
    """

    _structure = (
        'GAAGAC'  # BbsI
        'NN'
        '(NNNN)'  # Cassette overhang (start)
        '(N*?)'   # Target sequence
        '(NNNN)'  # Type specific overhang (end)
        'NN'
        'GTCTTC'  # BbsI
    )


class CIDARDevice(modules.Device):
    """A CIDAR MoClo device.
    """

    _structure = (
        'GGTCTC'   # BsaI
        'N'
        '(NNNN)'   # Upstream overhang
        '(N*)'     # Target sequence
        '(NNNN)'   # Downstream overhang
        'N'
        'GAGACC'   # BsaI
    )


### PARTS

_ent = 'GGTCTCN({start})(N*?)({end})NGAGACC'


@six.add_metaclass(abc.ABCMeta)
class CIDARPart(object):
    """A CIDAR MoClo standard part.

    A part is a plasmid with standardized flanking overhang sequences
    that allows immediate type recognition.
    """

    @property
    @abc.abstractmethod
    def _structure(self):
        return NotImplemented


class CIDARPromoter(CIDARPart, CIDAREntry):
    """A CIDAR Promoter part.

    .. image:: promoter.svg
       :align: center

    Parts of this type contain contain a promoter. The upstream overhangs can
    be changed to amend the order of assembly of a circuit from different
    cassettes. The upstream overhang can vary between four differents sequences,
    but the downstream overhang is unique.
    """

    # FIXME: enforce official upstream overhangs or not ?
    _structure = _ent.format(start='GGAG|GCTT|CGCT|TGCC', end='TACT')


class CIDARibosomeBindingSite(CIDARPart, CIDAREntry):
    """A CIDAR ribosome binding site.

    .. image:: rbs.svg
       :align: center

    Parts of this type contain a ribosome binding site (RBS). The downstream
    overhang doubles as the start codon for the subsequent coding sequence.

    """

    _structure = _ent.format(start='TACT', end='AATG')


class CIDARCodingSequence(CIDARPart, CIDAREntry):
    """A CIDAR coding sequence.

    .. image:: cds.svg
       :align: center

    Parts of this type contain a coding sequence, with the start codon located
    on the upstream overhang.

    Caution:
        Although the start codon is located on the upstream overhang, a STOP
        codon is expected to be found within this part target sequence before
        the downstream overhang.
    """

    _structure = _ent.format(start='AATG', end='AGGT')


class CIDARTerminator(CIDARPart, CIDAREntry):
    """A CIDAR terminator.

    .. image:: terminator.svg
       :align: center

    Parts of this type contain a terminator. The upstream overhang is always
    the same for the terminator to directly follow the coding sequence, but
    the downstream overhang can vary to specify an order for a following
    multigenic assembly within a device.
    """

    _structure = _ent.format(start='AGGT', end='GCTT|CGCT|TGCC|ACTA')


del _ent
