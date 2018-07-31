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

from Bio.Restriction import BsaI, BbsI

from ..core import parts, modules, vectors

__author__ = 'Martin Larralde <martin.larralde@ens-paris-saclay.fr>'
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

    cutter = BbsI

    @staticmethod
    def structure():  # noqa: D105
        return (
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

    cutter = BsaI

    @staticmethod
    def structure():  # noqa: D105
        return (
            'GAAGAC'  # BbsI
            'NN'
            '(NNNN)'  # Downstream overhang
            '(N'
            'GAGACC'  # BsaI
            'N*'      # Placeholder sequence
            'GGTCTC'  # BsaI
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

    cutter = BbsI

    @staticmethod
    def structure():  # noqa: D105
        return (
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

    cutter = BbsI


class CIDAREntry(modules.Entry):
    """A CIDAR MoClo entry.
    """

    cutter = BsaI


class CIDARCassette(modules.Cassette):
    """A CIDAR MoClo cassette.
    """

    cutter = BbsI


class CIDARDevice(modules.Device):
    """A CIDAR MoClo device.
    """

    cutter = BsaI


### PARTS

class CIDARPart(parts.AbstractPart):
    """A CIDAR MoClo standard part.

    A part is a plasmid with standardized flanking overhang sequences
    that allows immediate type recognition.
    """

    cutter = BsaI
    signature = NotImplemented


class CIDARPromoter(CIDARPart, CIDAREntry):
    """A CIDAR Promoter part.

    .. image:: promoter.svg
       :align: center

    Parts of this type contain contain a promoter. The upstream overhangs can
    be changed to amend the order of assembly of a circuit from different
    cassettes.

    Note:
        The CIDAR toolkit parts provide 4 different upstream overhangs:
        *GGAG*, *GCTT*, *CGCT*, and *TGCC*. These are not enforced in this
        module, and any upstream sequence will be accepted. The downstream
        sequence however is always *TACT*.

    """

    # FIXME: enforce official upstream overhangs or not ?
    signature = ('NNNN', 'TACT')


class CIDARRibosomeBindingSite(CIDARPart, CIDAREntry):
    """A CIDAR ribosome binding site.

    .. image:: rbs.svg
       :align: center

    Parts of this type contain a ribosome binding site (RBS). The downstream
    overhang doubles as the start codon for the subsequent coding sequence.

    """

    signature = ('TACT', 'AATG')


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

    signature = ('AATG', 'AGGT')


class CIDARTerminator(CIDARPart, CIDAREntry):
    """A CIDAR terminator.

    .. image:: terminator.svg
       :align: center

    Parts of this type contain a terminator. The upstream overhang is always
    the same for the terminator to directly follow the coding sequence, but
    the downstream overhang can vary to specify an order for a following
    multigenic assembly within a device.

    Note:
        The CIDAR toolkit parts provide 4 different downstream overhangs:
        *GCTT*, *CGCT*, *TGCC*, and *ACTA*. These are not enforced in this
        module, and any downstream sequence will be accepted. The upstream
        sequence however is always *AGGT*.

    """

    signature = ('AGGT', 'NNNN')
