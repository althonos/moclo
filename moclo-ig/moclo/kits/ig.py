# coding: utf-8
"""An implementation of the Icon Genetics ToolKit for the Python MoClo library.

References:
    1. `Weber, E., Engler, C., Gruetzner, R., Werner, S., Marillonnet, S. (2011).
       A Modular Cloning System for Standardized Assembly of Multigene Constructs.
       PLOS ONE, 6(2), e16765.
       <https://doi.org/10.1371/journal.pone.0016765>`_
    2. `Werner, S., Engler, C., Weber, E., Gruetzner, R., & Marillonnet, S. (2012).
       Fast track assembly of multigene constructs using Golden Gate cloning
       and the MoClo system. Bioengineered, 3(1), 38–43.
       <https://doi.org/10.4161/bbug.3.1.18223>`_
"""
from __future__ import absolute_import
from __future__ import unicode_literals

from Bio.Restriction import BsaI, BpiI

from ..core import parts, modules, vectors

__author__ = "Martin Larralde <martin.larralde@ens-paris-saclay.fr>"
__version__ = (
    __import__("pkg_resources")
    .resource_string(__name__, "ig.version")
    .strip()
    .decode("ascii")
)


# VECTORS ####################################################################


class IGEntryVector(vectors.EntryVector):
    """An Icon Genetics entry vector.

    References:
        *Weber et al.*, Figure 2A.

    """

    cutter = BpiI

    @classmethod
    def structure(cls):
        return (
            "GGTCTC"  # BsaI
            "N"
            "(NNNN)"  # Downstream overhang
            "(NN"
            "GTCTTC"  # BpiI
            "N*"  # Placeholder sequence
            "GAAGAC"  # BpiI
            "NN)"
            "(NNNN)"  # Upstream overhang
            "N"
            "GAGACC"  # BsaI
        )


class IGCassetteVector(vectors.CassetteVector):
    """An Icon Genetics cassette vector.

    References:
        *Weber et al.*, Figure 4A.

    """

    cutter = BsaI

    @classmethod
    def structure(cls):  # noqa: D105
        return (
            "GAAGAC"  # BpiI
            "NN"
            "NNNN"  # Cassette upstream overhang
            "(NNNN)"  # Cassette vector downstream overhang
            "(N"
            "GAGACC"  # BsaI
            "N*"  # Placeholder sequence
            "GGTCTC"  # BsaI
            "N)"
            "(NNNN)"  # Cassette vector upstream overhang
            "NNNN"  # Cassette downstream overhang
            "NN"
            "GTCTTC"  # BbsI
        )


class IGSingleCassetteVector(vectors.CassetteVector):
    """An Icon Genetics single cassette vector.

    These plasmids allow obtaining plasmids ready for plant transformation
    containing only a single transcription unit, assembled directly from
    entries, without using an intermediate *E. coli* plasmid.

    Danger:
        The resulting plasmid will not be a valid `IGCassette`!

    References:
        `AddGene MoClo <https://www.addgene.org/cloning/moclo/marillonnet/>`_
        construct reference.

    """

    cutter = BsaI


class IGDeviceVector(vectors.DeviceVector):
    """An Icon Genetics device vector.

    References:
        *Weber et al.*, Figure 4A.

    """

    cutter = BpiI


# MODULES ####################################################################


class IGProduct(modules.Product):
    """An Icon Genetics MoClo product.
    """

    cutter = BpiI


class IGEntry(modules.Entry):
    """An Icon Genetics MoClo entry.
    """

    cutter = BsaI


class IGCassette(modules.Cassette):
    """An Icon Genetics MoClo cassette.
    """

    cutter = BpiI


# PARTS ######################################################################

# Abstract ###################################################################


class IGPart(parts.AbstractPart):
    """An Icon Genetics MoClo standard part.

    A part is a plasmid with standardized flanking overhang sequences
    that allows immediate type recognition.
    """

    cutter = BsaI
    signature = NotImplemented


# Level 0 ####################################################################


class IGPromoter(IGPart, IGEntry):
    """An Icon Genetics promoter part.
    """

    signature = ("GGAG", "TACT")


class IGUntranslatedRegion(IGPart, IGEntry):
    """An Icon Genetics 5' UTR part.
    """

    signature = ("TACT", "AATG")


class IGSignalPeptide(IGPart, IGEntry):
    """An Icon Genetics signal peptide part.
    """

    signature = ("AATG", "AGGT")


class IGCodingSequence(IGPart, IGEntry):
    """An Icon Genetics CDS part.
    """

    signature = ("AGGT", "GCTT")


class IGTerminator(IGPart, IGEntry):
    """An Icon Genetics terminator part.
    """

    signature = ("GCTT", "CGCT")


# Level 0 ####################################################################


class IGEndLinker(IGPart, IGCassette):
    """An Icon Genetic end linker part.

    References:
        *Weber et al.*, Figure 5.

    """

    signature = ("NNNN", "GGGA")


# Level M ####################################################################


class IGLevelMEndLinker(IGPart, IGCassette):

    signature = ("NNNN", "GGGA")

    # FIXME: add prefix BsaI
    @classmethod
    def structure(cls):
        return super(IGLevelMEndLinker, cls).structure()


class IGLevelMVector(IGPart, IGDeviceVector):

    cutter = BpiI
    signature = ("GGGA", "NNNN")

    # FIXME: add prefix BsaI
    @classmethod
    def structure(cls):
        return super(IGLevelMVector, cls).structure()


# Level P ####################################################################


class IGLevelPEndLinker(IGPart, IGEntry):  # FIXME: hierarchy ?

    cutter = BsaI
    signature = ("NNNN", "GGGA")

    # FIXME: add prefix BpiI
    @classmethod
    def structure(cls):
        return super(IGLevelPEndLinker, cls).structure()


class IGLevelPVector(IGPart, IGCassetteVector):  # FIXME: hierarchy ?

    cutter = BsaI
    signature = ("GGGA", "NNNN")

    # FIXME: add prefix BpiI
    @classmethod
    def structure(cls):
        return super(IGLevelPVector, cls).structure()
