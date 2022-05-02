# coding: utf-8
"""An implementation of the original MoClo ToolKit for the Python MoClo library.

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


class MoCloVector(vectors.EntryVector):
    """A MoClo entry vector.

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


class MoCloCassetteVector(vectors.CassetteVector):
    """A MoClo cassette vector.

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


class MoCloSingleCassetteVector(vectors.CassetteVector):
    """An original MoClo single cassette vector.

    These plasmids allow obtaining plasmids ready for plant transformation
    containing only a single transcription unit, assembled directly from
    entries, without using an intermediate *E. coli* plasmid.

    Danger:
        The resulting plasmid will not be a valid `MoCloCassette`!

    References:
        `AddGene MoClo <https://www.addgene.org/cloning/moclo/marillonnet/>`_
        construct reference.

    """

    cutter = BsaI


class MoCloDeviceVector(vectors.DeviceVector):
    """An original MoClo device vector.

    References:
        *Weber et al.*, Figure 4A.

    """

    cutter = BpiI


# MODULES ####################################################################


class MoCloProduct(modules.Product):
    """An original MoClo product.
    """

    cutter = BpiI


class MoCloEntry(modules.Entry):
    """An original MoClo entry.
    """

    cutter = BsaI


class MoCloCassette(modules.Cassette):
    """An original MoClo cassette.
    """

    cutter = BpiI


# PARTS ######################################################################

# Abstract ###################################################################


class MoCloPart(parts.AbstractPart):
    """An original MoClo standard part.

    A part is a plasmid with standardized flanking overhang sequences
    that allows immediate type recognition.
    """

    cutter = BsaI
    signature = NotImplemented


# Level 0 ####################################################################


class MoCloPromoter(MoCloPart, MoCloEntry):
    """An original MoClo promoter part.
    """

    signature = ("GGAG", "TACT")


class MoCloUntranslatedRegion(MoCloPart, MoCloEntry):
    """An original MoClo 5' UTR part.
    """

    signature = ("TACT", "AATG")


class MoCloSignalPeptide(MoCloPart, MoCloEntry):
    """An original MoClo signal peptide part.
    """

    signature = ("AATG", "AGGT")


class MoCloCodingSequence(MoCloPart, MoCloEntry):
    """An original MoClo CDS part.
    """

    signature = ("AGGT", "GCTT")


class MoCloTerminator(MoCloPart, MoCloEntry):
    """An original MoClo terminator part.
    """

    signature = ("GCTT", "CGCT")


# Level 0 ####################################################################


class MoCloEndLinker(MoCloPart, MoCloCassette):
    """An Icon Genetic end linker part.

    References:
        *Weber et al.*, Figure 5.

    """

    signature = ("NNNN", "GGGA")


# Level M ####################################################################


class MoCloLevelMEndLinker(MoCloPart, MoCloCassette):

    signature = ("NNNN", "GGGA")

    # FIXME: add prefix BsaI
    @classmethod
    def structure(cls):
        return super(MoCloLevelMEndLinker, cls).structure()


class MoCloLevelMVector(MoCloPart, MoCloDeviceVector):

    cutter = BpiI
    signature = ("GGGA", "NNNN")

    # FIXME: add prefix BsaI
    @classmethod
    def structure(cls):
        return super(MoCloLevelMVector, cls).structure()


# Level P ####################################################################


class MoCloLevelPEndLinker(MoCloPart, MoCloEntry):  # FIXME: hierarchy ?

    cutter = BsaI
    signature = ("NNNN", "GGGA")

    # FIXME: add prefix BpiI
    @classmethod
    def structure(cls):
        return super(MoCloLevelPEndLinker, cls).structure()


class MoCloLevelPVector(MoCloPart, MoCloCassetteVector):  # FIXME: hierarchy ?

    cutter = BsaI
    signature = ("GGGA", "NNNN")

    # FIXME: add prefix BpiI
    @classmethod
    def structure(cls):
        return super(MoCloLevelPVector, cls).structure()
