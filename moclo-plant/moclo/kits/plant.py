# coding: utf-8
"""An implementation of the Plant Parts Kit for the Python MoClo library.

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
from .moclo import MoCloPart, MoCloEntry

__author__ = "Martin Larralde <martin.larralde@ens-paris-saclay.fr>"
__version__ = (
    __import__("pkg_resources")
    .resource_string(__name__, "plant.version")
    .strip()
    .decode("ascii")
)



class Plant5U(MoCloPart, MoCloEntry):
    """A Plant Part 5' UTR part.
    """

    signature = ("TACT", "CCAT")


class Plant5Uf(MoCloPart, MoCloEntry):
    """A Plant Part 5' UTR part to use with a fusion tag.
    """

    signature = ("TACT", "AATG")


class PlantPro5Uf(MoCloPart, MoCloEntry):
    """A Plant Part containing promoter and a 5'UTR region for fusion tag.
    """

    signature = ("GGAG", "CCAT")


class PlantPro5U(MoCloPart, MoCloEntry):
    """A Plant Part containing a promoter and 5'UTR region.
    """

    signature = ("GGAG", "AATG")


class PlantNSignal(MoCloPart, MoCloEntry):
    signature = ("CCAT", "AATG")


class PlantFullCDS(MoCloPart, MoCloEntry):
    signature = ("AATG", "GCTT")


class PlantCDS(MoCloPart, MoCloEntry):
    signature = ("AATG", "TTCG")


class PlantCDSNonStop(MoCloPart, MoCloEntry):
    signature = ("AGGT", "TTCG")


class PlantCSignal(MoCloPart, MoCloEntry):
    signature = ("TTCG", "GCTT")


class Plant3U(MoCloPart, MoCloEntry):
    signature = ("GCTT", "GGTA")


class PlantTer(MoCloPart, MoCloEntry):
    signature = ("GGTA", "CGCT")
