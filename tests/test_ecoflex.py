# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from moclo.kits import ecoflex
from moclo.registry.ecoflex import EcoFlexRegistry

from ._utils import PartsMetaCase


# --- Test Suite Metaclass ---------------------------------------------------

_Meta = PartsMetaCase("EcoFlex", EcoFlexRegistry, __name__)

exclude_tu2 = lambda item: item.id.startswith("pTU2")
exclude_prom = lambda item: item.id.startswith(("pTU2-A", "pTU2-a", "pTU2-2"))
exclude_tu2d = lambda item: item.id.startswith("pTU2-D")
exclude_cassette = lambda item: item.id.startswith(("pTU1", "pTU3"))


# --- Test EcoFlex Parts -----------------------------------------------------

TestEcoFlexPromoter = _Meta(ecoflex.EcoFlexPromoter, "Promoter", exclude_prom)
TestEcoFlexRBS = _Meta(ecoflex.EcoFlexRBS, "RBS")
TestEcoFlexTag = _Meta(ecoflex.EcoFlexTag, "Tag")
TestEcoFlexTerminator = _Meta(ecoflex.EcoFlexTerminator, "Terminator", exclude_tu2d)
TestEcoFlexCodingSequence = _Meta(ecoflex.EcoFlexCodingSequence, "CDS")
TestEcoFlexPromoterRBS = _Meta(ecoflex.EcoFlexPromoterRBS, "PromoterRBS")
TestEcoFlexTagLinker = _Meta(ecoflex.EcoFlexTagLinker, "TagLinker")


# --- Test EcoFlex Vectors ---------------------------------------------------

TestEcoFlexCassetteVector = _Meta(ecoflex.EcoFlexCassetteVector, "CassetteVector", exclude_tu2)
TestEcoFlexDeviceVector = _Meta(ecoflex.EcoFlexDeviceVector, "DeviceVector", exclude_cassette)
