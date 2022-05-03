# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from moclo.kits import moclo
from moclo.registry.plant import PlantRegistry

from ._utils import AssemblyTestCase, PartsMetaCase, build_registries

# --- Test Suite Metaclass ---------------------------------------------------

_Meta = PartsMetaCase("Plant", PlantRegistry, __name__)


# --- Test Plant Parts -------------------------------------------------------

TestMoCloPro = _Meta(moclo.MoCloPro, "Pro")
TestMoClo5U = _Meta(moclo.MoClo5U, "5U")
TestMoClo5Uf = _Meta(moclo.MoClo5Uf, "5Uf")
TestMoCloNTag = _Meta(moclo.MoCloNTag, "NTag")
TestMoCloPro5U = _Meta(moclo.MoCloPro5U, "Pro5U")
TestMoCloPro5Uf = _Meta(moclo.MoCloPro5Uf, "Pro5Uf")
TestMoCloCDS1 = _Meta(moclo.MoCloCDS1, "CDS1")
TestMoCloCDS1ns = _Meta(moclo.MoCloCDS1ns, "CDS1ns")
TestMoCloSP = _Meta(moclo.MoCloSP, "SP")
TestMoCloCDS2 = _Meta(moclo.MoCloCDS2, "CDS2")
TestMoCloCDS2ns = _Meta(moclo.MoCloCDS2ns, "CDS2ns")
TestMoCloCTag = _Meta(moclo.MoCloCTag, "CTag")
TestMoClo3U = _Meta(moclo.MoClo3U, "MoClo3U")
TestMoCloTer = _Meta(moclo.MoCloTer, "MoCloTer")
TestMoClo3UTer = _Meta(moclo.MoClo3UTer, "MoClo3UTer")
TestMoCloGene = _Meta(moclo.MoCloGene, "MoCloGene")
