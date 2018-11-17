# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from moclo.kits import ytk
from moclo.registry.base import CombinedRegistry
from moclo.registry.ytk import YTKRegistry, PTKRegistry

from ._utils import AssemblyTestCase, PartsMetaCase, build_registries

# --- Test Suite Metaclass ---------------------------------------------------

_Meta = PartsMetaCase(
    "YTK", lambda: CombinedRegistry() << YTKRegistry() << PTKRegistry(), __name__
)


def exclude_96(item):
    return item.id == "pYTK096"


# --- Test YTK Parts ---------------------------------------------------------

TestYTKPart1 = _Meta(ytk.YTKPart1, "1")
TestYTKPart2 = _Meta(ytk.YTKPart2, "2")
TestYTKPart3 = _Meta(ytk.YTKPart3, "3")
TestYTKPart3a = _Meta(ytk.YTKPart3a, "3a")
TestYTKPart3b = _Meta(ytk.YTKPart3b, "3b")
TestYTKPart4 = _Meta(ytk.YTKPart4, "4")
TestYTKPart4a = _Meta(ytk.YTKPart4a, "4a")
TestYTKPart4b = _Meta(ytk.YTKPart4b, "4b")
TestYTKPart234 = _Meta(ytk.YTKPart234, "234")
TestYTKPart234r = _Meta(ytk.YTKPart234r, "234r", exclude_96)
TestYTKPart5 = _Meta(ytk.YTKPart5, "5")
TestYTKPart6 = _Meta(ytk.YTKPart6, "6")
TestYTKPart7 = _Meta(ytk.YTKPart7, "7")
TestYTKPart8 = _Meta(ytk.YTKPart8, "8")
TestYTKPart8a = _Meta(ytk.YTKPart8a, "8a")
TestYTKPart8b = _Meta(ytk.YTKPart8b, "8b")
TestYTKPart678 = _Meta(ytk.YTKPart678, "678")

# --- Test YTK Vectors -------------------------------------------------------

# TestYTKEntryVector = _Meta(ytk.YTKEntryVector, "EntryVector")
# TestYTKCassetteVector = _Meta(ytk.YTKCassetteVector, "CassetteVector")
# TestYTKDeviceVector = _Meta(ytk.YTKDeviceVector, "DeviceVector")

# --- Test YTK Assemblies ----------------------------------------------------

# Generate test cases based on test assembly + reference assembly


class TestYTKAssembly(AssemblyTestCase, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build_registries("ytk")

    def test_ytk_device(self):
        """Check a YTK device assembly using several cassettes.

        DNA sequences courtesy of InBio's member `Sebastián Sosa Carrillo
        <https://research.pasteur.fr/en/member/sebastian-sosa-carrillo/>`_.
        """
        result, vector, modules = self.load_data("ytk_device")
        self.assertAssembly(
            ytk.YTKDeviceVector(vector), map(ytk.YTKCassette, modules.values()), result
        )

    def test_ytk_integration_vector(self):
        """Check the YTK integration vector assembly from the reference paper.
        """
        result, vector, modules = self.load_data("ytk_integration_vector")

        modules["pYTK008.gb"] = ytk.YTKPart1(modules["pYTK008.gb"])
        modules["pYTK047.gb"] = ytk.YTKPart234r(modules["pYTK047.gb"])
        modules["pYTK073.gb"] = ytk.YTKPart5(modules["pYTK073.gb"])
        modules["pYTK074.gb"] = ytk.YTKPart6(modules["pYTK074.gb"])
        modules["pYTK086.gb"] = ytk.YTKPart7(modules["pYTK086.gb"])
        modules["pYTK092.gb"] = ytk.YTKPart8b(modules["pYTK092.gb"])

        self.assertAssembly(ytk.YTKPart8a(vector), modules.values(), result)
