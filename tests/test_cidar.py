# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import gzip
import io
import itertools
import operator
import os
import textwrap
import sys
import unittest

import six
from Bio.Seq import Seq

from moclo.record import CircularRecord
from moclo.kits import cidar
from moclo.registry.cidar import CIDARRegistry

from ._utils import PartsMetaCase, AssemblyTestCase, build_registries


### Test CIDAR plasmids

# test suite metaclass
_Meta = PartsMetaCase("CIDAR", CIDARRegistry, __name__)
exclude_dva = lambda item: item.id.startswith("DVA_")
exclude_dvk = lambda item: item.id.startswith("DVK_")

# Generate test cases for each parts
TestCIDARPromoter = _Meta(cidar.CIDARPromoter, "Promoter", exclude_dva)
TestCIDARibosomeBindingSite = _Meta(cidar.CIDARRibosomeBindingSite, "RBS", exclude_dva)
TestCIDARCodingSequence = _Meta(cidar.CIDARCodingSequence, "CDS", exclude_dva)
TestCIDARTerminator = _Meta(cidar.CIDARTerminator, "Terminator", exclude_dva)
# ERRORS BELOW CAUSED BY 3 BsaI SITES IN R0063_AB SEQUENCES
TestCIDARPromoter.test_R0063_AB_is_Promoter = unittest.expectedFailure(
    TestCIDARPromoter.test_R0063_AB_is_Promoter
)

# Generate test cases for vectors
TestCIDAREntryVector = _Meta(cidar.CIDAREntryVector, "EntryVector", exclude_dvk)
TestCIDARCassetteVector = _Meta(
    cidar.CIDARCassetteVector, "CassetteVector", exclude_dva
)

# Generate test cases based on test assemblies
class TestCIDARAssembly(AssemblyTestCase):
    @classmethod
    def setUpClass(cls):
        build_registries("cidar")
        cls.registry = CIDARRegistry()

    def test_pJ02B2Rm_EF(self):
        expected = self.registry["pJ02B2Rm_EF"].entity.record
        vector = self.registry["DVK_EF"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_EB", "BCD2_BC", "E1010m_CD", "B0015_DF")
        ]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2Rm_AE(self):
        expected = self.registry["pJ02B2Rm_AE"].entity.record
        vector = self.registry["DVK_AE"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_AB", "BCD2_BC", "E1010m_CD", "B0015_DE")
        ]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2RmA_EF(self):
        expected = self.registry["pJ02B2Rm_EF(A)"].entity.record
        vector = self.registry["DVA_EF"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_EB", "BCD2_BC", "E1010m_CD", "B0015_DF")
        ]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2RmA_AE(self):
        expected = self.registry["pJ02B2Rm_AE(A)"].entity.record
        vector = self.registry["DVA_AE"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_AB", "BCD2_BC", "E1010m_CD", "B0015_DE")
        ]
        self.assertAssembly(vector, mods, expected)

    # ERRORS BELOW CAUSED BY NON-CONSISTENT GFP SEQUENCE BETWEEN:
    # - E0040 iGEM part sequence
    # - E0040m plasmid
    # - pJ02B2Gm* plasmid

    @unittest.expectedFailure
    def test_pJ02B2Gm_AE(self):
        expected = self.registry["pJ02B2Gm_AE"].entity.record
        vector = self.registry["DVK_AE"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_AB", "BCD2_BC", "E0040m_CD", "B0015_DE")
        ]
        self.assertAssembly(vector, mods, expected)

    @unittest.expectedFailure
    def test_pJ02B2Gm_EF(self):
        expected = self.registry["pJ02B2Gm_EF"].entity.record
        vector = self.registry["DVK_EF"].entity
        mods = [
            self.registry[x].entity
            for x in ("J23102_EB", "BCD2_BC", "E0040m_CD", "B0015_DF")
        ]
        self.assertAssembly(vector, mods, expected)

    @unittest.expectedFailure
    def test_pJ02B2RmGm_AF(self):
        expected = self.registry["pJ02B2Rm:Gm_AF"].entity.record
        vector = self.registry["DVA_AF"].entity
        mods = [self.registry[x].entity for x in ("pJ02B2Rm_AE", "pJ02B2Gm_EF")]
        from Bio.SeqIO import write

        write(expected, "/tmp/expected.gb", "gb")
        actual = vector.assemble(*mods)
        bbp = next(
            f
            for f in actual.features
            if "BioBrick prefix" in f.qualifiers.get("label", [])
        )
        write(actual << bbp.location.start - 1, "/tmp/actual.gb", "gb")
        self.assertAssembly(vector, mods, expected)
