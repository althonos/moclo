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


def setUpModule():
    build_registries('cidar')


### Test CIDAR plasmids

# test suite metaclass
_Meta = PartsMetaCase('CIDAR', 'cidar.tsv.xz', __name__)

# Generate test cases for each parts
TestCIDARPromoter = _Meta(cidar.CIDARPromoter, 'Promoter')
TestCIDARibosomeBindingSite = _Meta(cidar.CIDARRibosomeBindingSite, 'RBS')
TestCIDARCodingSequence = _Meta(cidar.CIDARCodingSequence, 'CDS')
TestCIDARTerminator = _Meta(cidar.CIDARTerminator, 'Terminator')

# Generate test cases based on test assemblies
class TestCIDARAssembly(AssemblyTestCase):

    @classmethod
    def setUpClass(cls):
        cls.registry = CIDARRegistry()

    def test_pJ02B2RmK_EF(self):
        expected = self.registry['pJ02B2RmK_EF'].entity.record
        vector = self.registry['DVK_EF'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_EB', 'BCD2_BC', 'E1010m_CD', 'B0015_DF')]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2RmK_AE(self):
        expected = self.registry['pJ02B2RmK_AE'].entity.record
        vector = self.registry['DVK_AE'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_AB', 'BCD2_BC', 'E1010m_CD', 'B0015_DE')]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2RmA_EF(self):
        expected = self.registry['pJ02B2RmA_EF'].entity.record
        vector = self.registry['DVA_EF'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_EB', 'BCD2_BC', 'E1010m_CD', 'B0015_DF')]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2RmA_AE(self):
        expected = self.registry['pJ02B2RmA_AE'].entity.record
        vector = self.registry['DVA_AE'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_AB', 'BCD2_BC', 'E1010m_CD', 'B0015_DE')]
        self.assertAssembly(vector, mods, expected)

    @unittest.expectedFailure # bad AddGene pJ02B2GmK_AE sequence
    def test_pJ02B2GmK_AE(self):
        expected = self.registry['pJ02B2GmK_AE'].entity.record
        vector = self.registry['DVK_AE'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_AB', 'BCD2_BC', 'E0040m_CD', 'B0015_DE')]
        self.assertAssembly(vector, mods, expected)

    def test_pJ02B2GmK_EF(self):
        expected = self.registry['pJ02B2GmK_EF'].entity.record
        vector = self.registry['DVK_EF'].entity
        mods = [self.registry[x].entity
                for x in ('J23102_EB', 'BCD2_BC', 'E0040m_CD', 'B0015_DF')]
        self.assertAssembly(vector, mods, expected)

    # FIXME !
    def test_pJ02B2RmGmA_AF(self):
        expected = self.registry['pJ02B2Rm:GmA_AF'].entity.record
        vector =  self.registry['DVA_AF'].entity
        mods = [self.registry[x].entity for x in ('pJ02B2RmK_AE', 'pJ02B2GmK_EF')]
        self.assertAssembly(vector, mods, expected)
