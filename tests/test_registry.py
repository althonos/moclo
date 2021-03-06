# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

import fs
import six
from Bio.SeqIO import write

from moclo.kits import ytk, cidar
from moclo.record import CircularRecord
from moclo.registry import base
from moclo.registry.ytk import YTKRegistry, PTKRegistry
from moclo.registry.cidar import CIDARRegistry

from ._utils import build_registries


class TestEmbeddedRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build_registries("ytk")
        build_registries("cidar")

    def test_circular_record(self):
        r = YTKRegistry()
        self.assertIsInstance(r["pYTK002"].entity.record, CircularRecord)

    def test_ytk_registry(self):
        r = YTKRegistry()

        # typecheck
        self.assertIsInstance(r["pYTK002"].entity, ytk.YTKPart1)
        self.assertIsInstance(r["pYTK047"].entity, ytk.YTKPart234r)
        self.assertIsInstance(r["pYTK096"].entity, ytk.YTKCassetteVector)

        # resistance check
        self.assertEqual(r["pYTK005"].resistance, "Chloramphenicol")
        self.assertEqual(r["pYTK085"].resistance, "Spectinomycin")

        # missing key check
        self.assertRaises(KeyError, r.__getitem__, "pYTK200")

    def test_ptk_registry(self):
        r = PTKRegistry()

        # typecheck
        self.assertIsInstance(r["pPTK004"].entity, ytk.YTKPart2)
        self.assertIsInstance(r["pPTK019"].entity, ytk.YTKPart4)
        self.assertIsInstance(r["pPTK013"].entity, ytk.YTKPart3a)

        # resistance check
        self.assertEqual(r["pPTK005"].resistance, "Chloramphenicol")

    def test_cidar_registry(self):
        r = CIDARRegistry()

        # typecheck
        self.assertIsInstance(r["C0062_CD"].entity, cidar.CIDARCodingSequence)
        self.assertIsInstance(r["BCD8_BC"].entity, cidar.CIDARRibosomeBindingSite)
        self.assertIsInstance(r["DVA_GB"].entity, cidar.CIDAREntryVector)
        self.assertIsInstance(r["DVK_GH"].entity, cidar.CIDARCassetteVector)

        # resistance check
        self.assertEqual(r["DVA_BC"].resistance, "Ampicillin")
        self.assertEqual(r["DVK_EF"].resistance, "Kanamycin")


class TestFilesystemRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build_registries("ytk")
        r = YTKRegistry()

        cls.memfs = fs.open_fs("mem://")

        buff = six.StringIO()
        write([r["pYTK002"].entity.record], buff, "genbank")
        with cls.memfs.open("pYTK002.gb", "w") as f:
            f.write(buff.getvalue())

        buff = six.StringIO()
        write([r["pYTK038"].entity.record], buff, "genbank")
        with cls.memfs.open("pYTK038.gb", "w") as f:
            f.write(buff.getvalue())

    @classmethod
    def tearDownClass(cls):
        cls.memfs.close()

    def test_mem_registry(self):
        r = base.FilesystemRegistry(self.memfs, ytk.YTKPart)
        self.assertIn("pYTK002", r)
        self.assertNotIn("pYTK003", r)
        self.assertIsInstance(r["pYTK002"].entity, ytk.YTKPart1)
        self.assertIsInstance(r["pYTK038"].entity, ytk.YTKPart3a)
        self.assertEqual(r["pYTK002"].resistance, "Chloramphenicol")

    def test_invalid_base(self):
        r = base.FilesystemRegistry(self.memfs, ytk.YTKPart8)
        self.assertRaises(RuntimeError, r.__getitem__, "pYTK002")

    def test_len(self):
        r = base.FilesystemRegistry(self.memfs, ytk.YTKPart)
        self.assertEqual(len(r), len(list(r)))
        self.assertEqual(len(r), len(self.memfs.listdir("/")))

    def test_iter(self):
        r = base.FilesystemRegistry(self.memfs, ytk.YTKPart)
        self.assertTrue(all(x in r for x in r))
        self.assertEqual(sorted(r), ["pYTK002", "pYTK038"])

    def test_circular_record(self):
        r = base.FilesystemRegistry(self.memfs, ytk.YTKPart)
        self.assertIsInstance(r["pYTK002"].entity.record, CircularRecord)


class TestCombinedRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build_registries("ytk")
        cls.ptk = PTKRegistry()
        cls.ytk = YTKRegistry()
        cls.registry = base.CombinedRegistry()
        cls.registry << cls.ytk << cls.ptk

    def test_contains(self):
        self.assertIn("pYTK001", self.registry)
        self.assertIn("pPTK001", self.registry)
        self.assertNotIn("B0030_AF", self.registry)

    def test_length(self):
        self.assertEqual(len(self.registry), len(YTKRegistry()) + len(PTKRegistry()))

    def test_getitem(self):
        self.assertIs(self.registry["pYTK001"], self.ytk["pYTK001"])
        self.assertIs(self.registry["pPTK003"], self.ptk["pPTK003"])
        self.assertRaises(KeyError, self.registry.__getitem__, "B0030_AF")

    def test_iter(self):
        self.assertEqual(sorted(self.registry), sorted(self.ptk) + sorted(self.ytk))
