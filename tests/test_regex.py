# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from moclo.regex import DNARegex
from moclo.record import CircularRecord


class TestDNARegex(unittest.TestCase):
    def test_type(self):
        dr = DNARegex("NN")
        self.assertRaises(TypeError, dr.search, "ATGC")

    def test_sequence_match(self):
        dr = DNARegex("AA(NN)")
        match = dr.search(Seq("ATGCAAGCAATA"))
        self.assertIsNotNone(match)
        self.assertEqual(match.start(), 4)
        self.assertEqual(match.end(), 8)
        self.assertEqual(match.span(1), (6, 8))
        self.assertEqual(match.group(1), Seq("GC"))

    def test_plasmid_match(self):
        dr = DNARegex("AA(NN)")
        match = dr.search(Seq("ATGCAGCATA"), linear=False)
        self.assertIsNotNone(match)
        self.assertEqual(match.span(0), (9, 13))
        self.assertEqual(match.span(1), (11, 13))
        self.assertEqual(match.group(1), Seq("TG"))
