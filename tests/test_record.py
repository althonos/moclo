# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from moclo.record import CircularRecord


class TestCircularRecord(unittest.TestCase):

    def test_init(self):
        """Assert a `CircularRecord` can be created from a `SeqRecord`.
        """
        sr = SeqRecord(seq=Seq("ATGCATGCATGC"), id="test_init")
        cr = CircularRecord(sr)

        self.assertIsInstance(cr, CircularRecord)
        self.assertNotIsInstance(cr.seq, SeqRecord)
        self.assertEqual(cr.seq, sr.seq)
        self.assertEqual(cr.id, sr.id)

    def test_shift_seq(self):
        """Assert a `CircularRecord` shifts its sequence as intended.
        """
        cr = CircularRecord(seq=Seq("ATGCATGCATGC"), id="test_shift_seq")

        self.assertEqual((cr >> 2).seq, Seq("GCATGCATGCAT"))
        self.assertEqual((cr >> 27).seq, Seq("TGCATGCATGCA"))
        self.assertEqual((cr >> len(cr)).seq, cr.seq)
        self.assertEqual((cr >> 0).seq, cr.seq)
        self.assertEqual((cr >> -1).seq, (cr << 1).seq)

        self.assertEqual((cr << 1).seq, "TGCATGCATGCA")
        self.assertEqual((cr << 14).seq, "GCATGCATGCAT")
        self.assertEqual((cr << 0).seq, cr.seq)
        self.assertEqual((cr << len(cr)).seq, cr.seq)
        self.assertEqual((cr << -5).seq, (cr >> 5).seq)

    def test_add(self):
        """Assert adding to a `CircularRecord` raises a type error.
        """
        cr = CircularRecord(seq=Seq("ATGCATGCATGC"), id="test_shift_seq")
        with self.assertRaises(TypeError):
            _ = cr + cr

    def test_radd(self):
        """Assert right-adding to a `CircularRecord` raises a type error.
        """
        cr = CircularRecord(seq=Seq("ATGCATGCATGC"), id="test_shift_seq")
        with self.assertRaises(TypeError):
            cr += cr
