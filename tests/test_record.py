# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

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

    def test_linear_error(self):
        """Assert a `CircularRecord` does not accept a linear record.
        """
        an = {"topology": "linear"}
        sr = SeqRecord(seq=Seq("ATGC"), id="linear", annotations=an)
        self.assertRaises(ValueError, CircularRecord, sr)

    def test_contains(self):
        """Assert `_ in CircularRecord` works as expected.
        """
        sr = SeqRecord(seq=Seq("ATGC"), id="test_init")
        cr = CircularRecord(sr)
        self.assertIn("ATGC", cr)
        self.assertIn("GCAT", cr)
        self.assertNotIn("ATGCAT", cr)

    def test_shift_features(self):
        """Assert a `CircularRecord` shifts its features as intended.
        """
        ft = [
            SeqFeature(
                FeatureLocation(ExactPosition(0), ExactPosition(2), strand=+1),
                type="promoter",
            ),
            SeqFeature(
                FeatureLocation(ExactPosition(2), ExactPosition(4), strand=+1),
                type="promoter",
            ),
        ]
        sr = SeqRecord(seq=Seq("ATGC"), id="feats", features=ft)
        cr = CircularRecord(sr)

        cr_1 = cr >> 1
        self.assertEqual(
            cr_1.features[0].location,
            FeatureLocation(ExactPosition(1), ExactPosition(3), strand=+1),
        )
        self.assertEqual(
            cr_1.features[1].location,
            FeatureLocation(ExactPosition(3), ExactPosition(5), strand=+1),
        )

        cr_2 = cr_1 >> 1
        self.assertEqual(
            cr_2.features[0].location,
            FeatureLocation(ExactPosition(2), ExactPosition(4), strand=+1),
        )
        self.assertEqual(
            cr_2.features[1].location,
            FeatureLocation(ExactPosition(0), ExactPosition(2), strand=+1),
        )

        cr_3 = cr_2 >> 1
        self.assertEqual(
            cr_3.features[0].location,
            FeatureLocation(ExactPosition(3), ExactPosition(5), strand=+1),
        )
        self.assertEqual(
            cr_3.features[1].location,
            FeatureLocation(ExactPosition(1), ExactPosition(3), strand=+1),
        )

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

        self.assertEqual((cr << -3).seq, (cr >> 3).seq)
        self.assertEqual((cr >> -3).seq, (cr << 3).seq)

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
