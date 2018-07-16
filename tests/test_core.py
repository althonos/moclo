# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BpiI

from moclo.core import EntryVector, Product
from moclo.record import CircularRecord


class TestEntryVector(unittest.TestCase):

    class MockEntryVector(EntryVector):
        cutter = BpiI

    class MockProduct(Product):
        cutter = BpiI

    def test_module_structure(self):
        self.assertEqual(self.MockProduct.structure(), (
            "GAAGAC"  # BpiI
            "NN"
            "(NNNN)"  # Product overhangs (start)
            "(NN*N)"  # Target
            "(NNNN)"  # Product overhangs (end)
            "NN"
            "GTCTTC"  # BpiI
        ))

    def test_vector_structure(self):
        self.assertEqual(self.MockEntryVector.structure(), (
            "N"
            "(NNNN)"
            "(NN"
            "GTCTTC"  # BpiI
            "N*"      # Placeholder
            "GAAGAC"  # BpiI
            "NN)"
            "(NNNN)"
            "N"
        ))

    def test_valid(self):
        """Assert a valid product is considered valid.
        """
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
        vector = self.MockEntryVector(SeqRecord(seqv, "vector"))
        self.assertTrue(vector.is_valid())

    def test_valid_rotated(self):
        """Assert a valid plasmid is considered valid, even after a rotation.
        """
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
        vector = self.MockEntryVector(CircularRecord(seqv, "vector"))
        self.assertTrue(vector.is_valid())
        vector = self.MockEntryVector(vector.record >> 10)
        self.assertTrue(vector.is_valid())
        vector = self.MockEntryVector(vector.record >> 10)
        self.assertTrue(vector.is_valid())

    def test_invalid(self):
        """Assert an invalid product is considered invalid.
        """
        seqv = Seq("ATG")
        vector = self.MockEntryVector(CircularRecord(seqv, "vector"))
        self.assertFalse(vector.is_valid())


    # def test_insert_linear(self):
    #     # Non-circular sequence
    #     seqp = Seq("TTTTGAAGACTTATGCAAAAAAAACGTATTGTCTTCTTTT")
    #     product = self.MockProduct(SeqRecord(seqp, "product"))
    #     seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
    #     vector = self.MockEntryVector(CircularRecord(seqv, "vector"))
    #
    #     self.assertTrue(product.is_valid())
    #     self.assertTrue(vector.is_valid())
    #
    #     self.assertEqual(vector.placeholder_sequence().seq, 'TTGTCTTCCACAGAAGACTT')
    #
    # def test_insert_plasmid(self):
    #     # Sequence with full structure in frame
    #     seqp = Seq("TTTTGAAGACTTATGCAAAAAAAACGTATTGTCTTCTTTT")
    #     product = self.MockProduct(CircularRecord(seqp, "product"))
    #     seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
    #     vector = self.MockEntryVector(CircularRecord(seqv, "vector"))
    #
    #     self.assertTrue(product.is_valid())
    #     self.assertTrue(vector.is_valid())
    #
    #     self.assertEqual(vector.placeholder_sequence().seq, 'TTGTCTTCCACAGAAGACTT')
    #     # self.assertEqual(vector.insert(product), 'CCATGCAAAAAAAACGTAGG')
    #
    #
    # def test_insert_plasmid_rotated(self):
    #     # Sequence needing to be rotated
    #     seqp = Seq("CAAAAAAAACGTATTGTCTTCTTTTTTTTGAAGACTTATG")
    #     product = self.MockProduct(CircularRecord(seqp, "product"))
    #     seqv = Seq("TCGTAGGCCATGCTTGTCTTCCACAGAAGACT")
    #     vector = self.MockEntryVector(CircularRecord(seqv, "vector"))
    #
    #     self.assertTrue(product.is_valid())
    #     self.assertTrue(vector.is_valid())
    #
    #     self.assertEqual(vector.placeholder_sequence(), 'TTGTCTTCCACAGAAGACTT')
    #     self.assertEqual(vector.insert(product), 'CCATGCAAAAAAAACGTAGG')  # FIXME: DnaRegex
