# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest
import warnings

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BpiI

from moclo import errors
from moclo.record import CircularRecord
from moclo.core.vectors import AbstractVector
from moclo.core.modules import AbstractModule


class TestAssembly(unittest.TestCase):

    class MockVector(AbstractVector):
        cutter = BpiI

    class MockModule(AbstractModule):
        cutter = BpiI

    def test_invalid_vector(self):
        """Assert an error is raised on a vector with invalid overhangs.
        """
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTATGCGG")
        vector = self.MockVector(CircularRecord(seqv, "vector"))

        seqm = Seq("GAAGACTTATGCCACAATGCTTGTCTTC")
        module = self.MockModule(CircularRecord(seqm, "module"))

        # vector is invalid because both overhangs are the same
        self.assertRaises(errors.InvalidSequence, vector.assemble, module)

    def test_duplicate_modules(self):
        """Assert an error is raised when assembling with duplicate modules.
        """
        # ATGC ---- CGTA
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
        vector = self.MockVector(CircularRecord(seqv, "vector"))

        # ATGC --- CGTA
        seqm1 = Seq("GAAGACTTATGCCACACGTATTGTCTTC")
        mod1 = self.MockModule(CircularRecord(seqm1, "mod1"))

        # CGTA --- ATGC
        seqm2 = Seq("GAAGACTTATGCTATACGTATTGTCTTC")
        mod2 = self.MockModule(CircularRecord(seqm2, "mod2"))

        with self.assertRaises(errors.DuplicateModules) as ctx:
            vector.assemble(mod1, mod2)
        self.assertEqual(set(ctx.exception.duplicates), {mod1, mod2})
        self.assertEqual(ctx.exception.details, "same start overhang: 'ATGC'")
        msg = "duplicate modules: mod1, mod2 (same start overhang: 'ATGC')"
        self.assertEqual(str(ctx.exception), msg)

    def test_missing_module(self):
        """Assert an error is raised when a module is missing.
        """
        # ATGC ---- CGTA
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
        vector = self.MockVector(CircularRecord(seqv, "vector"))

        # CGTA --- ATGA
        seqm1 = Seq("GAAGACTTATGACACACGTATTGTCTTC")
        mod1 = self.MockModule(CircularRecord(seqm1, "mod1"))

        with self.assertRaises(errors.MissingModule) as ctx:
            vector.assemble(mod1)
        msg = "no module with 'ATGC' start overhang"
        self.assertEqual(str(ctx.exception), msg)


    def test_unused_modules(self):
        """Assert an error is raised on unused modules during assembly.
        """
        # ATGC ---- CGTA
        seqv = Seq("CCATGCTTGTCTTCCACAGAAGACTTCGTAGG")
        vector = self.MockVector(CircularRecord(seqv, "vector"))

        # CGTA --- ATGC
        seqm1 = Seq("GAAGACTTATGCTATACGTATTGTCTTC")
        mod1 = self.MockModule(CircularRecord(seqm1, "mod1"))

        # AAAA --- CCCC
        seqm2 = Seq("GAAGACTTAAAACACACCCCTTGTCTTC")
        mod2 = self.MockModule(CircularRecord(seqm2, "mod2"))

        with warnings.catch_warnings(record=True) as captured:
            vector.assemble(mod1, mod2)

        self.assertEqual(len(captured), 1)
        self.assertIsInstance(captured[0].message, errors.UnusedModules)
        self.assertEqual(captured[0].message.remaining, (mod2,))
        self.assertEqual(str(captured[0].message), 'unused: mod2')
