# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import gzip
import io
import itertools
import os
import textwrap
import unittest

import six
from Bio.Seq import Seq

from moclo.record import CircularRecord
from moclo.kits import ytk

from .utils import plasmids



### Test Yeast ToolKit plasmids

# Metaclass for test suites: create a test suite per type with a method per
# plasmid that will check if that plasmid is recognized as the test suite
# part or not rightfully. The `_TestYTK` instance acts as a single test case.
class _TestYTK(unittest.TestCase):

    _plasmids = None

    @classmethod
    def make_suite(cls, part_cls, part_name, exclude=frozenset()):
        if cls._plasmids is None:
            cls._plasmids = list(plasmids('ytk.tsv.xz'))
        case_name = str('Test{}'.format(part_cls.__name__))
        tests = {
            test.__name__: test
            for plasmid in cls._plasmids
            if plasmid[0] not in exclude
            for test in (cls(*plasmid, part_cls=part_cls, part_name=part_name),)
        }
        return type(case_name, (unittest.TestCase,), tests)

    def __init__(self, id_, type_, name, desc, seq, part_cls, part_name):
        self.id = id_
        self.type = type_
        self.name = name
        self.desc = desc
        self.seq = seq
        self.rec = CircularRecord(Seq(seq), name=name, id=id_)
        self._part_cls = part_cls
        self._part_type = part_name

    def _assert_valid(self):
        err = '{} is not a valid Type {} but should be!'
        self.assertTrue(
            self._part_cls(self.rec).is_valid(),
            err.format(self.id, self._part_type)
        )

    def _assert_invalid(self):
        err = '{} is a valid Type {} but should not be!'
        self.assertFalse(
            self._part_cls(self.rec).is_valid(),
            err.format(self.id, self._part_type)
        )

    def __call__(self):
        if self.type == self._part_type:
            self._assert_valid()
        else:
            self._assert_invalid()

    @property
    def __name__(self):
        if self.type == self._part_type:
            msg = "test_{}_is_type{}"
        else:
            msg = "test_{}_is_not_type{}"
        return str(msg.format(self.id, self._part_type))

    @property
    def __doc__(self):
        if self.name == self._part_type:
            doc = 'Check that {} ({}) is a YTK Type {} part.\n'
        else:
            doc = 'Check that {} ({}) is not a YTK Type {} part.\n'
        return doc.format(self.id, self.name, self._part_type)


# Generate test cases
TestYTKPart1 = _TestYTK.make_suite(ytk.YTKPart1, '1')
TestYTKPart2 = _TestYTK.make_suite(ytk.YTKPart2, '2')
TestYTKPart3 = _TestYTK.make_suite(ytk.YTKPart3, '3')
TestYTKPart3a = _TestYTK.make_suite(ytk.YTKPart3a, '3a')
TestYTKPart3b = _TestYTK.make_suite(ytk.YTKPart3b, '3b')
TestYTKPart4 = _TestYTK.make_suite(ytk.YTKPart4, '4')
TestYTKPart4a = _TestYTK.make_suite(ytk.YTKPart4a, '4a')
TestYTKPart4b = _TestYTK.make_suite(ytk.YTKPart4b, '4b')
TestYTKPart234 = _TestYTK.make_suite(ytk.YTKPart234, '234')
TestYTKPart234r = _TestYTK.make_suite(ytk.YTKPart234r, '234r', exclude={'pYTK096'})
TestYTKPart5 = _TestYTK.make_suite(ytk.YTKPart5, '5')
TestYTKPart6 = _TestYTK.make_suite(ytk.YTKPart6, '6')
TestYTKPart7 = _TestYTK.make_suite(ytk.YTKPart7, '7')
TestYTKPart8 = _TestYTK.make_suite(ytk.YTKPart8, '8')
TestYTKPart8a = _TestYTK.make_suite(ytk.YTKPart8a, '8a')
TestYTKPart8b = _TestYTK.make_suite(ytk.YTKPart8b, '8b')
TestYTKPart678 = _TestYTK.make_suite(ytk.YTKPart678, '678')


### Test Yeast ToolKit reference paper construct

class TestYTKConstruct(unittest.TestCase):

    def test_cassette_vector(self):
        """Check the assembled YTK integration vector has the right sequence.
        """
        records = {
            p[0]: CircularRecord(Seq(p[4]), id=p[0], name=p[0])
            for p in plasmids('ytk.tsv.xz')
        }

        part1 = ytk.YTKPart1(records['pYTK008'])
        part234 = ytk.YTKPart234r(records['pYTK047'])
        part5 = ytk.YTKPart5(records['pYTK073'])
        part6 = ytk.YTKPart6(records['pYTK074'])
        part7 = ytk.YTKPart7(records['pYTK086'])
        part8a = ytk.YTKPart8a(records['pYTK090'])
        part8b = ytk.YTKPart8b(records['pYTK092'])

        assembly = part8a.assemble(part1, part234, part5, part6, part7, part8b)
        expected = records['pYTK096']

        self.assertEqual(len(assembly), len(expected), 'lengths differ')
        self.assertIn(assembly.seq, expected.seq + expected.seq, 'sequences differ')


### Test Yeast ToolKit multigene assembly

class TestYTKDevice(unittest.TestCase):

    def test_Device_assembly(self):
        """Check a YTK multigene assembly gives the expected result.
        """
        records = {
            p[0]: CircularRecord(Seq(p[4]), id=p[0], name=p[0])
            for p in plasmids('ytk-multigene.tsv.xz')
        }

        tu1 = ytk.YTKCassette(records['mCerulean'])
        tu2 = ytk.YTKCassette(records['mNeonGreen'])
        tu3 = ytk.YTKCassette(records['mScarlet'])
        vec = ytk.YTKDeviceVector(records['multigeneVec'])

        assembly = vec.assemble(tu1, tu2, tu3)
        expected = records['mCmNGmS']

        self.assertEqual(len(assembly), len(expected), 'lengths differ')
        self.assertIn(assembly.seq, expected.seq + expected.seq, 'sequences differ')
