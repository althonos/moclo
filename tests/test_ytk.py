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
from parameterized import parameterized

from moclo.record import CircularRecord
from moclo.kits import ytk

try:
    import lzma
except ImportError:
    from backports import lzma



### Plasmids from CSV loader

def plasmids(name):
    plasmids_file = os.path.join(__file__, '..', 'data', name)
    with io.TextIOWrapper(lzma.open(os.path.abspath(plasmids_file))) as f:
        for line in f:
            if not line.startswith('Plasmid Name'):
                yield line.strip().split('\t')


### Test Yeast ToolKit plasmids

# Metaclass for test cases: create a test case per type with a method per
# plasmid that will check if that plasmid is recognized as the test case
# part or not rightfully.
def test_ytk_part(_cls, _name, exclude=set()):

    def func(func, name, params):
        if params.args[1] == _name:
            msg = "test_{}_is_type{}"
        else:
            msg = "test_{}_is_not_type{}"
        return str(msg.format(params.args[0], _name))

    def doc(func, name, params):
        if params.args[1] == _name:
            doc = 'Check that {} ({}) is a YTK Type {} part.\n'
        else:
            doc = 'Check that {} ({}) is not a YTK Type {} part.\n'
        return doc.format(params.args[0], params.args[2], _name)

    test_plasmids = (p for p in plasmids('ytk.csv.xz') if not p[0] in exclude)

    class Test(unittest.TestCase):
        _part_cls = _cls
        _part_type = _name

        @parameterized.expand(test_plasmids, func, doc)
        def _test_plasmid(self, id_, type_, name, desc, seq):
            record = CircularRecord(Seq(seq), name=name, id=id_)
            if type_ == self._part_type:
                self.assertTrue(
                    self._part_cls(record).is_valid(),
                    '{} is not a valid Type {} but should be!'
                )
            else:
                self.assertFalse(self._part_cls(record).is_valid())

    Test.__name__ = str('Test{}'.format(_cls.__name__))
    return Test

# Generated test cases
TestYTKPart1 = test_ytk_part(ytk.YTKPart1, '1')
TestYTKPart2 = test_ytk_part(ytk.YTKPart2, '2')
TestYTKPart3 = test_ytk_part(ytk.YTKPart3, '3')
TestYTKPart3a = test_ytk_part(ytk.YTKPart3a, '3a')
TestYTKPart3b = test_ytk_part(ytk.YTKPart3b, '3b')
TestYTKPart4 = test_ytk_part(ytk.YTKPart4, '4')
TestYTKPart4a = test_ytk_part(ytk.YTKPart4a, '4a')
TestYTKPart4b = test_ytk_part(ytk.YTKPart4b, '4b')
TestYTKPart234 = test_ytk_part(ytk.YTKPart234, '234')
TestYTKPart234r = test_ytk_part(ytk.YTKPart234r, '234r', exclude={'pYTK096'})
TestYTKPart5 = test_ytk_part(ytk.YTKPart5, '5')
TestYTKPart6 = test_ytk_part(ytk.YTKPart6, '6')
TestYTKPart7 = test_ytk_part(ytk.YTKPart7, '7')
TestYTKPart8 = test_ytk_part(ytk.YTKPart8, '8')
TestYTKPart8a = test_ytk_part(ytk.YTKPart8a, '8a')
TestYTKPart8b = test_ytk_part(ytk.YTKPart8b, '8b')
TestYTKPart678 = test_ytk_part(ytk.YTKPart678, '678')


### Test Yeast ToolKit reference paper construct

class TestYTKConstruct(unittest.TestCase):

    def test_cassette_vector(self):
        """Check the assembled YTK integration vector has the right sequence.
        """
        records = {
            p[0]: CircularRecord(Seq(p[4]), id=p[0], name=p[0])
            for p in plasmids('ytk.csv.xz')
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

class TestYTKMultigene(unittest.TestCase):

    def test_multigene_assembly(self):
        """Check a YTK multigene assembly gives the expected result.
        """
        records = {
            p[0]: CircularRecord(Seq(p[4]), id=p[0], name=p[0])
            for p in plasmids('ytk-multigene.csv.xz')
        }

        print(records)

        tu1 = ytk.YTKCassette(records['mCerulean'])
        tu2 = ytk.YTKCassette(records['mNeonGreen'])
        tu3 = ytk.YTKCassette(records['mScarlet'])
        vec = ytk.YTKMultigeneVector(records['multigeneVec'])

        assembly = vec.assemble(tu1, tu2, tu3)
        expected = records['mCmNGmS']

        self.assertEqual(len(assembly), len(expected), 'lengths differ')
        self.assertIn(assembly.seq, expected.seq + expected.seq, 'sequences differ')
