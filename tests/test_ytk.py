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

def plasmids():
    plasmids_file = os.path.abspath(os.path.join(__file__, '..', 'ytk.csv.xz'))
    with io.TextIOWrapper(lzma.open(plasmids_file)) as f:
        for line in f:
            if not line.startswith('Plasmid Name'):
                yield line.strip().split('\t')


### Parameterized format functions



    # if params.args[1] not in ('cassette', 'entry vector'):
    #     doc =
    # else:
    #     doc =
    # return doc.format(params.args[0], params.args[1])


### Test Yeast ToolKit plasmids

# Metaclass for test cases: create a test case per type with a method per
# plasmid that will check if that plasmid is recognized as the test case
# part or not.
def test_ytk_part(_cls, _name, exclude=set()):

    def func(func, name, params):
        if params.args[1] == _name:
            msg = "test_{}_is_type{}"
        else:
            msg = "test_{}_is_not_type{}"
        return str(msg.format(params.args[0], _name))

    def doc(func, name, params):
        if params.args[1] == _name:
            doc = "Check that {} ('{}') is a YTK Type {} part.\n"
        else:
            doc = "Check that {} ('{}') is not a YTK Type {} part.\n"
        return doc.format(params.args[0], params.args[2], _name)

    class Test(unittest.TestCase):
        _part_cls = _cls
        _part_type = _name

        @parameterized.expand(plasmids(), func, doc)
        def _test_plasmid(self, id_, type_, name, desc, seq):
            record = CircularRecord(Seq(seq), name=name, id=id_)
            if type_ == self._part_type:
                self.assertTrue(self._part_cls(record).is_valid())
            else:
                self.assertFalse(self._part_cls(record).is_valid())
    Test.__name__ = str('Test{}'.format(_cls.__name__))
    return Test


TestYTKPart1 = test_ytk_part(ytk.YTKPart1, '1')
TestYTKPart2 = test_ytk_part(ytk.YTKPart2, '2')
TestYTKPart3 = test_ytk_part(ytk.YTKPart3, '3')
TestYTKPart3a = test_ytk_part(ytk.YTKPart3a, '3a')
TestYTKPart3b = test_ytk_part(ytk.YTKPart3b, '3b')
TestYTKPart4 = test_ytk_part(ytk.YTKPart4, '4')
TestYTKPart4a = test_ytk_part(ytk.YTKPart4a, '4a')
TestYTKPart4b = test_ytk_part(ytk.YTKPart4b, '4b')
TestYTKPart5 = test_ytk_part(ytk.YTKPart5, '5')
TestYTKPart6 = test_ytk_part(ytk.YTKPart6, '6')
TestYTKPart7 = test_ytk_part(ytk.YTKPart7, '7')
TestYTKPart8 = test_ytk_part(ytk.YTKPart8, '8')
TestYTKPart8a = test_ytk_part(ytk.YTKPart8a, '8a')
TestYTKPart8b = test_ytk_part(ytk.YTKPart8b, '8b')

# class TestYTKPlasmids(unittest.TestCase):
#
#     parts = {
#         '1': ytk.YTKPart1,
#         '2': ytk.YTKPart2,
#         '3': ytk.YTKPart3,
#         '3a': ytk.YTKPart3a,
#         '3b': ytk.YTKPart3b,
#         '4': ytk.YTKPart4,
#         '4a': ytk.YTKPart4a,
#         '4b': ytk.YTKPart4b,
#         '234': ytk.YTKPart234,
#         '234r': ytk.YTKPart234r,
#         '5': ytk.YTKPart5,
#         '6': ytk.YTKPart6,
#         '7': ytk.YTKPart7,
#         '8': ytk.YTKPart8,
#         '8a': ytk.YTKPart8a,
#         '8b': ytk.YTKPart8b,
#         '678': ytk.YTKPart678,
#         'entry vector': ytk.YTKEntryVector,
#         'cassette': ytk.YTKCassetteVector,
#     }
#
#     @parameterized.expand(plasmids(), _part_func_name, _part_doc_name)
#     def test_plasmid(self, id_, type_, name, desc, seq):
#
#         record = CircularRecord(seq=Seq(seq), name=name, id=id_)
#         part_cls = self.parts[type_]
#
#         for part in six.itervalues(self.parts):
#             # Check the part is valid
#             if part is self.parts[type_]:
#                 self.assertTrue(
#                     part(record).is_valid(),
#                     "{} is not a valid {} but should be!".format(type_, part)
#                 )
#             # A cassette can be mistaken for a YTKPart234r
#             elif type_ == 'cassette' and part in (ytk.YTKPart234r, ytk.YTKEntryVector):
#                 continue
#             # Part 234 can be mistaken for an entry vector
#             elif type_ == '234' and part is ytk.YTKEntryVector:
#                 continue
#             # All parts are also cassette vectors
#             elif part is not ytk.YTKCassetteVector:
#                 self.assertFalse(
#                     part(record).is_valid(),
#                     "{} is a valid {} but should not be!".format(type_, part)
#                 )



###

class TestYTKConstruct(unittest.TestCase):

    def test_cassette_vector(self):
        records = {
            p[0]: CircularRecord(Seq(p[4]), id=p[0], name=p[0])
            for p in plasmids()
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
        self.assertIn(assembly.seq, expected.seq + expected.seq, 'sequences are not rotation of each other')
