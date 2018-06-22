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


PLASMIDS_FILE = \
    os.path.abspath(os.path.join(__file__, '..', 'ytk.csv.xz'))
PLASMIDS = (
    line.split('\t')
    for line in io.TextIOWrapper(lzma.open(PLASMIDS_FILE))
    if not line.startswith('Plasmid Name')
)

### Parameterized format functions

def _part_func_name(func, name, params):
    if params.args[1] not in ('cassette', 'entry vector'):
        return "test_type_{}_is_type{}".format(params.args[0], params.args[1])
    else:
        return "test_type_{}_is_{}".format(params.args[0], params.args[1])

def _part_doc_name(func, name, params):
    if params[1] not in ('cassette', 'entry vector'):
        doc = "Check that {} is identified as a YTK Type {} part.\n"
    else:
        doc = "Check that {} is identified as a YTK {} part.\n"
    return doc.format(params.args[0], params.args[1])


### Test Yeast ToolKit plasmids

class TestYTKPlasmids(unittest.TestCase):

    parts = {
        '1': ytk.YTKPart1,
        '2': ytk.YTKPart2,
        '3': ytk.YTKPart3,
        '3a': ytk.YTKPart3a,
        '3b': ytk.YTKPart3b,
        '4': ytk.YTKPart4,
        '4a': ytk.YTKPart4a,
        '4b': ytk.YTKPart4b,
        '234': ytk.YTKPart234,
        '234r': ytk.YTKPart234r,
        '5': ytk.YTKPart5,
        '6': ytk.YTKPart6,
        '7': ytk.YTKPart7,
        '8': ytk.YTKPart8,
        '8a': ytk.YTKPart8a,
        '8b': ytk.YTKPart8b,
        '678': ytk.YTKPart678,
        'entry vector': ytk.YTKEntryVector,
        'cassette': ytk.YTKCassetteVector,
    }

    @parameterized.expand(PLASMIDS, _part_func_name, _part_doc_name)
    def test(self, id_, type_, name, desc, seq):

        record = CircularRecord(seq=Seq(seq), name=name, id=id_)
        part_cls = self.parts[type_]
        self.assertTrue(part_cls(record).is_valid())

        # for part in six.itervalues(self.parts):
            # if part is self.parts[type_]:
                # self.assertTrue(part(record).is_valid())
            # else:
                # self.assertFalse(part(record).is_valid())

        # info = self.parts[id_]
        # record = CircularRecord(seq=info['seq'], name=info['name'], id=id_)
        # part = ytk.YTKPart1(record)
        # if info['type'] == '1':
        #     self.assertTrue(part.is_valid())
        # else:
        #     self.assertFalse(part.is_valid())
