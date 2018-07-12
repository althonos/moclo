# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import io
import os
import unittest
import warnings

import Bio.SeqIO
import contexter
import fs.path
import fs.archive.tarfs
from Bio.Seq import Seq

try:
    import lzma
except ImportError:
    from backports import lzma

from moclo.record import CircularRecord


# fs.osfs.OSFS: FS where test data is located
DATAFS = fs.open_fs(os.path.join(__file__, '..', 'data'))


class PartsMetaCase(object):

    _plasmids = None

    def __init__(self, kit_name, archive, module):
        self.kit_name = kit_name
        self.archive = archive
        self.module = module

    def plasmids(self, datafs=DATAFS.opendir('parts')):
        with io.TextIOWrapper(lzma.open(datafs.openbin(self.archive))) as f:
            for line in f:
                if not line.startswith('Plasmid Name'):
                    yield line.strip().split('\t')

    def __call__(self, part_cls, part_name, exclude=frozenset()):
        if self._plasmids is None:
            self._plasmids = list(self.plasmids())
        attrs = {
            test.__name__: test
            for test in (
                self.make(*plasmid, part_cls=part_cls, part_name=part_name)
                for plasmid in self._plasmids
                if plasmid[0] not in exclude
            )
        }
        attrs['__name__'] = name = str('Test{}'.format(part_cls.__name__))
        attrs['__module__'] = self.module
        return type(name, (unittest.TestCase,), attrs)

    def make(self, id_, type_, name, desc, seq, part_cls, part_name):
        rec = CircularRecord(Seq(seq), name=name, id=id_)
        if type_ == part_name:
            def test(self_):
                err = '{} is not a valid {} {} but should be!'
                self_.assertTrue(
                    part_cls(rec).is_valid(),
                    err.format(id_, self.kit_name, part_name)
                )
            name = 'test_{}_is_{}'
            doc = 'Check that {} ({}) is a {} {}.\n'
        else:
            def test(self_):
                err = '{} is a valid {} {} but should not be!'
                self_.assertFalse(
                    part_cls(rec).is_valid(),
                    err.format(id_, self.kit_name, part_name)
                )
            name = "test_{}_is_not_{}"
            doc = 'Check that {} ({} - {}) is not a {} {}.\n'
        test.__name__ = str(name.format(id_, part_name))
        test.__doc__ = str(doc.format(id_, type_, name, self.kit_name, part_name))
        return test


class AssemblyTestCase(unittest.TestCase):

    def assertAssembly(self, vector, modules, result):
        assembly = vector.assemble(*modules)
        self.assertEqual(len(assembly), len(result), 'lengths differ')
        self.assertIn(assembly.seq, result.seq + result.seq, 'sequences differ')

    def load_data(self, name):
        archive_path = 'cases/{}.tar.xz'.format(name)
        with contexter.Contexter() as ctx:
            # open FASTA files
            casefs = ctx << fs.archive.open_archive(DATAFS, archive_path)
            result_fa = ctx << casefs.open('result.fa')
            vector_fa = ctx << casefs.open('vector.fa')
            modules_fa = ctx << casefs.open('modules.fa')
            # load records from FASTA handles
            res = CircularRecord(Bio.SeqIO.read(result_fa, 'fasta'))
            vec = CircularRecord(Bio.SeqIO.read(vector_fa, 'fasta'))
            mods = {
                record.id: CircularRecord(record)
                for record in Bio.SeqIO.parse(modules_fa, 'fasta')
            }
        return res, vec, mods
