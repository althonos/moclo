# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import subprocess
import sys

import Bio.SeqIO
import contexter
import fs.path
import fs.archive.tarfs

from moclo.record import CircularRecord

# --- Constants --------------------------------------------------------------

# fs.osfs.OSFS: FS located at the root of the project
PROJFS = fs.open_fs(os.path.join(__file__, "..", ".."))

# fs.osfs.OSFS: FS where test data is located
DATAFS = PROJFS.opendir("tests/data")


# --- Setup Helper -----------------------------------------------------------

class Registries(object):
    def __init__(self):
        self._built = set()

    def __call__(self, name):
        if name not in self._built:
            subprocess.Popen(
                args=[sys.executable, "setup.py", "build_ext", "-i"],
                cwd=PROJFS.getsyspath("moclo-{}".format(name)),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            ).communicate()
            self._built.add(name)


build_registries = Registries()


# --- Tests helpers ----------------------------------------------------------

def expectFailure(cls, method):
    setattr(cls, method, unittest.expectedFailure(getattr(cls, method)))


class PartsMetaCase(object):
    def __init__(self, kit_name, registry_factory, module):
        build_registries(kit_name.lower())
        self.kit_name = kit_name
        self.registry = registry_factory()
        self.module = module

    def __call__(self, part_cls, part_name, exclude=lambda item: False):
        """Create a whole test case for the given class.
        """
        tests = (
            self.make_test(item, part_cls=part_cls, part_name=part_name)
            for item in self.registry.values()
            if not exclude(item)
        )
        attrs = {test.__name__: test for test in tests}
        attrs["__name__"] = name = str("Test{}".format(part_cls.__name__))
        attrs["__module__"] = self.module
        return type(name, (unittest.TestCase,), attrs)

    def make_test(self, item, part_cls, part_name):
        """Create a single test for the given item.

        If ``ìtem.entity`` is a ``part_cls`` instance, it will check that
        ``part_cls(item.entity.record).is_valid()`` is `True` (`False`
        otherwise), i.e. it checks that the record is only identified as
        the actual type it should be.
        """
        rec = item.entity.record
        if isinstance(item.entity, part_cls):

            def test(self_):
                err = "{} is not a valid {} {} but should be!"
                self_.assertTrue(
                    part_cls(rec).is_valid(),
                    err.format(item.id, self.kit_name, part_name),
                )

            name = "test_{}_is_{}"
            doc = "Check that {} ({}) is a {} {}.\n"
        else:

            def test(self_):
                err = "{} is a valid {} {} but should not be!"
                self_.assertFalse(
                    part_cls(rec).is_valid(),
                    err.format(item.id, self.kit_name, part_name),
                )

            name = "test_{}_is_not_{}"
            doc = "Check that {} ({}) is not a {} {}.\n"
        test.__name__ = str(name.format(item.id, part_name))
        test.__doc__ = doc.format(
            item.id, item.entity.record.name, self.kit_name, part_name
        )
        return test


class AssemblyTestCase(unittest.TestCase):
    def assertAssembly(self, vector, modules, result):
        assembly = vector.assemble(*modules)
        self.assertEqual(len(assembly), len(result), "lengths differ")
        self.assertIn(assembly.seq, result.seq + result.seq, "sequences differ")

    def load_data(self, name):
        archive_path = "cases/{}.tar.xz".format(name)

        if not DATAFS.exists(archive_path):
            raise unittest.SkipTest("no test case found")

        with contexter.Contexter() as ctx:
            # open FASTA files
            casefs = ctx << fs.archive.open_archive(DATAFS, archive_path)
            result_fa = ctx << casefs.open("result.fa")
            vector_fa = ctx << casefs.open("vector.fa")
            modules_fa = ctx << casefs.open("modules.fa")
            # load records from FASTA handles
            res = CircularRecord(Bio.SeqIO.read(result_fa, "fasta"))
            vec = CircularRecord(Bio.SeqIO.read(vector_fa, "fasta"))
            mods = {
                record.id: CircularRecord(record)
                for record in Bio.SeqIO.parse(modules_fa, "fasta")
            }
        return res, vec, mods
