# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import doctest
import sys
import re
import types
import unittest
import warnings

import Bio.Seq
import Bio.SeqRecord
import Bio.Restriction

import moclo
import moclo.core
import moclo.kits


class IgnoreUnicodeChecker(doctest.OutputChecker):
    def check_output(self, want, got, optionflags):
        if sys.version_info[0] > 2:
            want = re.sub("u'(.*?)'", "'\\1'", want)
            want = re.sub('u"(.*?)"', '"\\1"', want)
        return doctest.OutputChecker.check_output(self, want, got, optionflags)


def _load_tests_from_module(tests, module, globs, setUp, tearDown):
    """Load tests from module, iterating through submodules.
    """
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            tests.addTests(
                doctest.DocTestSuite(
                    attr,
                    globs=globs,
                    setUp=setUp,
                    tearDown=tearDown,
                    checker=IgnoreUnicodeChecker(),
                )
            )
    return tests


def load_tests(loader=None, tests=unittest.TestSuite(), ignore=False):
    """load_test function used by unittest to find the doctests"""

    def _setUp(self):
        return

    def _tearDown(self):
        return

    globs = {
        # Biopython
        "Seq": Bio.Seq.Seq,
        "SeqRecord": Bio.SeqRecord.SeqRecord,
    }

    # Load all abstract base classes in the namespace
    for name in moclo.core.__all__:
        globs[name] = getattr(moclo.core, name)

    # Load all restriction enzymes in the namespace
    RestrictionType = Bio.Restriction.Restriction.RestrictionType
    for name in dir(Bio.Restriction):
        attr = getattr(Bio.Restriction, name)
        if isinstance(attr, type) and issubclass(attr, RestrictionType):
            globs[name] = attr

    if not sys.argv[0].endswith("green"):
        _load_tests_from_module(tests, moclo, globs, _setUp, _tearDown)
        _load_tests_from_module(tests, moclo.core, globs, _setUp, _tearDown)
        _load_tests_from_module(tests, moclo.kits, globs, _setUp, _tearDown)

    return tests


def setUpModule():
    warnings.simplefilter("ignore")


def tearDownModule():
    warnings.simplefilter(warnings.defaultaction)
