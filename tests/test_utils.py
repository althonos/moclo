# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import unittest

from moclo._utils import isabstract


class TestIsAbstract(unittest.TestCase):
    def test_abstract_method(self):
        self.assertTrue(isabstract(collections.Iterable))
        self.assertFalse(isabstract(int))

    def test_abstract_attribute(self):
        class Abstract(object):
            thingy = NotImplemented

        self.assertTrue(isabstract(Abstract))
