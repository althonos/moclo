# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import gzip
import io
import itertools
import os
import textwrap
import unittest
import warnings

import six
from Bio.Seq import Seq

from moclo.record import CircularRecord
from moclo.kits import ytk
from moclo.registry.base import CombinedRegistry
from moclo.registry.ytk import YTKRegistry, PTKRegistry

from ._utils import AssemblyTestCase, PartsMetaCase, build_registries


if six.PY3:

    def setUpModule():
        build_registries('ytk')
        build_registries('cidar')
        warnings.simplefilter('ignore', category=ResourceWarning)

    def tearDownModule():
        warnings.simplefilter(warnings.defaultaction)

else:

    def setUpModule():
        build_registries('ytk')
        build_registries('cidar')


### Test Yeast ToolKit plasmids

# metaclass for test suites
_Meta = PartsMetaCase(
    'YTK', lambda: CombinedRegistry() << YTKRegistry() << PTKRegistry(), __name__,
)

# Generate test cases
TestYTKPart1 = _Meta(ytk.YTKPart1, '1')
TestYTKPart2 = _Meta(ytk.YTKPart2, '2')
TestYTKPart3 = _Meta(ytk.YTKPart3, '3')
TestYTKPart3a = _Meta(ytk.YTKPart3a, '3a')
TestYTKPart3b = _Meta(ytk.YTKPart3b, '3b')
TestYTKPart4 = _Meta(ytk.YTKPart4, '4')
TestYTKPart4a = _Meta(ytk.YTKPart4a, '4a')
TestYTKPart4b = _Meta(ytk.YTKPart4b, '4b')
TestYTKPart234 = _Meta(ytk.YTKPart234, '234')
TestYTKPart234r = _Meta(
    ytk.YTKPart234r, '234r', exclude=lambda item: item.id == 'pYTK096'
)
TestYTKPart5 = _Meta(ytk.YTKPart5, '5')
TestYTKPart6 = _Meta(ytk.YTKPart6, '6')
TestYTKPart7 = _Meta(ytk.YTKPart7, '7')
TestYTKPart8 = _Meta(ytk.YTKPart8, '8')
TestYTKPart8a = _Meta(ytk.YTKPart8a, '8a')
TestYTKPart8b = _Meta(ytk.YTKPart8b, '8b')
TestYTKPart678 = _Meta(ytk.YTKPart678, '678')


### Test Yeast ToolKit multigene assembly

class TestYTKAssembly(AssemblyTestCase, unittest.TestCase):

    def test_ytk_device(self):
        """Check a YTK device assembly using several cassettes.

        DNA sequences courtesy of InBio's member `Sebastián Sosa Carrillo
        <https://research.pasteur.fr/en/member/sebastian-sosa-carrillo/>`_.
        """
        result, vector, modules = self.load_data('ytk_device')
        self.assertAssembly(
            ytk.YTKDeviceVector(vector),
            map(ytk.YTKCassette, modules.values()),
            result,
        )

    def test_ytk_integration_vector(self):
        """Check the YTK integration vector assembly from the reference paper.
        """
        result, vector, modules = self.load_data('ytk_integration_vector')

        modules['pYTK008.gb'] = ytk.YTKPart1(modules['pYTK008.gb'])
        modules['pYTK047.gb'] = ytk.YTKPart234r(modules['pYTK047.gb'])
        modules['pYTK073.gb'] = ytk.YTKPart5(modules['pYTK073.gb'])
        modules['pYTK074.gb'] = ytk.YTKPart6(modules['pYTK074.gb'])
        modules['pYTK086.gb'] = ytk.YTKPart7(modules['pYTK086.gb'])
        modules['pYTK092.gb'] = ytk.YTKPart8b(modules['pYTK092.gb'])

        self.assertAssembly(
            ytk.YTKPart8a(vector),
            modules.values(),
            result,
        )
