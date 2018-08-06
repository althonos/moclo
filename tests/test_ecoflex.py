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
from moclo.kits import ecoflex
from moclo.registry.ecoflex import EcoFlexRegistry

from ._utils import PartsMetaCase, build_registries


### Test EcoFlex MoClo plasmids

# metaclass for test suites
_Meta = PartsMetaCase('EcoFlex', EcoFlexRegistry, __name__)
exclude_tu2 = lambda item: item.id.startswith('pTU2')
exclude_tu2a = lambda item: item.id.startswith('pTU2-A')
exclude_tu2d = lambda item: item.id.startswith('pTU2-D')
exclude_tu3 = lambda item: item.id.startswith('pTU3')
exclude_cassette = lambda item: item.id.startswith(('pTU1', 'pTU3'))

# Generate test cases
TestEcoFlexPromoter = _Meta(ecoflex.EcoFlexPromoter, 'Promoter', exclude_tu2a)
TestEcoFlexRBS = _Meta(ecoflex.EcoFlexRBS, 'RBS')
TestEcoFlexTag = _Meta(ecoflex.EcoFlexTag, 'Tag')
TestEcoFlexTerminator = _Meta(ecoflex.EcoFlexTerminator, 'Terminator', exclude_tu2d)
TestEcoFlexCodingSequence = _Meta(ecoflex.EcoFlexCodingSequence, 'CDS')
TestEcoFlexPromoterRBS = _Meta(ecoflex.EcoFlexPromoterRBS, 'PromoterRBS')
TestEcoFlexCassetteVector = _Meta(ecoflex.EcoFlexCassetteVector, 'CassetteVector', exclude_tu2)
TestEcoFlexDeviceVector = _Meta(ecoflex.EcoFlexDeviceVector, 'DeviceVector', exclude_cassette)
