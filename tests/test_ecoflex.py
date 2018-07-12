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

from ._utils import PartsMetaCase



### Test EcoFlex MoClo plasmids

# metaclass for test suites
_Meta = PartsMetaCase('EcoFlex', 'ecoflex.tsv.xz', __name__)

# Generate test cases
TestEcoFlexPromoter = _Meta(ecoflex.EcoFlexPromoter, 'promoter')
TestEcoFlexRBS = _Meta(ecoflex.EcoFlexRBS, 'RBS')
TestEcoFlexTag = _Meta(ecoflex.EcoFlexTag, 'tag')
TestEcoFlexTerminator = _Meta(ecoflex.EcoFlexTerminator, 'terminator')
TestEcoFlexCodingSequence = _Meta(ecoflex.EcoFlexCodingSequence, 'CDS', exclude=('pBP-ORF'))
TestEcoFlexPromoterRBS = _Meta(ecoflex.EcoFlexPromoterRBS, 'promoterRBS')
