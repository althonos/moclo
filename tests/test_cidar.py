# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import gzip
import io
import itertools
import os
import textwrap
import sys
import unittest

import six
from Bio.Seq import Seq

from moclo.record import CircularRecord
from moclo.kits import cidar

from ._utils import PartsMetaCase


### Test CIDAR plasmids

# test suite metaclass
_Meta = PartsMetaCase('CIDAR', 'cidar.tsv.xz', __name__)

# Generate test cases for each parts
TestCIDARPromoter = _Meta(cidar.CIDARPromoter, 'Promoter')
TestCIDARibosomeBindingSite = _Meta(cidar.CIDARibosomeBindingSite, 'RBS')
TestCIDARCodingSequence = _Meta(cidar.CIDARCodingSequence, 'CDS')
TestCIDARTerminator = _Meta(cidar.CIDARTerminator, 'Terminator')
