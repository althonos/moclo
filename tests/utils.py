# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import io
import os

try:
    import lzma
except ImportError:
    from backports import lzma


### Plasmids from CSV loader

def plasmids(name):
    plasmids_file = os.path.join(__file__, '..', 'data', name)
    with io.TextIOWrapper(lzma.open(os.path.abspath(plasmids_file))) as f:
        for line in f:
            if not line.startswith('Plasmid Name'):
                yield line.strip().split('\t')
