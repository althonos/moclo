# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import io
import os
import fs

try:
    import lzma
except ImportError:
    from backports import lzma


# fs.osfs.OSFS: FS where test data is located
DATAFS = fs.open_fs(os.path.join(__file__, '..', 'data'))


def plasmids(name, datafs=DATAFS):
    """Load plasmids inventory archive.
    """
    with io.TextIOWrapper(lzma.open(DATAFS.openbin(name))) as f:
        for line in f:
            if not line.startswith('Plasmid Name'):
                yield line.strip().split('\t')
