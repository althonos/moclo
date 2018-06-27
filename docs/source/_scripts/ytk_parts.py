# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import gzip
import os

import Bio.Seq

# Path to doc files
DOCSRC_DIR = os.path.abspath(os.path.join(__file__, '..', '..'))
YTKKIT_DIR = os.path.join(DOCSRC_DIR, 'kits', 'ytk')

# Image configuration for each part
PARTS_CONF = {
    'type1':  ('CCCT', '84c4dd', 'TTGC', 'Type 1'),
    'type2':  ('AACG', 'a9d55e', 'ATAC', 'Type 2'),
    'type3':  ('TATG', 'f3eb91', 'TAGG', 'Type 3'),
    'type3a': ('TATG', 'f3eb91', 'AAGA', 'Type 3a'),
    'type3b': ('TTCT', 'f3eb91', 'TAGG', 'Type 3b'),
    'type4':  ('ATCC', 'ee8d9c', 'CGAC', 'Type 4'),
    'type4a': ('ATCC', 'ee8d9c', 'ACCG', 'Type 4a'),
    'type4b': ('TGGC', 'ee8d9c', 'CGAC', 'Type 4b'),
    'type5':  ('GCTG', '8c69d5', 'ATGT', 'Type 5'),
    'type6':  ('TACA', 'f4d091', 'CTCA', 'Type 6'),
    'type7':  ('GAGT', '8c6239', 'GGCT', 'Type 7'),
    'type8':  ('CCGA', '7f7f7f', 'GGGA', 'Type 8'),
    'type8a': ('CCGA', '7f7f7f', 'GTTA', 'Type 8a'),
    'type8b': ('CAAT', '7f7f7f', 'GGGA', 'Type 8b'),
}


def generate_svg(workdir=YTKKIT_DIR,
                 template=os.path.join(YTKKIT_DIR, 'part.svg.gz'),
                 conf=PARTS_CONF):

    with gzip.open(template, 'rt') as handle:
        template = handle.read()

    for id_, (up, color, down_r, name) in conf.items():
        filename = os.path.join(YTKKIT_DIR, '{}.svg'.format(id_))
        up_seq = Bio.Seq.Seq(up)
        down_r_seq = Bio.Seq.Seq(down_r)

        with open(filename, 'w') as handle:
            handle.write(template.format(
                upstream=up_seq,
                upstream_r=up_seq.complement(),
                downstream=down_r_seq.complement(),
                downstream_r=down_r,
                type=name,
                color=color
            ))
