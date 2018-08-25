# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from Bio.SeqFeature import SeqFeature, FeatureLocation


def cutter_check(cutter, name):
    if cutter is NotImplemented:
        raise NotImplementedError("{} does not declare a cutter".format(name))
    elif cutter.is_blunt():
        raise ValueError("cannot use a blunt cutter for Golden Gate")
    elif cutter.is_unknown():
        raise ValueError("cannot use an unknown cutter for Golden Gate")


def add_as_source(src_record, dst_record, location=None):
    quals = {
        "organism": "synthetic DNA construct",
        "mol_type": "other DNA",
        "plasmid": src_record.id,
        "label": "source: {}".format(src_record.id),
    }
    location = location or FeatureLocation(0, len(dst_record))
    feat = SeqFeature(location, type="source", qualifiers=quals)
    dst_record.features.append(feat)
    return dst_record
