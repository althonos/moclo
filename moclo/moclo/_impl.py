# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import importlib


def _import_from(*names):
    for name in names:
        if name is None:
            return name
        try:
            return importlib.import_module(name)
        except ImportError:
            pass
    raise ImportError("no module found among: {}".format(", ".join(names)))


bz2 = _import_from("bz2file", "bz2")
json = _import_from("hyperjson", "ujson", "yajl", "rapidjson", "simplejson", "json")
ssl = _import_from("ssl", None)
