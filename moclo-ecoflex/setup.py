#!/usr/bin/env python3
# coding: utf-8

import bz2
import json
import glob
import os
import setuptools
import sys

from setuptools.command.build_ext import build_ext as _build_ext


if sys.version_info[0] == 2:

    def bz2_open_write(filename):
        return bz2.BZ2File(filename, "w")


else:

    def bz2_open_write(filename):
        return bz2.open(filename, "wt")


class Registry(setuptools.Extension):
    def __init__(self, name):
        directory = os.path.join("registry", name)
        sources = glob.glob(os.path.join(directory, "*.gb"))
        setuptools.Extension.__init__(self, name, sources)


class build_ext(_build_ext):
    def get_ext_filename(self, ext_name):
        return os.path.join(*ext_name.split(".")) + ".json.bz2"

    def build_extension(self, ext):
        # find directories
        registry = []
        gb_dir = os.path.dirname(ext.sources[0])
        dst_file = self.get_ext_fullpath(ext.name)
        dst_dir = os.path.dirname(dst_file)
        # read all records
        self.announce("collecting records from {}".format(gb_dir), 2)
        for gb_file in sorted(ext.sources):
            id_, _ = os.path.splitext(os.path.basename(gb_file))
            with open(gb_file) as gb_rec:
                registry.append({"id": id_, "gb": gb_rec.read()})
        # write the compressed registry
        self.mkpath(dst_dir)
        self.announce("writing {} records to {}".format(len(registry), dst_file), 2)
        with bz2_open_write(dst_file) as reg:
            json.dump(registry, reg)


if __name__ == "__main__":
    setuptools.setup(
        ext_package="moclo.registry",
        ext_modules=[Registry("ecoflex")],
        cmdclass={"build_ext": build_ext},
    )
