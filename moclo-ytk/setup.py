#!/usr/bin/env python3
# coding: utf-8

import bz2
import json
import glob
import os
import setuptools
import sys
import tarfile

from setuptools.command.build_ext import build_ext as _build_ext


class Registry(setuptools.Extension):
    def __init__(self, name):
        directory = os.path.join("registry", name)
        sources = glob.glob(os.path.join(directory, "*.gb"))
        setuptools.Extension.__init__(self, name, sources)


class build_ext(_build_ext):
    def get_ext_filename(self, ext_name):
        return os.path.join(*ext_name.split(".")) + ".tar.gz"

    def build_extension(self, ext):
        # find directories
        registry = []
        gb_dir = os.path.dirname(ext.sources[0])
        dst_file = self.get_ext_fullpath(ext.name)
        dst_dir = os.path.dirname(dst_file)
        # copy sequences
        self.mkpath(dst_dir)
        with tarfile.open(dst_file, mode="w:gz") as tar:
            for gb_file in sorted(ext.sources):
                arcname, _ = os.path.splitext(os.path.basename(gb_file))
                tar.add(gb_file, arcname=arcname)


if __name__ == "__main__":
    setuptools.setup(
        ext_package="moclo.registry",
        ext_modules=[Registry("ytk"), Registry("ptk")],
        cmdclass={"build_ext": build_ext},
    )
