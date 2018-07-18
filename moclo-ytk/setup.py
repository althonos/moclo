#!/usr/bin/env python3
# coding: utf-8

import bz2
# import distutils.cmd
import distutils.core
import json
import glob
import os
import setuptools

from setuptools.command.build_ext import build_ext as _build_ext


class Registry(distutils.core.Extension):

    def __init__(self, name):
        directory = os.path.join('registry', name)
        sources = glob.glob(os.path.join(directory, '*.gb'))
        super(Registry, self).__init__(name, sources)


class build_ext(_build_ext):

    def build_extension(self, ext):

        if not isinstance(ext, Registry):
            return _build_ext.build_extension(ext)

        registry = []
        gb_dir = os.path.dirname(ext.sources[0])
        dst_dir = os.path.dirname(self.get_ext_fullpath(ext.name))
        dst_file = os.path.join(dst_dir, '{}.json.bz2').format(ext.name)

        self.announce('collecting records from {}'.format(gb_dir), 2)
        for gb_file in sorted(ext.sources):
            id_, _ = os.path.splitext(os.path.basename(gb_file))
            with open(gb_file) as gb_rec:
                registry.append({'id': id_, 'gb': gb_rec.read()})

        self.announce('writing {} records to {}'.format(len(registry), dst_file), 2)
        with bz2.open(dst_file, 'wt') as reg:
            json.dump(registry, reg)



if __name__ == '__main__':
    setuptools.setup(
        ext_package='moclo.registry',
        ext_modules=[Registry('ytk')],
        cmdclass={'build_ext': build_ext},
    )
