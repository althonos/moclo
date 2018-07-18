#!/usr/bin/env python3
# coding: utf-8

import bz2
import distutils.cmd
import json
import glob
import os
import setuptools


class BuildRegistry(distutils.cmd.Command):

    description = 'compile registry of records to a single compressed file'
    user_options = []

    def initialize_options(self):
        self.project_dir = None
        self.src_dir = None
        self.dst_dir = None

    def finalize_options(self):
        self.project_dir = os.path.dirname(__file__)
        self.src_dir = os.path.join(self.project_dir, 'registry')
        self.dst_dir = os.path.join(self.project_dir, 'moclo', 'registry')

    def run(self):
        self._build_registry('ytk')
        #self._build_registry('ptk')

    def _build_registry(self, name, ):
        registry = []
        gb_dir = os.path.join(self.src_dir, name)
        dst_file = os.path.join(self.dst_dir, '{}.json.bz2').format(name)

        self.announce('collecting records from {}'.format(gb_dir), 2)
        for gb_file in sorted(glob.glob(os.path.join(gb_dir, '*.gb'))):
            id_, _ = os.path.splitext(os.path.basename(gb_file))
            with open(gb_file) as gb_rec:
                registry.append({'id': id_, 'gb': gb_rec.read()})

        self.announce('writing {} records to {}'.format(len(registry), dst_file), 2)
        with bz2.open(dst_file, 'wt') as reg:
            json.dump(registry, reg)



if __name__ == '__main__':
    setuptools.setup(cmdclass={'build_registry': BuildRegistry})
