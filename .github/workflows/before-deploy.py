#!/usr/bin/env python3
#coding: utf-8

import glob
import os
import re
import subprocess
import sys

ref = os.getenv('GITHUB_REF')
if not ref.startswith("refs/tags/"):
    raise RuntimeError("not a tagged commit: {!r}".format(ref))

if ref.count("/") == 2:
    libdir = 'moclo'
elif ref.count("/") == 3:
    lib = tag.split('/')[-2]
    libdir = 'moclo-{}'.format(lib)
else:
    raise RuntimeError("Could not recognize ref: {!r}".format(ref))

args = [
    sys.executable,
    'setup.py',
    'check',
    '-rms',
    'sdist',
    '-d',
    os.path.join(os.pardir, 'dist'),
    'bdist_wheel',
    '-d',
    os.path.join(os.pardir, 'dist'),
]

print("Deploying", libdir)
subprocess.Popen(args, cwd=libdir).communicate()

wheel = next(glob.iglob(os.path.join('dist', '*.whl')))
new_wheel = re.sub('cp38-cp38m-linux_x86_64', 'py2.py3-none-any', wheel)

os.rename(wheel, new_wheel)
