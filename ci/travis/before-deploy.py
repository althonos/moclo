#!/usr/bin/env python3
#coding: utf-8

import glob
import os
import re
import subprocess
import sys

tag = os.getenv('TRAVIS_TAG')
lib = tag.split('/')[0]

if lib == 'ytk':
    libdir = 'moclo-ytk'
elif lib == 'cidar':
    libdir = 'moclo-cidar'
elif lib == 'eclofex':
    libdir = 'moclo-ecoflex'
else:
    libdir = 'moclo'

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
new_wheel = re.sub('cp36-cp36m-linux_x86_64', 'py2.py3-none-any', wheel)

os.rename(wheel, new_wheel)
