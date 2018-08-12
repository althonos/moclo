# coding: utf-8
import os
import sys
import subprocess

DOCSRC_DIR = os.path.abspath(os.path.join(__file__, '..', '..'))
_EXT_DIR = os.path.abspath(os.path.join(DOCSRC_DIR, '..', '..', 'moclo-{}'))

def build_registries(name):
    cmd = subprocess.Popen(
        [sys.executable, 'setup.py', 'build_ext', '--inplace'],
        cwd=_EXT_DIR.format(name),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    cmd.communicate()
