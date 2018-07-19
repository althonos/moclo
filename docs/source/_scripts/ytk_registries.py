# coding: utf-8
import os
import sys
import subprocess

DOCSRC_DIR = os.path.abspath(os.path.join(__file__, '..', '..'))
YTK_DIR = os.path.abspath(os.path.join(DOCSRC_DIR, '..', '..', 'moclo-ytk'))

def build_registries():
    cmd = subprocess.Popen(
        [sys.executable, 'setup.py', 'build_ext', '--inplace'],
        cwd=YTK_DIR,
    )
    cmd.communicate()
