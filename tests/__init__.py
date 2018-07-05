# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import importlib
import os
import sys

# Patch the PYTHONPATH to find the base moclo package
proj = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, os.path.join(proj, 'moclo'))

# Load the kits namespace and add additional plugins packages
import moclo.kits

# CIDAR
moclo.kits.__path__.append(os.path.join(proj, 'moclo-cidar', 'moclo', 'kits'))
import moclo.kits.cidar

# Yeast ToolKit
moclo.kits.__path__.append(os.path.join(proj, 'moclo-ytk', 'moclo', 'kits'))
import moclo.kits.ytk
