# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys

# Patch the PYTHONPATH to find the base moclo package
proj = os.path.abspath(os.path.join(__file__, "..", ".."))
sys.path.insert(0, os.path.join(proj, "moclo"))

# Load the kits namespace and add additional plugins packages
import moclo.kits  # noqa: E402
import moclo.registry  # noqa: E402

# Load extensions
for extension in ["cidar", "ytk", "ecoflex", "moclo", "plant"]:
    extension_dir = os.path.join(proj, "moclo-{}".format(extension))
    moclo.kits.__path__.append(os.path.join(extension_dir, "moclo", "kits"))
    moclo.registry.__path__.append(os.path.join(extension_dir, "moclo", "registry"))
