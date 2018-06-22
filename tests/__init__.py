# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

try:
    # Try to load the system version
    import moclo.kits.ytk

except ImportError:

    import importlib
    import os
    import sys

    # Patch the PYTHONPATH to find the base moclo package
    proj = os.path.abspath(os.path.join(__file__, '..', '..'))
    sys.path.insert(0, os.path.join(proj, 'moclo'))
    importlib.reload(sys.modules['moclo'])

    # Load the kits namespace and add additional plugins packages
    import moclo.kits

    # Yeast ToolKit
    moclo.kits.__path__.append(os.path.join(proj, 'moclo-ytk', 'moclo', 'kits'))
    import moclo.kits.ytk
