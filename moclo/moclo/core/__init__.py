# coding: utf-8
# noqa: F401
"""MoClo object model and abstract classes.

This module contains base classes that implement the Modular Cloning logic.
All of these classes are abstract, as they need a DNA sequence structure to
be specified. Implementations can be found in the ``moclo.kits`` namespace
once installed. The MoClo system relies on the Golden Gate assembly combined
to clever sequence design to create genetic constructs in a simple and
deterministic way.
"""
from __future__ import absolute_import

from .parts import AbstractPart
from .modules import AbstractModule, Cassette, Entry, Device, Product
from .vectors import AbstractVector, CassetteVector, EntryVector, DeviceVector

__all__ = [
    "AbstractPart",
    "AbstractModule",
    "AbstractVector",
    "Cassette",
    "CassetteVector",
    "Entry",
    "EntryVector",
    "Device",
    "DeviceVector",
    "Product",
]
