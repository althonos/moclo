# coding: utf-8
"""Definition of the MoClo base classes.

This module contains base classes that implement the Modular Cloning logic.
All of these classes are abstract, as they need a DNA sequence structure to
be specified. Implementations can be found in the ``moclo.kits`` namespace
once installed.

The MoClo system relies on the Golden Gate assembly combined to clever
overhang design to create genetic constructs in a simple and deterministic
way. These concepts are transposed in the object model:

**Module**
    A module is a sequence of DNA that contains a sequence of interest, such
    as a promoter, a CDS, a protein binding site, etc., organised in a way
    it can be combined to other modules to create an assembly. This involves
    flanking that target sequence with Type IIS restriction sites, which
    depend on the level of the module.

**Vector**
    A vector is a plasmidic DNA sequence that is used to hold combination of
    several modules. Several modules of the same level can be assembled into
    a single vector of the same level to produce a module of the next level.
"""
from __future__ import absolute_import
from __future__ import unicode_literals

from .modules import AbstractModule, Cassette, Entry, Multigene, Product
from .vectors import AbstractVector, CassetteVector, EntryVector, MultigeneVector

__all__ = list(map(str, (
    'AbstractModule',
    'AbstractVector',
    'Cassette',
    'CassetteVector',
    'Entry',
    'EntryVector',
    'Multigene',
    'MultigeneVector',
    'Product',
)))
