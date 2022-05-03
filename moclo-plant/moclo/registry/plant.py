# coding: utf-8
"""Icon Genetics ToolKit sequences registry.

Sequences were obtained from the MoClo plasmid files distributed with the kit
in a zip archive, available under the *Protocol & Resources* tab of the
Icon Genetics MoClo repository.

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import six

from ..kits import moclo, plant
from .base import EmbeddedRegistry
from ._utils import find_resistance


class PlantRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "plant.tar.gz"

    _types = {}

    def _load_entity(self, record):
        if record.id in self._types:
            return self._types[record.id](record)
        return moclo.MoCloPart.characterize(record)
