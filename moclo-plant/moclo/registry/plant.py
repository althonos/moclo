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
    _file = "plant.json.bz2"

    _types = {}

    def _load_name(self, raw, index):
        return raw["record"].name

    def _load_id(self, raw, index):
        return raw["record"].id

    def _load_resistance(self, raw, index):
        try:
            return find_resistance(raw["record"])
        except RuntimeError:
            msg = "could not find antibiotics resistance of '{}'"
            six.raise_from(RuntimeError(msg.format(raw["record"].id)), None)

    def _load_entity(self, raw, index):
        if raw["record"].id in self._types:
            return self._types[raw["record"].id](raw["record"])
        return moclo.MoCloPart.characterize(raw["record"])
