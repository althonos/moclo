# coding: utf-8
"""Icon Genetics ToolKit sequences registry.

Sequences were obtained from the MoClo plasmid files distributed with the kit
in a zip archive, available under the *Protocol & Resources* tab of the
Icon Genetics MoClo repository.

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import six

from ..kits import ig
from .base import EmbeddedRegistry
from ._utils import find_resistance


class IGRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ig.json.bz2"

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

    _VECTORS = {
        "pTU1": ig.EcoFlexCassetteVector,
        "pTU2": ecoflex.EcoFlexDeviceVector,
        "pTU3": ecoflex.EcoFlexCassetteVector,
    }

    def _load_entity(self, raw, index):
        for prefix, cls in six.iteritems(self._VECTORS):
            if raw["record"].id.startswith(prefix):
                return cls(raw["record"])
        return ig.IGPart.characterize(raw["record"])
