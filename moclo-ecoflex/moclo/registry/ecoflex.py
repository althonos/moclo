# coding: utf-8
"""EcoFlex ToolKit sequences registry.

See Also:
    The annotation script running on Python 3 in the `GitHub repository
    <https://github.com/althonos/moclo/blob/master/scripts/ecoflex_registry.py>_`

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import collections

import six

from ..kits import ecoflex
from .base import EmbeddedRegistry
from ._utils import find_resistance


class EcoFlexRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ecoflex.json.bz2"

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
        return raw["resistance"]

    _VECTORS = collections.OrderedDict(
        [
            ("pTU1", ecoflex.EcoFlexCassetteVector),
            ("pTU2", ecoflex.EcoFlexDeviceVector),
            ("pTU3", ecoflex.EcoFlexCassetteVector),
        ]
    )

    def _load_entity(self, raw, index):
        for prefix, cls in six.iteritems(self._VECTORS):
            if raw["record"].id.startswith(prefix):
                return cls(raw["record"])
        return ecoflex.EcoFlexPart.characterize(raw["record"])
