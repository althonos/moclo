# coding: utf-8
"""EcoFlex ToolKit sequences registry.

Sequences were obtained from two sources:
* The EcoFlex plasmid files distributed with the kit in a zip archive,
  available under the *Protocol & Resources* tab of the CIDAR repository.
* The individual plasmid files from their dedicated AddGene webpages, using
  the *full depositor* sequences. Missing sequences (*DVA_AE* and *DVA_AF*)
  were obtained by editing the overhangs of *DVA_EF*.

Plasmids were rotated to share the same origin, using the start of the
*BioBrick* prefix as a reference location. This ensures no feature overlaps
the zero coordinate, to ensure a complete ``biopython`` compatibility.

Common features were colored using the same palette as in the *Yeast ToolKit*.
*AmpR* and *CmR* received additional cross-references from Swiss-Prot, as well
as *GFP* and *mRFP1*.

See Also:
    The annotation script running on Python 3 in the `GitHub repository
    <https://github.com/althonos/moclo/blob/master/scripts/ecoflex_registry.py>_`.

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import six

from ..kits import ecoflex
from .base import EmbeddedRegistry
from ._utils import find_resistance


class EcoFlexRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ecoflex.tar.gz"

    _VECTORS = {
        "pTU1": ecoflex.EcoFlexCassetteVector,
        "pTU2": ecoflex.EcoFlexDeviceVector,
        "pTU3": ecoflex.EcoFlexCassetteVector,
    }

    def _load_entity(self, record):
        for prefix, cls in six.iteritems(self._VECTORS):
            if record.id.startswith(prefix):
                return cls(record)
        return ecoflex.EcoFlexPart.characterize(record)
