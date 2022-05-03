# coding: utf-8
"""CIDAR ToolKit sequences registry.

Sequences were obtained from two sources:
* The CIDAR plasmid files distributed with the kit in a zip archive,
  available under the *Protocol & Resources* tab of the CIDAR repository.
* The individual plasmid files from their dedicated AddGene webpages, using
  the *full depositor* sequences. Missing sequences (*DVA_AE* and *DVA_AF*)
  were obtained by editing the overhangs of *DVA_EF*.

In case of mismatch between the two sources, the Zip sequences were used
preferably.

Plasmids were rotated to share the same origin, using the start of the
*BioBrick* prefix as a reference location. This ensures no feature overlaps
the zero coordinate, which was the case beforehand, to ensure a complete
Biopython compatibility.

Common features were colored using the same palette as in the `Yeast ToolKit`.
*AmpR* and *KanR* received additional cross-references from Swiss-Prot, and
a *AmpR* and *KanR* promoter features were added, based on the sequence of
YTK parts. Promoters, RBS and terminators feature types were changed from
*misc_feature* to the expected type.

See Also:
    The annotation script running on Python 3 in the `GitHub repository
    <https://github.com/althonos/moclo/blob/master/scripts/cidar_registry.py>_`

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import re

import six

from ..kits import cidar
from .base import EmbeddedRegistry
from ._utils import find_resistance


class CIDARRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "cidar.tar.gz"

    _ENTITY_RX = re.compile(r"MoClo (.*): ([^\-\[\(]*)")

    _CLASSES = {
        "Cassette Vector": cidar.CIDARCassetteVector,
        "Entry Vector": cidar.CIDAREntryVector,
        "Device": cidar.CIDARDevice,
        "Transcriptional Unit": cidar.CIDARCassette,
        "Basic Part": cidar.CIDARPart,
    }

    _TYPES = {
        "Double terminator": cidar.CIDARTerminator,
        "RBS": cidar.CIDARRibosomeBindingSite,
        "CDS": cidar.CIDARCodingSequence,
        "Controllable promoter": cidar.CIDARPromoter,
        "Constitutive promoter": cidar.CIDARPromoter,
    }

    def _load_entity(self, record):
        match = self._ENTITY_RX.match(record.description)
        if match is not None:
            class_, type_ = map(six.text_type.strip, match.groups())
            if class_ == "Destination Vector":
                if record.id.startswith("DVA"):
                    class_ = "Entry Vector"
                elif record.id.startswith("DVK"):
                    class_ = "Cassette Vector"
            if class_ != "Basic Part":
                return self._CLASSES[class_](record)
            if type_ in self._TYPES:
                return self._TYPES[type_](record)
        raise RuntimeError("could not find type of '{}'".format(record.id))
