# coding: utf-8
"""CIDAR ToolKit sequences registry.

Sequences were obtained from the AddGene depositor full sources. Missing
sequences (*DVA_AE* and *DVA_AF*) were obtained by editing the overhangs of
*DVA_EF*.

Plasmids were rotated to share the same origin, using the the start of the
*BioBrick* prefix as a reference location. This ensures no feature overlaps
the zero coordinate, which was the case beforehand, to ensure a complete
``biopython`` compatibility.

Common features were colored using the same palette as in the `Yeast ToolKit`.
*AmpR* and *KanR* received additional cross-references from Swiss-Prot, and
a *AmpR* and *KanR* promoter features were added, based on the sequence of
YTK parts. Some modules lack annotations of their target sequence; this was
fixed for some promoters (*R0010*, *R0040* and *R0060*).

See Also:
    The annotation script running in Python 3 in the `GitHub repository
    <https://github.com/althonos/moclo/blob/master/scripts/cidar_registry.py>_`

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import re

import six

from ..kits import cidar
from .base import EmbeddedRegistry
from .utils import find_resistance


class CIDARRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "cidar.json.bz2"

    # _types = {
    #     '1': ytk.YTKPart1,
    #     '2': ytk.YTKPart2,
    #     '3': ytk.YTKPart3,
    #     '3a': ytk.YTKPart3a,
    #     '3b': ytk.YTKPart3b,
    #     '4': ytk.YTKPart4,
    #     '4a': ytk.YTKPart4a,
    #     '4b': ytk.YTKPart4b,
    #     '234': ytk.YTKPart234,
    #     '234r': ytk.YTKPart234r,
    #     '5': ytk.YTKPart5,
    #     '6': ytk.YTKPart6,
    #     '7': ytk.YTKPart7,
    #     '8': ytk.YTKPart8,
    #     '8a': ytk.YTKPart8a,
    #     '8b': ytk.YTKPart8b,
    #     '678': ytk.YTKPart678,
    #     'cassette vector': ytk.YTKCassetteVector,
    #     'entry vector': ytk.YTKEntryVector,
    # }

    def __init__(self, location='CIDAR Plate'):
        super(CIDARRegistry, self).__init__()
        self.location = location

    def _load_name(self, raw, index):
        return raw['record'].name

    def _load_id(self, raw, index):
        return raw['record'].id

    def _load_resistance(self, raw, index):
        try:
            return find_resistance(raw['record'])
        except RuntimeError:
            msg = "could not find antibiotics resistance of '{}'"
            six.raise_from(RuntimeError(msg.format(raw['record'].id)), None)
        return raw['resistance']

    _ENTITY_RX = re.compile(r'MoClo (.*): ([^\-\[\(]*)')

    _CLASSES = {
        'Cassette Vector': cidar.CIDARCassetteVector,
        'Entry Vector': cidar.CIDAREntryVector,
        'Device': cidar.CIDARDevice,
        'Transcriptional Unit': cidar.CIDARCassette,
        'Basic Part': cidar.CIDARPart,
    }

    _TYPES = {
        'Double terminator': cidar.CIDARTerminator,
        'RBS': cidar.CIDARRibosomeBindingSite,
        'CDS': cidar.CIDARCodingSequence,
        'Controllable promoter': cidar.CIDARPromoter,
        'Constitutive promoter': cidar.CIDARPromoter,
    }

    def _load_entity(self, raw, index):
        match = self._ENTITY_RX.match(raw['record'].description)
        if match is not None:
            class_, type_ = map(str.strip, match.groups())
            if class_ == 'Destination Vector':
                if raw['id'].startswith('DVA'):
                    class_ = 'Entry Vector'
                elif raw['id'].startswith('DVK'):
                    class_ = 'Cassette Vector'
            if class_ != 'Basic Part':
                return self._CLASSES[class_](raw['record'])
            if type_ in self._TYPES:
                return self._TYPES[type_](raw['record'])
        raise RuntimeError("could not find type of '{}'".format(raw['id']))

    def _load_location(self, raw, index):
        x, y = index % 12, index // 12
        return '{}, {}{}'.format(self.location, chr(ord('A') + y), x+1)
