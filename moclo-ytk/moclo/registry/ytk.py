# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six

from ..kits import ytk
from .base import EmbeddedRegistry


class YTKRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ytk.json.bz2"

    _antibio = {
        'CmR': 'Chloramphenicol',
        'AmpR': 'Ampicillin',
        'KanR': 'Kanamycin',
        'SmR': 'Spectinomycin',
    }

    _types = {
        '1': ytk.YTKPart1,
        '2': ytk.YTKPart2,
        '3': ytk.YTKPart3,
        '3a': ytk.YTKPart3a,
        '3b': ytk.YTKPart3b,
        '4': ytk.YTKPart4,
        '4a': ytk.YTKPart4a,
        '4b': ytk.YTKPart4b,
        '234': ytk.YTKPart234,
        '234r': ytk.YTKPart234r,
        '5': ytk.YTKPart5,
        '6': ytk.YTKPart6,
        '7': ytk.YTKPart7,
        '8': ytk.YTKPart8,
        '8a': ytk.YTKPart8a,
        '8b': ytk.YTKPart8b,
        '678': ytk.YTKPart678,
        'cassette vector': ytk.YTKCassetteVector,
        'entry vector': ytk.YTKEntryVector,
    }

    def __init__(self, location='YTK Plate'):
        super(YTKRegistry, self).__init__()
        self.location = location

    def _load_name(self, raw, index):
        return raw['record'].description

    def _load_id(self, raw, index):
        return raw['record'].id

    def _load_resistance(self, raw, index):
        try:
            gene = next(
                f.qualifiers['label'][0] for f in raw['record'].features
                if any(a in f.qualifiers.get('label', []) for a in self._antibio)
            )
            return self._antibio[gene]
        except StopIteration:
            msg = "could not find antibiotics resistance of '{}'"
            six.raise_from(RuntimeError(msg.format(raw['record'].id)), None)
        return raw['resistance']

    def _load_entity(self, raw, index):
        comments = raw['record'].annotations['comment'].splitlines()
        hint = next(c for c in comments if c.startswith('YTK:'))
        comments.remove(hint)
        raw['record'].annotations['comment'] = "\n".join(comments)
        _, type_ = hint.strip().split(':', 1)
        return self._types[type_](raw['record'])

    def _load_location(self, raw, index):
        x, y = index % 12, index // 12
        return '{}, {}{}'.format(self.location, chr(ord('A') + y), x+1)


class PTKRegistry(YTKRegistry):

    _file = 'ptk.inv'

    def __init__(self, location='PTK Plate'):
        super(PTKRegistry, self).__init__(location)
