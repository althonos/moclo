# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from ..kits import ytk
from .base import EmbeddedRegistry


class YTKRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ytk.inv"
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
        return 'pYTK{:03}'.format(index+1)

    def _load_location(self, raw, index):
        x, y = index % 12, index // 12
        return '{}, {}{}'.format(self.location, chr(ord('A') + y), x+1)
