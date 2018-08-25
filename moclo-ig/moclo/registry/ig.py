# coding: utf-8
"""Icon Genetics ToolKit sequences registry.
"""
from __future__ import absolute_import
from __future__ import unicode_literals

import re

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
        # TODO: sanitize sequences in `moclo-ig/registry/ig` to use the
        #       right sequence labels instead of faulty ones
        #
        # try:
        #     return find_resistance(raw['record'])
        # except RuntimeError:
        #     msg = "could not find antibiotics resistance of '{}'"
        #     six.raise_from(RuntimeError(msg.format(raw['record'].id)), None)
        # return raw['resistance']
        #
        labels = [
            label
            for f in raw["record"].features
            for label in f.qualifiers.get("label", [])
        ]
        if "AP r" in labels or "AP\\r" in labels or "AP(R)" in labels:
            return "Ampicillin"
        elif (
            "Sm/Sp\\no\\DraIII" in labels
            or "Sm/Sp" in labels
            or "spec" in labels
            or "spec\orf?" in labels
        ):
            return "Spectinomycin"
        elif "NPTII" in labels or "Kan\(no\BpiI)" in labels:
            return "Kanamycin"
        else:
            raise RuntimeError("antibio of " + raw["record"].name)

    # _ENTITY_RX = re.compile(r'MoClo (.*): ([^\-\[\(]*)')
    #
    # _CLASSES = {
    #     'Cassette Vector': cidar.CIDARCassetteVector,
    #     'Entry Vector': cidar.CIDAREntryVector,
    #     'Device': cidar.CIDARDevice,
    #     'Transcriptional Unit': cidar.CIDARCassette,
    #     'Basic Part': cidar.CIDARPart,
    # }
    #
    # _TYPES = {
    #     'Double terminator': cidar.CIDARTerminator,
    #     'RBS': cidar.CIDARRibosomeBindingSite,
    #     'CDS': cidar.CIDARCodingSequence,
    #     'Controllable promoter': cidar.CIDARPromoter,
    #     'Constitutive promoter': cidar.CIDARPromoter,
    # }

    def _load_entity(self, raw, index):

        types = [
            ig.IGLevelPEndLinker,
            ig.IGLevelMEndLinker,
            ig.IGLevelPVector,
            ig.IGLevelMVector,
            ig.IGEndLinker,
            ig.IGTerminator,
            ig.IGCodingSequence,
            ig.IGSignalPeptide,
            ig.IGUntranslatedRegion,
            ig.IGPromoter,
            ig.IGCassette,
            ig.IGCassetteVector,
            ig.IGSingleCassetteVector,
            ig.IGDeviceVector,
        ]

        for t in types:
            entity = t(raw["record"])
            if entity.is_valid():
                return entity
        raise RuntimeError("entity of " + raw["record"].name)

        # match = self._ENTITY_RX.match(raw['record'].description)
        # if match is not None:
        #     class_, type_ = map(six.text_type.strip, match.groups())
        #     if class_ == 'Destination Vector':
        #         if raw['id'].startswith('DVA'):
        #             class_ = 'Entry Vector'
        #         elif raw['id'].startswith('DVK'):
        #             class_ = 'Cassette Vector'
        #     if class_ != 'Basic Part':
        #         return self._CLASSES[class_](raw['record'])
        #     if type_ in self._TYPES:
        #         return self._TYPES[type_](raw['record'])
        # raise RuntimeError("could not find type of '{}'".format(raw['id']))
