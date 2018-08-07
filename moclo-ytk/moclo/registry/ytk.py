# coding: utf-8
"""Yeast ToolKit and Pichia ToolKit sequences registry.

Yeast ToolKit:
    Sequences were obtained from two different AddGene sources:
    * The YTK plasmid files distributed with the kit in a zip archive,
      available under the *Protocol & Resources* tab of the YTK repository.
    * The individual plasmid files from their dedicated AddGene webpages,
      using the *full depositor* sequences.

    Records were then merged using a Python script and ``biopython``; details
    can be found in the history of the ``git`` repository. Some annotations
    were also added to common genes and sequences (such as the Chloramphenicol
    resistance cassette `CmR <https://www.uniprot.org/uniprot/P62577>`_, or
    the `H4 ARS consensus <https://www.yeastgenome.org/locus/ARS209>`_).
    Duplicated features refering to the same location were merged directly,
    and duplicated features with unidentical locations were manually reviewed.

Pichia ToolKit:
    Sequences were obtained from two different AddGene sources:
    * The PTK plasmid files available on each AddGene plasmid webpage as
      *Supplemental Documents*
    * The indiviual plasmid files from their dedicated AddGene webpages,
      using the *full AddGene* sequences.

    Since both sequences did not share the same origins, plasmids were rotated
    for the upstream *BsaI* recognition site to span on nucleotides 3 to 11,
    in order for the PTK plasmids to be organized like the YTK plasmids. Using
    ``biopython``, sequences were annotated using a custom annotation script
    (see `source <https://github.com/althonos/moclo/blob/master/scripts>_`).

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import six

from ..kits import ytk
from .base import EmbeddedRegistry
from ._utils import find_resistance


class YTKRegistry(EmbeddedRegistry):

    _module = __name__
    _file = "ytk.json.bz2"

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

    def _load_name(self, raw, index):
        return raw['record'].description

    def _load_id(self, raw, index):
        return raw['record'].id

    def _load_resistance(self, raw, index):
        try:
            return find_resistance(raw['record'])
        except StopIteration:
            msg = "could not find antibiotics resistance of '{}'"
            six.raise_from(RuntimeError(msg.format(raw['record'].id)), None)

    def _load_entity(self, raw, index):
        comments = raw['record'].annotations['comment'].splitlines()
        hint = next(c for c in comments if c.startswith('YTK:'))
        comments.remove(hint)
        raw['record'].annotations['comment'] = "\n".join(comments)
        _, type_ = hint.strip().split(':', 1)
        return self._types[type_](raw['record'])


class PTKRegistry(YTKRegistry):

    _file = "ptk.json.bz2"
