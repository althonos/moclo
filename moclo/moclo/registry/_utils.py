# coding: utf-8

from .._utils import isabstract

_ANTIBIOTICS = {
    "KanR": "Kanamycin",
    "CamR": "Chloramphenicol",
    "CmR": "Chloramphenicol",
    "KnR": "Kanamycin",
    "AmpR": "Ampicillin",
    "SmR": "Spectinomycin",
    "SpecR": "Spectinomycin",
}


def find_resistance(record):
    """Infer the antibiotics resistance of the given record.

    Arguments:
        record (`~Bio.SeqRecord.SeqRecord`): an annotated sequence.

    Raises:
        RuntimeError: when there's not exactly one resistance cassette.

    """
    for feature in record.features:
        labels = set(feature.qualifiers.get("label", []))
        cassettes = labels.intersection(_ANTIBIOTICS)
        if len(cassettes) > 1:
            raise RuntimeError("multiple resistance cassettes detected")
        elif len(cassettes) == 1:
            return _ANTIBIOTICS.get(cassettes.pop())
    raise RuntimeError("could not find the resistance of '{}'".format(record.id))
