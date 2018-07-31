# coding: utf-8

from .._utils import isabstract

_ANTIBIOTICS = {
    'KanR': 'Kanamycin',
    'CamR': 'Chloramphenicol',
    'CmR': 'Chloramphenicol',
    'KnR': 'Kanamycin',
    'AmpR': 'Ampicillin',
    'SmR': 'Spectinomycin',
    'SpecR': 'Spectinomycin',
}


def find_resistance(record):
    """Infer the antibiotics resistance of the given record.

    Arguments:
        record (`~Bio.SeqRecord.SeqRecord`): an annotated sequence.

    Raises:
        RuntimeError: when there's not exactly one resistance cassette.

    """
    for feature in record.features:
        labels = set(feature.qualifiers.get('label', []))
        cassettes = labels.intersection(_ANTIBIOTICS)
        if len(cassettes) > 1:
            raise RuntimeError('multiple resistance cassettes detected')
        elif len(cassettes) == 1:
            return _ANTIBIOTICS.get(cassettes.pop())
    raise RuntimeError("could not find the resistance of '{}'".format(record.id))


def find_type(record, base):
    """Infer the type of the record from the given base type.
    """
    classes = list(base.__subclasses__())
    if not isabstract(base):
        classes.append(base)
    for cls in classes:
        entity = cls(record)
        if entity.is_valid():
            return entity
    raise RuntimeError("could not find the type for '{}'".format(record.id))
