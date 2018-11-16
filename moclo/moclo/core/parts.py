# coding: utf-8
"""Moclo part classes.
"""

import typing

from Bio.Seq import Seq

from .._utils import isabstract
from .modules import AbstractModule
from .vectors import AbstractVector
from ._structured import StructuredRecord
from ._utils import cutter_check

if typing.TYPE_CHECKING:
    from typing import Union  # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401


__all__ = ["AbstractPart"]


class AbstractPart(StructuredRecord):
    """An abstract modular cloning part.

    Parts can be either modules or vectors, but are determined by their
    flanking overhangs sequences, declared in the ``signature`` class
    attribute. The part structure is derived from the part class (module
    of vector), signature, and restriction enzyme.

    Example:
        >>> class ExamplePart(AbstractPart, Entry):
        ...     cutter = BsaI
        ...     signature = ('ATGC', 'ATTC')
        ...
        >>> ExamplePart.structure()
        'GGTCTCN(ATGC)(NN*N)(ATTC)NGAGACC'
        
    """

    cutter = NotImplemented
    signature = NotImplemented

    def __new__(cls, *args, **kwargs):
        cutter_check(cls.cutter, name=cls.__name__)
        return super(AbstractPart, cls).__new__(cls)

    @classmethod
    def structure(cls):
        # type: () -> Text
        """Get the part structure, as a DNA regex pattern.

        The structure of most parts can be obtained automatically from the
        part signature and the restriction enzyme used in the Golden Gate
        assembly.

        Warning:
            If overloading this method, the returned pattern must include 3
            capture groups to capture the following features:

            1. The upstream (5') overhang sequence
            2. The vector placeholder sequence
            3. The downstream (3') overhang sequence

        """

        if cls.signature is NotImplemented:
            raise NotImplementedError("no signature defined")

        up = cls.cutter.elucidate()
        down = str(Seq(up).reverse_complement())
        ovhg = cls.cutter.ovhgseq
        upsig, downsig = cls.signature

        if cls.cutter.is_5overhang():
            upsite = "^{}_".format(ovhg)
            downsite = "_{}^".format(Seq(ovhg).reverse_complement())
        else:
            upsite = "_{}^".format(ovhg)
            downsite = "^{}_".format(Seq(ovhg).reverse_complement())

        if issubclass(cls, AbstractModule):
            return "".join(
                [
                    up.replace(upsite, "({})(".format(upsig)),
                    "N*",
                    down.replace(downsite, ")({})".format(downsig)),
                ]
            )
        elif issubclass(cls, AbstractVector):
            return "".join(
                [
                    down.replace(downsite, "({})(".format(downsig)),
                    "N*",
                    up.replace(upsite, ")({})".format(upsig)),
                ]
            )
        else:
            raise RuntimeError("Part must be either a module or a vector!")

    @classmethod
    def characterize(cls, record):
        """Load the record in a concrete subclass of this type.
        """
        classes = list(cls.__subclasses__())
        if not isabstract(cls):
            classes.append(cls)
        for subclass in classes:
            entity = subclass(record)
            if entity.is_valid():
                return entity
        raise RuntimeError("could not find the type for '{}'".format(record.id))
