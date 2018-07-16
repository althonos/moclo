# coding: utf-8
"""Moclo part classes.
"""

import abc
import typing

import six
from Bio.Seq import Seq

from .._utils import classproperty
from .modules import AbstractModule
from .vectors import AbstractVector
from ._structured import StructuredRecord
from ._utils import cutter_check

if typing.TYPE_CHECKING:
    from typing import Union             # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401


@six.add_metaclass(abc.ABCMeta)
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

        if cls.signature is NotImplemented:
            raise NotImplementedError('no signature defined')

        up = cls.cutter.elucidate()
        down = str(Seq(up).reverse_complement())
        ovhg = cls.cutter.ovhgseq
        upsig, downsig = cls.signature

        if cls.cutter.is_5overhang():
            upsite = '^{}_'.format(ovhg)
            downsite = '_{}^'.format(Seq(ovhg).reverse_complement())
        else:
            upsite = '_{}^'.format(ovhg)
            downsite = '^{}_'.format(Seq(ovhg).reverse_complement())

        if issubclass(cls, AbstractModule):
            return ''.join([
                up.replace(upsite, '({})('.format(upsig)),
                'N*',
                down.replace(downsite, ')({})'.format(downsig))
            ])
        elif issubclass(cls, AbstractVector):
            return ''.join([
                down.replace(downsite, '({})('.format(downsig)),
                'N*',
                up.replace(upsite, ')({})'.format(upsig))
            ])
        else:
            raise RuntimeError("Part must be either a module or a vector!")
