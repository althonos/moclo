# coding: utf-8
"""Moclo part classes.
"""

import abc
import typing

import six
from Bio.Seq import Seq

from ..utils import classproperty
from .modules import AbstractModule
from .vectors import AbstractVector
from ._structured import StructuredRecord
from ._utils import cutter_check

if typing.TYPE_CHECKING:
    from typing import Union             # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401


@six.add_metaclass(abc.ABCMeta)
class AbstractPart(StructuredRecord):

    cutter = NotImplemented
    signature = NotImplemented

    def __new__(cls, *args, **kwargs):
        cutter_check(cls.cutter, name=cls.__name__)
        return super(AbstractPart, cls).__new__(cls)

    @classproperty
    def _structure(cls):
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
