# coding: utf-8

import abc
import warnings

import cached_property
import six
import typing
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .. import errors
from ..record import CircularRecord
from ..utils import catch_warnings
from ._structured import StructuredRecord

if typing.TYPE_CHECKING:
    from typing import Any, MutableMapping, Union
    from .modules import AbstractModule


class AbstractVector(StructuredRecord):
    """An abstract modular cloning vector.
    """
    _level = None  # type: Union[None, int]

    def overhang_start(self):
        # type: () -> Seq
        return self._match.group(3).seq

    def overhang_end(self):
        # type: () -> Seq
        return self._match.group(1).seq

    def placeholder_sequence(self):
        # type: () -> SeqRecord
        return self._match.group(2)

    @catch_warnings('ignore', category=BiopythonWarning)
    def assemble(self, module, *modules, **kwargs):
        # type: (AbstractModule, *AbstractModule, **Any) -> SeqRecord

        # If the start and end overhangs are the same, the assembly will
        # not be the only stable product in the bioreactor.
        if self.overhang_start() == self.overhang_end():
            details = 'vector is not suitable for assembly'
            raise errors.InvalidSequence(self, details=details)

        # Identify all modules by their respective overhangs, checking
        # for possible duplicates.
        modmap = {module.overhang_start(): module}  # type: MutableMapping[Seq, AbstractModule]
        for mod in modules:
            mod2 = modmap.setdefault(mod.overhang_start(), mod)
            if mod2 is not mod:
                details = "same start overhang: '{}'".format(mod.overhang_start())
                raise errors.DuplicateModules(mod2, mod, details=details)

        # Generate the complete inserted sequence
        try:
            overhang_next = self.overhang_end()
            assembly = SeqRecord(overhang_next, id='assembly')
            while overhang_next != self.overhang_start():
                module = modmap.pop(overhang_next)
                assembly += module.target_sequence()
                overhang_next = module.overhang_end()
                if overhang_next != self.overhang_end():
                    assembly += overhang_next
        except KeyError as ke:
            # Raise the MissingModule error without the KeyError traceback
            raise six.raise_from(errors.MissingModule(ke.args[0]), None)

        # Check all modules were used
        if modmap:
            warnings.warn(errors.UnusedModules(*modmap.values()))

        # Replace placeholder in the vector while keeping annotations
        ph_start, ph_end = self._match.span(0)
        rec = (self.record << ph_start)
        return CircularRecord(assembly + rec[ph_end - ph_start:])


class EntryVector(AbstractVector):
    """Level 0 vector.
    """

    _level = 0


class CassetteVector(AbstractVector):
    """Level 1 vector.
    """

    _level = 1


class MultigeneVector(AbstractVector):
    """Level 2 vector.
    """

    _level = 2
