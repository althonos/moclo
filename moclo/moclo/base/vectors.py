# coding: utf-8
"""MoClo vector classes.

A vector is a plasmidic DNA sequence that can hold a combination of modules of
the same level to create a single module of the following level. Vectors
contain a placeholder sequence that is replaced by the concatenation of the
modules during the Golden Gate assembly.
"""

import warnings

import six
import typing
from Bio import BiopythonWarning
from Bio.SeqRecord import SeqRecord

from .. import errors
from ..record import CircularRecord
from ..utils import catch_warnings
from ._structured import StructuredRecord

if typing.TYPE_CHECKING:
    from typing import Any, MutableMapping, Union   # noqa: F401
    from Bio.Seq import Seq                         # noqa: F401
    from .modules import AbstractModule             # noqa: F401


class AbstractVector(StructuredRecord):
    """An abstract modular cloning vector.
    """

    _level = None  # type: Union[None, int]

    def overhang_start(self):
        # type: () -> Seq
        """Get the upstream overhang of the vector sequence.
        """
        return self._match.group(3).seq

    def overhang_end(self):
        # type: () -> Seq
        """Get the downstream overhang of the vector sequence.
        """
        return self._match.group(1).seq

    def placeholder_sequence(self):
        # type: () -> SeqRecord
        """Get the placeholder sequence in the vector.

        The placeholder sequence is replaced by the concatenation of modules
        during the assembly. It often contains a dropout sequence, such as a
        GFP expression cassette that can be used to measure the progress of
        the assembly.
        """
        return self._match.group(2)

    @catch_warnings('ignore', category=BiopythonWarning)
    def assemble(self, module, *modules, **kwargs):
        # type: (AbstractModule, *AbstractModule, **Any) -> SeqRecord
        """Assemble the provided modules into the vector.

        Arguments:
            module (`~moclo.base.modules.AbstractModule`): a module to insert
                in the vector.
            modules (`~moclo.base.modules.AbstractModule`, optional): additional
                modules to insert in the vector. The order of the parameters
                is not important, since modules will be sorted by their start
                overhang in the function.

        Returns:
            `~Bio.SeqRecord.SeqRecord`: the assembled sequence with sequence
            annotations inherited from the vector and the modules.

        Raises:
            `~moclo.errors.DuplicateModules`: when two different modules share
                the same start overhang, leading in possibly non-deterministic
                constructs.
            `~moclo.errors.MissingModule`: when a module has an end overhang
                that is not shared by any other module, leading to a partial
                construct only
            `~moclo.errors.InvalidSequence`: when one of the modules does not
                match the required module structure (missing site, wrong
                overhang, etc.).
            `~moclo.errors.UnusedModules`: when some modules were not used
                during the assembly.

        """
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


class DeviceVector(AbstractVector):
    """Level 2 vector.
    """

    _level = 2
