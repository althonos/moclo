# coding: utf-8
"""MoClo vector classes.

A vector is a plasmidic DNA sequence that can hold a combination of modules of
the same level to create a single module of the following level. Vectors
contain a placeholder sequence that is replaced by the concatenation of the
modules during the Golden Gate assembly.
"""

import abc
import typing
import warnings

import six
from Bio import BiopythonWarning
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .. import errors
from .._utils import catch_warnings, classproperty
from ..record import CircularRecord
from ._utils import cutter_check, add_as_source
from ._structured import StructuredRecord

if typing.TYPE_CHECKING:
    from typing import Any, MutableMapping, Union           # noqa: F401
    from Bio.Restriction.Restriction import RestrictionType # noqa: F401
    from .modules import AbstractModule                     # noqa: F401


class AbstractVector(StructuredRecord):
    """An abstract modular cloning vector.
    """

    _level = None           # type: Union[None, int]
    cutter = NotImplemented # type: Union[NotImplemented, RestrictionType]

    def __new__(cls, *args, **kwargs):
        cutter_check(cls.cutter, name=cls.__name__)
        return super(AbstractVector, cls).__new__(cls)

    @classmethod
    def structure(cls):
        downstream = cls.cutter.elucidate()
        upstream = str(Seq(downstream).reverse_complement())
        return ''.join([
            upstream.replace('^', ')(').replace('_', '('),
            'N*',
            downstream.replace('^', ')(').replace('_', ')')
        ])

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
        if self.cutter.is_3overhang():
            return self._match.group(2) + self.overhang_end()
        else:
            return self.overhang_start() + self._match.group(2)

    def target_sequence(self):
        # type: () -> SeqRecord
        """Get the target sequence in the vector.

        The target sequence if the part of the plasmid that is not discarded
        during the assembly (everything except the placeholder sequence).
        """
        if self.cutter.is_3overhang():
            start, end = self._match.span(2)[0], self._match.span(3)[1]
        else:
            start, end = self._match.span(1)[0], self._match.span(2)[1]
        return add_as_source(self.record, (self.record << start)[end - start:])

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
                during the assembly (mostly caused by duplicate parts).

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

        # Get the alphabet using the largest alphabet from source records
        alphabets = [
            mod.seq.alphabet for mod in six.itervalues(modmap)
            if mod.seq.alphabet.letters is not None
        ]
        alphabets.append(IUPAC.unambiguous_dna)
        alphabet = max(alphabets, key=lambda a: len(a.letters))

        # Generate the complete inserted sequence
        try:
            overhang_next = self.overhang_end()
            assembly = SeqRecord(Seq('', alphabet), id='assembly')
            while overhang_next != self.overhang_start():
                module = modmap.pop(overhang_next)
                assembly += module.target_sequence()
                overhang_next = module.overhang_end()
        except KeyError as ke:
            # Raise the MissingModule error without the KeyError traceback
            raise six.raise_from(errors.MissingModule(ke.args[0]), None)

        # Check all modules were used
        if modmap:
            warnings.warn(errors.UnusedModules(*modmap.values()))

        # Replace placeholder in the vector
        return CircularRecord(assembly + self.target_sequence())


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
