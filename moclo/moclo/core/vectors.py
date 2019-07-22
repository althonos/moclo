# coding: utf-8
"""MoClo vector classes.

A vector is a plasmidic DNA sequence that can hold a combination of modules of
the same level to create a single module of the following level. Vectors
contain a placeholder sequence that is replaced by the concatenation of the
modules during the Golden Gate assembly.
"""

import typing

from Bio.Seq import Seq
from property_cached import cached_property

from .. import errors
from ._assembly import AssemblyManager
from ._utils import cutter_check, add_as_source
from ._structured import StructuredRecord

if typing.TYPE_CHECKING:
    from typing import Any, MutableMapping, Union  # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401
    from Bio.Restriction.Restriction import RestrictionType  # noqa: F401
    from .modules import AbstractModule  # noqa: F401


class AbstractVector(StructuredRecord):
    """An abstract modular cloning vector.
    """

    _level = None  # type: Union[None, int]
    cutter = NotImplemented  # type: Union[NotImplemented, RestrictionType]

    def __new__(cls, *args, **kwargs):
        cutter_check(cls.cutter, name=cls.__name__)
        return super(AbstractVector, cls).__new__(cls)

    @classmethod
    def structure(cls):
        # type: () -> Text
        """Get the vector structure, as a DNA regex pattern.

        Warning:
            If overloading this method, the returned pattern must include 3
            capture groups to capture the following features:

            1. The downstream (3') overhang sequence
            2. The vector placeholder sequence
            3. The upstream (5') overhang sequence

        """
        downstream = cls.cutter.elucidate()
        upstream = str(Seq(downstream).reverse_complement())
        return "".join(
            [
                upstream.replace("^", ")(").replace("_", "("),
                "N*",
                downstream.replace("^", ")(").replace("_", ")"),
            ]
        )

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
        return add_as_source(self.record, (self.record << start)[end - start :])

    @cached_property
    def _match(self):
        _match = super(AbstractVector, self)._match
        if len(self.cutter.catalyse(_match.group(0).seq)) > 3:
            raise errors.IllegalSite(self.seq)
        return _match

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

        mgr = AssemblyManager(
            vector=self,
            modules=[module] + list(modules),
            name=kwargs.get("name", "assembly"),
            id_=kwargs.get("id", "assembly"),
        )
        return mgr.assemble()


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
