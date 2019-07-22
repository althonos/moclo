# coding: utf-8
"""Moclo module classes.

A module is a sequence of DNA that contains a sequence of interest, such as a
promoter, a CDS, a protein binding site, etc., organised in a way it can be
combined to other modules to create an assembly. This involves flanking that
target sequence with Type IIS restriction sites, which depend on the level of
the module, as well as the chosen MoClo protocol.
"""

import typing

from Bio.Seq import Seq
from property_cached import cached_property

from .. import errors
from ._structured import StructuredRecord
from ._utils import cutter_check, add_as_source

if typing.TYPE_CHECKING:
    from typing import Union, Text  # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401


class AbstractModule(StructuredRecord):
    """An abstract modular cloning module.

    Attributes:
        cutter (`~Bio.Restriction.Restriction.RestrictionType`): the enzyme
            used to cut the target sequence from the backbone plasmid during
            Golden Gate assembly.

    """

    _level = None  # type: Union[None, int]
    cutter = NotImplemented  # type: Union[NotImplemented, RestrictionType]

    def __new__(cls, *args, **kwargs):
        cutter_check(cls.cutter, name=cls.__name__)
        return super(AbstractModule, cls).__new__(cls)

    @classmethod
    def structure(cls):
        # type: () -> Text
        """Get the module structure, as a DNA regex pattern.

        Warning:
            If overloading this method, the returned pattern must include 3
            capture groups to capture the following features:

            1. The upstream (5') overhang sequence
            2. The module target sequence
            3. The downstream (3') overhang sequence

        """
        upstream = cls.cutter.elucidate()
        downstream = str(Seq(upstream).reverse_complement())
        return "".join(
            [
                upstream.replace("^", "(").replace("_", ")("),
                "N*",
                downstream.replace("^", ")").replace("_", ")("),
            ]
        )

    def overhang_start(self):
        # type: () -> Seq
        """Get the upstream overhang of the target sequence.

        Returns:
            `~Bio.Seq.Seq`: the downstream overhang.

        """
        return self._match.group(1).seq

    def overhang_end(self):
        # type: () -> Seq
        """Get the downstream overhang of the target sequence.

        Returns:
            `~Bio.Seq.Seq`: the downstream overhang.

        """
        return self._match.group(3).seq

    def target_sequence(self):
        # type: () -> SeqRecord
        """Get the target sequence of the module.

        Modules are often stored in a standardized way, and contain more than
        the sequence of interest: for instance they can contain an antibiotic
        marker, that will not be part of the assembly when that module is
        assembled into a vector; only the target sequence is inserted.

        Returns:
            `~Bio.SeqRecord.SeqRecord`: the target sequence with annotations.

        Note:
            Depending on the cutting direction of the restriction enzyme used
            during assembly, the overhang will be left at the beginning or at
            the end, so the obtained record is exactly the sequence the enzyme
            created during restriction.

        """
        if self.cutter.is_3overhang():
            start, end = self._match.span(2)[0], self._match.span(3)[1]
        else:
            start, end = self._match.span(1)[0], self._match.span(2)[1]
        return add_as_source(self.record, (self.record << start)[: end - start])

    @cached_property
    def _match(self):
        _match = super(AbstractModule, self)._match
        if len(self.cutter.catalyse(_match.group(0).seq)) > 3:
            raise errors.IllegalSite(self.seq)
        return _match


class Product(AbstractModule):
    """A level -1 module, often obtained as a PCR product.

    Modules of this level are the lowest components of the MoClo system, but
    are not practical to work with until they are assembled in a standard
    vector to obtain *entries*.

    """

    _level = -1


class Entry(AbstractModule):
    """A level 0 module, often obtained from the official toolkits plamisds.

    Entries are assembled from products into a standard vector suitable for
    selection and storage.

    """

    _level = 0


class Cassette(AbstractModule):
    """A level 1 module, also refered as a Transcriptional Unit.

    Cassettes can either express genes in their target organism, or be
    assembled into *multigene* modules for expressing many genes at once,
    depending on the chosen cassette vector during level 0 assembly.

    """

    _level = 1


class Device(AbstractModule):
    """A level 2 module, also refered as a Multigene plasmid.

    Modules of this level are assembled from several transcriptional units so
    that they contain several genes that can be expressed all at once. Most of
    the MoClo implementations are designed so that multiple devices can be
    assembled into a module that is also a valid level 1 module, as does the
    **Golden Braid** system with its **α** and **Ω** plasmids.

    """

    _level = 2
