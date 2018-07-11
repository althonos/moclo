# coding: utf-8
"""Moclo module classes.

A module is a sequence of DNA that contains a sequence of interest, such as a
promoter, a CDS, a protein binding site, etc., organised in a way it can be
combined to other modules to create an assembly. This involves flanking that
target sequence with Type IIS restriction sites, which depend on the level of
the module, as well as the chosen MoClo protocol.
"""

import typing

from ._structured import StructuredRecord

if typing.TYPE_CHECKING:
    from typing import Union             # noqa: F401
    from Bio.Seq import Seq              # noqa: F401
    from Bio.SeqRecord import SeqRecord  # noqa: F401


class AbstractModule(StructuredRecord):
    """An abstract modular cloning module.
    """

    _level = None  # type: Union[None, int]

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

        Danger:
            The start and end overhangs are not included in the returned
            sequence.

        """
        return self._match.group(2)


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
    the MoClo implementations are desinged so that multiple devices can be
    assembled into a module that is also a valid level 1 module, as does the
    **Golden Braid** system with its **α** and **Ω** plasmids.

    """

    _level = 2
