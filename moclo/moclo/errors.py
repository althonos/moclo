# coding: utf-8

import typing

import six

if typing.TYPE_CHECKING:
    from typing import Any  # noqa: F401
    from Bio.Seq import Seq  # noqa: F401
    from .base.modules import AbstractModule  # noqa: F401


class MocloError(Exception):
    """Base class for all MoClo-related exceptions.
    """

    pass


@six.python_2_unicode_compatible
class InvalidSequence(ValueError, MocloError):
    """Invalid sequence provided.
    """

    _msg = "invalid sequence: {}"

    def __init__(self, sequence, exc=None, details=None):
        self.sequence = sequence
        self.exc = exc
        self.details = details

    def __str__(self):
        s = self._msg
        if self.details is not None:
            s = "".join([s, " ", "(", self.details, ")"])
        return s.format(self.sequence)


class IllegalSite(InvalidSequence):
    """Sequence with illegal site provided.
    """

    _msg = "illegal site in sequence: {}"


class AssemblyError(MocloError, RuntimeError):
    """Assembly-specific run-time error.
    """


@six.python_2_unicode_compatible
class DuplicateModules(AssemblyError):
    """Several modules share the same overhangs.
    """

    def __init__(self, *duplicates, **options):
        # type: (*AbstractModule, **Any) -> None
        self.duplicates = duplicates
        self.details = options.pop("details", None)

    def __str__(self):
        s = "duplicate modules: {}"
        if self.details is not None:
            s = "".join([s, " ", "(", self.details, ")"])
        return s.format(", ".join(d.record.id for d in self.duplicates))


@six.python_2_unicode_compatible
class MissingModule(AssemblyError):
    """A module is missing in the assembly.
    """

    def __init__(self, start_overhang, **options):
        self.start_overhang = start_overhang
        self.details = options.pop("details", None)

    def __str__(self):
        s = "no module with '{}' start overhang"
        if self.details is not None:
            s = "".join([s, " ", "(", self.details, ")"])
        return s.format(self.start_overhang)


class AssemblyWarning(MocloError, Warning):
    """Assembly-specific run-time warning.

    Warnings can be turned into errors using the `warnings.catch_warnings`
    decorator combined to `warnings.simplefilter` with `action` set to
    `"error"`.
    """


@six.python_2_unicode_compatible
class UnusedModules(AssemblyWarning):
    """Not all modules were used during assembly.
    """

    def __init__(self, *remaining, **options):
        self.remaining = remaining
        self.details = options.pop("details", None)

    def __str__(self):
        s = "unused: {}"
        if self.details is not None:
            s = "".join([s, " ", "(", str(self.details), ")"])
        return s.format(", ".join(r.record.id for r in self.remaining))
