# coding: utf-8

import typing

if typing.TYPE_CHECKING:
    from typing import Any
    from Bio.Seq import Seq
    from .base.modules import AbstractModule


class MocloError(Exception):
    pass

class InvalidSequence(ValueError, MocloError):
    """Invalid sequence provided.
    """

    def __init__(self, sequence, exc=None, details=None):
        self.sequence = sequence
        self.exc = exc
        self.details = details

    def __str__(self):
        s = 'invalid sequence: {}'
        if self.details is not None:
            s = ''.join([s, ' ', '(', self.details, ')'])
        return s.format(self.sequence)


class AssemblyError(MocloError, RuntimeError):
    """Assembly-specific run-time error.
    """


class DuplicateModules(AssemblyError):
    """Several modules share the same overhangs.
    """

    def __init__(self, *duplicates, **options):
        # type: (*AbstractModule, **Any) -> None
        self.duplicates = duplicates
        self.details = options.pop('details', None)

    def __str__(self):
        s = 'duplicate modules: {}'
        if self.details is not None:
            s = ''.join([s, ' ', '(', self.details, ')'])
        return s.format(', '.join(self.duplicates))


class MissingModule(AssemblyError):
    """A module is missing in the assembly.
    """

    def __init__(self, start_overhang, **options):
        self.start_overhang = start_overhang
        self.details = options.pop('details', None)

    def __str__(self):
        s = "no module with '{}' start overhang"
        if self.details is not None:
            s = ''.join([s, ' ', '(', self.details, ')'])
        return s.format(self.start_overhang)


class AssemblyWarning(MocloError, Warning):
    """Assembly-specific run-time warning.
    """


class UnusedModules(AssemblyWarning):
    """Not all modules were used during assembly.
    """

    def __init__(self, *remaining, **options):
        self.remaining = remaining
        self.details = options.pop('details', None)

    def __str__(self):
        s = 'unused: {}'
        if self.details is not None:
            s = ''.join([s, ' ', '(', self.details, ')'])
        return s.format(', '.join(map(str, self.remaining)))
