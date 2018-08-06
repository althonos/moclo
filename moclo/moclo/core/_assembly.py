# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import copy
import re
import itertools
import warnings

import six
from Bio import BiopythonWarning
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .. import errors, __version__
from ..record import CircularRecord
from .._utils import catch_warnings


class AssemblyManager(object):

    def __init__(self, vector, modules, id_="assembly", name="assembly"):

        if vector.overhang_start() == vector.overhang_end():
            details = 'vector is not suitable for assembly'
            raise errors.InvalidSequence(vector, details=details)

        self.vector = vector
        self.modules = modules
        self.elements = modules + [vector]
        self.name = name
        self.id = id_

    @catch_warnings('ignore', category=BiopythonWarning)
    def assemble(self):
        modmap = self._generate_modules_map()
        alphabet = self._get_alphabet(modmap)

        for elem in self.elements:
            self._deref_citations(elem.record)

        assembly = self._generate_assembly(modmap, alphabet)

        self._annotate_assembly(assembly)
        self._ref_citations(assembly)
        for elem in self.elements:
            self._ref_citations(elem.record)

        return assembly

    def _generate_modules_map(self):
        modmap = {}
        for mod in self.modules:
            m = modmap.setdefault(mod.overhang_start(), mod)
            if m is not mod:
                details = "same start overhang: '{}'".format(mod.overhang_start())
                raise errors.DuplicateModules(m, mod, details=details)
        return modmap

    def _get_alphabet(self, modmap):
        alphabets = [
            mod.seq.alphabet for mod in six.itervalues(modmap)
            if mod.seq.alphabet.letters is not None
        ]
        alphabets.append(IUPAC.unambiguous_dna)
        return max(alphabets, key=lambda a: len(a.letters))

    def _generate_assembly(self, modmap, alphabet):
        try:
            overhang_next = self.vector.overhang_end()
            assembly = SeqRecord(Seq('', alphabet))
            while overhang_next != self.vector.overhang_start():
                module = modmap.pop(overhang_next)
                assembly += module.target_sequence()
                overhang_next = module.overhang_end()
        except KeyError as ke:
            raise six.raise_from(errors.MissingModule(ke.args[0]), None)
        if modmap:
            warnings.warn(errors.UnusedModules(*modmap.values()))
        return CircularRecord(assembly + self.vector.target_sequence())

    _CITATION_RX = re.compile(r'\[(\d*)\]')

    def _deref_citations(self, record):
        references = record.annotations.get('references', [])
        for feature in record.features:
            for i, ref in enumerate(feature.qualifiers.get('citation', [])):
                match = self._CITATION_RX.match(ref)
                if match is None:
                    raise ValueError("invalid citation: '{}'".format(ref))
                ref_index = int(match.group(1)) - 1
                feature.qualifiers['citation'][i] = references[ref_index]

    def _ref_citations(self, record):
        references = record.annotations.setdefault('references', [])
        for feature in record.features:
            for i, ref in enumerate(feature.qualifiers.get('citation', [])):
                if ref not in references:
                    references.append(ref)
                ref_index = references.find(ref) + 1
                feature.qualifiers['citation'][i] = "{}".format(ref_index)

    def _annotate_assembly(self, assembly):
        assembly.id = self.id
        assembly.name = self.name
        ants = assembly.annotations

        ants['topology'] = 'circular'
        ants['organism'] = ants['source'] = 'synthetic DNA construct'
        ants['source'] = 'synthetic DNA construct'
        ants['molecule_type'] = 'ds-DNA'
        ants['data_file_division'] = 'SYN'
        ants['comment'] = [
            'Generated with moclo v{}'.format(__version__),
            'Vector: {}'.format(self.vector.record.id),
            'Modules: {}'.format(', '.join(
                mod.record.id for mod in self.modules
            )),
        ]
