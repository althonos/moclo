# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import abc
import io
import itertools
import typing

import Bio.SeqIO
import cached_property
import pkg_resources
import six

from .._impl import bz2, json
from ..record import CircularRecord
from ..core import AbstractModule, AbstractVector


class Item(typing.NamedTuple('Item', [
            ('id', typing.Text),
            ('name', typing.Text),
            ('entity', typing.Union[AbstractModule, AbstractVector]),
            ('location', typing.Text),
            ('resistance', typing.Text),
        ])):
    """A registry item.
    """


class AbstractRegistry(typing.Mapping[typing.Text, Item]):
    """An abstract registry holding MoClo plasmids.
    """


class CombinedRegistry(AbstractRegistry):
    """A registry combining several registries into a single collection.
    """

    def __init__(self):
        self._data = {}

    def __lshift__(self, registry):
        self.add_registry(registry)

    def add_registry(self, registry):
        for item in six.itervalues(registry):
            self._data.setdefault(item.id, item)

    def __contains__(self, item):
        return item in self._data

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)


class EmbeddedRegistry(AbstractRegistry):

    _module = NotImplemented
    _file = NotImplemented
    _types = NotImplemented

    def _load_id(self, raw, index):
        return raw['id']

    def _load_name(self, raw, index):
        return raw['name']

    def _load_location(self, raw, index):
        return raw['location']

    def _load_resistance(self, raw, index):
        return raw['resistance']

    def _load_entity(self, raw, index):
        return self._types[raw['type']](raw['record'])

    @cached_property.cached_property
    def _data(self):
        with pkg_resources.resource_stream(self._module, self._file) as rs:
            with io.TextIOWrapper(bz2.BZ2File(rs)) as decomp:
                raw_data = json.load(decomp)

        for raw in raw_data:
            record = Bio.SeqIO.read(io.StringIO(raw['gb']), 'gb')
            raw['record'] = CircularRecord(record)

        return {
            item.id: item
            for item in (
                Item(
                    id=self._load_id(raw, index),
                    name=self._load_name(raw, index),
                    location=self._load_location(raw, index),
                    resistance=self._load_resistance(raw, index),
                    entity=self._load_entity(raw, index),
                )
                for index, raw in enumerate(raw_data)
            )
        }

    def __len__(self):
        return len(self._data)

    def __getitem__(self, item):
        return self._data[item]

    def __iter__(self):
        return iter(self._data)
