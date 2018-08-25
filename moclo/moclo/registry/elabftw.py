# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import typing

import Bio.SeqIO
import six

from ..core import AbstractModule, AbstractVector, AbstractPart
from ..record import CircularRecord
from .._impl import json, ssl
from .base import AbstractRegistry, Item
from ._utils import find_resistance

if typing.TYPE_CHECKING:
    from typing import Collection, Optional, Text  # noqa: F401


class ELabFTWRegistry(AbstractRegistry):
    """A registry in a running ``eLabFTW`` server.
    """

    def __init__(
        self,
        server,  # type: Text
        token,  # type: Text
        base,  # type: Type[AbstractPart]
        category="Plasmids",  # type: Text
        include_tags=None,  # type: Optional[Collection[Text]]
        exclude_tags=None,  # type: Optional[Collection[Text]]
        strict_ssl=False,  # type: bool
        ignore_unknown=True,  # type: bool
    ):
        # type: (...) -> None
        """Open the registry located on the given ``eLabFTW`` instance.

        Arguments:
            server (str): the URL to the server, with a possible non-standard
                port (e.g. ``https://elab.myserver.org:3418``).
            token (str): an ``eLabFTW`` session token. Tokens can be obtained
                from your profile page in the ``eLabFTW`` Web UI.
            base (type): the base Python type of all inventory items (e.g.
                `~moclo.kits.ytk.YTKPart` if your inventory contains only
                valid Yeast ToolKit parts).

        Keyword Arguments:
            category (str): the item category to filter items with in the
                ``eLabFTW`` item inventory. [default: ``"Plasmids"``]
            strict_ssl (bool): set to `True` to enforce verification of the SSL
                certificate when connecting to the server [default: `False`].
            ignore_unknown (bool): ignore any record which type cannot be
                identified. Set to `False` to raise a `RuntimeError` instead.
                This will also silence entries without attached GenBank files.
                [default: `True`]
            include_tags (list, optional): a list of tags to include items
                from the online inventory with. [default: `None`]
            exclude_tags (list, optional): a list of tags to exclude items
                from the online inventory with. [default: `None`]

        """

        bases = (AbstractPart, AbstractModule, AbstractVector)
        if not isinstance(base, type) or not issubclass(base, (bases)):
            raise TypeError("base cannot be '{}'".format(base))

        if not isinstance(server, six.string_types):
            msg = "server must be a string, not '{}'"
            raise TypeError(msg.format(type(server).__name__))
        elif not server.startswith("http"):
            raise ValueError("bad server URL: '{}'".format(server))

        self.base = base
        self.server = server
        self.token = token
        self.category = category
        self._strict = strict_ssl
        self._ignore_unknown = ignore_unknown

        self._include = set(include_tags) if include_tags is not None else None
        self._exclude = set(exclude_tags) if exclude_tags is not None else None

    def _request(self, url):

        req = six.moves.urllib.request.Request(url)  # noqa: B310
        req.add_header("Authorization", self.token)

        if ssl is None:
            ctx = None
        elif self._strict:
            ctx = ssl.create_default_context()
            ctx.check_hostname = True
            ctx.verify_mode = ssl.CERT_OPTIONAL
        else:
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE

        return six.moves.urllib.request.urlopen(req, context=ctx)

    def _get_all_items(self):
        url = "{}/api/v1/items/".format(self.server)
        with self._request(url) as res:
            data = json.load(res)
        # TODO: error check
        return data

    def _get_item(self, item_id):
        url = "{}/api/v1/items/{}".format(self.server, item_id)
        with self._request(url) as res:
            data = json.load(res)
        # TODO: error check
        return data

    def _item_from_id(self, item_id):

        plasmid = self._get_item(item_id)

        for attachment in plasmid.get("uploads", []):
            url = "{}/uploads/{}".format(self.server, attachment["long_name"])
            try:
                with self._request(url) as res:
                    data = six.StringIO(res.read().decode("utf-8"))
                    record = CircularRecord(Bio.SeqIO.read(data, "genbank"))
            except (ValueError, UnicodeDecodeError):
                continue
            else:
                break
        else:
            raise RuntimeError("no GenBank file for: '{}'".format(item_id))

        entity = self.base.characterize(record)
        resistance = find_resistance(record)

        # When exported with Snapgene, the record does not have any name / id
        if record.name == "Exported":
            record.id = record.name = plasmid["title"]

        return Item(
            entity=entity,
            id=plasmid["title"],
            name=plasmid["title"],
            resistance=resistance,
        )

    def __len__(self):
        return sum(1 for item in self)

    def __length_hint__(self):
        return sum(
            1 for item in self._get_all_items() if item["category"] == self.category
        )

    def __getitem__(self, key):
        for item in self._get_all_items():
            if item["category"] == self.category and item["title"] == key:
                return self._item_from_id(item["id"])
        raise KeyError(key)

    def __iter__(self):
        for item in self._get_all_items():
            tags = item["tags"].split("|") if item["tags"] else []

            if item["category"] != self.category:
                continue
            if self._include is not None and self._include.isdisjoint(tags):
                continue
            if self._exclude is not None and not self._exclude.isdisjoint(tags):
                continue

            try:
                self._item_from_id(item["id"])
                yield item["title"]
            except RuntimeError:
                if not self._ignore_unknown:
                    raise
