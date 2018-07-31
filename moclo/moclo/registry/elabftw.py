# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import Bio.SeqIO
import six

from ..core import AbstractModule, AbstractVector, AbstractPart
from .._utils import catch_warnings
from .._impl import json, ssl
from .base import AbstractRegistry, Item
from ._utils import find_resistance, find_type


class ELabFTWRegistry(AbstractRegistry):
    """A registry in a running ``eLabFTW`` server.
    """

    def __init__(self,
                 server,
                 token,
                 base,
                 category="Plasmids",      # type: Text
                 strict_ssl=False,         # type: bool
                ):
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
                ``eLabFTW`` item inventory. [default: ``Plasmids``]
            strict_ssl (bool): set to `True` to enforce verification of the SSL
                certificate when connecting to the server [default: `False`].

        """

        bases = (AbstractPart, AbstractModule, AbstractVector)
        if not isinstance(base, type) or not issubclass(base, (bases)):
            raise TypeError("base cannot be '{}'".format(base))

        self.base = base
        self.server = server
        self.token = token
        self.category = category
        self._strict = strict_ssl

    def _request(self, url):

        req = six.moves.urllib.request.Request(url)
        req.add_header('Authorization', self.token)

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

        for attachment in plasmid.get('uploads', []):
            url = '{}/uploads/{}'.format(self.server, attachment['long_name'])
            try:
                with self._request(url) as res:
                    data = six.StringIO(res.read().decode('utf-8'))
                    record = Bio.SeqIO.read(data, 'genbank')
            except (ValueError, UnicodeDecodeError) as err:
                continue
            else:
                break
        else:
            raise RuntimeError("no GenBank file for: '{}'".format(item_id))

        entity = find_type(record, self.base)
        resistance = find_resistance(record)

        return Item(
            entity=entity,
            id=plasmid['title'],
            name=plasmid['title'],
            resistance=resistance,
        )

    def __len__(self):
        return sum(
            1
            for item in self._get_all_items()
            if item['category'] == self.category
        )

    def __getitem__(self, key):
        for item in self._get_all_items():
            if item['category'] == self.category and item['title'] == key:
                return self._item_from_id(item['id'])
        else:
            raise KeyError(key)

    def __iter__(self):
        for item in self._get_all_items():
            if item['category'] == self.category:
                yield item['title']
