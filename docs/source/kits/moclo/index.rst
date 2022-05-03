MoClo Kit
=========

.. currentmodule:: moclo.kits.moclo

.. automodule:: moclo.kits.moclo

Level -1
--------

Module
^^^^^^
.. autoclass:: MoCloProduct(Product)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: MoCloEntryVector(EntryVector)
   :members:
   :inherited-members:
   :special-members: __init__


Level 0
-------

Module
^^^^^^
.. autoclass:: MoCloEntry(Entry)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: MoCloCassetteVector(CassetteVector)
   :members:
   :inherited-members:
   :special-members: __init__


Parts
^^^^^

.. autoclass:: MoCloPro(MoCloPart, MoCloEntry)
.. autoclass:: MoClo5U(MoCloPart, MoCloEntry)
.. autoclass:: MoClo5Uf(MoCloPart, MoCloEntry)
.. autoclass:: MoCloNTag(MoCloPart, MoCloEntry)
.. autoclass:: MoCloPro5U(MoCloPart, MoCloEntry)
.. autoclass:: MoCloPro5Uf(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCDS1(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCDS1ns(MoCloPart, MoCloEntry)
.. autoclass:: MoCloSP(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCDS2(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCDS2ns(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCTag(MoCloPart, MoCloEntry)
.. autoclass:: MoClo3U(MoCloPart, MoCloEntry)
.. autoclass:: MoCloTer(MoCloPart, MoCloEntry)
.. autoclass:: MoClo3UTer(MoCloPart, MoCloEntry)
.. autoclass:: MoCloGene(MoCloPart, MoCloEntry)


Level 1
-------

Module
^^^^^^
.. autoclass:: MoCloCassette(Cassette)
   :members:

Vector
^^^^^^
.. autoclass:: MoCloDeviceVector(DeviceVector)
   :members:

Parts
^^^^^
.. autoclass:: MoCloEndLinker(MoCloPart, MoCloCassette)
   :members:


Level M
-------

Parts
^^^^^
.. autoclass:: MoCloLevelMVector(MoCloPart, MoCloDeviceVector)
   :members:

.. autoclass:: MoCloLevelMEndLinker(MoCloPart, MoCloCassette)
   :members:


Level P
-------

Parts
^^^^^
.. autoclass:: MoCloLevelPVector(MoCloPart, MoCloCassetteVector)
   :members:

.. autoclass:: MoCloLevelPEndLinker(MoCloPart, MoCloEntry)
   :members:
