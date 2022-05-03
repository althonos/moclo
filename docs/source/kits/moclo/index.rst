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

.. autoclass:: MoCloPromoter(MoCloPart, MoCloEntry)
.. autoclass:: MoCloUntranslatedRegion(MoCloPart, MoCloEntry)
.. autoclass:: MoCloSignalPeptide(MoCloPart, MoCloEntry)
.. autoclass:: MoCloCodingSequence(MoCloPart, MoCloEntry)
.. autoclass:: MoCloTerminator(MoCloPart, MoCloEntry)


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
