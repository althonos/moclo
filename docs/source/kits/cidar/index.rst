CIDAR Kit
=========

.. currentmodule:: moclo.kits.cidar

.. automodule:: moclo.kits.cidar

Level -1
--------

Module
^^^^^^
.. autoclass:: CIDARProduct(Product)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: CIDAREntryVector(EntryVector)
   :members:
   :inherited-members:
   :special-members: __init__


Level 0
-------

Module
^^^^^^
.. autoclass:: CIDAREntry(Entry)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: CIDARCassetteVector(CassetteVector)
   :members:
   :inherited-members:
   :special-members: __init__


Parts
^^^^^

.. autoclass:: CIDARPromoter(CIDARPart, CIDAREntry)
.. autoclass:: CIDARRibosomeBindingSite(CIDARPart, CIDAREntry)
.. autoclass:: CIDARCodingSequence(CIDARPart, CIDAREntry)
.. autoclass:: CIDARTerminator(CIDARPart, CIDAREntry)



Level 1
-------

Module
^^^^^^
.. autoclass:: CIDARCassette(Cassette)
   :members:

Vector
^^^^^^
.. autoclass:: CIDARDeviceVector(DeviceVector)
   :members:


Level 2
-------

Module
^^^^^^

.. autoclass:: CIDARDevice(Device)
   :members:
