Icon Genetics Kit
=================

.. currentmodule:: moclo.kits.ig

.. automodule:: moclo.kits.ig

Level -1
--------

Module
^^^^^^
.. autoclass:: IGProduct(Product)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: IGEntryVector(EntryVector)
   :members:
   :inherited-members:
   :special-members: __init__


Level 0
-------

Module
^^^^^^
.. autoclass:: IGEntry(Entry)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: IGCassetteVector(CassetteVector)
   :members:
   :inherited-members:
   :special-members: __init__


Parts
^^^^^

.. autoclass:: IGPromoter(IGPart, IGEntry)
.. autoclass:: IGUntranslatedRegion(IGPart, IGEntry)
.. autoclass:: IGSignalPeptide(IGPart, IGEntry)
.. autoclass:: IGCodingSequence(IGPart, IGEntry)
.. autoclass:: IGTerminator(IGPart, IGEntry)


Level 1
-------

Module
^^^^^^
.. autoclass:: IGCassette(Cassette)
   :members:

Vector
^^^^^^
.. autoclass:: IGDeviceVector(DeviceVector)
   :members:

Parts
^^^^^
.. autoclass:: IGEndLinker(IGPart, IGCassette)
   :members:


Level M
-------

Parts
^^^^^
.. autoclass:: IGLevelMVector(IGPart, IGDeviceVector)
   :members:

.. autoclass:: IGLevelMEndLinker(IGPart, IGCassette)
   :members:


Level P
-------

Parts
^^^^^
.. autoclass:: IGLevelPVector(IGPart, IGCassetteVector)
   :members:

.. autoclass:: IGLevelPEndLinker(IGPart, IGEntry)
   :members:
