EcoFlex Kit
===========

.. currentmodule:: moclo.kits.ecoflex

.. automodule:: moclo.kits.ecoflex

.. Level -1
.. --------

.. Module
.. ^^^^^^
.. .. autoclass:: EcoFlexProduct(Product)
..    :members:
..    :inherited-members:
..    :special-members: __init__


.. Vector
.. ^^^^^^
.. .. autoclass:: EcoFlexEntryVector(EntryVector)
..    :members:
..    :inherited-members:
..    :special-members: __init__


Level 0
-------

Module
^^^^^^
.. autoclass:: EcoFlexEntry(Entry)
   :members:
   :inherited-members:
   :special-members: __init__


Vector
^^^^^^
.. autoclass:: EcoFlexCassetteVector(CassetteVector)
   :members:
   :inherited-members:
   :special-members: __init__


Parts
^^^^^

.. autoclass:: EcoFlexPromoter(EcoFlexPart, EcoFlexEntry)
.. autoclass:: EcoFlexRBS(EcoFlexPart, EcoFlexEntry)
.. autoclass:: EcoFlexTagLinker(EcoFlexPart, EcoFlexEntry)
.. autoclass:: EcoFlexTag(EcoFlexPart, EcoFlexEntry)
.. autoclass:: EcoFlexCodingSequence(EcoFlexPart, EcoFlexEntry)
.. autoclass:: EcoFlexTerminator(EcoFlexPart, EcoFlexEntry)



Level 1
-------

Module
^^^^^^
.. autoclass:: EcoFlexCassette(Cassette)
   :members:

Vector
^^^^^^
.. autoclass:: EcoFlexDeviceVector(DeviceVector)
   :members:


Level 2
-------

Module
^^^^^^

.. autoclass:: EcoFlexDevice(Device)
   :members:
