Library Reference
=================

.. toctree::
   :hidden:

   record
   registry
   modules
   vectors
   parts
   errors


Record (`moclo.record`)
-----------------------

.. currentmodule:: moclo.record
.. autosummary::
   :nosignatures:

   CircularRecord


Registry (`moclo.registry.base`)
--------------------------------

.. currentmodule:: moclo.registry.base
.. autosummary::
   :nosignatures:

   Item
   AbstractRegistry
   CombinedRegistry
   EmbeddedRegistry


.. currentmodule:: moclo.core
.. automodule:: moclo.core


Modules (`moclo.core.modules`)
------------------------------

.. currentmodule:: moclo.core.modules
.. autosummary::
   :nosignatures:

   AbstractModule
   Entry
   Cassette
   Device


Vectors (`moclo.core.vectors`)
------------------------------

.. currentmodule:: moclo.core.vectors
.. autosummary::
   :nosignatures:

   AbstractVector
   EntryVector
   CassetteVector
   DeviceVector

Parts (`moclo.core.parts`)
--------------------------

.. currentmodule:: moclo.core.parts
.. autosummary::
   :nosignatures:

   AbstractPart


Errors (`moclo.errors`)
-----------------------

.. currentmodule:: moclo.errors
.. rubric:: Base classes
.. autosummary::
   :nosignatures:

   MocloError
   AssemblyError
   AssemblyWarning

.. rubric:: Errors
.. autosummary::
   :nosignatures:

   DuplicateModules
   InvalidSequence
   MissingModule

.. rubric:: Warnings
.. autosummary::
   :nosignatures:

   UnusedModules
