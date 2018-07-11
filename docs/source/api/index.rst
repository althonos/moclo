Base API
========

.. toctree::
   :hidden:

   parts
   modules
   vectors
   errors

.. currentmodule:: moclo.base
.. automodule:: moclo.base


Modules (`moclo.base.modules`)
------------------------------

.. currentmodule:: moclo.base.modules
.. autosummary::
   :nosignatures:

   AbstractModule
   Entry
   Cassette
   Device


Vectors (`moclo.base.vectors`)
------------------------------

.. currentmodule:: moclo.base.vectors
.. autosummary::
   :nosignatures:

   AbstractVector
   EntryVector
   CassetteVector
   DeviceVector

Parts (`moclo.base.parts`)
--------------------------

.. currentmodule:: moclo.base.parts
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


Record (`moclo.record`)
-----------------------

.. currentmodule:: moclo.record
.. autosummary::
   :nosignatures:

   CircularRecord
