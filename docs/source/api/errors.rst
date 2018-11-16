Errors
======

.. automodule:: moclo.errors
.. currentmodule:: moclo.errors


Base classes
------------

.. autoclass:: MocloError(Exception)
.. autoclass:: AssemblyError(MocloError, RuntimeError)
.. autoclass:: AssemblyWarning(MocloError, Warning)


Errors
------

.. autoclass:: DuplicateModules(AssemblyError)
.. autoclass:: InvalidSequence(MocloError, ValueError)
.. autoclass:: IllegalSite(InvalidSequence)
.. autoclass:: MissingModule(AssemblyError)


Warnings
--------

.. autoclass:: UnusedModules(AssemblyWarning)
