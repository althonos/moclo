``moclo``
=========

*A Python implementation of the* `MoClo <https://www.addgene.org/cloning/moclo/>`__ *system logic.*

|Source| |Travis| |Docs| |Codecov| |Codacy| |License|

.. |Codacy| image:: https://img.shields.io/codacy/grade/5b29a9c0d91f4e82944a46997bd9a480/master.svg?style=flat-square&maxAge=300
   :target: https://www.codacy.com/app/althonos/moclo

.. |Codecov| image:: https://img.shields.io/codecov/c/github/althonos/moclo/master.svg?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/moclo

.. |Travis| image:: https://img.shields.io/travis/althonos/moclo.svg?style=flat-square&maxAge=3600
   :target: https://travis-ci.org/althonos/moclo/branches

.. |License| image:: https://img.shields.io/github/license/althonos/moclo.svg?style=flat-square&maxAge=300
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/moclo

.. |Docs| image:: https://img.shields.io/readthedocs/moclo.svg?maxAge=3600&style=flat-square
   :target: https://moclo.readthedocs.io/


Documentation
-------------

The documentation is hosted on ``readthedocs.org``, and built against the latest
commit of the development repository. It contains a comprehensive API reference
as well as examples compiled from Jupyter notebooks at each build.


Base module
-----------

The base logic is handled by the core `moclo <https://github.com/althonos/moclo/tree/master/moclo-ytk>`_
module. It embbeds an object model of the MoClo system logic, but does not enforce
any specific sequence structure, and is not usable alone. You must install a kit
(listed below) to be able to validate and compute assemblies.


Kits
----

Additional kits can be installed separately depending on what's needed. The
following implementations are available:

* `Yeast ToolKit (moclo-ytk) <https://github.com/althonos/moclo/tree/master/moclo-ytk>`_
* `CIDAR Kit (moclo-cidar) <https://github.com/althonos/moclo/tree/master/moclo-cidar>`_
* `EcoFlex Kit (moclo-ecoflex) <https://github.com/althonos/moclo/tree/master/moclo-ecoflex>`_

Once installed, kits are available in the ``moclo.kits`` package namespace.
`Kit-specific documentation <https://moclo.readthedocs.io/en/latest/#kits>`_ is
available as well.


About
-----

This project is licensed under the `MIT License <http://choosealicense.com/licenses/mit/>`_.
It was developed during a placement at the
`InBio team <https://research.pasteur.fr/en/team/experimental-and-computational-methods-for-modeling-cellular-processes/>`_
at the Institut Pasteur of Paris during the summer of 2018.
