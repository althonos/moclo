``moclo``
=========

*A Python implementation of the* `MoClo <https://www.addgene.org/cloning/moclo/>`__ *system logic.*

|Source| |PyPI| |Travis| |Docs| |Codecov| |Codacy| |Format| |License|

.. |Codacy| image:: https://img.shields.io/codacy/grade/5b29a9c0d91f4e82944a46997bd9a480/master.svg?style=flat-square&maxAge=300
   :target: https://www.codacy.com/app/althonos/moclo

.. |Codecov| image:: https://img.shields.io/codecov/c/github/althonos/moclo/master.svg?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/moclo

.. |PyPI| image:: https://img.shields.io/pypi/v/moclo.svg?style=flat-square&maxAge=300
   :target: https://pypi.python.org/pypi/moclo

.. |Travis| image:: https://img.shields.io/travis/althonos/moclo.svg?style=flat-square&maxAge=3600
   :target: https://travis-ci.org/althonos/moclo/branches

.. |Format| image:: https://img.shields.io/pypi/format/moclo.svg?style=flat-square&maxAge=300
   :target: https://pypi.python.org/pypi/moclo

.. |Versions| image:: https://img.shields.io/pypi/pyversions/moclo.svg?style=flat-square&maxAge=300
   :target: https://travis-ci.org/althonos/moclo/

.. |License| image:: https://img.shields.io/pypi/l/moclo.svg?style=flat-square&maxAge=300
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/moclo/tree/master/moclo

.. |Docs| image:: https://img.shields.io/readthedocs/moclo.svg?maxAge=3600&style=flat-square
   :target: https://moclo.readthedocs.io/


Requirements
------------

+---------------------+----------------------------+------------------------+--------------------------+---------------------------+
| **biopython**       |  Sequence handling         | |PyPI biopython|       | |Source biopython|       | |License biopython|       |
+---------------------+----------------------------+------------------------+--------------------------+---------------------------+
| **cached-property** |  Lazy regex evaluation     | |PyPI cached-property| | |Source cached-property| | |License cached-property| |
+---------------------+----------------------------+------------------------+--------------------------+---------------------------+
| **six**             | Python 2/3 compatibility   | |PyPI six|             | |Source six|             | |License six|             |
+---------------------+----------------------------+------------------------+--------------------------+---------------------------+

.. |PyPI cached-property| image:: https://img.shields.io/pypi/v/cached-property.svg?style=flat-square&maxAge=600
   :target: https://pypi.python.org/pypi/cached-property

.. |PyPI biopython| image:: https://img.shields.io/pypi/v/biopython.svg?style=flat-square&maxAge=600
   :target: https://pypi.org/project/biopython/

.. |PyPI six| image:: https://img.shields.io/pypi/v/six.svg?style=flat-square&maxAge=600
   :target: https://pypi.org/project/six/

.. |Source cached-property| image:: https://img.shields.io/badge/source-GitHub-303030.svg?style=flat-square&maxAge=600
   :target: https://github.com/pydanny/cached-property

.. |Source biopython| image:: https://img.shields.io/badge/source-GitHub-303030.svg?style=flat-square&maxAge=600
   :target: https://github.com/biopython/biopython

.. |Source six| image:: https://img.shields.io/badge/source-GitHub-303030.svg?style=flat-square&maxAge=600
   :target: https://github.com/benjaminp/six

.. |License cached-property| image:: https://img.shields.io/pypi/l/cached-property.svg?style=flat-square&maxAge=600
   :target: https://choosealicense.com/licenses/bsd-3-clause/

.. |License biopython| image:: https://img.shields.io/badge/license-BSD%2FBioPython-blue.svg?style=flat-square&maxAge=600
   :target: https://choosealicense.com/licenses/bsd-3-clause/

.. |License six| image:: https://img.shields.io/pypi/l/six.svg?style=flat-square&maxAge=600
   :target: https://choosealicense.com/licenses/mit/


Installation
------------

This package is available as a *wheel*, and can be installed with ``pip``::

  $ pip install --user moclo

To see more ways of installing, head over to the `Installation <https://moclo.readthedocs.io/en/latest/install.html>`__
page of the online documentation.


Kits
----

By itself, ``moclo`` is not very useful. To be able to simulate MoClo assemblies
you can install some of the following toolkits:

- `moclo-ytk <https://pypi.org/project/moclo-ytk>`_: MoClo Yeast ToolKit,
  *John Dueber Lab*, and Pichia ToolKit, *Volker Sieber Lab*
- `moclo-cidar <https://pypi.org/project/moclo-cidar>`_: MoClo CIDAR kit,
  *Douglas Densmore Lab*
- `moclo-ecoflex <https://pypi.org/project/moclo-ecoflex>`_: MoClo EcoFlex,
  *Paul Freemont Lab*
- `moclo-ig <https://pypi.org/project/moclo-ig>`_: Icon Genetics MoClo,
  *Sylvestre Marillonnet Lab*
- `moclo-gb3 <https://pypi.org/project/moclo-gb3>`_: Golden Braid 3.0,
  *Diego Orzaez Lab*

Toolkits ship with concrete implementation of the MoClo logic (using the DNA
signatures and restriction enzymes from the reference paper), as well as official
sequences obtained from `AddGene <https://www.addgene.org>`_ and manually
annotated with higher-quality features. These sequences can be accessed through
the ``moclo.registry`` module, using the *registry* interface.


About
-----

This library is licensed under the `MIT License <http://choosealicense.com/licenses/mit/>`_.
It was developed during a placement at the
`InBio team <https://research.pasteur.fr/en/team/experimental-and-computational-methods-for-modeling-cellular-processes/>`_
at the Institut Pasteur of Paris during the summer of 2018.
