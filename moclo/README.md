# `moclo`

*A Python implementation of the [MoClo](https://www.addgene.org/cloning/moclo/) system logic.*

[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square)](https://github.com/althonos/moclo/tree/master/moclo)
[![PyPI](https://img.shields.io/pypi/v/moclo.svg?style=flat-square&maxAge=300)](https://pypi.python.org/pypi/moclo)
[![Build](https://img.shields.io/github/workflow/status/althonos/moclo/Test?style=flat-square&maxAge=3600)](https://github.com/althonos/moclo/actions)
[![Docs](https://img.shields.io/readthedocs/moclo.svg?maxAge=3600&style=flat-square)](https://moclo.readthedocs.io/)
[![Codecov](https://img.shields.io/codecov/c/github/althonos/moclo/master.svg?style=flat-square&maxAge=600)](https://codecov.io/gh/althonos/moclo)
[![Codacy](https://img.shields.io/codacy/grade/5b29a9c0d91f4e82944a46997bd9a480/master.svg?style=flat-square&maxAge=300)](https://www.codacy.com/app/althonos/moclo)
[![Format](https://img.shields.io/pypi/format/moclo.svg?style=flat-square&maxAge=300)](https://pypi.python.org/pypi/moclo)
[![Versions](https://img.shields.io/pypi/pyversions/moclo.svg?style=flat-square&maxAge=300)](https://pypi.python.org/pypi/moclo)
[![License](https://img.shields.io/pypi/l/moclo.svg?style=flat-square&maxAge=300)](https://choosealicense.com/licenses/mit/)


## 📚 Documentation

The documentation is hosted on [ReadTheDocs](https://moclo.readthedocs.org),
and built against the latest commit of the development repository. It contains
a comprehensive API reference as well as examples compiled from Jupyter
notebooks at each build.


## 🛠️ Installation

This package is available as a *wheel*, and can be installed with ``pip``:
```console
$ pip install --user moclo
```

To see more ways of installing, head over to the
[Installation](https://moclo.readthedocs.io/en/latest/install.html)
page of the online documentation.


## 🧰 Kits

By itself, `moclo` is not very useful. To be able to simulate MoClo assemblies
you can install some of the following toolkits:

- [`moclo-ytk`](https://pypi.org/project/moclo-ytk): MoClo Yeast ToolKit,
  *John Dueber Lab*, and Pichia ToolKit, *Volker Sieber Lab*
- [`moclo-cidar`](https://pypi.org/project/moclo-cidar): MoClo CIDAR kit,
  *Douglas Densmore Lab*
- [`moclo-ecoflex`](https://pypi.org/project/moclo-ecoflex): MoClo EcoFlex,
  *Paul Freemont Lab*
- [`moclo-ig`](https://pypi.org/project/moclo-ig): Icon Genetics MoClo,
  *Sylvestre Marillonnet Lab*
- [`moclo-gb3`](https://pypi.org/project/moclo-gb3): Golden Braid 3.0,
  *Diego Orzaez Lab*

Toolkits ship with concrete implementation of the MoClo logic (using the DNA
signatures and restriction enzymes from the reference paper), as well as official
sequences obtained from [Addgene](https://www.addgene.org) and manually
annotated with higher-quality features. These sequences can be accessed through
the `moclo.registry` module, using the *registry* interface.


## 📜 License

This project is licensed under the [MIT License](http://choosealicense.com/licenses/mit/).

*This project is in no way affiliated, sponsored, or otherwise endorsed by [Addgene](https://www.addgene.org) or any of the MoClo toolkit creators.
It was developed by [Martin Larralde](https://github.com/althonos/pyhmmer)
during a placement at the [InBio team](https://research.pasteur.fr/en/team/experimental-and-computational-methods-for-modeling-cellular-processes/)
at the Institut Pasteur of Paris during the summer of 2018.*
