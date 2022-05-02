# `moclo` [![Stars](https://img.shields.io/github/stars/althonos/moclo.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/moclo/stargazers)

*A Python implementation of the [MoClo](https://www.addgene.org/cloning/moclo/) system logic.*

[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square)](https://github.com/althonos/moclo)
[![Build](https://img.shields.io/github/workflow/status/althonos/moclo/Test?style=flat-square&maxAge=3600)](https://github.com/althonos/moclo/actions)
[![Docs](https://img.shields.io/readthedocs/moclo.svg?maxAge=3600&style=flat-square)](https://moclo.readthedocs.io/)
[![Codecov](https://img.shields.io/codecov/c/github/althonos/moclo/master.svg?style=flat-square&maxAge=600)](https://codecov.io/gh/althonos/moclo)
[![Codacy](https://img.shields.io/codacy/grade/5b29a9c0d91f4e82944a46997bd9a480/master.svg?style=flat-square&maxAge=300)](https://www.codacy.com/app/althonos/moclo)
[![License](https://img.shields.io/pypi/l/moclo.svg?style=flat-square&maxAge=300)](https://choosealicense.com/licenses/mit/)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.1401815-blue.svg?style=flat-square&maxAge=31536000)](https://zenodo.org/badge/latestdoi/138012703)

## 📚 Documentation

The documentation is hosted on [ReadTheDocs](https://moclo.readthedocs.org),
and built against the latest commit of the development repository. It contains
a comprehensive API reference as well as examples compiled from Jupyter
notebooks at each build.


## 🔩 Base module

The base logic is handled by the core [`moclo`](https://github.com/althonos/moclo/tree/master/moclo)
module. It embeds an object model of the MoClo system logic, but does not enforce
any specific sequence structure, and is not usable alone. You must install a kit
(listed below) to be able to validate and compute assemblies.


## 🧰 Kits

Additional kits can be installed separately depending on what's needed. The
following implementations are available:

* [Yeast ToolKit (`moclo-ytk`)](https://github.com/althonos/moclo/tree/master/moclo-ytk)
* [CIDAR Kit (`moclo-cidar`)](https://github.com/althonos/moclo/tree/master/moclo-cidar)
* [EcoFlex Kit (`moclo-ecoflex`)](https://github.com/althonos/moclo/tree/master/moclo-ecoflex)
* [Icon Genetics / original MoClo (`moclo-ig`)](https://github.com/althonos/moclo/tree/master/moclo-ig)
* [GoldenBraid 3.0 (`moclo-gb3`)](https://github.com/althonos/moclo/tree/master/moclo-gb3)

Once installed, kits are available in the `moclo.kits` namespace module.
[Kit-specific documentation](https://moclo.readthedocs.io/en/latest/#kits) is
available as well.


## 🗂️ Registries

Kit-specific modules and vectors are distributed with the library files, so that
each library provides the base parts needed to create an assembly. They can be
found in the `moclo.registry` namespace. See also the documentation of each
`moclo.registry` submodule for a detail of how sequences were obtained. The
embedded sequences are distributed in GenBank format with the source distributions
of each plugin.


## 🗒️ Notebook

[![Docker Build Status](https://img.shields.io/docker/build/althonos/moclo.svg?style=flat-square&maxAge=3600)](https://hub.docker.com/r/althonos/moclo/builds/) [![Docker Pulls](https://img.shields.io/docker/pulls/althonos/moclo.svg?style=flat-square&maxAge=3600)](https://hub.docker.com/r/althonos/moclo/)

This repository provides a YTK-specific Jupyter notebook as a Docker image,
which can be used to generate a protocol for YTK MoClo assembly. Run it locally
using the following command:
```console
docker run --rm -it -p 8888:8888 althonos/moclo
```
and visit [https://localhost:8888/](https://localhost:8888/) to start interacting
with the notebook.


## ⚖️ License

This project is licensed under the [MIT License](http://choosealicense.com/licenses/mit/).

*This project is in no way affiliated, sponsored, or otherwise endorsed by [Addgene](https://www.addgene.org) or any of the MoClo toolkit creators.
It was developed by [Martin Larralde](https://github.com/althonos/pyhmmer)
during a placement at the [InBio team](https://research.pasteur.fr/en/team/experimental-and-computational-methods-for-modeling-cellular-processes/)
at the Institut Pasteur of Paris during the summer of 2018.*
