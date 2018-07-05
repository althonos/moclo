Installation
============

The ``moclo`` module is designed to be modular, and as such, you only need to
install whatever functionalities you are willing to use. Packages are distributed
on PyPI, and it is advised to use ``pip`` to install them. See the
`pip documentation <https://pip.pypa.io/en/stable/installing/>`_ to get pip if
it is not installed on your system.

Commands below use ``pip`` in user mode: the packages will be installed in a
user-dependent location, and no additional permissions are needed. If for some
reason you need a system-wide setup, remove the ``--user`` flag. Installing in
user-mode should be prefered to avoid dependency issues, in particular when on
an OS which provides a package manager (such as ``aptitude`` on Debian, or even
``homebrew`` on Mac OSX).


PyPI + ``pip`` |PyPI|
---------------------

.. |PyPI| image:: https://img.shields.io/pypi/v/moclo.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/moclo

To download the latest release from the Python Package Index:

.. code-block:: console

  $ pip install --user moclo moclo-ytk moclo-cidar


GitHub + ``pip`` |Travis|
-------------------------

.. |Travis| image:: https://img.shields.io/travis/althonos/moclo.svg?style=flat-square&maxAge=3600
   :target: https://travis-ci.org/althonos/moclo

To download the development version from the source repository, you can specify
a subfolder in the installation command and directly install it:

.. code-block:: console

  $ pip install --user git+https://github.com/althonos/moclo#subdirectory=moclo
  $ pip install --user git+https://github.com/althonos/moclo#subdirectory=moclo-ytk
  $ pip install --user git+https://github.com/althonos/moclo#subdirectory=moclo-cidar

Check the CI build is passing, or else you may be installing a broken version of
the library !
