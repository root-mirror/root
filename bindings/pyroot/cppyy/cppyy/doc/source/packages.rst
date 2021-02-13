.. _packages:

PyPI Packages
=============

Cppyy
-----

The ``cppyy`` module is a frontend (see :ref:`Package Structure
<package-structure>`), and most of the code is elsewhere. However, it does
contain the docs for all of the modules, which are built using
Sphinx: http://www.sphinx-doc.org/en/stable/ and published to
http://cppyy.readthedocs.io/en/latest/index.html using a webhook. To create
the docs::

    $ pip install sphinx_rtd_theme
    Collecting sphinx_rtd_theme
    ...
    Successfully installed sphinx-rtd-theme-0.2.4
    $ cd docs
    $ make html

The Python code in this module supports:

* Interfacing to the correct backend for CPython or PyPy.
* Pythonizations (TBD)

Cppyy-backend
-------------

The ``cppyy-backend`` module contains two areas:

* A patched copy of cling
* Wrapper code


Package structure
-----------------
.. _package-structure:

There are four PyPA packages involved in a full installation, with the
following structure::

               (A) _cppyy (PyPy)
           /                        \
 (1) cppyy                            (3) cling-backend -- (4) cppyy-cling
           \                        /
             (2) CPyCppyy (CPython)

The user-facing package is always ``cppyy`` (1).
It is used to select the other (versioned) required packages, based on the
python interpreter for which it is being installed.

Below (1) follows a bifurcation based on interpreter.
This is needed for functionality and performance: for CPython, there is the
CPyCppyy package (2).
It is written in C++, makes use of the Python C-API, and installs as a Python
extension module.
For PyPy, there is the builtin module ``_cppyy`` (A).
This is not a PyPA package.
It is written in RPython as it needs access to low-level pointers, JIT hints,
and the ``_cffi_backend`` backend module (itself builtin).

Shared again across interpreters is the backend, which is split in a small
wrapper (3) and a large package that contains Cling/LLVM (4).
The former is still under development and expected to be updated frequently.
It is small enough to download and build very quickly.
The latter, however, takes a long time to build, but since it is very stable,
splitting it off allows the creation of binary wheels that need updating
only infrequently (expected about twice a year).

All code is publicly available; see the
:doc:`section on repositories <repositories>`.
