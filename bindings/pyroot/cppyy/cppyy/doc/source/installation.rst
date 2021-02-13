.. _installation:

Installation
============

cppyy requires a (modern) C++ compiler.
When installing through `conda-forge`_, ``conda`` will install the compiler
for you, to match the other conda-forge packages.
When using ``pip`` and the wheels from `PyPI`_, you minimally need gcc5,
clang5, or MSVC'17.
When installing from source, the only requirement is full support for C++11
(e.g. minimum gcc 4.8.1 on GNU/Linux), but older compilers than the ones
listed for the wheels have not been tested.

With CPython on Linux or MacOS, probably by far the easiest way to install
cppyy, is through conda-forge on `Anaconda`_ (or `miniconda`_).
A Windows recipe for ``conda`` is not available yet, but is forthcoming, so
use ``pip`` for that platform for now (see below).
PyPI always has the authoritative releases (conda-forge pulls the sources
from there), so conda-forge may sometimes lag PyPI.
If you absolutely need the latest release, use PyPI or consider
:ref:`building from source <building_from_source>`.

To install using ``conda``, create and/or activate your (new) work environment
and install from the conda-forge channel::

  $ conda create -n WORK
  $ conda activate WORK
  (WORK) $ conda install -c conda-forge cppyy
  (WORK) [current compiler] $

To install with ``pip`` through `PyPI`_, it is recommend to use
`virtualenv`_ (or module `venv`_ for modern pythons).
The use of virtualenv prevents pollution of any system directories and allows
you to wipe out the full installation simply by removing the virtualenv
created directory ("WORK" in this example)::

  $ virtualenv WORK
  $ source WORK/bin/activate
  (WORK) $ python -m pip install cppyy
  (WORK) $

If you use the ``--user`` option to ``pip`` and use ``pip`` directly on the
command line, instead of through ``python``, make sure that the ``PATH``
envar points to the bin directory that will contain the installed entry
points during the installation, as the build process needs them.
You may also need to install ``wheel`` first if you have an older version of
``pip`` and/or do not use virtualenv (which installs wheel by default).
Example::

 $ python -m pip install wheel --user
 $ PATH=$HOME/.local/bin:$PATH python -m pip install cppyy --user


Wheels on PyPI
--------------

Wheels for the backend (``cppyy-cling``) are available on PyPI for GNU/Linux,
MacOS-X, and MS Windows (both 32b and 64b).

The Linux wheels are built on manylinux, but with gcc 5.5, not the 4.8.2 that
ships with manylinux1, since cppyy exposes C++ APIs and g++ introduced
ABI incompatibilities starting with its 5 series forward.
Using 4.8.2 would have meant that any software using cppyy would have to
be (re)compiled for the older gcc ABI, which the odds don't favor.
Note that building cppyy fully with 4.8.2 (and requiring the old ABI across
the board) does work.

The wheels for MS Windows were build with MSVC Community Edition 2017.

There are no wheels for the ``CPyCppyy`` and ``cppyy`` packages, to allow
the C++ standard chosen to match the local compiler.


pip with conda
--------------

Although installing ``cppyy`` through `conda-forge`_ is recommended, it is
possible to build/install with ``pip`` under Anaconda/miniconda.

Typical Python extensions only expose a C interface for use through the
Python C-API, requiring only calling conventions (and the Python C-API
version, of course) to match to be binary compatible.
Here, cppyy differs because it exposes C++ APIs: it thus requires a C++
run-time that is ABI compatible with the C++ compiler that was used during
build-time.

A set of modern compilers is available through conda-forge, but are only
intended for use with ``conda-build``.
In particular, the corresponding run-time is installed (for use through rpath
when building), but not set up.
That is, the conda compilers are added to ``PATH`` but not their libraries
to ``LD_LIBRARY_PATH`` (Mac, Linux; ``PATH`` for both on MS Windows).
Thus, you get the conda compilers and your system libraries mixed in the same
build environment, unless you set ``LD_LIBRARY_PATH`` (``PATH`` on Windows)
explicitly, e.g. by adding ``$CONDA_PREFIX/lib``.
Note that the conda documentation recommends against this.
Furthermore, the compilers from conda-forge are not vanilla distributions:
header files have been modified, which can can lead to parsing problems if
your system C library does not support C11, for example.

Nevertheless, with the above caveats, if your system C/C++ run-times are new
enough, the following can be made to work::

 $ conda create -n WORK
 $ conda activate WORK
 (WORK) $ conda install python
 (WORK) $ conda install -c conda-forge compilers
 (WORK) [current compiler] $ python -m pip install cppyy


C++ standard with pip
---------------------

The C++17 standard is the default for Mac and Linux (both PyPI and
conda-forge); but it is C++14 for MS Windows (compiler limitation).
When installing from PyPI using ``pip``, you can control the standard
selection by setting the ``STDCXX`` envar to '17', '14', or '11' (for Linux,
the backend does not need to be recompiled).
Note that the build will lower your choice if the compiler used does not
support a newer standard.


Install from source
-------------------
.. _installation_from_source:

To build an existing release from source, tell ``pip`` to not download any
binary wheels.
Build-time only dependencies are ``cmake`` (for general build), ``python``
(obviously, but also for LLVM), and a modern C++ compiler (one that supports
at least C++11).
Use the envar ``STDCXX`` to control the C++ standard version; ``MAKE`` to
change the ``make`` command, ``MAKE_NPROCS`` to control the maximum number of
parallel jobs allowed, and ``VERBOSE=1`` to see full build/compile commands.
Example (using ``--verbose`` to see ``pip`` progress)::

 $ STDCXX=17 MAKE_NPROCS=32 pip install --verbose cppyy --no-binary=cppyy-cling

Compilation of the backend, which contains a customized version of
Clang/LLVM, can take a long time, so by default the setup script will use all
cores (x2 if hyperthreading is enabled).
Once built, however, the wheel of ``cppyy-cling`` is reused by pip for all
versions of CPython and for PyPy, thus the long compilation is needed only
once for all different versions of Python on the same machine.

See the :ref:`section on repos <building_from_source>` for more
details/options.


PyPy
----

PyPy 5.7 and 5.8 have a built-in module ``cppyy``.
You can still install the cppyy package, but the built-in module takes
precedence.
To use cppyy, first import a compatibility module::

 $ pypy
 [PyPy 5.8.0 with GCC 5.4.0] on linux2
 >>>> import cppyy_compat, cppyy
 >>>>

You may have to set ``LD_LIBRARY_PATH`` appropriately if you get an
``EnvironmentError`` (it will indicate the needed directory).

Note that your python interpreter (whether CPython or ``pypy-c``) may not have
been linked by the C++ compiler.
This can lead to problems during loading of C++ libraries and program shutdown.
In that case, re-linking is highly recommended.

Very old versions of PyPy (5.6.0 and earlier) have a built-in ``cppyy`` based
on `Reflex`_, which is less feature-rich and no longer supported.
However, both the :doc:`distribution utilities <utilities>` and user-facing
Python codes are very backwards compatible, making migration straightforward.


Precompiled header
------------------

For performance reasons (reduced memory and CPU usage), a precompiled header
(PCH) of the system and compiler header files will be installed or, failing
that, generated on startup.
Obviously, this PCH is not portable and should not be part of any wheel.

Some compiler features, such as AVX, OpenMP, fast math, etc. need to be
active during compilation of the PCH, as they depend both on compiler flags
and system headers (for intrinsics, or API calls).
You can control compiler flags through the ``EXTRA_CLING_ARGS`` envar and thus
what is active in the PCH.
In principle, you can also change the C++ language standard by setting the
appropriate flag on ``EXTRA_CLING_ARGS`` and rebuilding the PCH.
However, if done at this stage, that disables some automatic conversion for
C++ types that were introduced after C++11 (such as ``string_view`` and
``optional``).

If you want multiple PCHs living side-by-side, you can generate them
yourself (note that the given path must be absolute)::

 >>> import cppyy_backend.loader as l
 >>> l.set_cling_compile_options(True)         # adds defaults to EXTRA_CLING_ARGS
 >>> install_path = '/full/path/to/target/location/for/PCH'
 >>> l.ensure_precompiled_header(install_path)

You can then select the appropriate PCH with the ``CLING_STANDARD_PCH`` envar::

 $ export CLING_STANDARD_PCH=/full/path/to/target/location/for/PCH/allDict.cxx.pch

Or disable it completely by setting that envar to "none".


.. _`conda-forge`: https://anaconda.org/conda-forge/cppyy
.. _`Anaconda`: https://www.anaconda.com/distribution/
.. _`miniconda`: https://docs.conda.io/en/latest/miniconda.html
.. _`PyPI`: https://pypi.python.org/pypi/cppyy/
.. _`virtualenv`: https://pypi.python.org/pypi/virtualenv
.. _`venv`: https://docs.python.org/3/library/venv.html
.. _`Reflex`: https://root.cern.ch/how/how-use-reflex
