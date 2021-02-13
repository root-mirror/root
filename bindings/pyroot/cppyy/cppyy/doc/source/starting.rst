.. _starting:

Trying it out
=============

This is a basic guide to try cppyy and see whether it works for you.
Large code bases will benefit from more advanced features such as
:doc:`pythonizations <pythonizations>` for a cleaner interface to clients;
precompiled modules for faster parsing and reduced memory usage;
":ref:`dictionaries <dictionaries>`" to package locations and manage
dependencies; and mapping files for automatic, lazy, loading.
You can, however, get very far with just the basics and it may even be
completely sufficient for small packages with fewer classes.

cppyy works by parsing C++ definitions through ``cling``, generating tiny
wrapper codes to honor compile-time features and create standardized
interfaces, then compiling/linking those wrappers with the ``clang`` JIT.
It thus requires only those two ingredients: *C++ definitions* and
*linker symbols*.
All cppyy uses, the basic and the more advanced, are variations on the
theme of bringing these two together at the point of use.

Definitions typically live in header files and symbols in libraries.
Headers can be loaded with ``cppyy.include`` and libraries with the
``cppyy.load_library`` call.
Loading the header is sufficient to start exploring, with ``cppyy.gbl`` the
starting point of all things C++, while the linker symbols are only needed at 
the point of first use.

Here is an example using the `zlib`_ library, which is likely available on
your system:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('zlib.h')        # bring in C++ definitions
    >>> cppyy.load_library('libz')     # load linker symbols
    >>> cppyy.gbl.zlibVersion()        # use a zlib API
    '1.2.11'
    >>>

Since header files can include other header files, it is easy to aggregate
all relevant ones into a single header to include.
If there are project-specific include paths, you can add those paths through
``cppyy.add_include_path``.
If a header is C-only and not set for use with C++, use ``cppyy.c_include``,
which adds ``extern "C"`` around the header.

Library files can be aggregated by linking all relevant ones to a single
library to load.
Using the linker for this purpose allows regular system features such as
``rpath`` and envars such as ``LD_LIBRARY_PATH`` to be applied as usual.
Note that any mechanism that exposes the library symbols will work.
For example, you could also use the standard module ``ctypes`` through
``ctypes.CDLL`` with the ``ctypes.RTLD_GLOBAL`` option.

To explore, start from ``cppyy.gbl`` to access your namespaces, classes,
functions, etc., etc. directly; or use python's ``dir`` (or tab-completion)
to see what is available.
Use python's ``help`` to see list the methods and data members of classes and
see the interfaces of functions.

Now try this out for some of your own headers, libraries, and APIs!

.. _`zlib`: https://en.wikipedia.org/wiki/Zlib
