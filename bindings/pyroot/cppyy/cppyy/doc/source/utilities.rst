.. _utilities:

Utilities
=========

The ``cppyy-backend`` package brings in the following utilities to help
with repackaging and redistribution:

  * cling-config: for compile time flags
  * rootcling and genreflex: for dictionary generation
  * cppyy-generator: part of the :doc:`CMake interface <cmake_interface>`


Compiler/linker flags
---------------------

``cling-config`` is a small utility to provide access to the as-installed
configuration, such as compiler/linker flags and installation directories, of
other components.
Usage examples::

    $ cling-config --help
    Usage: cling-config [--cflags] [--cppflags] [--cmake]
    $ cling-config --cmake
    /usr/local/lib/python2.7/dist-packages/cppyy_backend/cmake


.. _dictionaries:

Dictionaries
------------

Loading header files or code directly into ``cling`` is fine for interactive
work and smaller packages, but large scale applications benefit from
pre-compiling code, using the automatic class loader, and packaging
dependencies in so-called "dictionaries."

A `dictionary` is a generated C++ source file containing references to the
header locations used when building (and any additional locations provided),
a set of forward declarations to reduce the need of loading header files, and
a few I/O helper functions.
The name "dictionary" is historic: before ``cling`` was used, it contained
the complete generated C++ reflection information, whereas now that is
derived at run-time from the header files.
It is still possible to fully embed header files rather than only storing
their names and search locations, to make the dictionary more self-contained.

After generating the dictionary, it should be compiled into a shared library.
This provides additional dependency control: by linking it directly with any
further libraries needed, you can use standard mechanisms such as ``rpath``
to locate those library dependencies.
Alternatively, you can add the additional libraries to load to the mapping
files of the class loader (see below).

In tandem with any dictionary, a pre-compiled module (.pcm) file will be
generated.
C++ modules are still on track for inclusion in the C++20 standard and most
modern C++ compilers, ``clang`` among them, already have implementations.
The benefits for cppyy include faster bindings generation, lower memory
footprint, and isolation from preprocessor macros and compiler flags.
The use of modules is transparent, other than the requirement that they
need to be co-located with the compiled dictionary shared library.

Optionally, the dictionary generation process also produces a mapping file,
which lists the libraries needed to load C++ classes on request (for details,
see the section on the class loader below).

Structurally, you could have a single dictionary for a project as a whole,
but more likely a large project will have a pre-existing functional
decomposition that can be followed, with a dictionary per functional unit.


Generation
^^^^^^^^^^

There are two interfaces onto the same underlying dictionary generator:
``rootcling`` and ``genreflex``.
The reason for having two is historic and they are not complete duplicates,
so one or the other may suit your preference better.
It is foreseen that both will be replaced once C++ modules become more
mainstream, as that will allow simplification and improved robustness.

rootcling
"""""""""

The first interface is called ``rootcling``::

    $ rootcling
    Usage: rootcling [-v][-v0-4] [-f] [out.cxx] [opts] file1.h[+][-][!] file2.h[+][-][!] ...[Linkdef.h]
    For more extensive help type: /usr/local/lib/python2.7/dist-packages/cppyy_backend/bin/rootcling -h

Rather than providing command line options, the main steering of
``rootcling`` behavior is done through
`#pragmas in a Linkdef.h <https://root.cern.ch/root/html/guides/users-guide/AddingaClass.html#the-linkdef.h-file>`_
file, with most pragmas dedicated to selecting/excluding (parts of) classes
and functions.
Additionally, the Linkdef.h file may contain preprocessor macros.

The output consists of a dictionary file (to be compiled into a shared
library), a C++ module, and an optional mapping file, as described above.

genreflex
"""""""""

The second interface is called ``genreflex``::

    $ genreflex
    Generates dictionary sources and related ROOT pcm starting from an header.
    Usage: genreflex headerfile.h [opts] [preproc. opts]
    ...

``genreflex`` has a richer command line interface than ``rootcling`` as can
be seen from the full help message.

.. _selection-files:

Selection/exclusion is driven through a `selection file`_ using an XML format
that allows both exact and pattern matching to namespace, class, enum,
function, and variable names.

.. _`selection file`: https://linux.die.net/man/1/genreflex


Example
"""""""

Consider the following basic example code, living in a header "MyClass.h":

  .. code-block:: C++

    class MyClass {
    public:
        MyClass(int i) : fInt(i) {}
        int get_int() { return fInt; }

    private:
        int fInt;
    };

and a corresponding "Linkdef.h" file, selecting only ``MyClass``::

    #ifdef __ROOTCLING__
    #pragma link off all classes;
    #pragma link off all functions;
    #pragma link off all globals;
    #pragma link off all typedef;

    #pragma link C++ class MyClass;

    #endif

For more pragmas, see the `rootcling manual`_.
E.g., a commonly useful pragma is one that selects all C++ entities that are
declared in a specific header file::

    #pragma link C++ defined_in "MyClass.h";

Next, use ``rootcling`` to generate the dictionary (here:
``MyClass_rflx.cxx``) and module files::

    $ rootcling -f MyClass_rflx.cxx MyClass.h Linkdef.h

Alternatively, define a "myclass_selection.xml" file::

    <lcgdict>
        <class name="MyClass" />
    </lcgdict>

serving the same purpose as the Linkdef.h file above (in fact, ``rootcling``
accepts a "selection.xml" file in lieu of a "Linkdef.h").
For more tags, see the `selection file`_ documentation.
Commonly used are ``namespace``, ``function``, ``enum``, or ``variable``
instead of the ``class`` tag, and ``pattern`` instead of ``name`` with
wildcarding in the value string.

Next, use ``genreflex`` to generate the dictionary (here:
``MyClass_rflx.cxx``) and module files::

    $ genreflex MyClass.h --selection=myclass_selection.xml -o MyClass_rflx.cxx

From here, compile and link the generated dictionary file with the project
and/or system specific options and libraries into a shared library, using
``cling-config`` for the relevant cppyy compiler/linker flags.
(For work on MS Windows, this `helper script`_ may be useful.)
To continue the example, assuming Linux::

    $ g++ `cling-config --cppflags` -fPIC -O2 -shared MyClass_rflx.cxx -o MyClassDict.so

Instead of loading the header text into ``cling``, you can now load the
dictionary:

.. code-block:: python

    >>> import cppyy
    >>> cppyy.load_reflection_info('MyClassDict')
    >>> cppyy.gbl.MyClass(42)
    <cppyy.gbl.MyClass object at 0x7ffb9f230950>
    >>> print(_.get_int())
    42
    >>>

and use the selected C++ entities as if the header was loaded.

The dictionary shared library can be relocated, as long as it can be found
by the dynamic loader (e.g. through ``LD_LIBRARY_PATH``) and the header file
is fully embedded or still accessible (e.g. through a path added to
``cppyy.add_include_path`` at run-time, or with ``-I`` to
``rootcling``/``genreflex`` during build time).
When relocating the shared library, move the .pcm with it.
Once support for C++ modules is fully fleshed out, access to the header file
will no longer be needed.

.. _`rootcling manual`: https://root.cern.ch/root/html/guides/users-guide/AddingaClass.html#the-linkdef.h-file
.. _`helper script`: https://bitbucket.org/wlav/cppyy/src/master/test/make_dict_win32.py


Class loader
^^^^^^^^^^^^

Explicitly loading dictionaries is fine if this is hidden under the hood of
a Python package and thus transparently done on ``import``.
Otherwise, the automatic class loader is more convenient, as it allows direct
use without having to manually find and load dictionaries (assuming these are
locatable by the dynamic loader).

The class loader utilizes so-called rootmap files, which by convention should
live alongside the dictionary shared library (and C++ module file).
These are simple text files, which map C++ entities (such as classes) to the
dictionaries and other libraries that need to be loaded for their use.

With ``genreflex``, the mapping file can be automatically created with
``--rootmap-lib=MyClassDict``, where "MyClassDict" is the name of the shared
library (without the extension) build from the dictionary file.
With ``rootcling``, create the same mapping file with
``-rmf MyClassDict.rootmap -rml MyClassDict``.
It is necessary to provide the final library name explicitly, since it is
only in the separate linking step where these names are fixed and those names
may not match the default choice.

With the mapping file in place, the above example can be rerun without
explicit loading of the dictionary:

.. code-block:: python

    >>> import cppyy
    >>> from cppyy.gbl import MyClass
    >>> MyClass(42).get_int()
    42
    >>>


.. _cppyy-generator:

Bindings collection
-------------------

``cppyy-generator`` is a clang-based utility program which takes a set of C++
header files and generates a JSON output file describing the objects found in
them.
This output is intended to support more convenient access to a set of
cppyy-supported bindings::

    $ cppyy-generator --help
    usage: cppyy-generator [-h] [-v] [--flags FLAGS] [--libclang LIBCLANG]
                           output sources [sources ...]
    ...

This utility is mainly used as part of the
:doc:`CMake interface <cmake_interface>`.
