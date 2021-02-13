.. _dictionaries:

Dictionaries
============

Loading code directly into Cling is fine for interactive work and small
scripts, but large scale applications should take advantage of pre-packaging
code, linking in libraries, and describing other dependencies.
The necessary tools are installed as part of the backend.


Dictionary generation
---------------------

A "reflection dictionary" makes it simple to combine the necessary headers and
libraries into a single package for use and distribution.
The relevant headers are read by a tool called `genreflex`_ which generates
C++ files that are to be compiled into a shared library.
That library can further be linked with any relevant project libraries that
contain the implementation of the functionality declared in the headers.
For example, given a file called ``project_header.h`` and an implementation
residing in ``libproject.so``, the following will generate a
``libProjectDict.so`` reflection dictionary::

    $ genreflex project_header.h
    $ g++ -std=c++17 -fPIC -rdynamic -O2 -shared `genreflex --cppflags` project_header_rflx.cpp -o libProjectDict.so -L$PROJECTHOME/lib -lproject

Instead of loading the header text into Cling, you can now load the
dictionary:

.. code-block:: python

    >>> import cppyy
    >>> cppyy.load_reflection_info('libProjectDict.so')
    <CPPLibrary object at 0xb6fd7c4c>
    >>> from cppyy.gbl import SomeClassFromProject
    >>>

and use the C++ entities from the header as before.

.. _`genreflex`: https://linux.die.net/man/1/genreflex


Automatic class loader
----------------------

Explicitly loading dictionaries is fine if this is hidden under the hood of
a Python package and thus simply done on import.
Otherwise, the automatic class loader is more convenient, as it allows direct
use without having to manually find and load dictionaries.

The class loader utilizes so-called rootmap files, which by convention should
live alongside the dictionaries in places reachable by LD_LIBRARY_PATH.
These are simple text files, which map C++ entities (such as classes) to the
dictionaries and other libraries that need to be loaded for their use.

The ``genreflex`` tool can produce rootmap files automatically.
For example::

    $ genreflex project_header.h --rootmap=libProjectDict.rootmap --rootmap-lib=libProjectDict.so
    $ g++ -std=c++17 -fPIC -rdynamic -O2 -shared `genreflex --cppflags` project_header_rflx.cpp -o libProjectDict.so -L$CPPYYHOME/lib -lCling -L$PROJECTHOME/lib -lproject

where the first option (``--rootmap``) specifies the output file name, and the
second option (``--rootmap-lib``) the name of the reflection library.
It is necessary to provide that name explicitly, since it is only in the
separate linking step where these names are fixed (if the second option is not
given, the library is assumed to be libproject_header.so).

With the rootmap file in place, the above example can be rerun without explicit
loading of the reflection info library:

.. code-block:: python

    >>> import cppyy
    >>> from cppyy.gbl import SomeClassFromProject
    >>>


Selection files
---------------
.. _selection-files:

Sometimes it is necessary to restrict or expand what genreflex will pick up
from the header files.
For example, to add or remove standard classes or to hide implementation
details.
This is where `selection files`_ come in.
These are XML specifications that allow exact or pattern matching to classes,
functions, etc.
See ``genreflex --help`` for a detailed specification and add
``--selection=project_selection.xml`` to the ``genreflex`` command line.

With the aid of a selection file, a large project can be easily managed:
simply ``#include`` all relevant headers into a single header file that is
handed to ``genreflex``.

.. _`selection files`: https://linux.die.net/man/1/genreflex
