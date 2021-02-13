.. _exceptions:


Exceptions
==========

All C++ exceptions are converted to Python exceptions and all Python
exceptions are converted to C++ exceptions, to allow exception propagation
through multiple levels of callbacks, while retaining the option to handle
the outstanding exception as needed in either language.
To preserve an exception across the language boundaries, it must derive from
``std::exception``.
If preserving the exception (or its type) is not possible, generic exceptions
are used to propagate the exception: ``Exception`` in Python or
``CPyCppyy::PyException`` in C++.

In the most common case of an instance of a C++ exception class derived from
``std::exception`` that is thrown from a compiled library and which is
copyable, the exception can be caught and handled like any other bound C++
object (or with ``Exception`` on the Python and ``std::exception`` on the
C++ side).
If the exception is not copyable, but derived from ``std::exception``, the
result of its ``what()`` reported with an instance of Python's ``Exception``.
In all other cases, including exceptions thrown from interpreted code (due to
limitations of the Clang JIT), the exception will turn into an instance of
``Exception`` with a generic message.

The standard C++ exceptions are explicitly not mapped onto standard Python
exceptions, since other than a few simple cases, the mapping is too crude to
be useful as the typical usage in each standard library is too different.
Thus, for example, a thrown ``std::runtime_error`` instance will become a
``cppyy.gbl.std.runtime_error`` instance on the Python side (with Python's
``Exception`` as its base class), not a ``RuntimeError`` instance.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>

In addition, the examples require the ``throw`` to be in compiled code.
Save the following and build it into a shared library ``libfeatures.so`` (or
``libfeatures.dll`` on MS Windows):

  .. code-block:: C++

    #include "features.h"

    void throw_an_error(int i) {
        if (i) throw SomeError{"this is an error"};
        throw SomeOtherError{"this is another error"};
    }

And load the resulting library:

  .. code-block:: python

    >>> cppyy.load_library('libfeatures')
    >>>

Then try it out:

  .. code-block:: python

    >>> cppyy.gbl.throw_an_error(1)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    cppyy.gbl.SomeError: void ::throw_an_error(int i) =>
        SomeError: this is an error
    >>> 

Note how the full type is preserved and how the result of ``what()`` is used
for printing the exception.
By preserving the full C++ type, it is possible to call any other member
functions the exception may provide beyond ``what`` or access any additional
data it carries.

To catch the exception, you can either use the full type, or any of its base
classes, including ``Exception`` and ``cppyy.gbl.std.exception``:

  .. code-block:: python

    >>> try:
    ...     cppyy.gbl.throw_an_error(0)
    ... except cppyy.gbl.SomeOtherError as e:  # catch by exact type
    ...     print("received:", e)
    ... 
    received: <cppyy.gbl.SomeOtherError object at 0x7f9e11d3db10>
    >>> try:
    ...     cppyy.gbl.throw_an_error(0)
    ... except Exception as e:                 # catch through base class
    ...     print("received:", e)
    ... 
    received: <cppyy.gbl.SomeOtherError object at 0x7f9e11e00310>
    >>> 

