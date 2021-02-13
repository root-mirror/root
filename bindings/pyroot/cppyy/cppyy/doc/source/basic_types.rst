.. _basic_types:


Basic types
===========

C++ has a far richer set of builtin types than Python.
Most Python code can remain relatively agnostic to that, and ``cppyy``
provides automatic conversions as appropriate.
On the other hand, Python builtin types such as lists and maps are far
richer than any builtin types in C++.
These are mapped to their Standard Template Library equivalents instead.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>


`Builtins`
""""""""""

The selection of builtin data types varies greatly between Python and C++.
Where possible, builtin data types map onto the expected equivalent Python
types, with the caveats that there may be size differences, different
precision or rounding, etc.
For example, a C++ ``float`` is returned as a Python ``float``, which is in
fact a C++ ``double``.
If sizes allow, conversions are automatic.
For example, a C++ ``unsigned int`` becomes a Python2 ``long`` or Python3
``int``, but unsigned-ness is still honored:

  .. code-block:: python

    >>> cppyy.gbl.gUint
    0L
    >>> type(cppyy.gbl.gUint)
    <type 'long'>
    >>> cppyy.gbl.gUint = -1
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ValueError: cannot convert negative integer to unsigned
    >>>

Some types are builtin in Python, but (STL) classes in C++.
Examples are ``str`` vs. ``std::string`` and ``complex`` vs. ``std::complex``.
These classes have been pythonized to behave the same wherever possible.
For example, string comparison work directly, and ``std::complex`` has
``real`` and ``imag`` properties:

  .. code-block:: python

    >>> c = cppyy.gbl.std.complex['double'](1, 2)
    >>> c
    (1+2j)
    >>> c.real, c.imag
    (1.0, 2.0)
    >>> s = cppyy.gbl.std.string("aap")
    >>> type(s)
    <class cppyy.gbl.std.string at 0x7fa75edbf8a0>
    >>> s == "aap"
    True
    >>>

To pass an argument through a C++ ``char`` (signed or unsigned) use a Python
string of size 1.
In many cases, the explicit C types from module ``ctypes`` can also be used,
but that module does not have a public API (for type conversion or otherwise),
so support is somewhat limited.

There are automatic conversions between C++'s ``std::vector`` and Python's
``list`` and ``tuple``, where possible, as they are often used in a similar
manner.
These datatypes have completely different memory layouts, however, and the
``std::vector`` requires that all elements are of the same type and laid
out consecutively in memory.
Conversion thus requires type checks, memory allocation, and copies.
This can be rather expensive.
See the section on :ref:`STL <stl>`.


`Arrays`
""""""""

Builtin arrays are supported through arrays from module ``array`` (or any
other builtin-type array that implements the Python buffer interface, such
as numpy arrays) and a low-level view type from ``cppyy`` for returns and
variable access (that implements the buffer interface as well).
Out-of-bounds checking is limited to those cases where the size is known at
compile time.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> from array import array
    >>> c = Concrete()
    >>> c.array_method(array('d', [1., 2., 3., 4.]), 4)
    1 2 3 4
    >>> c.m_data[4] # static size is 4, so out of bounds
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    IndexError: buffer index out of range
    >>>

Arrays of arrays are supported through the C++ low-level view objects.
This only works well if sizes are known at compile time or can be inferred.
If sizes are not known, the size is set to a large integer (depending on the
array element size) to allow access.
It is then up to the developer not to access the array out-of-bounds.
There is limited support for arrays of instances, but those should be avoided
in C++ anyway:

  .. code-block:: python

    >>> cppyy.cppdef('std::string str_array[3][2] = {{"aa", "bb"}, {"cc", "dd"}, {"ee", "ff"}};')
    True
    >>> type(cppyy.gbl.str_array[0][1])
    <class cppyy.gbl.std.string at 0x7fd650ccb650>
    >>> cppyy.gbl.str_array[0][1]
    'bb'
    >>> cppyy.gbl.str_array[4][0]
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    IndexError: tuple index out of range
    >>>


`Pointers`
""""""""""

When the C++ code takes a pointer or reference type to a specific builtin
type (such as an ``unsigned int`` for example), then types need to match
exactly.
``cppyy`` supports the types provided by the standard modules ``ctypes`` and
``array`` for those cases.
Example of using a reference to builtin:

  .. code-block:: python

    >>> from ctypes import c_uint
    >>> u = c_uint(0)
    >>> c.uint_ref_assign(u, 42)
    >>> u.value
    42
    >>>

For objects, an object, a pointer to an object, and a smart pointer to an
object are represented the same way, with the necessary (de)referencing
applied automatically.
Pointer variables are also bound by reference, so that updates on either the
C++ or Python side are reflected on the other side as well.


`Enums`
"""""""

Named, anonymous, and class enums are supported.
The Python-underlying type of an enum is implementation dependent and may even
be different for different enums on the same compiler.
Typically, however, the types are ``int`` or ``unsigned int``, which
translates to Python's ``int`` or ``long`` on Python2 or class ``int`` on
Python3.
Separate from the underlying, all enums have their own Python type to allow
them to be used in template instantiations:

  .. code-block:: python

    >>> from cppyy.gbl import kBanana   # classic enum, globally available
    >>> print(kBanana)
    29
    >>> cppyy.gbl.EFruit
    <class '__main__.EFruit'>
    >>> print(cppyy.gbl.EFruit.kApple)
    78
    >>> cppyy.gbl.E1                    # C++11 class enum, scoped
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    AttributeError: <namespace cppyy.gbl at 0x7ff2766a4af0> has no attribute 'E1'.
    >>> cppyy.gbl.NamedClassEnum.E1
    42
    >>>

