Low-level code
==============

.. toctree::
   :hidden:

C code and older C++ code sometimes makes use of low-level features such as
pointers to builtin types, some of which do not have any Python equivalent
(e.g. ``unsigned short*``).
Furthermore, such codes tend to be ambiguous: the information from header
file is not sufficient to determine the full purpose.
For example, an ``int*`` type may refer to the address of a single ``int``
(an out-parameter, say) or it may refer to an array of ``int``, the ownership
of which is not clear either.
cppyy provides a few low-level helpers and integration with the Python
`ctypes module`_ to cover these cases.

Use of these low-level helpers will obviously lead to very "C-like" code and
it is recommended to :doc:`pythonize <pythonizations>` the code, perhaps
using the Cling JIT and embedded C++.

Note: the low-level module is not loaded by default (since its use is, or
should be, uncommon).
It needs to be imported explicitly:

  .. code-block:: python

    >>> import cppyy.ll
    >>>


`C/C++ casts`
-------------

C++ instances are auto-casted to the most derived available type, so do not
require explicit casts even when a function returns a pointer to a base
class or interface.
However, when given only a ``void*`` or ``intptr_t`` type on return, a cast
is required to turn it into something usable.

* **bind_object**: This is the preferred method to proxy a C++ address,
  and lives in ``cppyy``, not ``cppyy.ll``, as it is not a low-level C++
  cast, but a ``cppyy`` API that is also used internally.
  It thus plays well with object identity, references, etc.
  Example:

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... struct MyStruct { int fInt; };
    ... void* create_mystruct() { return new MyStruct{42}; }
    ... """)
    ... 
    >>> s = cppyy.gbl.create_mystruct()
    >>> print(s)
    <cppyy.LowLevelView object at 0x10559d430>
    >>> sobj = cppyy.bind_object(s, 'MyStruct')
    >>> print(sobj)
    <cppyy.gbl.MyStruct object at 0x7ff25e28eb20>
    >>> print(sobj.fInt)
    42
    >>>

  Instead of the type name as a string, ``bind_object`` can also take the
  actual class (here: ``cppyy.gbl.MyStruct``).

* **Typed nullptr**: A Python side proxy can pass through a pointer to
  pointer function argument, but if the C++ side allocates memory and
  stores it in the pointer, the result is a memory leak.
  In that case, use ``bind_object`` to bind ``cppyy.nullptr`` instead, to
  get a typed nullptr to pass to the function.
  Example (continuing from the example above):

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... void create_mystruct(MyStruct** ptr) { *ptr = new MyStruct{42}; }
    ... """)
    ...
    >>> s = cppyy.bind_object(cppyy.nullptr, 'MyStruct')
    >>> print(s)
    <cppyy.gbl.MyStruct object at 0x0>
    >>> cppyy.gbl.create_mystruct(s)
    >>> print(s)
    <cppyy.gbl.MyStruct object at 0x7fc7d85b91c0>
    >>> print(s.fInt)
    42
    >>>

* **C-style cast**: This is the simplest option for builtin types.
  The syntax is "template-style", example:

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... void* get_data(int sz) {
    ... int* iptr = (int*)malloc(sizeof(int)*sz);
    ... for (int i=0; i<sz; ++i) iptr[i] = i;
    ... return iptr;
    ... }""")
    ...
    >>> NDATA = 4
    >>> d = cppyy.gbl.get_data(NDATA)
    >>> print(d)
    <cppyy.LowLevelView object at 0x1068cba30>
    >>> d = cppyy.ll.cast['int*'](d)
    >>> d.reshape((NDATA,))
    >>> print(list(d))
    [0, 1, 2, 3]
    >>>

* **C++-style casts**: Similar to the C-style cast, there are
  ``ll.static_cast`` and ``ll.reinterpret_cast``.
  There should never be a reason for a ``dynamic_cast``, since that only
  applies to objects, for which auto-casting will work.
  The syntax is "template-style", just like for the C-style cast above.


`NumPy casts`
-------------

The ``cppyy.LowLevelView`` type returned for pointers to basic types,
including for ``void*``, is a simple and light-weight view on memory, given a
pointer, type, and number of elements (or unchecked, if unknown).
It only supports basic operations such as indexing and iterations, but also
the buffer protocol for integration with full-fledged functional arrays such
as NumPy`s ``ndarray``.

In addition, specifically when dealing with ``void*`` returns, you can use
NumPy's low-level ``frombuffer`` interface to perform the cast.
Example:

  .. code-block:: python

     >>> cppyy.cppdef("""
     ... void* create_float_array(int sz) {
     ...     float* pf = (float*)malloc(sizeof(float)*sz);
     ...     for (int i = 0; i < sz; ++i) pf[i] = 2*i;
     ...     return pf;
     ... }""")
     ...
     >>> import numpy as np
     >>> NDATA = 8
     >>> arr = cppyy.gbl.create_float_array(NDATA)
     >>> print(arr)
     <cppyy.LowLevelView object at 0x109f15230>
     >>> arr.reshape((NDATA,))   # adjust the llv's size
     >>> v = np.frombuffer(arr, dtype=np.float32, count=NDATA)  # cast to float
     >>> print(len(v))
     8
     >>> print(v)
     array([ 0.,  2.,  4.,  6.,  8., 10., 12., 14.], dtype=float32)
     >>>

Note that NumPy will internally check the total buffer size, so if the size
you are casting *to* is larger than the size you are casting *from*, then
the number of elements set in the ``reshape`` call needs to be adjusted
accordingly.


`Capsules`
----------

It is not possible to pass proxies from cppyy through function arguments of
another binder (and vice versa, with the exception of ``ctypes``, see below),
because each will use a different internal representation, including for type
checking and extracting the C++ object address.
However, all Python binders are able to rebind (just like ``bind_object``
above for cppyy) the result of at least one of the following:

* **ll.addressof**: Takes a cppyy bound C++ object and returns its address as
  an integer value.
  Takes an optional ``byref`` parameter and if set to true, returns a pointer
  to the address instead.

* **ll.as_capsule**: Takes a cppyy bound C++ object and returns its address as
  a PyCapsule object.
  Takes an optional ``byref`` parameter and if set to true, returns a pointer
  to the address instead.

* **ll.as_cobject**: Takes a cppyy bound C++ object and returns its address as
  a PyCObject object for Python2 and a PyCapsule object for Python3.
  Takes an optional ``byref`` parameter and if set to true, returns a pointer
  to the address instead.

* **as_ctypes**: Takes a cppyy bound C++ object and returns its address as
  a ``ctypes.c_void_p`` object.
  Takes an optional ``byref`` parameter and if set to true, returns a pointer
  to the address instead.


`ctypes`
--------

The `ctypes module`_ has been part of Python since version 2.5 and provides a
Python-side foreign function interface.
It is clunky to use and has very bad performance, but it is guaranteed to be
available.
It does not have a public C interface, only the Python one, but its internals
have been stable since its introduction, making it safe to use for tight and
efficient integration at the C level (with a few Python helpers to assure
lazy lookup).

Objects from ``ctypes`` can be passed through arguments of functions that
take a pointer to a single C++ builtin, and ``ctypes`` pointers can be passed 
when a pointer-to-pointer is expected, e.g. for array out-parameters.
This leads to the following set of possible mappings:

========================================  ========================================
C++                                       ctypes
========================================  ========================================
by value (ex.: ``int``)                   ``.value`` (ex.: ``c_int(0).value``)
by const reference (ex.: ``const int&``)  ``.value`` (ex.: ``c_int(0).value``)
by reference (ex.: ``int&``)              direct (ex.: ``c_int(0)``)
by pointer (ex.: ``int*``)                direct (ex.: ``c_int(0)``)
by ptr-ref (ex.: ``int*&``)               ``pointer`` (ex.: ``pointer(c_int(0))``)
by ptr-ptr **in** (ex.: ``int**``)        ``pointer`` (ex.: ``pointer(c_int(0))``)
by ptr-ptr **out** (ex.: ``int**``)       ``POINTER`` (ex.: ``POINTER(c_int)()``)
========================================  ========================================

The ``ctypes`` pointer objects (from ``POINTER``, ``pointer``, or ``byref``)
can also be used for pass by reference or pointer, instead of the direct
object, and ``ctypes.c_void_p`` can pass through all pointer types.
The addresses will be adjusted internally by cppyy.

Note that ``ctypes.c_char_p`` is expected to be a NULL-terminated C string,
not a character array (see the `ctypes module`_ documentation), and that
``ctypes.c_bool`` is a C ``_Bool`` type, not C++ ``bool``.


`Memory`
--------

C++ has three ways of allocating heap memory (``malloc``, ``new``, and
``new[]``) and three corresponding ways of deallocation (``free``,
``delete``, and ``delete[]``).
Direct use of ``malloc`` and ``new`` should be avoided for C++ classes, as
these may override ``operator new`` to control their allocation own.
However these low-level allocators can be necessary for builtin types on
occassion if the C++ side takes ownership (otherwise, prefer either
``array`` from the builtin module ``array`` or ``ndarray`` from Numpy).

The low-level module adds the following functions:

* **ll.malloc**: an interface on top of C's malloc.
  Use it as a template with the number of elements (not the number types) to
  be allocated.
  The result is a ``cppyy.LowLevelView`` with the proper type and size:

  .. code-block:: python

    >>> arr = cppyy.ll.malloc[int](4)   # allocates memory for 4 C ints
    >>> print(len(arr))
    4
    >>> print(type(arr[0]))
    <type 'int'>
    >>>

  The actual C malloc can also be used directly, through ``cppyy.gbl.malloc``,
  taking the number of *bytes* to be allocated and returning a ``void*``.

* **ll.free**: an interface to C's free, to deallocate memory allocated by
  C's malloc.
  To continue to example above:

  .. code-block:: python

    >>> cppyy.ll.free(arr)
    >>>

  The actual C free can also be used directly, through ``cppyy.gbl.free``.

* **ll.array_new**: an interface on top of C++'s ``new[]``.
  Use it as a template; the result is a ``cppyy.LowLevelView`` with the
  proper type and size:

  .. code-block:: python

    >>> arr = cppyy.ll.array_new[int](4)   # allocates memory for 4 C ints
    >>> print(len(arr))
    4
    >>> print(type(arr[0]))
    <type 'int'>
    >>>

* **ll.array_delete**: an interface on top of C++'s ``delete[]``.
  To continue to example above:

  .. code-block:: python

    >>> cppyy.ll.array_delete(arr)
    >>>


.. _`ctypes module`: https://docs.python.org/3/library/ctypes.html
