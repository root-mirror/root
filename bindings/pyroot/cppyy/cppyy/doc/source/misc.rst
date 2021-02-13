.. _features:

Miscellaneous
=============

.. toctree::
   :hidden:

   cppyy_features_header

This is a collection of a few more features listed that do not have a proper
place yet in the rest of the documentation.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>


`Special variables`
-------------------

There are several conventional "special variables" that control behavior of
functions or provide (internal) information.
Often, these can be set/used in pythonizations to handle memory management or
Global Interpreter Lock (GIL) release.

* ``__python_owns__``: a flag that every bound instance carries and determines
  whether Python or C++ owns the C++ instance (and associated memory).
  If Python owns the instance, it will be destructed when the last Python
  reference to the proxy disappears.
  You can check/change the ownership with the __python_owns__ flag that every
  bound instance carries.
  Example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> c = Concrete()
    >>> c.__python_owns__         # True: object created in Python
    True
    >>>

* ``__creates__``: a flag that every C++ overload carries and determines
  whether the return value is owned by C++ or Python: if ``True``, Python owns
  the return value, otherwise C++.

* ``__set_lifeline__``: a flag that every C++ overload carries and determines
  whether the return value should place a back-reference on ``self``, to
  prevent the latter from going out of scope before the return value does.
  The default is ``False``, but will be automatically set at run-time if a
  return value's address is a C++ object pointing into the memory of ``this``,
  or if ``self`` is a by-value return.

* ``__release_gil__``: a flag that every C++ overload carries and determines
  whether the Global Interpreter Lock (GIL) should be released during the C++
  call to allow multi-threading.
  The default is ``False``.

* ``__useffi__``: a flag that every C++ overload carries and determines
  whether generated wrappers or direct foreign functions should be used.
  This is for PyPy only; the flag has no effect on CPython.

* ``__sig2exc__``: a flag that every C++ overload carries and determines
  whether C++ signals (such as SIGABRT) should be converted into Python
  exceptions.

* ``__cppname__``: a string that every C++ bound class carries and contains
  the actual C++ name (as opposed to ``__name__`` which has the Python name).
  This can be useful for template instantiations, documentation, etc.


`STL algorithms`
----------------

It is usually easier to use a Python equivalent or code up the effect of an
STL algorithm directly, but when operating on a large container, calling an
STL algorithm may offer better performance.
It is important to note that all STL algorithms are templates and need the
correct types to be properly instantiated.
STL containers offer typedefs to obtain those exact types and these should
be used rather than relying on the usual implicit conversions of Python types
to C++ ones.
For example, as there is no ``char`` type in Python, the ``std::remove`` call
below can not be instantiated using a Python string, but the
``std::string::value_type`` must be used instead:

  .. code-block:: python

    >>> cppstr = cppyy.gbl.std.string
    >>> n = cppstr('this is a C++ string')
    >>> print(n)
    this is a C++ string
    >>> n.erase(cppyy.gbl.std.remove(n.begin(), n.end(), cppstr.value_type(' ')))
    <cppyy.gbl.__wrap_iter<char*> object at 0x7fba35d1af50>
    >>> print(n)
    thisisaC++stringing
    >>>


`Reduced typing`
----------------

Typing ``cppyy.gbl`` all the time gets old rather quickly, but the dynamic
nature of ``cppyy`` makes something like ``from cppyy.gbl import *``
impossible.
For example, classes can be defined dynamically after that statement and then
they would be missed by the import.
In scripts, it is easy enough to rebind names to achieve a good amount of
reduction in typing (and a modest performance improvement to boot, because of
fewer dictionary lookups), e.g.:

  .. code-block:: python

    import cppyy
    std = cppyy.gbl.std
    v = std.vector[int](range(10))

But even such rebinding becomes annoying for (brief) interactive sessions.

For CPython only (and not with tools such as IPython or in IDEs that replace
the interactive prompt), there is a fix, using
``from cppyy.interactive import *``.
This makes lookups in the global dictionary of the current frame also
consider everything under ``cppyy.gbl``.
This feature comes with a performance `penalty` and is not meant for
production code.
Example usage:

  .. code-block:: python

    >>> from cppyy.interactive import *
    >>> v = std.vector[int](range(10))
    >>> print(list(v))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>>
    >>> cppdef("struct SomeStruct {};")
    True
    >>> s = SomeStruct()          # <- dynamically made available
    >>> s
    <cppyy.gbl.SomeStruct object at 0x7fa9b8624320>
    >>>

For PyPy, IPython, etc. ``cppyy.gbl`` is simply rebound as ``g`` and
``cppyy.gbl.std`` is made available as ``std``.
Not as convenient as full lookup, and missing any other namespaces that may be
available, but still saves some typing in may cases.


`Odds and ends`
---------------

* **namespaces**: Are represented as python classes.
  Namespaces are more open-ended than classes, so sometimes initial access may
  result in updates as data and functions are looked up and constructed
  lazily.
  Thus the result of ``dir()`` on a namespace shows the classes and functions
  available for binding, even if these may not have been created yet.
  Once created, namespaces are registered as modules, to allow importing from
  them.
  The global namespace is ``cppyy.gbl``.

* **NULL**: Is represented as ``cppyy.nullptr``.
  Starting C++11, the keyword ``nullptr`` is used to represent ``NULL``.
  For clarity of intent, it is recommended to use this instead of ``None``
  (or the integer ``0``, which can serve in some cases), as ``None`` is better
  understood as ``void`` in C++.
