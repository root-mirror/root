.. _python:


Python
======

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>


`PyObject`
----------

Arguments and return types of ``PyObject*`` can be used, and passed on to
CPython API calls (or through ``cpyext`` in PyPy).


`Doc strings`
-------------

The documentation string of a method or function contains the C++
arguments and return types of all overloads of that name, as applicable.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> print Concrete.array_method.__doc__
    void Concrete::array_method(int* ad, int size)
    void Concrete::array_method(double* ad, int size)
    >>>


`Help`
------

Bound C++ class is first-class Python and can thus be inspected like any
Python objects can.
For example, we can ask for ``help()``:

  .. code-block:: python

    >>> help(Concrete)
    Help on class Concrete in module gbl:

    class Concrete(Abstract)
     |  Method resolution order:
     |      Concrete
     |      Abstract
     |      CPPInstance
     |      __builtin__.object
     |
     |  Methods defined here:
     |
     |  __assign__(self, const Concrete&)
     |      Concrete& Concrete::operator=(const Concrete&)
     |
     |  __init__(self, *args)
     |      Concrete::Concrete(int n = 42)
     |      Concrete::Concrete(const Concrete&)
     |
     etc. ....

