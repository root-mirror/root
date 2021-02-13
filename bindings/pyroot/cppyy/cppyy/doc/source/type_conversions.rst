.. _type_conversions:


Type conversions
================

Most type conversions are done automatically, e.g. between Python ``str``
and C++ ``std::string`` and ``const char*``, but low-level APIs exist to
perform explicit conversions.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>


 .. _sec-auto-casting-label:

`Auto-casting`
--------------

Object pointer returns from functions provide the most derived class known
(i.e. exposed in header files) in the hierarchy of the object being returned.
This is important to preserve object identity as well as to make casting,
a pure C++ feature after all, superfluous.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Abstract, Concrete
    >>> c = Concrete()
    >>> Concrete.show_autocast.__doc__
    'Abstract* Concrete::show_autocast()'
    >>> d = c.show_autocast()
    >>> type(d)
    <class '__main__.Concrete'>
    >>>

As a consequence, if your C++ classes should only be used through their
interfaces, then no bindings should be provided to the concrete classes
(e.g. by excluding them using a :ref:`selection file <selection-files>`).
Otherwise, more functionality will be available in Python than in C++.

Sometimes, however, full control over a cast is needed.
For example, if the instance is bound by another tool or even a 3rd party,
hand-written, extension library.
Assuming the object supports the ``PyCapsule`` or ``CObject`` abstraction,
then a C++-style reinterpret_cast (i.e. without implicitly taking offsets
into account), can be done by taking and rebinding the address of an
object:

  .. code-block:: python

    >>> from cppyy import addressof, bind_object
    >>> e = bind_object(addressof(d), Abstract)
    >>> type(e)
    <class '__main__.Abstract'>
    >>>


`Operators`
-----------

If conversion operators are defined in the C++ class and a Python equivalent
exists (i.e. all builtin integer and floating point types, as well as
``bool``), then these will map onto those Python conversions.
Note that ``char*`` is mapped onto ``__str__``.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> print(Concrete())
    Hello operator const char*!
    >>>

C++ code can overload conversion operators by providing methods in a class or
global functions.
Special care needs to be taken for the latter: first, make sure that they are
actually available in some header file.
Second, make sure that headers are loaded in the desired order.
I.e. that these global overloads are available before use.

