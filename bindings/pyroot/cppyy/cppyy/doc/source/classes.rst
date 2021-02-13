.. _classes:


Classes
=======

Both Python and C++ support object-oriented code through classes and thus
it is logical to expose C++ classes as Python ones, including the full
inheritance hierarchy.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>


`Basics`
--------

All bound C++ code starts off from the global C++ namespace, represented in
Python by ``gbl``.
This namespace, as any other namespace, is treated as a module after it has
been loaded.
Thus, we can import C++ classes that live underneath it:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> Concrete
    <class cppyy.gbl.Concrete at 0x2058e30>
    >>>

Placing classes in the same structure as imposed by C++ guarantees identity,
even if multiple Python modules bind the same class.
There is, however, no necessity to expose that structure to end-users: when
developing a Python package that exposes C++ classes through ``cppyy``,
consider ``cppyy.gbl`` an "internal" module, and expose the classes in any
structure you see fit.
The C++ names will continue to follow the C++ structure, however, as is needed
for e.g. pickling:

  .. code-block:: python

    >>> from cppyy.gbl import Namespace
    >>> Concrete == Namespace.Concrete
    False
    >>> n = Namespace.Concrete.NestedClass()
    >>> type(n)
    <class cppyy.gbl.Namespace.Concrete.NestedClass at 0x22114c0>
    >>> type(n).__name__
    NestedClass
    >>> type(n).__module__
    cppyy.gbl.Namespace.Concrete
    >>> type(n).__cpp_name__
    Namespace::Concrete::NestedClass
    >>>


`Constructors`
--------------

Python and C++ both make a distinction between allocation (``__new__`` in
Python, ``operator new`` in C++) and initialization (``__init__`` in Python,
the constructor call in C++).
When binding, however, there comes a subtle semantic difference: the Python
``__new__`` allocates memory for the proxy object only, and ``__init__``
initializes the proxy by creating or binding the C++ object.
Thus, no C++ memory is allocated until ``__init__``.
The advantages are simple: the proxy can now check whether it is initialized,
because the pointer to C++ memory will be NULL if not; it can be a reference
to another proxy holding the actual C++ memory; and it can now transparently
implement a C++ smart pointer.
If ``__init__`` is never called, eg. when a call to the base class
``__init__`` is missing in a derived class override, then accessing the proxy
will result in a Python ``ReferenceError`` exception.


`Destructors`
-------------

There should not be a reason to call a destructor directly in CPython, but
PyPy uses a garbage collector and that makes it sometimes useful to destruct
a C++ object where you want it destroyed.
Destructors are accessible through the conventional ``__destruct__`` method.
Accessing an object after it has been destroyed will result in a Python
``ReferenceError`` exception.


`Inheritance`
-------------

The output of help shows the inheritance hierarchy, constructors, public
methods, and public data.
For example, ``Concrete`` inherits from ``Abstract`` and it has
a constructor that takes an ``int`` argument, with a default value of 42.
Consider:

  .. code-block:: python

    >>> from cppyy.gbl import Abstract
    >>> issubclass(Concrete, Abstract)
    True
    >>> a = Abstract()
    Traceback (most recent call last):
      File "<console>", line 1, in <module>
    TypeError: cannot instantiate abstract class 'Abstract'
    >>> c = Concrete()
    >>> isinstance(c, Concrete)
    True
    >>> isinstance(c, Abstract)
    True
    >>> d = Concrete(13)
    >>>

Just like in C++, interface classes that define pure virtual methods, such
as ``Abstract`` does, can not be instantiated, but their concrete
implementations can.
As the output of ``help`` showed, the ``Concrete`` constructor takes
an integer argument, that by default is 42.


`Cross-inheritance`
-------------------

Python classes that derive from C++ classes can override virtual methods as
long as those methods are declared on class instantiation (adding methods to
the Python class after the fact will not provide overrides on the C++ side,
only on the Python side).
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Abstract, call_abstract_method
    >>> class PyConcrete(Abstract):
    ...     def abstract_method(self):
    ...         print("Hello, Python World!\n")
    ...     def concrete_method(self):
    ...         pass
    ...
    >>> pc = PyConcrete()
    >>> call_abstract_method(pc)
    Hello, Python World!
    >>>

Note that it is not necessary to provide a constructor (``__init__``), but
if you do, you *must* call the base class constructor through the ``super``
mechanism.


 .. _sec-methods-label:

`Methods`
---------

C++ methods are represented as Python ones: these are first-class objects and
can be bound to an instance.
If a method is virtual in C++, the proper concrete method is called, whether
or not the concrete class is bound.
Similarly, if all classes are bound, the normal Python rules apply:

  .. code-block:: python

    >>> c.abstract_method()
    called Concrete::abstract_method
    >>> c.concrete_method()
    called Concrete::concrete_method
    >>> m = c.abstract_method
    >>> m()
    called Concrete::abstract_method
    >>>


`Data members`
--------------

Data members are implemented as properties, using descriptors.
For example, The ``Concrete`` instances have a public data member ``m_int``:

  .. code-block:: python

    >>> c.m_int, d.m_int
    (42, 13)
    >>>

Note however, that the data members are typed: setting them results in a
memory write on the C++ side.
This is different in Python, where references are replaced, and thus any
type will do:

  .. code-block:: python

    >>> c.m_int = 3.14   # a float does not fit in an int
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: int/long conversion expects an integer object
    >>> c.m_int = int(3.14)
    >>> c.m_int, d.m_int
    (3, 13)
    >>>

Private and protected data members are not accessible, contrary to Python
data members, and C++ const-ness is respected:

  .. code-block:: python

    >>> c.m_const_int = 71    # declared 'const int' in class definition
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: assignment to const data not allowed
    >>>

Static C++ data members act like Python class-level data members.
They are also represented by property objects and both read and write access
behave as expected:

  .. code-block:: python

    >>> Concrete.s_int       # access through class
    321
    >>> c.s_int = 123        # access through instance
    >>> Concrete.s_int
    123


 .. _sec-operators-label:

`Operators`
-----------

Many C++ operators can be mapped to their Python equivalent.
When the operators are part of the C++ class definition, this is done
directly.
If they are defined globally, the lookup is done lazily (ie. can resolve
after the class definition by loading the global definition or by defining
them interactively).
Some operators have no Python equivalent and are instead made available by
mapping them onto the following conventional functions:

===================  ===================
C++                  Python
===================  ===================
``operator=``        ``__assign__``
``operator++(int)``  ``__postinc__``
``operator++()``     ``__preinc__``
``operator--(int)``  ``__postdec__``
``operator--()``     ``__predec__``
``unary operator*``  ``__deref__``
``operator->``       ``__follow__``
``operator&&``       ``__dand__``
``operator||``       ``__dor__``
``operator,``        ``__comma__``
===================  ===================

Here is an example of operator usage, using STL iterators directly (note that
this is not necessary in practice as STL and STL-like containers work
transparently in Python for-loops):

  .. code-block:: python

    >>> v = cppyy.gbl.std.vector[int](range(3))
    >>> i = v.begin()
    >>> while (i != v.end()):
    ...    print(i.__deref__())
    ...    _ = i.__preinc__()
    ...
    0
    1
    2
    >>>

Overridden ``operator new`` and ``operator delete``, as well as their array
equivalents, are not accessible but will be called as appropriate.


`Templates`
-----------

Templated classes are instantiated using square brackets.
(For backwards compatibility reasons, parentheses work as well.)
The instantiation of a templated class yields a class, which can then
be used to create instances.

Templated classes need not pre-exist in the bound code, just their
declaration needs to be available.
This is true for e.g. all of STL:

  .. code-block:: python

    >>> cppyy.gbl.std.vector                # template metatype
    <cppyy.Template 'std::vector' object at 0x7fffed2674d0>
    >>> cppyy.gbl.std.vector(int)           # instantiates template -> class
    <class cppyy.gbl.std.vector<int> at 0x1532190>
    cppyy.gbl.std.vector[int]()             # instantiates class -> object
    <cppyy.gbl.std.vector<int> object at 0x2341ec0>
    >>>

The template arguments may be actual types or their names as a string,
whichever is more convenient.
Thus, the following are equivalent:

  .. code-block:: python

     >>> from cppyy.gbl.std import vector
     >>> type1 = vector[Concrete]
     >>> type2 = vector['Concrete']
     >>> type1 == type2
     True
     >>>


`Typedefs`
----------

Typedefs are simple python references to the actual classes to which
they refer.

  .. code-block:: python

    >>> from cppyy.gbl import Concrete_t
    >>> Concrete is Concrete_t
    True
    >>>

