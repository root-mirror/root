.. _functions:


Functions
=========

C++ functions are first-class objects in Python and can be used wherever
Python functions can be used, including for dynamically constructing
classes.

The C++ code used for the examples below can be found
:doc:`here <cppyy_features_header>`, and it is assumed that that code is
loaded at the start of any session.
Download it, save it under the name ``features.h``, and load it:

  .. code-block:: python

    >>> import cppyy
    >>> cppyy.include('features.h')
    >>>

Function argument type conversions follow the expected rules, with implicit
conversions allowed, including between Python builtin types and STL types,
but it is rather more efficient to make conversions explicit.


`Free functions`
----------------

All bound C++ code starts off from the global C++ namespace, represented in
Python by ``gbl``.
This namespace, as any other namespace, is treated as a module after it has
been loaded.
Thus, we can directly import C++ functions from it and other namespaces that
themselves may contain more functions.
All lookups on namespaces are done lazily, thus if loading more headers bring
in more functions (incl. new overloads), these become available dynamically.

  .. code-block:: python

    >>> from cppyy.gbl import global_function, Namespace
    >>> global_function == Namespace.global_function
    False
    >>> from cppyy.gbl.Namespace import global_function
    >>> global_function == Namespace.global_function
    True
    >>> from cppyy.gbl import global_function
    >>>

Free functions can be bound to a class, following the same rules as apply to
Python functions: unless marked as static, they will turn into member
functions when bound to an instance, but act as static functions when called
through the class.
Consider this example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete, call_abstract_method
    >>> c = Concrete()
    >>> Concrete.callit = call_abstract_method
    >>> Concrete.callit(c)
    called Concrete::abstract_method
    >>> c.callit()
    called Concrete::abstract_method
    >>> Concrete.callit = staticmethod(call_abstract_method)
    >>> c.callit()
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: void ::call_abstract_method(Abstract* a) =>
        TypeError: takes at least 1 arguments (0 given)
    >>> c.callit(c)
    called Concrete::abstract_method
    >>>


`Static methods`
----------------

Class static functions are treated the same way as free functions, except
that they are accessible either through the class or through an instance,
just like Python's ``staticmethod``.


`Methods`
---------

For class methods, see the :ref:`methods section <sec-methods-label>` under
the :doc:`classes heading<classes>`.


`Operators`
-----------

Globally defined operators are found lazily (ie. can resolve after the class
definition by loading the global definition or by defining them interactively)
and are mapped onto a Python equivalent when possible.
See the :ref:`operators section <sec-operators-label>` under the
:doc:`classes heading<classes>` for more details.


`Templates`
-----------

Templated functions (and class methods) can either be called using square
brackets (``[]``) to provide the template arguments explicitly, or called
directly, through automatic lookup.
The template arguments may either be a string of type names (this results
in faster code, as it needs no further lookup/verification) or a list of
the actual types to use (which tends to be more convenient).

**Note**: the Python type ``float`` maps to the C++ type ``float``, even
as Python uses a C ``double`` as its internal representation.
The motivation is that doing so makes the Python code more readable (and
Python may anyway change its internal representation in the future).
The same has been true for Python ``int``, which used to be a C ``long``
internally.

Examples, using multiply from :doc:`features.h <cppyy_features_header>`:

  .. code-block:: python

   >>> mul = cppyy.gbl.multiply
   >>> mul(1, 2)
   2
   >>> mul(1., 5)
   5.0
   >>> mul[int](1, 1)
   1
   >>> mul[int, int](1, 1)
   1
   >>> mul[int, int, float](1, 1)
   1.0
   >>> mul[int, int](1, 'a')
    TypeError: Template method resolution failed:
    none of the 6 overloaded methods succeeded. Full details:
    int ::multiply(int a, int b) =>
      TypeError: could not convert argument 2 (int/long conversion expects an integer object)
    ...
    Failed to instantiate "multiply(int,std::string)"
   >>> mul['double, double, double'](1., 5)
   5.0
   >>>


`Overloading`
-------------

C++ supports overloading, whereas Python supports "duck typing", thus C++
overloads have to be selected dynamically in response to the available
"ducks".
This may lead to additional lookups or template instantiations.
However, pre-existing methods (incl. auto-instantiated methods) are always
preferred over new template instantiations:

  .. code-block:: python

    >>> global_function(1.)        # selects 'double' overload
    2.718281828459045
    >>> global_function(1)         # selects 'int' overload
    42
    >>>

C++ does a static dispatch at compile time based on the argument types.
The dispatch is a selection among overloads (incl. templates) visible at that
point in the *translation unit*.
Bound C++ in Python does a dynamic dispatch: it considers all overloads
visible *globally* at that point in the execution.
Because the dispatch is fundamentally different (albeit in line with the
expectation of the respective languages), differences can occur.
Especially if overloads live in different header files and are only an
implicit conversion apart, or if types that have no direct equivalent in
Python, such as e.g. ``unsigned short``, are used.

There are two rounds to finding an overload.
If all overloads fail argument conversion during the first round, where
implicit conversions are not allowed, _and_ at least one converter has
indicated that it can do implicit conversions, a second round is tried.
In this second round, implicit conversions are allowed, including class
instantiation of temporaries.
During some template calls, implicit conversions are not allowed at all, to
make sure new instantiations happen instead.

In the rare occasion where the automatic overload selection fails, the
``__overload__`` function can be called to access a specific overload
matching a specific function signature:

  .. code-block:: python

     >>> global_function.__overload__('double')(1)   # int implicitly converted
     2.718281828459045
     >>>

Note that ``__overload__`` only does a lookup; it performs no (implicit)
conversions.
To see all available overloads, use ``help()`` or look at the ``__doc__``
string of the function:

  .. code-block:: python

     >>> print(global_function.__doc__)
     int ::global_function(int)
     double ::global_function(double)
     >>>


`Return values`
---------------

Most return types are readily amenable to automatic memory management: builtin
returns, by-value returns, (const-)reference returns to internal data, smart
pointers, etc.
The important exception is pointer returns.
 
A function that returns a pointer to an object over which Python should claim
ownership, should have its ``__creates__`` flag set through its
:doc:`pythonization <pythonizations>`.
Well-written APIs will have clear clues in their naming convention about the
ownership rules.
For example, functions called ``New...``, ``Clone...``, etc.  can be expected
to return freshly allocated objects.
A simple name-matching in the pythonization then makes it simple to mark all
these functions as creators.

The return values are :ref:`auto-casted <sec-auto-casting-label>`.


`\*args and \*\*kwds`
---------------------

C++ default arguments work as expected.
Keywords, however, are a Python language feature that does not exist in C++.
Many C++ function declarations do have formal arguments, but these are not
part of the C++ interface (the argument names are repeated in the definition,
making the names in the declaration irrelevant: they do not even need to be
provided).
Thus, although ``cppyy`` will map keyword argument names to formal argument
names from the C++ declaration, use of this feature is not recommended unless
you have a guarantee that the names in C++ the interface are maintained.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import Concrete
    >>> c = Concrete()       # uses default argument
    >>> c.m_int
    42
    >>> c = Concrete(13)     # uses provided argument
    >>> c.m_int
    13
    >>> args = (27,)
    >>> c = Concrete(*args)  # argument pack
    >>> c.m_int
    27
    >>> c = Concrete(n=17)
    >>> c.m_int
    17
    >>> kwds = {'n' : 18}
    >>> c = Concrete(**kwds)
    >>> c.m_int
    18
    >>>


`Callbacks`
-----------

Python callables (functions/lambdas/instances) can be passed to C++ through
function pointers and/or ``std::function``.
This involves creation of a temporary wrapper, which has the same life time as
the Python callable it wraps, so the callable needs to be kept alive on the
Python side if the C++ side stores the callback.
Example:

  .. code-block:: python

    >>> from cppyy.gbl import call_int_int
    >>> print(call_int_int.__doc__)
    int ::call_int_int(int(*)(int,int) f, int i1, int i2)
    >>> def add(a, b):
    ...    return a+b
    ...
    >>> call_int_int(add, 3, 7)
    7
    >>> call_int_int(lambda x, y: x*y, 3, 7)
    21
    >>>

