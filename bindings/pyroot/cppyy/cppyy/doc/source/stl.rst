.. _stl:


STL
===

Parts of the Standard Template Library (STL), in particular its container
types, are the de facto equivalent of Python's builtin types.
STL is written in C++ and Python bindings of it are fully functional as-is,
but are much more useful when pluggable into idiomatic expressions where
Python builtin containers are expected (e.g. in list contractions).

There are two extremes to achieve such drop-in behavior: copy into Python
builtins, so that the Python-side always deals with true Python objects; or
adjust the C++ interfaces to be the same as their Python equivalents.
Neither is very satisfactory: the former is not because of the existence of
global/static variables and return-by-reference.
If only a copy is available, then expected modifications do not propagate.
Copying is also either slow (when copying every time) or memory intensive (if
the results are cached).
Filling out the interfaces may look more appealing, but all operations then
involve C++ function calls, which can be slower than the Python equivalents,
and C++-style error handling.

Given that neither choice will satisfy all cases, ``cppyy`` aims to maximize
functionality and minimum surprises based on common use.
Thus, for example, ``std::vector`` grows a pythonistic ``__len__`` method,
but does not lose its C++ ``size`` method.
Passing a Python container through a const reference to a ``std::vector``
will trigger automatic conversion, but such an attempt through a non-const
reference will fail since a non-temporary C++ object is required [#f1]_ to
return any updates/changes.

``std::string`` is almost always converted to Python's ``str`` on function
returns (the exception is return-by-reference when assigning), but not when
its direct use is more likely such as in the case of (global) variables or
when iterating over a ``std::vector<std::string>``.

The rest of this section shows examples of how STL containers can be used in
a natural, pythonistic, way.


`vector`
--------

A ``std::vector`` is the most commonly used C++ container type because it is
more efficient and performant than specialized types such as ``list`` and
``map``, unless the number of elements gets very large.
Python has several similar types, from the builtin ``tuple`` and ``list``,
the ``array`` from builtin module ``array``, to "as-good-as-builtin"
``numpy.ndarray``.
A vector is more like the latter two in that it can contain only one type,
but more like the former two in that it can contain objects.
In practice, it can interplay well with all these containers, but e.g.
efficiency and performance can differ significantly.

A vector can be instantiated from any sequence, including generators, and
vectors of objects can be recursively constructed:

  .. code-block:: python

    >>> from cppyy.gbl.std import vector, pair
    >>> v = vector[int](range(10))
    >>> len(v)
    10
    >>> vp = vector[pair[int, int]](((1, 2), (3, 4)))
    >>> len(vp)
    2
    >>> vp[1][0]
    3
    >>>

To extend a vector in-place with another sequence object, use ``+=``, just as
would work for Python's list:

  .. code-block:: python

    >>> v += range(10, 20)
    >>> len(v)
    20
    >>>
    
The easiest way to print the full contents of a vector, is by using a list
and printing that instead.
Indexing and slicing of a vector follows the normal Python slicing rules:

  .. code-block:: python

    >>> v[1]
    1
    >>> v[-1]
    19
    >>> v[-4:]
    <cppyy.gbl.std.vector<int> object at 0x7f9051057650>
    >>> list(v[-4:])
    [16, 17, 18, 19]
    >>>

The usual iteration operations work on vector, but the C++ rules still apply,
so a vector that is being iterated over can `not` be modified in the loop
body.
(On the plus side, this makes it much faster to iterate over a vector than,
say, a numpy ndarray.)

  .. code-block:: python

    >>> for i in v[2:5]:
    ...     print(i)
    ...
    2
    3
    4
    >>> 2 in v
    True
    >>> sum(v)
    190
    >>>

When a function takes a non-l-value (const-ref, move, or by-value) vector as
a parameter, another sequence can be used and cppyy will automatically
generate a temporary.
Typically, this will be faster than coding up such a temporary on the Python
side, but if the same sequence is used multiple times, creating a temporary
once and re-using it will be the most efficient approach.o

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... int sumit1(const std::vector<int>& data) {
    ...   return std::accumulate(data.begin(), data.end(), 0);
    ... }
    ... int sumit2(std::vector<int> data) {
    ...   return std::accumulate(data.begin(), data.end(), 0);
    ... }
    ... int sumit3(const std::vector<int>&& data) {
    ...   return std::accumulate(data.begin(), data.end(), 0);
    ... }""")
    ...
    True
    >>> cppyy.gbl.sumit1(range(5))
    10
    >>> cppyy.gbl.sumit2(range(6))
    16
    >>> cppyy.gbl.sumit3(range(7))
    21
    >>>

The temporary vector is created using the vector constructor taking an
``std::initializer_list``, which is more flexible than constructing a
temporary vector and filling it: it allows the data in the container to be
implicitly converted (e.g. from ``int`` to ``double`` type, or from
pointer to derived to pointer to base class).
As a consequence, however, with STL containers being allowed where Python
containers are, this in turn means that you can pass e.g. an
``std::vector<int>`` (or ``std::list<int>``) where a ``std::vector<double>``
is expected and a temporary is allowed:

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... double sumit4(const std::vector<double>& data) {
    ...   return std::accumulate(data.begin(), data.end(), 0);
    ... }""")
    ...
    True
    >>> cppyy.gbl.sumit4(vector[int](range(7)))
    21.0
    >>>

Normal overload resolution rules continue to apply, however, thus if an
overload were available that takes an ``const std::vector<int>&``, it would
be preferred.

When templates are involved, overload resolution is stricter, to ensure that
a better matching instantiation is preferred over an implicit conversion.
However, that does mean that as-is, C++ is actually more flexible: it has the
curly braces initializer syntax to explicitly infer an
``std::initializer_list``, with no such equivalent in Python.

Although in general this approach guarantees the intended result, it does put
some strictures on the Python side, requiring careful use of types.
However, an easily fixable error is preferable over an implicitly wrong
result.
Note the type of the init argument in the call resulting in an (attempted)
implicit instantiation in the following example:

  .. code-block:: python

    >>> cppyy.cppdef("""
    ... template<class T>
    ... T sumit_T(const std::vector<T>& data, T init) {
    ...  return std::accumulate(data.begin(), data.end(), init);
    ... }""")
    ...
    True
    >>> cppyy.gbl.sumit_T(vector['double'](range(7)), 0)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: Template method resolution failed:
      Failed to instantiate "sumit_T(std::vector<double>&,int)"
      Failed to instantiate "sumit_T(std::vector<double>*,int)"
      Failed to instantiate "sumit_T(std::vector<double>,int)"
    >>> cppyy.gbl.sumit_T(vector['double'](range(7)), 0.)
    21.0
    >>>

To be sure, the code is `too` strict in the simplistic example above, and
with a future version of Cling it should be possible to lift some of these
restrictions without causing incorrect results.

.. rubric:: Footnotes

.. [#f1] The meaning of "temporary" differs between Python and C++: in a statement such as ``func(std.vector[int]((1, 2, 3)))``, there is no temporary as far as Python is concerned, even as there clearly is in the case of a similar statement in C++. Thus that call will succeed even if ``func`` takes a non-const reference.
