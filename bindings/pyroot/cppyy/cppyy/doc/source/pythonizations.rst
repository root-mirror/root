.. _pythonizations:

Pythonizations
==============

Automatic bindings generation mostly gets the job done, but unless a C++
library was designed with expressiveness and interactivity in mind, using it
will feel stilted.
Thus, if you are not the end-user of a set of bindings, it is beneficial to
implement *pythonizations*.
Some of these are already provided by default, e.g. for STL containers.
Consider the following code, iterating over an STL map, using naked bindings
(i.e. "the C++ way"):

.. code-block:: python

   >>> from cppyy.gbl import std
   >>> m = std.map[int, int]()
   >>> for i in range(10):
   ...     m[i] = i*2
   ...
   >>> b = m.begin()
   >>> while b != m.end():
   ...     print(b.__deref__().second, end=' ')
   ...     b.__preinc__()
   ...
   0 2 4 6 8 10 12 14 16 18 
   >>>   

Yes, that is perfectly functional, but it is also very clunky.
Contrast this to the (automatic) pythonization:

.. code-block:: python

   >>> for key, value in m:
   ...    print(value, end=' ')
   ...
   0 2 4 6 8 10 12 14 16 18
   >>>

Such a pythonization can be written completely in Python using the bound C++
methods, with no intermediate language necessary.
Since it is written on abstract features, there is also only one such
pythonization that works for all STL map instantiations.


Python callbacks
----------------

Since bound C++ entities are fully functional Python ones, pythonization can
be done explicitly in an end-user facing Python module.
However, that would prevent lazy installation of pythonizations, so instead a
callback mechanism is provided.

A callback is a function or callable object taking two arguments: the Python
proxy class to be pythonized and its C++ name.
The latter is provided to allow easy filtering.
This callback is then installed through ``cppyy.py.add_pythonization`` and
ideally only for the relevant namespace (installing callbacks for classes in
the global namespace is supported, but beware of name clashes).

Pythonization is most effective of well-structured C++ libraries that have
idiomatic behaviors.
It is then straightforward to use Python reflection to write rules.
For example, consider this callback that looks for the conventional C++
function ``GetLength`` and replaces it with Python's ``__len__``:

.. code-block:: python

    import cppyy

    def replace_getlength(klass, name):
        try:
            klass.__len__ = klass.__dict__['GetLength']
        except KeyError:
            pass

    cppyy.py.add_pythonization(replace_getlength, 'MyNamespace')

    cppyy.cppdef("""
    namespace MyNamespace {
    class MyClass {
    public:
        MyClass(int i) : fInt(i) {}
        int GetLength() { return fInt; }

    private:
        int fInt;
    };
    }""")

    m = cppyy.gbl.MyNamespace.MyClass(42)
    assert len(m) == 42


C++ callbacks
-------------

If you are familiar with the Python C-API, it may sometimes be beneficial to
add unique optimizations to your C++ classes to be picked up by the
pythonization layer.
There are two conventional function that cppyy will look for (no registration
of callbacks needed):

.. code-block:: C++

    static void __cppyy_explicit_pythonize__(PyObject* klass, const std::string&);

which is called *only* for the class that declares it.
And:

.. code-block:: C++

    static void __cppyy_pythonize__(PyObject* klass, const std::string&);

which is also called for all derived classes.

Just as with the Python callbacks, the first argument will be the Python
class proxy, the second the C++ name, for easy filtering.
When called, cppyy will be completely finished with the class proxy, so any
and all changes, including such low-level ones such as the replacement of
iteration or buffer protocols, are fair game.
