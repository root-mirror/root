.. _debugging:
   
Debugging
=========

By default, the ``clang`` JIT as used by cppyy does not generate debugging
information.
This is first of all because it has proven to be not reliable in all cases,
but also because in a production setting this information, being internal to
the wrapper generation, goes unused.
However, that does mean that a debugger that starts from python will not be
able to step through JITed code into the C++ function that needs debugging,
even when such information is available for that C++ function.

To enable debugging information in JITed code, set the ``EXTRA_CLING_ARGS``
envar to ``-g`` (and any further compiler options you need, e.g. add ``-O2``
to debug optimized code).

On a crash in C++, the backend will attempt to provide a stack trace.
This works quite well on Linux (through ``gdb``) and decently on MacOS
(through ``unwind``), but is currently unreliable on MS Windows.
To prevent printing of this trace, which can be slow to produce, set the
envar ``CPPYY_CRASH_QUIET`` to '1'.

It is even more useful to obtain a traceback through the Python code that led
up to the problem in C++.
Many modern debuggers allow mixed-mode C++/Python debugging (for example
`gdb`_ and `MSVC`_), but cppyy can also turn abortive C++ signals (such as a
segmentation violation) into Python exceptions, yielding a normal traceback.
This is particularly useful when working with cross-inheritance and other
cross-language callbacks.

To enable the signals to exceptions conversion, import the lowlevel module
``cppyy.ll`` and use:

  .. code-block:: python

    import cppyy.ll
    cppyy.ll.set_signals_as_exception(True)

Call ``set_signals_as_exception(False)`` to disable the conversion again.
It is recommended to only have the conversion enabled around the problematic
code, as it comes with a performance penalty.
If the problem can be localized to a specific function, you can use its
``__sig2exc__`` flag to only have the conversion active in that function.
Finally, for convenient scoping, you can also use:

  .. code-block:: python

    with cppyy.ll.signals_as_exception():
        # crashing code goes here

The translation of signals to exceptions is as follows (all of the exceptions
are subclasses of ``cppyy.ll.FatalError``):

========================================  ========================================
C++ signal                                Python exception
========================================  ========================================
``SIGSEGV``                               ``cppyy.ll.SegmentationViolation``
``SIGBUS``                                ``cppyy.ll.BusError``
``SIGABRT``                               ``cppyy.ll.AbortSignal``
``SIGILL``                                ``cppyy.ll.IllegalInstruction``
========================================  ========================================

As an example, consider the following cross-inheritance code that crashes
with a segmentation violation in C++, because a ``nullptr`` is dereferenced:

  .. code-block:: python

    import cppyy
    import cppyy.ll

    cppyy.cppdef("""
       class Base {
       public:
          virtual ~Base() {}
          virtual int runit() = 0;
       };

       int callback(Base* b) {
           return b->runit();
       }

       void segfault(int* i) { *i = 42; }
    """)

    class Derived(cppyy.gbl.Base):
        def runit(self):
            print("Hi, from Python!")
            cppyy.gbl.segfault(cppyy.nullptr)

If now used with ``signals_as_exception``, e.g. like so:

  .. code-block:: python

    d = Derived()
    with cppyy.ll.signals_as_exception():
        cppyy.gbl.callback(d)

it produces the following, very informative, Python-side trace::

    Traceback (most recent call last):
      File "crashit.py", line 25, in <module>
        cppyy.gbl.callback(d)
    cppyy.ll.SegmentationViolation: int ::callback(Base* b) =>
        SegmentationViolation: void ::segfault(int* i) =>
        SegmentationViolation: segfault in C++; program state was reset

whereas without, there would be no Python-side information at all.


.. _`gdb`: https://wiki.python.org/moin/DebuggingWithGdb
.. _`MSVC`: https://docs.microsoft.com/en-us/visualstudio/python/debugging-mixed-mode-c-cpp-python-in-visual-studio
