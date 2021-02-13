""" Low-level utilities, to be used for "emergencies only".
"""

import cppyy

try:
    import __pypy__
    del __pypy__
    ispypy = True
except ImportError:
    ispypy = False

__all__ = [
    'cast',
    'static_cast',
    'reinterpret_cast',
    'dynamic_cast',
    'malloc',
    'free',
    'array_new',
    'array_detele',
    'signals_as_exception',
    'set_signals_as_exception'
    'FatalError'
    'BusError',
    'SegmentationViolation',
    'IllegalInstruction',
    'AbortSignal',
    ]


# import low-level python converters
for _name in ['addressof', 'as_cobject', 'as_capsule', 'as_ctypes']:
    try:
        exec('%s = cppyy._backend.%s' % (_name, _name))
        __all__.append(_name)
    except AttributeError:
        pass
del _name


# create low-level helpers
cppyy.cppdef("""namespace __cppyy_internal {
// type casting
    template<typename T, typename U>
    T cppyy_cast(U val) { return (T)val; }

    template<typename T, typename U>
    T cppyy_static_cast(U val) { return static_cast<T>(val); }

    template<typename T, typename U>
    T cppyy_reinterpret_cast(U val) { return reinterpret_cast<T>(val); }

    template<typename T, typename S>
    T* cppyy_dynamic_cast(S* obj) { return dynamic_cast<T*>(obj); }

// memory allocation/free-ing
    template<typename T>
    T* cppyy_malloc(size_t count=1) { return (T*)malloc(sizeof(T*)*count); }

    template<typename T>
    T* cppyy_array_new(size_t count) { return new T[count]; }

    template<typename T>
    void cppyy_array_delete(T* ptr) { delete[] ptr; }
}""")


# helper for sizing arrays
class ArraySizer(object):
    def __init__(self, func):
        self.func = func
    def __getitem__(self, t):
        self.array_type = t
        return self
    def __call__(self, size):
        res = self.func[self.array_type](size)
        res.reshape((size,))
        return res

# import casting helpers
cast             = cppyy.gbl.__cppyy_internal.cppyy_cast
static_cast      = cppyy.gbl.__cppyy_internal.cppyy_static_cast
reinterpret_cast = cppyy.gbl.__cppyy_internal.cppyy_reinterpret_cast
dynamic_cast     = cppyy.gbl.__cppyy_internal.cppyy_dynamic_cast

# import memory allocation/free-ing helpers
malloc           = ArraySizer(cppyy.gbl.__cppyy_internal.cppyy_malloc)
free             = cppyy.gbl.free      # for symmetry
array_new        = ArraySizer(cppyy.gbl.__cppyy_internal.cppyy_array_new)
array_delete     = cppyy.gbl.__cppyy_internal.cppyy_array_delete

# signals as exceptions
if not ispypy:
    FatalError            = cppyy._backend.FatalError
    BusError              = cppyy._backend.BusError
    SegmentationViolation = cppyy._backend.SegmentationViolation
    IllegalInstruction    = cppyy._backend.IllegalInstruction
    AbortSignal           = cppyy._backend.AbortSignal

    class signals_as_exception:
        def __enter__(self):
            cppyy._backend.SetGlobalSignalPolicy(1)

        def __exit__(self, type, value, traceback):
            cppyy._backend.SetGlobalSignalPolicy(0)

    set_signals_as_exception = cppyy._backend.SetGlobalSignalPolicy

else:
    class FatalError(Exception):
        pass
    class BusError(FatalError):
        pass
    class SegmentationViolation(FatalError):
        pass
    class IllegalInstruction(FatalError):
        pass
    class AbortSignal(FatalError):
        pass

    class signals_as_exception:
        def __enter__(self):
            pass   # not yet implemented

        def __exit__(self, type, value, traceback):
            pass   # not yet implemented

    def set_signals_as_exception(seton):
        return False

del ispypy

