"""Dynamic C++ bindings generator.

This module provides dynamic bindings to C++ through Cling, the LLVM-based C++
interpreter, allowing interactive mixing of Python and C++. Example:

    >>> import cppyy
    >>> cppyy.cppdef(\"\"\"
    ... class MyClass {
    ... public:
    ...     MyClass(int i) : m_data(i) {}
    ...     int m_data;
    ... };\"\"\")
    True
    >>> from cppyy.gbl import MyClass
    >>> m = MyClass(42)
    >>> cppyy.cppdef(\"\"\"
    ... void say_hello(MyClass* m) {
    ...     std::cout << "Hello, the number is: " << m->m_data << std::endl;
    ... }\"\"\")
    True
    >>> MyClass.say_hello = cppyy.gbl.say_hello
    >>> m.say_hello()
    Hello, the number is: 42
    >>> m.m_data = 13
    >>> m.say_hello()
    Hello, the number is: 13
    >>>

For full documentation, see:
   https://cppyy.readthedocs.io/

"""

__author__ = 'Wim Lavrijsen <WLavrijsen@lbl.gov>'

__all__ = [
    'cppdef',                 # declare C++ source to Cling
    'include',                # load and jit a header file
    'c_include',              # load and jit a C header file
    'load_library',           # load a shared library
    'nullptr',                # unique pointer representing NULL
    'sizeof',                 # size of a C++ type
    'typeid',                 # typeid of a C++ type
    'add_include_path',       # add a path to search for headers
    'add_library_path',       # add a path to search for headers
    'add_autoload_map',       # explicitly include an autoload map
    ]

from ._version import __version__

import os, sys, sysconfig, warnings
import importlib

# import libcppyy with Python version number
major, minor = sys.version_info[0:2]
py_version_str = '{}_{}'.format(major, minor)
libcppyy_mod_name = 'libcppyy' + py_version_str

try:
    importlib.import_module(libcppyy_mod_name)
except ImportError:
    raise ImportError(
            'Failed to import {}. Please check that ROOT has been built for Python {}.{}'.format(
                libcppyy_mod_name, major, minor))

# ensure 'import libcppyy' will find the versioned module
sys.modules['libcppyy'] = sys.modules[libcppyy_mod_name]

# tell cppyy that libcppyy_backend is versioned
def _check_py_version(lib_name, cbl_var):
    import re
    if re.match('^libcppyy_backend\d_\d$', lib_name):
       # library name already has version
       if lib_name.endswith(py_version_str):
           return ''
       else:
           raise RuntimeError('CPPYY_BACKEND_LIBRARY variable ({})'
                              ' does not match Python version {}.{}'
                              .format(cbl_var, major, minor))
    else:
       return py_version_str

if 'CPPYY_BACKEND_LIBRARY' in os.environ:
    clean_cbl = False
    cbl_var = os.environ['CPPYY_BACKEND_LIBRARY']
    start = 0
    last_sep = cbl_var.rfind(os.path.sep)
    if last_sep >= 0:
        start = last_sep + 1
    first_dot = cbl_var.find('.', start)
    if first_dot >= 0:
        # lib_name = [/path/to/]libcppyy_backend[py_version_str]
        # suffix = so | ...
        lib_name = cbl_var[:first_dot]
        suff = cbl_var[first_dot+1:]
        ver = _check_py_version(lib_name[start:], cbl_var)
        os.environ['CPPYY_BACKEND_LIBRARY'] = '.'.join([
            lib_name + ver, suff])
    else:
        ver = _check_py_version(cbl_var[start:], cbl_var)
        os.environ['CPPYY_BACKEND_LIBRARY'] += ver
else:
    clean_cbl = True
    if not any(var in os.environ for var in ('LD_LIBRARY_PATH','DYLD_LIBRARY_PATH')):
        # macOS SIP can prevent DYLD_LIBRARY_PATH from having any effect.
        # Set cppyy env variable here to make sure libraries are found.
        _lib_dir = os.path.dirname(os.path.dirname(__file__))
        _lcb_path = os.path.join(_lib_dir, 'libcppyy_backend')
        os.environ['CPPYY_BACKEND_LIBRARY'] = _lcb_path + py_version_str
    else:
        os.environ['CPPYY_BACKEND_LIBRARY'] = 'libcppyy_backend' + py_version_str

if not 'CLING_STANDARD_PCH' in os.environ:
    local_pch = os.path.join(os.path.dirname(__file__), 'allDict.cxx.pch')
    if os.path.exists(local_pch):
        os.putenv('CLING_STANDARD_PCH', local_pch)
        os.environ['CLING_STANDARD_PCH'] = local_pch

try:
    import __pypy__
    del __pypy__
    ispypy = True
except ImportError:
    ispypy = False

# import separately instead of in the above try/except block for easier to
# understand tracebacks
if ispypy:
    from ._pypy_cppyy import *
else:
    from ._cpython_cppyy import *
if clean_cbl:
    del os.environ['CPPYY_BACKEND_LIBRARY'] # we set it, so clean it


#- allow importing from gbl --------------------------------------------------
sys.modules['cppyy.gbl'] = gbl
sys.modules['cppyy.gbl.std'] = gbl.std


#- external typemap ----------------------------------------------------------
from . import _typemap
_typemap.initialize(_backend)


#- pythonization factories ---------------------------------------------------
from . import _pythonization as py
py._set_backend(_backend)

def _standard_pythonizations(pyclass, name):
  # pythonization of tuple; TODO: placed here for convenience, but verify that decision
    if name.find('std::tuple<', 0, 11) == 0 or name.find('tuple<', 0, 6) == 0:
        import cppyy
        pyclass._tuple_len = cppyy.gbl.std.tuple_size(pyclass).value
        def tuple_len(self):
            return self.__class__._tuple_len
        pyclass.__len__ = tuple_len
        def tuple_getitem(self, idx, get=cppyy.gbl.std.get):
            if idx < self.__class__._tuple_len:
                return get[idx](self)
            raise IndexError(idx)
        pyclass.__getitem__ = tuple_getitem

if not ispypy:
    py.add_pythonization(_standard_pythonizations)   # should live on std only
# TODO: PyPy still has the old-style pythonizations, which require the full
# class name (not possible for std::tuple ...)

# std::make_shared creates needless templates: rely on Python's introspection
# instead. This also allows Python derived classes to be handled correctly.
class py_make_shared(object):
    def __init__(self, cls):
        self.cls = cls
    def __call__(self, *args):
        if len(args) == 1 and type(args[0]) == self.cls:
            obj = args[0]
        else:
            obj = self.cls(*args)
        obj.__python_owns__ = False     # C++ to take ownership
        return gbl.std.shared_ptr[self.cls](obj)

class make_shared(object):
    def __getitem__(self, cls):
        return py_make_shared(cls)

gbl.std.make_shared = make_shared()
del make_shared


#--- CFFI style interface ----------------------------------------------------
def cppdef(src):
    """Declare C++ source <src> to Cling."""
    if not gbl.gInterpreter.Declare(src):
        return False
    return True

def load_library(name):
    """Explicitly load a shared library."""
    gSystem = gbl.gSystem
    if name[:3] != 'lib':
        if not gSystem.FindDynamicLibrary(gbl.TString(name), True) and\
               gSystem.FindDynamicLibrary(gbl.TString('lib'+name), True):
            name = 'lib'+name
    sc = gSystem.Load(name)
    if sc == -1:
        raise RuntimeError("Unable to load library "+name)

def include(header):
    """Load (and JIT) header file <header> into Cling."""
    gbl.gInterpreter.Declare('#include "%s"' % header)

def c_include(header):
    """Load (and JIT) header file <header> into Cling."""
    gbl.gInterpreter.Declare("""extern "C" {
#include "%s"
}""" % header)

def add_include_path(path):
    """Add a path to the include paths available to Cling."""
    if not os.path.isdir(path):
        raise OSError("no such directory: %s" % path)
    gbl.gInterpreter.AddIncludePath(path)

def add_library_path(path):
    """Add a path to the library search paths available to Cling."""
    if not os.path.isdir(path):
        raise OSError("no such directory: %s" % path)
    gbl.gSystem.AddDynamicPath(path)

# add access to Python C-API headers
apipath = sysconfig.get_path('include', 'posix_prefix' if os.name == 'posix' else os.name)
if os.path.exists(apipath):
    add_include_path(apipath)
elif ispypy:
  # possibly structured without 'pythonx.y' in path
    apipath = os.path.dirname(apipath)
    if os.path.exists(apipath) and os.path.exists(os.path.join(apipath, 'Python.h')):
        add_include_path(apipath)

# add access to extra headers for dispatcher (CPyCppyy only (?))
if not ispypy:
    if 'CPPYY_API_PATH' in os.environ:
        apipath_extra = os.environ['CPPYY_API_PATH']
    else:
        apipath_extra = os.path.join(os.path.dirname(apipath), 'site', 'python'+sys.version[:3])
        if not os.path.exists(os.path.join(apipath_extra, 'CPyCppyy')):
            import glob, libcppyy
            apipath_extra = os.path.dirname(libcppyy.__file__)
          # a "normal" structure finds the include directory 3 levels up,
          # ie. from lib/pythonx.y/site-packages
            for i in range(3):
                if not os.path.exists(os.path.join(apipath_extra, 'include')):
                    apipath_extra = os.path.dirname(apipath_extra)

            apipath_extra = os.path.join(apipath_extra, 'include')
          # add back pythonx.y or site/pythonx.y if available
            for p in glob.glob(os.path.join(apipath_extra, 'python'+sys.version[:3]+'*'))+\
                     glob.glob(os.path.join(apipath_extra, '*', 'python'+sys.version[:3]+'*')):
                if os.path.exists(os.path.join(p, 'CPyCppyy')):
                    apipath_extra = p
                    break

    cpycppyy_path = os.path.join(apipath_extra, 'CPyCppyy')
    if apipath_extra.lower() != 'none':
        if not os.path.exists(cpycppyy_path):
            warnings.warn("CPyCppyy API path not found (tried: %s); set CPPYY_API_PATH to fix" % os.path.dirname(cpycppyy_path))
        else:
            add_include_path(apipath_extra)

if os.getenv('CONDA_PREFIX'):
  # MacOS, Linux
    include_path = os.path.join(os.getenv('CONDA_PREFIX'), 'include')
    if os.path.exists(include_path): add_include_path(include_path)

  # Windows
    include_path = os.path.join(os.getenv('CONDA_PREFIX'), 'Library', 'include')
    if os.path.exists(include_path): add_include_path(include_path)

# assuming that we are in PREFIX/lib/python/site-packages/cppyy, add PREFIX/include to the search path
include_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, os.path.pardir, os.path.pardir, 'include'))
if os.path.exists(include_path): add_include_path(include_path)

del include_path, apipath, ispypy

def add_autoload_map(fname):
    """Add the entries from a autoload (.rootmap) file to Cling."""
    if not os.path.isfile(fname):
        raise OSError("no such file: %s" % fname)
    gbl.gInterpreter.LoadLibraryMap(fname)

def _get_name(tt):
    if type(tt) == str:
        return tt
    try:
        ttname = tt.__cpp_name__
    except AttributeError:
        ttname = tt.__name__
    return ttname

_sizes = {}
def sizeof(tt):
    """Returns the storage size (in chars) of C++ type <tt>."""
    if not isinstance(tt, type) and not type(tt) == str:
        tt = type(tt)
    try:
        return _sizes[tt]
    except KeyError:
        sz = gbl.gInterpreter.ProcessLine("sizeof(%s);" % (_get_name(tt),))
        _sizes[tt] = sz
        return sz

_typeids = {}
def typeid(tt):
    """Returns the C++ runtime type information for type <tt>."""
    if not isinstance(tt, type):
        tt = type(tt)
    try:
        return _typeids[tt]
    except KeyError:
        tidname = 'typeid_'+str(len(_typeids))
        gbl.gInterpreter.ProcessLine(
            "namespace _cppyy_internal { auto* %s = &typeid(%s); }" %\
            (tidname, _get_name(tt),))
        tid = getattr(gbl._cppyy_internal, tidname)
        _typeids[tt] = tid
        return tid
