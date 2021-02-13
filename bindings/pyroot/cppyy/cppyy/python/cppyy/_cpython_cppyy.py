""" CPython-specific touch-ups
"""

from . import _stdcpp_fix
from cppyy_backend import loader

__all__ = [
    'gbl',
    'load_reflection_info',
    'addressof',
    'bind_object',
    'nullptr',
    '_backend',
    ]

# first load the dependency libraries of the backend, then pull in the
# libcppyy extension module
c = loader.load_cpp_backend()
import libcppyy as _backend
_backend._cpp_backend = c

# explicitly expose APIs from libcppyy
import ctypes
_w = ctypes.CDLL(_backend.__file__, ctypes.RTLD_GLOBAL)


# some beautification for inspect (only on p2)
import sys
if sys.hexversion < 0x3000000:
  # TODO: this reliese on CPPOverload cooking up a func_code object, which atm
  # is simply not implemented for p3 :/

  # convince inspect that PyROOT method proxies are possible drop-ins for python
  # methods and classes for pydoc
    import inspect

    inspect._old_isfunction = inspect.isfunction
    def isfunction(object):
        if type(object) == _backend.CPPOverload and not object.im_class:
            return True
        return inspect._old_isfunction( object )
    inspect.isfunction = isfunction

    inspect._old_ismethod = inspect.ismethod
    def ismethod(object):
        if type(object) == _backend.CPPOverload:
            return True
        return inspect._old_ismethod(object)
    inspect.ismethod = ismethod
    del isfunction, ismethod


### template support ---------------------------------------------------------
class Template(object):  # expected/used by ProxyWrappers.cxx in CPyCppyy
    def __init__(self, name):
        self.__name__ = name

    def __repr__(self):
        return "<cppyy.Template '%s' object at %s>" % (self.__name__, hex(id(self)))

    def __call__(self, *args):
        newargs = [self.__name__]
        for arg in args:
            if type(arg) == str:
                arg = ','.join(map(lambda x: x.strip(), arg.split(',')))
            newargs.append(arg)
        pyclass = _backend.MakeCppTemplateClass(*newargs)

      # special case pythonization (builtin_map is not available from the C-API)
        if 'push_back' in pyclass.__dict__ and not '__iadd__' in pyclass.__dict__:
            if 'reserve' in pyclass.__dict__:
                def iadd(self, ll):
                    self.reserve(len(ll))
                    for x in ll: self.push_back(x)
                    return self
            else:
                def iadd(self, ll):
                    for x in ll: self.push_back(x)
                    return self
            pyclass.__iadd__ = iadd

        return pyclass

    def __getitem__(self, *args):
        if args and type(args[0]) == tuple:
            return self.__call__(*(args[0]))
        return self.__call__(*args)

_backend.Template = Template


#- :: and std:: namespaces ---------------------------------------------------
gbl = _backend.CreateScopeProxy('')
gbl.__class__.__repr__ = lambda cls : '<namespace cppyy.gbl at 0x%x>' % id(cls)
gbl.std =  _backend.CreateScopeProxy('std')
# for move, we want our "pythonized" one, not the C++ template
gbl.std.move  = _backend.move


#- add to the dynamic path as needed -----------------------------------------
import os
def add_default_paths():
    gSystem = gbl.gSystem
    if os.getenv('CONDA_PREFIX'):
      # MacOS, Linux
        lib_path = os.path.join(os.getenv('CONDA_PREFIX'), 'lib')
        if os.path.exists(lib_path): gSystem.AddDynamicPath(lib_path)

      # Windows
        lib_path = os.path.join(os.getenv('CONDA_PREFIX'), 'Library', 'lib')
        if os.path.exists(lib_path): gSystem.AddDynamicPath(lib_path)

  # assuming that we are in PREFIX/lib/python/site-packages/cppyy, add PREFIX/lib to the search path
    lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, os.path.pardir))
    if os.path.exists(lib_path): gSystem.AddDynamicPath(lib_path)

    try:
        with open('/etc/ld.so.conf') as ldconf:
            for line in ldconf:
                f = line.strip()
                if (os.path.exists(f)):
                    gSystem.AddDynamicPath(f)
    except IOError:
        pass
add_default_paths()
del add_default_paths


#- exports -------------------------------------------------------------------
addressof     = _backend.addressof
bind_object   = _backend.bind_object
nullptr       = _backend.nullptr

def load_reflection_info(name):
    sc = gbl.gSystem.Load(name)
    if sc == -1:
        raise RuntimeError("Unable to load reflection library "+name)
