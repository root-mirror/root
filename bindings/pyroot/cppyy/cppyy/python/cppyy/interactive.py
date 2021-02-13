import sys
try:
    import __pypy__
    del __pypy__
    ispypy = True
except ImportError:
    ispypy = False


#- fake namespace for interactive lazy lookups -------------------------------
class InteractiveLazy(object):
    def __init__(self, hook_okay):
        self._hook_okay = hook_okay

    def __getattr__(self, attr):
        import cppyy, sys
        if attr == '__all__':
          # copy all exported items from cppyy itself
            for v in cppyy.__all__:
                self.__dict__[v] = getattr(cppyy, v)

          # add the lookup hook into cppyy.gbl if legal, or put it under 'g'
          # if not (PyPy and IPython for now)
            if self._hook_okay:
                caller = sys.modules[sys._getframe(1).f_globals['__name__']]
                cppyy._backend._set_cpp_lazy_lookup(caller.__dict__)
                return cppyy.__all__
            else:
                self.__dict__['g']   = cppyy.gbl
                self.__dict__['std'] = cppyy.gbl.std
                return ['g', 'std']+cppyy.__all__
        return getattr(cppyy, attr)

sys.modules['cppyy.interactive'] = InteractiveLazy(\
    not ispypy and not (hasattr(__builtins__, '__IPYTHON__') or 'IPython' in sys.modules))
del InteractiveLazy, ispypy
