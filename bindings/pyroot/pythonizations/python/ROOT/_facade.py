import types
import sys
import os
from functools import partial

import libcppyy as cppyy_backend
from cppyy import gbl as gbl_namespace
from cppyy import cppdef, include
from libROOTPythonizations import gROOT, CreateBufferFromAddress

from ._application import PyROOTApplication
_numba_pyversion = (2, 7, 5)
if sys.version_info[:3] > _numba_pyversion:
    # Python <= 2.7.5 cannot use exec in an inner function
    from ._numbadeclare import _NumbaDeclareDecorator


class PyROOTConfiguration(object):
    """Class for configuring PyROOT"""

    def __init__(self):
        self.IgnoreCommandLineOptions = True
        self.ShutDown = True
        self.DisableRootLogon = False
        self.StartGUIThread = True


class _gROOTWrapper(object):
    """Internal class to manage lookups of gROOT in the facade.
       This wrapper calls _finalSetup on the facade when it
       receives a lookup, unless the lookup is for SetBatch.
       This allows to evaluate the command line parameters
       before checking if batch mode is on in _finalSetup
    """

    def __init__(self, facade):
        self.__dict__['_facade'] = facade
        self.__dict__['_gROOT'] = gROOT

    def __getattr__( self, name ):
        if name != 'SetBatch' and self._facade.__dict__['gROOT'] != self._gROOT:
            self._facade._finalSetup()
        return getattr(self._gROOT, name)

    def __setattr__(self, name, value):
        return setattr(self._gROOT, name, value)


class ROOTFacade(types.ModuleType):
    """Facade class for ROOT module"""

    def __init__(self, module, is_ipython):
        types.ModuleType.__init__(self, module.__name__)

        self.module = module
        # Importing all will be customised later
        self.module.__all__ = []

        self.__doc__  = module.__doc__
        self.__name__ = module.__name__
        self.__file__ = module.__file__

        # Inject gROOT global
        self.gROOT = _gROOTWrapper(self)

        # Expose some functionality from CPyCppyy extension module
        self._cppyy_exports = [ 'nullptr', 'bind_object', 'as_cobject', 'addressof',
                                'SetMemoryPolicy', 'kMemoryHeuristics', 'kMemoryStrict',
                                'SetOwnership' ]
        for name in self._cppyy_exports:
            setattr(self, name, getattr(cppyy_backend, name))
        # For backwards compatibility
        self.MakeNullPointer = partial(self.bind_object, 0)
        self.BindObject = self.bind_object
        self.AsCObject = self.as_cobject

        # Initialize configuration
        self.PyConfig = PyROOTConfiguration()

        self._is_ipython = is_ipython

        # Redirect lookups to temporary helper methods
        # This lets the user do some actions before all the machinery is in place:
        # - Set batch mode in gROOT
        # - Set options in PyConfig
        self.__class__.__getattr__ = self._getattr
        self.__class__.__setattr__ = self._setattr

        # Setup import hook
        self._set_import_hook()

    def AddressOf(self, obj):
        # Return an indexable buffer of length 1, whose only element
        # is the address of the object.
        # The address of the buffer is the same as the address of the
        # address of the object

        # addr is the address of the address of the object
        addr = self.addressof(instance = obj, byref = True)

        # Create a buffer (LowLevelView) from address
        return CreateBufferFromAddress(addr)

    def _set_import_hook(self):
        # This hook allows to write e.g:
        # from ROOT.A import a
        # instead of the longer:
        # from ROOT import A
        # from A import a
        try:
            import __builtin__
        except ImportError:
            import builtins as __builtin__  # name change in p3
        _orig_ihook = __builtin__.__import__
        def _importhook(name, *args, **kwds):
            if name[0:5] == 'ROOT.':
                try:
                    sys.modules[name] = getattr(self, name[5:])
                except Exception:
                    pass
            return _orig_ihook(name, *args, **kwds)
        __builtin__.__import__ = _importhook

    def _handle_import_all(self):
        # Called if "from ROOT import *" is executed in the app.
        # Customises lookup in Python's main module to also
        # check in C++'s global namespace

        # Get caller module (jump over the facade frames)
        num_frame = 2
        frame = sys._getframe(num_frame).f_globals['__name__']
        while frame == 'ROOT._facade':
            num_frame += 1
            frame = sys._getframe(num_frame).f_globals['__name__']
        caller = sys.modules[frame]

        # Install the hook
        cppyy_backend._set_cpp_lazy_lookup(caller.__dict__)

    def _fallback_getattr(self, name):
        # Try:
        # - in the global namespace
        # - in the ROOT namespace
        # - in gROOT (ROOT lists such as list of files,
        #   memory mapped files, functions, geometries ecc.)
        # The first two attempts allow to lookup
        # e.g. ROOT.ROOT.Math as ROOT.Math

        if name == '__all__':
            self._handle_import_all()
            # Make the attributes of the facade be injected in the
            # caller module
            raise AttributeError()
        # Note that hasattr caches the lookup for getattr
        elif hasattr(gbl_namespace, name):
            return getattr(gbl_namespace, name)
        elif hasattr(gbl_namespace.ROOT, name):
            return getattr(gbl_namespace.ROOT, name)
        else:
            res = gROOT.FindObject(name)
            if res:
                return res
        raise AttributeError("Failed to get attribute {} from ROOT".format(name))

    def _finalSetup(self):
        # Prevent this method from being re-entered through the gROOT wrapper
        self.__dict__['gROOT'] = gROOT

        # Setup interactive usage from Python
        self.__dict__['app'] = PyROOTApplication(self.PyConfig, self._is_ipython)
        if not self.gROOT.IsBatch() and self.PyConfig.StartGUIThread:
            self.app.init_graphics()

        # Set memory policy to kUseHeuristics.
        # This restores the default in PyROOT which was changed
        # by new Cppyy
        self.SetMemoryPolicy(self.kMemoryHeuristics)

        # Redirect lookups to cppyy's global namespace
        self.__class__.__getattr__ = self._fallback_getattr
        self.__class__.__setattr__ = lambda self, name, val: setattr(gbl_namespace, name, val)

        # Run rootlogon if exists
        self._run_rootlogon()

    def _getattr(self, name):
        # Special case, to allow "from ROOT import gROOT" w/o starting the graphics
        if name == '__path__':
            raise AttributeError(name)

        self._finalSetup()

        return getattr(self, name)

    def _setattr(self, name, val):
        self._finalSetup()

        return setattr(self, name, val)

    def _run_rootlogon(self):
        # Run custom logon file (must be after creation of ROOT globals)
        hasargv = hasattr(sys, 'argv')
        # -n disables the reading of the logon file, just like with root
        if hasargv and not '-n' in sys.argv and not self.PyConfig.DisableRootLogon:
            file_path = os.path.expanduser('~/.rootlogon.py')
            if os.path.exists(file_path):
                # Could also have used execfile, but import is likely to give fewer surprises
                module_name = 'rootlogon'
                if sys.version_info >= (3,5):
                    import importlib.util
                    spec = importlib.util.spec_from_file_location(module_name, file_path)
                    module = importlib.util.module_from_spec(spec)
                    sys.modules[module_name] = module
                    spec.loader.exec_module(module)
                else:
                    import imp
                    imp.load_module(module_name, open(file_path, 'r'), file_path, ('.py','r',1))
                    del imp
            else:
                # If the .py version of rootlogon exists, the .C is ignored (the user can
                # load the .C from the .py, if so desired).
                # System logon, user logon, and local logon (skip Rint.Logon)
                name = '.rootlogon.C'
                logons = [
                    os.path.join(str(self.TROOT.GetEtcDir()), 'system' + name),
                    os.path.expanduser(os.path.join('~', name))
                    ]
                if logons[-1] != os.path.join(os.getcwd(), name):
                    logons.append(name)
                for rootlogon in logons:
                    if os.path.exists(rootlogon):
                        self.TApplication.ExecuteFile(rootlogon)

    # Inject version as __version__ property in ROOT module
    @property
    def __version__(self):
        return self.gROOT.GetVersion()

    # Overload VecOps namespace
    # The property gets the C++ namespace, adds the pythonizations and
    # eventually deletes itself so that following calls go directly
    # to the C++ namespace. This mechanic ensures that we pythonize the
    # namespace lazily.
    @property
    def VecOps(self):
        ns = self._fallback_getattr('VecOps')
        try:
            from libROOTPythonizations import AsRVec
            ns.AsRVec = AsRVec
        except:
            raise Exception('Failed to pythonize the namespace VecOps')
        del type(self).VecOps
        return ns

    # Overload RDF namespace
    @property
    def RDF(self):
        ns = self._fallback_getattr('RDF')
        try:
            from libROOTPythonizations import MakeNumpyDataFrame
            ns.MakeNumpyDataFrame = MakeNumpyDataFrame
        except:
            raise Exception('Failed to pythonize the namespace RDF')
        del type(self).RDF
        return ns

    # Overload TMVA namespace
    @property
    def TMVA(self):
        ns = self._fallback_getattr('TMVA')
        try:
            from libROOTPythonizations import AsRTensor
            ns.Experimental.AsRTensor = AsRTensor
        except:
            raise Exception('Failed to pythonize the namespace TMVA')
        del type(self).TMVA
        return ns

    # Create and overload Numba namespace
    @property
    def Numba(self):
        if sys.version_info[:3] <= _numba_pyversion:
            raise Exception('ROOT.Numba requires Python above version {}.{}.{}'.format(*_numba_pyversion))
        cppdef('namespace Numba {}')
        ns = self._fallback_getattr('Numba')
        ns.Declare = staticmethod(_NumbaDeclareDecorator)
        del type(self).Numba
        return ns

    # Get TPyDispatcher for programming GUI callbacks
    @property
    def TPyDispatcher(self):
        include('ROOT/TPyDispatcher.h')
        tpd = gbl_namespace.TPyDispatcher
        type(self).TPyDispatcher = tpd
        return tpd
