# Author: Enric Tejedor, Danilo Piparo CERN  06/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

# macOS SIP can prevent DYLD_LIBRARY_PATH from having any effect.
# Set cppyy env variable here to make sure libraries are found.
from os import environ, path
if not any(var in environ for var in ('LD_LIBRARY_PATH','DYLD_LIBRARY_PATH')):
    _lib_dir = path.dirname(path.dirname(__file__))
    _lcb_path = path.join(_lib_dir, 'libcppyy_backend')
    environ['CPPYY_BACKEND_LIBRARY'] = _lcb_path

# Inform cppyy that we have already generated the PCH
_etc_dir = path.join(path.dirname(path.dirname(path.dirname(__file__))), 'etc')
environ['CLING_STANDARD_PCH'] = path.join(_etc_dir, 'allDict.cxx.pch')

import cppyy
import ROOT.pythonization as pyz

import functools
import importlib
import pkgutil

def pythonization(lazy = True):
    """
    Pythonizor decorator to be used in pythonization modules for pythonizations.
    These pythonizations functions are invoked upon usage of the class.
    Parameters
    ----------
    lazy : boolean
        If lazy is true, the class is pythonized upon first usage, otherwise
        upon import of the ROOT module.
    """
    def pythonization_impl(fn):
        """
        The real decorator. This structure is adopted to deal with parameters
        fn : function
            Function that implements some pythonization.
            The function must accept two parameters: the class
            to be pythonized and the name of that class.
        """
        if lazy:
            cppyy.py.add_pythonization(fn)
        else:
            fn()
    return pythonization_impl

# Trigger the addition of the pythonizations
for _, module_name, _ in  pkgutil.walk_packages(pyz.__path__):
    module = importlib.import_module(pyz.__name__ + '.' + module_name)

# Configure ROOT facade module
import sys
from ._facade import ROOTFacade
sys.modules[__name__] = ROOTFacade(sys.modules[__name__])
