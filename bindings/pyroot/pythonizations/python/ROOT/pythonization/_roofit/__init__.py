# Authors:
# * Harshal Shende  03/2021
# * Hinnerk C. Schmidt 02/2021
# * Jonas Rembser 03/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

import sys
import cppyy

from ROOT import pythonization

from ._rooabscollection import RooAbsCollection
from ._rooabsdata import RooAbsData
from ._rooabspdf import RooAbsPdf
from ._rooabsreal import RooAbsReal
from ._rooabsreallvalue import RooAbsRealLValue
from ._rooarglist import RooArgList
from ._rooargset import RooArgSet
from ._roocategory import RooCategory
from ._roochi2var import RooChi2Var
from ._roodatahist import RooDataHist
from ._roodataset import RooDataSet
from ._roodecays import RooDecay, RooBDecay, RooBCPGenDecay, RooBCPEffDecay, RooBMixDecay
from ._roogenfitstudy import RooGenFitStudy
from ._rooglobalfunc import (
    FitOptions,
    Format,
    Frame,
    MultiArg,
    YVar,
    ZVar,
    Slice,
    Import,
    Link,
    LineColor,
    FillColor,
    MarkerColor,
    LineStyle,
    FillStyle,
    MarkerStyle,
)
from ._roomcstudy import RooMCStudy
from ._roomsgservice import RooMsgService
from ._roonllvar import RooNLLVar
from ._rooprodpdf import RooProdPdf
from ._roosimultaneous import RooSimultaneous
from ._roosimwstool import RooSimWSTool
from ._rooworkspace import RooWorkspace


# list of python classes that are used to pythonize RooFit classes
python_classes = [
    RooAbsCollection,
    RooAbsData,
    RooAbsPdf,
    RooAbsReal,
    RooAbsRealLValue,
    RooArgList,
    RooArgSet,
    RooBCPGenDecay,
    RooBCPEffDecay,
    RooBDecay,
    RooBMixDecay,
    RooCategory,
    RooChi2Var,
    RooDataHist,
    RooDataSet,
    RooDecay,
    RooGenFitStudy,
    RooMCStudy,
    RooMsgService,
    RooNLLVar,
    RooProdPdf,
    RooSimultaneous,
    RooSimWSTool,
    RooWorkspace,
]

# list of python functions that are used to pythonize RooGlobalFunc function in RooFit
python_roofit_functions = [
    FitOptions,
    Format,
    Frame,
    MultiArg,
    YVar,
    ZVar,
    Slice,
    Import,
    Link,
    LineColor,
    FillColor,
    MarkerColor,
    LineStyle,
    FillStyle,
    MarkerStyle,
]

# create a dictionary for convenient access to python classes
python_classes_dict = dict()
for python_class in python_classes:
    python_classes_dict[python_class.__name__] = python_class


def get_defined_attributes(klass, consider_base_classes=False):
    """
    Get all class attributes that are defined in a given class or optionally in
    any of its base classes (except for `object`).
    """

    blacklist = ["__dict__", "__doc__", "__hash__", "__module__", "__weakref__"]

    if not consider_base_classes:
        return [attr for attr in klass.__dict__.keys() if attr not in blacklist]

    # get a list of this class and all its base classes, excluding `object`
    method_resolution_order = klass.mro()
    if object in method_resolution_order:
        method_resolution_order.remove(object)

    def is_defined(funcname):

        if funcname in blacklist:
            return False

        in_any_dict = False

        for mro_class in method_resolution_order:
            if funcname in mro_class.__dict__:
                in_any_dict = True

        return in_any_dict

    return [attr for attr in dir(klass) if is_defined(attr)]


def rebind_instancemethod(to_class, from_class, func_name):
    """
    Bind the instance method `from_class.func_name` also to class `to_class`.
    """

    from_method = getattr(from_class, func_name)

    import sys

    if sys.version_info > (3, 0):
        to_method = from_method
    else:
        import new

        to_method = new.instancemethod(from_method.__func__, None, to_class)

    setattr(to_class, func_name, to_method)


def make_func_name_orig(func_name):
    """Return the name that we will give to the original cppyy function."""
    # special treatment of magic functions, e.g.: __getitem__ > _getitem
    if func_name.startswith("__") and func_name.endswith("__"):
        func_name = func_name[2:-2]

    return "_" + func_name


@pythonization()
def pythonize_roofit_class(klass, name):
    # Parameters:
    # klass: class to pythonize
    # name: string containing the name of the class

    if not name.startswith("Roo"):
        return

    if not name in python_classes_dict:
        return

    python_klass = python_classes_dict[name]

    # list of functions to pythonize, which are assumed to be all functions in
    # that are manually defined in the Python classes or their superclasses
    func_names = get_defined_attributes(python_klass)

    for func_name in func_names:

        # if the RooFit class already has a function with the same name as our
        # pythonization, we rename it and prefix it with an underscore
        if hasattr(klass, func_name):
            # new name for original function
            func_name_orig = make_func_name_orig(func_name)

            setattr(klass, func_name_orig, getattr(klass, func_name))

        rebind_instancemethod(klass, python_klass, func_name)

    return


def pythonize_roofit_namespace(ns):

    for python_func in python_roofit_functions:
        func_name = python_func.__name__
        func_name_orig = "_" + func_name

        if sys.version_info <= (3, 0):
            # In Python 2 the RooFit is treated like a class and the global
            # functions in the namespace must be static methods.
            python_func = staticmethod(python_func)

        setattr(ns, func_name_orig, getattr(ns, func_name))
        setattr(ns, func_name, python_func)

    return ns
