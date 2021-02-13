# Author: Stefan Wunsch, Enric Tejedor CERN  06/2018

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from libROOTPythonizations import AddPrettyPrintingPyz
from ROOT import pythonization

def _add_getitem_checked(klass):
    # Parameters:
    # - klass: class where to add a __getitem__ method that raises
    # IndexError if index is out of range

    def getitem_checked(o, i):
        # Get item of `o` at `i` or raise IndexError if index is
        # out of range.
        # Assumes `o` has `__len__`.
        # Parameters:
        # - o: object
        # - i: index to be checked in object
        # Returns:
        # - o[i]
        if i >= 0 and i < len(o):
            return o._getitem__unchecked(i)
        else:
            raise IndexError('index out of range')

    klass._getitem__unchecked = klass.__getitem__
    klass.__getitem__ = getitem_checked


@pythonization()
def pythonizegeneric(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    # Add pretty printing via setting the __str__ special function

    # Exclude classes which have the method __str__ already defined in C++
    m = getattr(klass, '__str__', None)
    has_cpp_str = True if m is not None and type(m).__name__ == 'CPPOverload' else False

    # Exclude std::string with its own pythonization from cppyy
    exclude = [ 'std::string' ]

    if name not in exclude and not has_cpp_str:
        AddPrettyPrintingPyz(klass)

    return True
