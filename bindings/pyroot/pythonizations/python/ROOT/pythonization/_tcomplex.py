# Author: Massimiliano Galli CERN  05/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
import cppyy

def _rsub(self, other):
    # Parameters:
    # - self: complex number
    # - other: other (in general not complex) number
    return -self+other

def _perform_division(self, other):
    # Parameters:
    # - self: complex number
    # - other: int, float of long (Py2) number
    TComplex = cppyy.gbl.TComplex
    other_complex = TComplex.TComplex(other,0)
    return other_complex/self

def _rdiv(self, other):
    # Parameters:
    # - self: complex number
    # - other: other term
    if isinstance(other, (int, long, float)):
        return _perform_division(self, other)
    else:
        return NotImplemented

def _rtruediv(self, other):
    # Parameters:
    # - self: complex number
    # - other: other term
    if isinstance(other, (int, float)):
        return _perform_division(self, other)
    else:
        return NotImplemented

@pythonization()
def pythonize_tcomplex(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class
    if name == 'TComplex':

        # implements __radd__ as equal to __add__
        klass.__radd__ = klass.__add__

        # implements __rsub__ by assigning the function previously defined
        klass.__rsub__ = _rsub

        # implements __rmul__ as equal to __mul__
        klass.__rmul__ = klass.__mul__

        # implements __rtruediv__ by assigning the function previously defined
        # necessary for Python3
        klass.__rtruediv__ = _rtruediv

        # implements __rtruediv__ by assigning the function previously defined
        # necessary for Python2
        klass.__rdiv__ = _rdiv

    return True

