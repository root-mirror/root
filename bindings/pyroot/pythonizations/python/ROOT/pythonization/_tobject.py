# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
from libROOTPythonizations import AddTObjectEqNePyz
import cppyy

# Searching

def _contains(self, o):
    # Relies on TObject::FindObject
    # Parameters:
    # - self: object where to search
    # - o: object to be searched in self
    # Returns:
    # - True if self contains o
    return bool(self.FindObject(o))

# Comparison operators

def _lt(self, o):
    if isinstance(o, cppyy.gbl.TObject):
        return self.Compare(o) == -1
    else:
        return NotImplemented

def _le(self, o):
    if isinstance(o, cppyy.gbl.TObject):
        return self.Compare(o) <= 0
    else:
        return NotImplemented

def _gt(self, o):
    if isinstance(o, cppyy.gbl.TObject):
        return self.Compare(o) == 1
    else:
        return NotImplemented

def _ge(self, o):
    if isinstance(o, cppyy.gbl.TObject):
        return self.Compare(o) >= 0
    else:
        return NotImplemented


@pythonization(lazy = False)
def pythonize_tobject():
    klass = cppyy.gbl.TObject

    # Allow 'obj in container' syntax for searching
    klass.__contains__ = _contains

    # Inject comparison operators
    AddTObjectEqNePyz(klass)
    klass.__lt__ = _lt
    klass.__le__ = _le
    klass.__gt__ = _gt
    klass.__ge__ = _ge

    return True
