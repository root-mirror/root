# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization


@pythonization()
def pythonize_tstring(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: name of the class

    if name == 'TString':
        # Support `len(s)` as `s.Length()`
        klass.__len__ = klass.Length

        # Add string representation
        klass.__str__  = klass.Data
        klass.__repr__ = lambda self: "'{}'".format(self)

        # Add comparison operators
        klass.__eq__ = lambda self, o: str(self) == o
        klass.__ne__ = lambda self, o: str(self) != o
        klass.__lt__ = lambda self, o: str(self) <  o
        klass.__le__ = lambda self, o: str(self) <= o
        klass.__gt__ = lambda self, o: str(self) >  o
        klass.__ge__ = lambda self, o: str(self) >= o

    return True
