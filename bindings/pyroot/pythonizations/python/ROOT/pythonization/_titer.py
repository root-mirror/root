# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
import cppyy

def _next_pyz(self):
    # Parameters:
    # - self: iterator on a collection
    # Returns:
    # - next object in the collection, or raises StopIteration if none
	o = self.Next()
	if o:
		return o
	else:
		raise StopIteration()


@pythonization(lazy = False)
def pythonize_titer():
    klass = cppyy.gbl.TIter

    # Make TIter a Python iterable
    klass.__iter__ = lambda self: self

    # Make TIter a Python iterator
    klass.__next__ = _next_pyz  # Py3
    klass.next     = _next_pyz  # Py2

    return True
