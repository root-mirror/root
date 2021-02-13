# Author: Enric Tejedor CERN  03/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
import cppyy

def set_size(self, buf):
    # Parameters:
    # - self: graph object
    # - buf: buffer of doubles
    # Returns:
    # - buffer whose size has been set
    buf.reshape((self.GetN(),))
    return buf

# Create a composite pythonizor.
#
# A composite is a type of pythonizor, i.e. it is a callable that expects two
# parameters: a class proxy and a string with the name of that class.
# A composite is created with the following parameters:
# - A string to match the class/es to be pythonized
# - A string to match the method/s to be pythonized in the class/es
# - A callable that will post-process the return value of the matched method/s
#
# Here we create a composite that will match TGraph, TGraph2D and their error
# subclasses, and will pythonize their getter methods of the X,Y,Z coordinate
# and error arrays, which in C++ return a pointer to a double.
# The pythonization consists in setting the size of the array that the getter
# method returns, so that it is known in Python and the array is fully usable
# (its length can be obtained, it is iterable).
comp = cppyy.py.compose_method('^TGraph(2D)?$|^TGraph.*Errors$', # class to match
                               'GetE?[XYZ]$',                    # method to match
                               set_size)                         # post-process function

# Invoke the decorator with the composite as argument to add it to the list of pythonizors
pythonization()(comp)
