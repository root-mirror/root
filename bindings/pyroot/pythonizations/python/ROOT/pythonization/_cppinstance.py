# Author: Massimiliano Galli CERN  06/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
import cppyy
import libcppyy
from libROOTPythonizations import AddCPPInstancePickling

@pythonization(lazy = False)
def pythonize_cppinstance():
    klass = libcppyy.CPPInstance

    AddCPPInstancePickling(klass)

    return True
