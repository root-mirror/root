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
from libROOTPythonizations import AddTClassDynamicCastPyz


@pythonization(lazy = False)
def pythonize_tclass():
    klass = cppyy.gbl.TClass

    # DynamicCast
    AddTClassDynamicCastPyz(klass)

    return True
