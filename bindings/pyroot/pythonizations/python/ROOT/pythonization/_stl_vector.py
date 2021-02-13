# Author: Stefan Wunsch CERN  08/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
from ROOT.pythonization._rvec import add_array_interface_property


@pythonization()
def pythonize_stl_vector(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name.startswith("std::vector<"):
        # Add numpy array interface
        # NOTE: The pythonization is reused from ROOT::VecOps::RVec
        add_array_interface_property(klass, name)

    return True
