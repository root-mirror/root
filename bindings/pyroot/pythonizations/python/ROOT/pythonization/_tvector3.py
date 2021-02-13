# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization

from ._generic import _add_getitem_checked


@pythonization()
def pythonize_tvector3(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name == 'TVector3':
        # `len(v)` is always 3
        klass.__len__ = lambda _: 3

        # Add checked __getitem__.
        # Allows to throw pythonic IndexError when index is out of range
        # and to iterate over the vector.
        _add_getitem_checked(klass)

    return True
