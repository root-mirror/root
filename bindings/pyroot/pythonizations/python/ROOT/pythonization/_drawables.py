# Author: Enric Tejedor CERN  04/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
from cppyy.gbl import kCanDelete
from libcppyy import SetOwnership


def _Draw(self, *args):
    # Parameters:
    # self: Object being drawn
    # args: arguments for Draw

    self._OriginalDraw(*args)

    # When drawing a TPad, it gets added to the list of primititves of its
    # mother TPad (fMother) with kCanDelete == 1. This means that, when
    # fMother is destructed, it will attempt to destroy its child TPad too.
    # To prevent a double delete, here we instruct the Python proxy of the
    # child C++ TPad being drawn not to destroy the latter (ROOT-10060).
    # 
    # A similar principle is applied to TButton, TColorWheel, TPolyLine3D,
    # TPolyMarker and TPolyMarker3D, whose kCanDelete bit is set in one of
    # their constructors. Later, when being drawn, they are appended to
    # the list of primitives of gPad.
    if self.TestBit(kCanDelete):
        SetOwnership(self, False)

    self.Draw = self._OriginalDraw

def _init(self, *args):
    # Parameters:
    # self: Object being initialized
    # args: arguments for __init__

    self._original__init__(*args)

    # TSlider is a special case, since it is appended to gPad already
    # in one of its constructors, after setting kCanDelete.
    # Therefore, we need to set the ownership here and not in Draw
    # (TSlider does not need to be drawn). This is ROOT-10095.
    if self.TestBit(kCanDelete):
        SetOwnership(self, False)
        # We have already set the ownership while initializing,
        # so we do not need the custom Draw inherited from TPad to
        # do it again in case it is executed.
        self.Draw = self._OriginalDraw

@pythonization()
def pythonize_drawables(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: name of the class

    if name in {'TPad', 'TButton', 'TColorWheel',
                'TPolyLine3D', 'TPolyMarker', 'TPolyMarker3D'}:
        # Draw
        klass._OriginalDraw = klass.Draw
        klass.Draw = _Draw
    elif name == 'TSlider':
        # __init__
        klass._original__init__ = klass.__init__
        klass.__init__ = _init

    return True
