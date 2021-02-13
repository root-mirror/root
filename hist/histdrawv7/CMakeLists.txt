# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT hist/histdraw package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTHistDraw
  HEADERS
    ROOT/RHistDisplayItem.hxx
    ROOT/RHistDrawable.hxx
    ROOT/RHistStatBox.hxx
  SOURCES
    src/RHistDisplayItem.cxx
    src/RHistDrawable.cxx
    src/RHistStatBox.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    MathCore
    ROOTHist
    ROOTGpadv7
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
