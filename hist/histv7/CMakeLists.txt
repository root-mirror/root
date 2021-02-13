# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT hist/hist package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTHist
  HEADERS
    ROOT/RAxis.hxx
    ROOT/RAxisConfig.hxx
    ROOT/RHist.hxx
    ROOT/RHistBinIter.hxx
    ROOT/RHistBufferedFill.hxx
    ROOT/RHistConcurrentFill.hxx
    ROOT/RHistData.hxx
    ROOT/RHistImpl.hxx
    ROOT/RHistUtils.hxx
    ROOT/RHistView.hxx
  SOURCES
    src/RAxis.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    MathCore
    Matrix
    RIO
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
