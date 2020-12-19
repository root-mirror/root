# Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT R package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(RMVA
  HEADERS
    TMVA/MethodC50.h
    TMVA/MethodRSNNS.h
    TMVA/MethodRSVM.h
    TMVA/MethodRXGB.h
    TMVA/RMethodBase.h
  SOURCES
    src/MethodC50.cxx
    src/MethodRSNNS.cxx
    src/MethodRSVM.cxx
    src/MethodRXGB.cxx
    src/RMethodBase.cxx
  DEPENDENCIES
    Core
    Matrix
    RInterface
    RIO
    Thread
    TMVA
)
