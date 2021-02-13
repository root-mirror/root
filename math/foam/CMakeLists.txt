# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/foam package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(Foam
  HEADERS
    TFoam.h
    TFoamCell.h
    TFoamIntegrand.h
    TFoamMaxwt.h
    TFoamSampler.h
    TFoamVect.h
  SOURCES
    src/TFoam.cxx
    src/TFoamCell.cxx
    src/TFoamIntegrand.cxx
    src/TFoamMaxwt.cxx
    src/TFoamSampler.cxx
    src/TFoamVect.cxx
  DEPENDENCIES
    Hist
    MathCore
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
