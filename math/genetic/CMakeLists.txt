# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/genetic package
# @author Pere Mato, CERN
############################################################################
ROOT_STANDARD_LIBRARY_PACKAGE(Genetic
  HEADERS
    Math/GeneticMinimizer.h
  SOURCES
    src/GeneticMinimizer.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    MathCore
    TMVA
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
