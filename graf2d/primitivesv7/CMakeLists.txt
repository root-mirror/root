# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/primitives package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTGraphicsPrimitives
  HEADERS
    ROOT/RBox.hxx
    ROOT/RFrameTitle.hxx
    ROOT/RLegend.hxx
    ROOT/RLine.hxx
    ROOT/RMarker.hxx
    ROOT/RPaveText.hxx
    ROOT/RText.hxx
  SOURCES
    src/RBox.cxx
    src/RLegend.cxx
    src/RLine.cxx
    src/RMarker.cxx
    src/RText.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    ROOTGpadv7
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
