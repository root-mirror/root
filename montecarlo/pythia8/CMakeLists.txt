# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT montecarlo/pythia8 package
# @author Pere Mato, CERN
############################################################################

# required also for cling
include_directories(${PYTHIA8_INCLUDE_DIR})

ROOT_STANDARD_LIBRARY_PACKAGE(EGPythia8
  HEADERS
    TPythia8.h
    TPythia8Decayer.h
  SOURCES
    src/TPythia8.cxx
    src/TPythia8Decayer.cxx
  DEPENDENCIES
    Core
    EG
    Graf
    Physics
)

ROOT_ADD_CXX_FLAG(_EGPythia8_FLAGS -Wno-delete-non-virtual-dtor)
ROOT_ADD_CXX_FLAG(_EGPythia8_FLAGS -Wno-misleading-indentation)
ROOT_ADD_CXX_FLAG(_EGPythia8_FLAGS -Wno-overloaded-virtual)

separate_arguments(_EGPythia8_FLAGS)

target_compile_options(EGPythia8 PRIVATE ${_EGPythia8_FLAGS})
target_link_libraries(EGPythia8 PUBLIC ${PYTHIA8_LIBRARIES})
