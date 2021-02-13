# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

###########################################################
# CMakeLists.txt file for building TPython
###########################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTTPython
  HEADERS
    TPyArg.h
    TPyReturn.h
    TPython.h
  SOURCES
    src/TPyArg.cxx
    src/TPyClassGenerator.cxx
    src/TPyReturn.cxx
    src/TPython.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
  LIBRARIES
    cppyy${PYTHON_UNDER_VERSION_STRING_Development_Main}
    # We link libTPython against Python libraries to compensate for the fact that libcppyy
    # is built with unresolved symbols. If we didn't do this, invoking TPython from C++
    # would not work.
    ${PYTHON_LIBRARIES_Development_Main}
)

# Disables warnings originating from deprecated register keyword in Python
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_STANDARD GREATER_EQUAL 11)
    target_compile_options(ROOTTPython PRIVATE -Wno-register)
endif()
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" AND CMAKE_CXX_STANDARD GREATER_EQUAL 11)
    target_compile_options(ROOTTPython PRIVATE -Wno-register)
    target_compile_options(ROOTTPython PRIVATE -Wno-deprecated-register)
endif()
