# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT net/davix package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(RDAVIX
  HEADERS
    ROOT/RRawFileDavix.hxx
    TDavixFile.h
    TDavixSystem.h
  SOURCES
    src/RRawFileDavix.cxx
    src/TDavixFile.cxx
    src/TDavixSystem.cxx
  LIBRARIES
    Davix::Davix
  DEPENDENCIES
    Net
    RIO
)

target_compile_options(RDAVIX PRIVATE -Wno-deprecated-declarations)

ROOT_ADD_TEST_SUBDIRECTORY(test)
