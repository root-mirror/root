# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT net/netxng package
# @author Lukasz Janyst <ljanyst@cern.ch>
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(NetxNG
  HEADERS
    TNetXNGFile.h
    TNetXNGFileStager.h
    TNetXNGSystem.h
  SOURCES
    src/TNetXNGFile.cxx
    src/TNetXNGFileStager.cxx
    src/TNetXNGSystem.cxx
  LIBRARIES
    ${XROOTD_LIBRARIES}
  DEPENDENCIES
    Net
    RIO
    Thread
  BUILTINS
    XROOTD
)

target_include_directories(NetxNG PRIVATE ${XROOTD_INCLUDE_DIRS})

target_compile_options(NetxNG PRIVATE -Wno-shadow)
