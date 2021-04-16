# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/foundationiolite package
############################################################################

if(WIN32)
  set(rawfile_local_headers ROOT/RRawFileWin.hxx)
  set(rawfile_local_sources src/RRawFileWin.cxx)
else()
  set(rawfile_local_headers ROOT/RRawFileUnix.hxx)
  set(rawfile_local_sources src/RRawFileUnix.cxx)
endif()

if(uring)
  list(APPEND rawfile_local_headers ROOT/RIoUring.hxx)
endif ()

# Used for generation of Core dictionary
set(Foundationiolite_dict_headers
  ${rawfile_local_headers}
  ROOT/RRawFile.hxx
  PARENT_SCOPE
)

set(ROOTFOUNDATIONIOLITE_SOURCES
  src/RRawFile.cxx
  ${rawfile_local_sources}
)

ROOT_LINKER_LIBRARY(ROOTFoundationIOLite ${ROOTFOUNDATIONIOLITE_SOURCES} DEPENDENCIES ROOTFoundation)

ROOT_INSTALL_HEADERS()
ROOT_ADD_TEST_SUBDIRECTORY(test)
