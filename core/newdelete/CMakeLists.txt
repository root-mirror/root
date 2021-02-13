# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/newdelete package
############################################################################

ROOT_LINKER_LIBRARY(New
  src/MemCheck.cxx
  src/NewDelete.cxx
  DEPENDENCIES
    Core
)

target_include_directories(New PRIVATE ${CMAKE_SOURCE_DIR}/core/clib/res)

ROOT_INSTALL_HEADERS()
