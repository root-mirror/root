# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

find_package(LZ4 REQUIRED)
find_package(xxHash REQUIRED)

ROOT_OBJECT_LIBRARY(Lz4 src/ZipLZ4.cxx)
target_include_directories(Lz4 PRIVATE
   ${LZ4_INCLUDE_DIR}
   ${xxHash_INCLUDE_DIR}
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
