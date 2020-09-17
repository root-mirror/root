# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

find_package(LZ4 REQUIRED)
find_package(ZLIB REQUIRED)

ROOT_OBJECT_LIBRARY(Zip
  src/Bits.c
  src/ZDeflate.c
  src/ZTrees.c
  src/ZInflate.c
  src/Compression.cxx
  src/RZip.cxx
  src/Bitshuffle.c
  src/Bitshuffle_core.c
  src/Bitshuffle_iochain.c
)

target_include_directories(Zip PRIVATE 
   ${LZ4_INCLUDE_DIR} 
   ${ZLIB_INCLUDE_DIR}
   ${CMAKE_SOURCE_DIR}/core/lzma/inc
   ${CMAKE_SOURCE_DIR}/core/lz4/inc
   ${CMAKE_SOURCE_DIR}/core/zstd/inc
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_SOURCE_DIR}/core/meta/inc
   ${CMAKE_SOURCE_DIR}/core/thread/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
