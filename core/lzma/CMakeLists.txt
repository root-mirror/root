# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/lzma package
############################################################################

ROOT_OBJECT_LIBRARY(Lzma src/ZipLZMA.c BUILTINS LZMA)

target_include_directories(Lzma PRIVATE
   ${LIBLZMA_INCLUDE_DIR}
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
