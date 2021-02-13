# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/clib package
############################################################################

set(Clib_dict_headers strlcpy.h
                      snprintf.h
                      strtok.h PARENT_SCOPE)

ROOT_OBJECT_LIBRARY(Clib
  src/Demangle.c
  src/attach.c
  src/detach.c
  src/getpagesize.c
  src/keys.c
  src/mcalloc.c
  src/mfree.c
  src/mmalloc.c
  src/mmapsup.c
  src/mmcheck.c
  src/mmemalign.c
  src/mmstats.c
  src/mmtrace.c
  src/mrealloc.c
  src/mvalloc.c
  src/sbrksup.c
  src/snprintf.c
  src/strlcat.c
  src/strlcpy.c
)

target_include_directories(Clib PRIVATE
   ${CMAKE_SOURCE_DIR}/core/clib/res
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
