# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# - Locate gl2ps library
# Defines:
#
#  GL2PS_FOUND
#  GL2PS_INCLUDE_DIR
#  GL2PS_INCLUDE_DIRS (not cached)
#  GL2PS_LIBRARIES

find_path(GL2PS_INCLUDE_DIR NAMES gl2ps.h HINTS ${GL2PS_DIR}/include $ENV{GL2PS_DIR}/include /usr/include)
find_library(GL2PS_LIBRARY NAMES gl2ps HINTS ${GL2PS_DIR}/lib $ENV{GL2PS_DIR}/lib)

set(GL2PS_INCLUDE_DIRS ${GL2PS_INCLUDE_DIR})
if(GL2PS_LIBRARY)
  set(GL2PS_LIBRARIES ${GL2PS_LIBRARY})
endif()


# handle the QUIETLY and REQUIRED arguments and set GL2PS_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(gl2ps DEFAULT_MSG GL2PS_LIBRARY GL2PS_INCLUDE_DIR)

mark_as_advanced(GL2PS_FOUND GL2PS_INCLUDE_DIR GL2PS_LIBRARY)
