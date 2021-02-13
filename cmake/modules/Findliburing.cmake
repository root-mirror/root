# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# - Locate liburing library
#
# Defines:
#
# LIBURING_FOUND
# LIBURING_LIBRARY
# LIBURING_LIBRARY_PATH
# LIBURING_INCLUDE_DIR

find_library(LIBURING_LIBRARY NAMES uring)
find_path(LIBURING_INCLUDE_DIR NAMES liburing.h)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(liburing DEFAULT_MSG LIBURING_LIBRARY LIBURING_INCLUDE_DIR)

mark_as_advanced(LIBURING_FOUND LIBURING_LIBRARY LIBURING_INCLUDE_DIR)
get_filename_component(LIBURING_LIBRARY_PATH ${LIBURING_LIBRARY} DIRECTORY)
