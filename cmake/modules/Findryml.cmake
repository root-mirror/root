# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# - Locate ryml library
#
# Defines:
#
# RYML_FOUND
# RYML_LIBRARY
# RYML_LIBRARY_PATH
# RYML_INCLUDE_DIR

find_library(RYML_LIBRARY NAMES ryml HINTS ${RYML_DIR}/lib $ENV{RYML_DIR}/lib)
find_path(RYML_INCLUDE_DIR NAMES ryml.hpp HINTS ${RYML_DIR}/include $ENV{RYML_DIR}/include)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ryml DEFAULT_MSG RYML_LIBRARY RYML_INCLUDE_DIR)

mark_as_advanced(RYML_FOUND RYML_LIBRARY RYML_INCLUDE_DIR)
get_filename_component(RYML_LIBRARY_PATH ${RYML_LIBRARY} DIRECTORY)
