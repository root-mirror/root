# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# Try to find Graphviz.
# This will define:
# GRAPHVIZ_FOUND - system has Graphviz
# GRAPHVIZ_INCLUDE_DIR - the Graphviz include directory
# GRAPHVIZ_xxx_LIBRARY - Graphviz libraries
# GRAPHVIZ_LIBRARIES - Link these to use Graphviz (not cached)

find_path(GRAPHVIZ_INCLUDE_DIR graphviz/gvc.h HINTS ${GRAPHVIZ_DIR} ENV GRAPHVIZ_DIR PATH_SUFFIXES include)

find_library(GRAPHVIZ_cdt_LIBRARY NAMES cdt HINTS ${GRAPHVIZ_DIR} ENV GRAPHVIZ_DIR PATH_SUFFIXES lib)
find_library(GRAPHVIZ_gvc_LIBRARY NAMES gvc HINTS ${GRAPHVIZ_DIR} ENV GRAPHVIZ_DIR PATH_SUFFIXES lib)
find_library(GRAPHVIZ_graph_LIBRARY NAMES graph cgraph HINTS ${GRAPHVIZ_DIR} ENV GRAPHVIZ_DIR PATH_SUFFIXES lib)
find_library(GRAPHVIZ_pathplan_LIBRARY NAMES pathplan HINTS ${GRAPHVIZ_DIR} ENV GRAPHVIZ_DIR PATH_SUFFIXES lib)

set(GRAPHVIZ_LIBRARIES ${GRAPHVIZ_gvc_LIBRARY} ${GRAPHVIZ_graph_LIBRARY} ${GRAPHVIZ_cdt_LIBRARY} ${GRAPHVIZ_pathplan_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set GRAPHVIZ_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Graphviz DEFAULT_MSG GRAPHVIZ_INCLUDE_DIR
                                                       GRAPHVIZ_cdt_LIBRARY
                                                       GRAPHVIZ_gvc_LIBRARY
                                                       GRAPHVIZ_graph_LIBRARY
                                                       GRAPHVIZ_pathplan_LIBRARY)

mark_as_advanced(GRAPHVIZ_INCLUDE_DIR
                 GRAPHVIZ_cdt_LIBRARY
                 GRAPHVIZ_graph_LIBRARY
                 GRAPHVIZ_gvc_LIBRARY
                 GRAPHVIZ_pathplan_LIBRARY)
