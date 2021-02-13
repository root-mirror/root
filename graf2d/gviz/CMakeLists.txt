# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/gviz package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(Gviz
  HEADERS
    TGraphEdge.h
    TGraphNode.h
    TGraphStruct.h
  SOURCES
    src/TGraphEdge.cxx
    src/TGraphNode.cxx
    src/TGraphStruct.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    ${GRAPHVIZ_LIBRARIES}
  DEPENDENCIES
    Hist
    Graf
)

target_include_directories(Gviz PRIVATE ${GRAPHVIZ_INCLUDE_DIR}/graphviz)
