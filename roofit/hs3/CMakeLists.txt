# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT rootfit/histfactory package
# @author Carsten Burgard, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(HS3
  HEADERS
  JSONInterface.h
  RooJSONFactoryWSTool.h
  HistFactoryJSONTool.h
  JSONParser.h
  RYMLParser.h
  SOURCES
  src/RooJSONFactoryWSTool.cxx
  src/HistFactoryJSONTool.cxx
  src/JSONParser.cxx
  src/RYMLParser.cxx
  src/JSONFactories_RooFitCore.cxx
  src/JSONFactories_HistFactory.cxx
  DICTIONARY_OPTIONS
    "-writeEmptyRootPCM"
  DEPENDENCIES
    RooFit
    RooFitCore
    RooStats
    HistFactory    
)

if(ryml)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DINCLUDE_RYML")
  target_include_directories(HS3 PRIVATE
    ${RYML_INCLUDE_DIRS}
    )
  target_link_libraries(HS3 PRIVATE -lc4core -lryml)
  target_link_directories(HS3 PRIVATE
    ${RYML_LIB_DIR}
    )
else()
  target_include_directories(HS3 PRIVATE
    ${CMAKE_SOURCE_DIR}/io/io/res
    )
endif()
  
