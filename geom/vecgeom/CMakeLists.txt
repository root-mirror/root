# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT geom/vecgeom package
# @author Mihaela Gheata, CERN
############################################################################

include_directories(AFTER SYSTEM ${VECGEOM_INCLUDE_DIRS})
if ( Vc_FOUND )
   include_directories( AFTER SYSTEM ${Vc_INCLUDE_DIRS})
   if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
       AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=6")
   endif()
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(ConverterVG
  HEADERS
    TGeoVGConverter.h
    TGeoVGShape.h
  SOURCES
    src/TGeoVGConverter.cxx
    src/TGeoVGShape.cxx
  LIBRARIES
    ${VECGEOM_LIBRARIES}
    ${CMAKE_DL_LIBS}
  DEPENDENCIES
    Geom
)
