# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/x11ttf package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(GX11TTF
  HEADERS
    TGX11TTF.h
  SOURCES
    src/TGX11TTF.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    ${X11_LIBRARIES}
    ${X11_Xpm_LIB}
    ${X11_Xft_LIB}
    ${FREETYPE_LIBRARIES}
    ZLIB::ZLIB
    ${CMAKE_THREAD_LIBS_INIT}
  DEPENDENCIES
    GX11
    Graf
)

target_include_directories(GX11TTF PRIVATE ${FREETYPE_INCLUDE_DIRS} ${X11_INCLUDE_DIR})
