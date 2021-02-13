# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf3d/x3d package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(X3d
  HEADERS
    TViewerX3D.h
    TX3DFrame.h
  SOURCES
    src/x3d.c
    src/TViewerX3D.cxx
    src/TX3DFrame.cxx
  LIBRARIES
    ${X11_LIBRARIES}
  DEPENDENCIES
    Graf3d
    Gui
)

target_include_directories(X3d PRIVATE ${X11_INCLUDE_DIR})
