# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT gui/webgui6 package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(WebGui6
  HEADERS
    TWebCanvas.h
    TWebGuiFactory.h
    TWebMenuItem.h
    TWebPadOptions.h
    TWebPadPainter.h
    TWebPainting.h
    TWebPS.h
    TWebSnapshot.h
  SOURCES
    src/TWebCanvas.cxx
    src/TWebGuiFactory.cxx
    src/TWebMenuItem.cxx
    src/TWebPadPainter.cxx
    src/TWebPainting.cxx
    src/TWebPS.cxx
    src/TWebSnapshot.cxx
  DEPENDENCIES
    Core
    Gui
    ROOTWebDisplay
)
