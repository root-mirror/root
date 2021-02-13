# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/asimage package
# @author Pere Mato, CERN
############################################################################

if(cocoa)  # Not needed to locate PNG, JPG, TIFF libraries
  set(ASEXTRA_LIBRARIES)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(ASImage
  NO_INSTALL_HEADERS
  HEADERS
    TASImage.h
    TASImagePlugin.h
    TASPluginGS.h
  SOURCES
    TASImage.cxx
    TASPluginGS.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    ${AFTERIMAGE_LIBRARIES}
    ${ASEXTRA_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${X11_LIBRARIES}
    ZLIB::ZLIB
  DEPENDENCIES
    Core
    Graf
  BUILTINS
    AFTERIMAGE
)

target_include_directories(ASImage SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS} ${AFTERIMAGE_INCLUDE_DIR} ${X11_INCLUDE_DIR})

ROOT_STANDARD_LIBRARY_PACKAGE(ASImageGui
  NO_INSTALL_HEADERS
  HEADERS
    TASPaletteEditor.h
  SOURCES
    src/TASPaletteEditor.cxx
  LINKDEF
    LinkDefGui.h
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    ${AFTERIMAGE_LIBRARIES}
    ${ASEXTRA_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${X11_LIBRARIES}
    ZLIB::ZLIB
  DEPENDENCIES
    Gui
    ASImage
)

target_include_directories(ASImageGui SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS} ${AFTERIMAGE_INCLUDE_DIR} ${X11_INCLUDE_DIR})

ROOT_INSTALL_HEADERS()
