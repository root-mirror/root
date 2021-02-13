# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/meta package
############################################################################

set(GuiCore_dict_headers
  GuiTypes.h
  TApplicationImp.h
  TBrowser.h
  TBrowserImp.h
  TCanvasImp.h
  TClassMenuItem.h
  TContextMenu.h
  TContextMenuImp.h
  TControlBarImp.h
  TGuiFactory.h
  TInspectorImp.h
  TObjectSpy.h
  TToggleGroup.h
  TToggle.h
  PARENT_SCOPE
)

ROOT_OBJECT_LIBRARY(GuiCore
  src/InitGui.cxx
  src/TApplicationImp.cxx
  src/TBrowser.cxx
  src/TBrowserImp.cxx
  src/TCanvasImp.cxx
  src/TClassMenuItem.cxx
  src/TContextMenu.cxx
  src/TContextMenuImp.cxx
  src/TControlBarImp.cxx
  src/TGuiFactory.cxx
  src/TInspectorImp.cxx
  src/TObjectSpy.cxx
  src/TToggle.cxx
  src/TToggleGroup.cxx
)

target_include_directories(GuiCore PRIVATE 
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_SOURCE_DIR}/core/meta/inc
   ${CMAKE_SOURCE_DIR}/core/cont/inc
   ${CMAKE_SOURCE_DIR}/core/thread/inc
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
