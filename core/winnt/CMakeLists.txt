# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/winnt package
############################################################################

set(Winnt_dict_headers
  TWinNTSystem.h
  PARENT_SCOPE
)

ROOT_OBJECT_LIBRARY(Winnt
  src/TWin32SplashThread.cxx
  src/TWinNTSystem.cxx
  src/Win32Splash.cxx
)

target_include_directories(Winnt PRIVATE
   ${CMAKE_SOURCE_DIR}/core/foundation/res
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_SOURCE_DIR}/core/meta/inc
   ${CMAKE_SOURCE_DIR}/core/cont/inc
   ${CMAKE_SOURCE_DIR}/core/thread/inc
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_SOURCE_DIR}/core/textinput/inc
   ${CMAKE_SOURCE_DIR}/core/gui/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

ROOT_INSTALL_HEADERS()
