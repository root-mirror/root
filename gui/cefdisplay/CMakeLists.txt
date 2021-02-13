# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT gui/cefdisplay package
############################################################################

if(DEFINED ENV{CEF_ROOT})

set(libname ROOTCefDisplay)

set(CEF_sources src/simple_app.cxx src/gui_handler.cxx src/RCefWebDisplayHandle.cxx)

set(CEF_root $ENV{CEF_ROOT})

# ensure bin directory is exists where CEF files should be copied
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(ROOT_ARCHITECTURE MATCHES macosx)
   set(CEF_platform src/gui_handler_mac.mm)
   set(CEF_RELEASE_DIR ${CEF_root}/Release/Chromium\ Embedded\ Framework.framework)
   set(CEF_LIBRARY ${CEF_RELEASE_DIR}/Chromium\ Embedded\ Framework)
   set(CEF_LIB_DEPENDENCY "-framework Foundation")
   set(CEF_DLL_WRAPPER ${CEF_root}/build/libcef_dll_wrapper/libcef_dll_wrapper.a)
   file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Frameworks/Chromium\ Embedded\ Framework.framework)
   file(COPY ${CEF_RELEASE_DIR} DESTINATION ${CMAKE_BINARY_DIR}/Frameworks/)
elseif(MSVC)
   set(CEF_platform src/gui_handler_win.cc)
   set(CEF_RESOURCES ${CEF_root}/Resources)
   set(CEF_RELEASE_DIR ${CEF_root}/Release)
   set(CEF_LIBRARY ${CEF_RELEASE_DIR}/libcef.lib)
   set(CEF_DLL_WRAPPER ${CEF_root}/build/libcef_dll_wrapper/Release/libcef_dll_wrapper.lib)
   file(COPY ${CEF_RELEASE_DIR}/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
   file(COPY ${CEF_RESOURCES}/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
else()
   set(CEF_platform src/gui_handler_linux.cxx)
   set(CEF_RESOURCES ${CEF_root}/Resources)
   set(CEF_RELEASE_DIR ${CEF_root}/Release)
   set(CEF_LIBRARY ${CEF_RELEASE_DIR}/libcef.so)
   set(CEF_DLL_WRAPPER ${CEF_root}/build/libcef_dll_wrapper/libcef_dll_wrapper.a)
   set(CEF_LIB_DEPENDENCY ${X11_LIBRARIES})
   file(COPY ${CEF_RELEASE_DIR}/ DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
   file(COPY ${CEF_RESOURCES}/ DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
   file(COPY ${CEF_RELEASE_DIR}/swiftshader DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

include_directories($ENV{CEF_ROOT} inc)
set(CEF_MAIN src/cef_main.cxx)

ROOT_LINKER_LIBRARY(${libname} ${CEF_sources} ${CEF_platform}
                    LIBRARIES ${CMAKE_DL_LIBS} ${CEF_LIBRARY} ${CEF_DLL_WRAPPER} ${CEF_LIB_DEPENDENCY}
                    DEPENDENCIES RHTTP ROOTWebDisplay)

ROOT_EXECUTABLE(cef_main ${CEF_MAIN} LIBRARIES ${CEF_LIBRARY} ${CEF_DLL_WRAPPER})

endif()
