# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/imt package
############################################################################

if(NOT WIN32)
  set(MULTIPROC_LIB "MultiProc")
endif()

ROOT_LINKER_LIBRARY(Imt
    src/base.cxx
    src/TExecutor.cxx
    src/TTaskGroup.cxx
  DEPENDENCIES
    ${MULTIPROC_LIB}
  BUILTINS
    TBB
)

target_include_directories(Clib PRIVATE
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

target_link_libraries(Imt PRIVATE Thread INTERFACE Core)

if(imt)
  ROOT_GENERATE_DICTIONARY(G__Imt STAGE1
    ROOT/TFuture.hxx
    ROOT/TTaskGroup.hxx
    ROOT/RTaskArena.hxx
    ROOT/TExecutor.hxx
    ROOT/TThreadExecutor.hxx
    LINKDEF
      LinkDef.h
    MODULE
      Imt
    DEPENDENCIES
      Core
      Thread
      MultiProc
    BUILTINS
      TBB
  )

  # G__Imt.cxx is automatically added by ROOT_GENERATE_DICTIONARY()
  target_sources(Imt PRIVATE
    src/RTaskArena.cxx
    src/TImplicitMT.cxx
    src/TThreadExecutor.cxx
  )

  target_include_directories(Imt SYSTEM PRIVATE ${TBB_INCLUDE_DIRS})
  target_link_libraries(Imt PRIVATE ${TBB_LIBRARIES})
  set_target_properties(Imt PROPERTIES COMPILE_FLAGS "${TBB_CXXFLAGS}")

  ROOT_ADD_TEST_SUBDIRECTORY(test)
else()
  ROOT_GENERATE_DICTIONARY(G__Imt STAGE1
    ROOT/TExecutor.hxx
    LINKDEF
      LinkDef.h
    MODULE
      Imt
    DEPENDENCIES
      Core
      MultiProc
  )
endif()

ROOT_INSTALL_HEADERS()
