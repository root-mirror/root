# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/vecops package
############################################################################

if(builtin_vdt)
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTVecOps
  HEADERS
    ROOT/RAdoptAllocator.hxx
    ROOT/RVec.hxx
  SOURCES
    src/RAdoptAllocator.cxx
    src/RVec.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
)

if(builtin_vdt OR vdt)
   target_include_directories(ROOTVecOps PRIVATE ${VDT_INCLUDE_DIRS} INTERFACE $<BUILD_INTERFACE:${VDT_INCLUDE_DIRS}>)
endif()

if(builtin_vdt)
  target_link_libraries(ROOTVecOps PRIVATE ${VDT_LIBRARIES})
elseif(vdt)
  target_link_libraries(ROOTVecOps PUBLIC ${VDT_LIBRARIES})
endif()

if(MSVC)
  target_compile_options(ROOTVecOps PRIVATE -O2 /fp:fast)
  target_compile_definitions(ROOTVecOps PRIVATE _USE_MATH_DEFINES)
else()
  target_compile_options(ROOTVecOps PRIVATE -O3 -ffast-math)
endif()

include(CheckCXXSymbolExists)
check_symbol_exists(m __sqrt_finite HAVE_FINITE_MATH)
if(NOT HAVE_FINITE_MATH AND NOT MSVC)
  target_compile_options(ROOTVecOps PRIVATE -fno-finite-math-only)
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
