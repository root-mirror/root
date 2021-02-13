# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(periodic)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

file(GLOB DBsources RELATIVE ${CMAKE_SOURCE_DIR} Ndb*.cxx)
file(GLOB XSsources RELATIVE ${CMAKE_SOURCE_DIR} XS*.cxx)
file(GLOB headers RELATIVE ${CMAKE_SOURCE_DIR} *.h)
list(REMOVE_ITEM headers NdbLinkDef.h)

ROOT_GENERATE_DICTIONARY(G__NdbDict ${headers} LINKDEF NdbLinkDef.h)

#---Create a shared library with geneated dictionary

add_library(NdbDict SHARED ${DBsources} G__NdbDict.cxx)

target_include_directories(NdbDict PRIVATE ${CMAKE_SOURCE_DIR})

# Do not add -Dname_EXPORTS to the command-line when building files in this
# target. Doing so is actively harmful for the modules build because it
# creates extra module variants, and not useful because we don't use these
# macros.
set_target_properties(${NdbDict} PROPERTIES DEFINE_SYMBOL "")

if(MSVC)
  set_target_properties(NdbDict PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()
target_link_libraries(NdbDict ROOT::Core)

#---Create  a main program using the library
add_executable(XSGui ${XSsources})
target_link_libraries(XSGui NdbDict ROOT::Gui)
