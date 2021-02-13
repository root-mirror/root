# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#---------------------------------------------------------------------------------------------------
#  RootCPack.cmake
#   - basic setup for packaging ROOT using CTest
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Package up needed system libraries - only for WIN32?
#
include(InstallRequiredSystemLibraries)

#----------------------------------------------------------------------------------------------------
# General packaging setup - variable relavant to all package formats
#
set(CPACK_PACKAGE_DESCRIPTION "ROOT project")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ROOT project")
set(CPACK_PACKAGE_VENDOR "ROOT project")
set(CPACK_PACKAGE_VERSION ${ROOT_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${ROOT_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${ROOT_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${ROOT_PATCH_VERSION})

string(REGEX REPLACE "^([0-9]+).*$" "\\1" CXX_MAJOR ${CMAKE_CXX_COMPILER_VERSION})
string(REGEX REPLACE "^([0-9]+)\\.([0-9]+).*$" "\\2" CXX_MINOR ${CMAKE_CXX_COMPILER_VERSION})

#---Resource Files-----------------------------------------------------------------------------------
configure_file(LICENSE LICENSE.txt COPYONLY)
configure_file(LGPL2_1.txt LGPL2_1.txt COPYONLY)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")
if (APPLE)
  # Apple productbuild cannot handle .md files as CPACK_PACKAGE_DESCRIPTION_FILE;
  # convert to HTML instead.
  find_program(CONVERTER textutil)
  if (NOT CONVERTER)
    message(FATAL_ERROR "textutil executable not found")
  endif()
  execute_process(COMMAND ${CONVERTER} -convert html "${CMAKE_SOURCE_DIR}/README.md" -output "${CMAKE_BINARY_DIR}/README.html")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_BINARY_DIR}/README.html")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_BINARY_DIR}/README.html")
else()
  configure_file(README.md README.md COPYONLY)
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_BINARY_DIR}/README.md")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_BINARY_DIR}/README.md")
endif()

#---Source package settings--------------------------------------------------------------------------
set(CPACK_SOURCE_IGNORE_FILES
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/tests
    "~$"
    "/CVS/"
    "/.svn/"
    "/\\\\\\\\.svn/"
    "/.git/"
    "/\\\\\\\\.git/"
    "\\\\\\\\.swp$"
    "\\\\\\\\.swp$"
    "\\\\.swp"
    "\\\\\\\\.#"
    "/#"
)
set(CPACK_SOURCE_STRIP_FILES "")

#---Binary package setup-----------------------------------------------------------------------------
if(MSVC)
  if (MSVC_VERSION LESS 1900)
    math(EXPR VS_VERSION "${VC_MAJOR} - 6")
  elseif(MSVC_VERSION LESS 1910)
    math(EXPR VS_VERSION "${VC_MAJOR} - 5")
  elseif(MSVC_VERSION LESS 1919)
    math(EXPR VS_VERSION "${VC_MAJOR} - 4")
  elseif(MSVC_VERSION LESS 1925)
    math(EXPR VS_VERSION "${VC_MAJOR} - 3")
  endif()
  set(COMPILER_NAME_VERSION ".vc${VS_VERSION}")
else()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(COMPILER_NAME_VERSION "-gcc${CXX_MAJOR}.${CXX_MINOR}")
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(COMPILER_NAME_VERSION "-clang${CXX_MAJOR}${CXX_MINOR}")
  endif()
endif()

#---Processor architecture---------------------------------------------------------------------------

set(arch ${CMAKE_SYSTEM_PROCESSOR})

#---OS and version-----------------------------------------------------------------------------------
if(APPLE)
  execute_process(COMMAND sw_vers "-productVersion"
                  COMMAND cut -d . -f 1-2
                  OUTPUT_VARIABLE osvers OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(OS_NAME_VERSION macos-${osvers}-${CMAKE_SYSTEM_PROCESSOR})
elseif(WIN32)
  set(OS_NAME_VERSION win32)
else()
  execute_process(COMMAND lsb_release -is OUTPUT_VARIABLE osid OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND lsb_release -rs OUTPUT_VARIABLE osvers OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(osid MATCHES Ubuntu)
    string(REGEX REPLACE "([0-9]+)[.].*" "\\1" osvers "${osvers}")
    set(OS_NAME_VERSION Linux-ubuntu${osvers}-${arch})
  elseif(osid MATCHES Scientific)
    string(REGEX REPLACE "([0-9]+)[.].*" "\\1" osvers "${osvers}")
    set(OS_NAME_VERSION Linux-slc${osvers}-${arch})
  elseif(osid MATCHES Fedora)
    string(REGEX REPLACE "([0-9]+)" "\\1" osvers "${osvers}")
    set(OS_NAME_VERSION Linux-fedora${osvers}-${arch})
  elseif(osid MATCHES CentOS)
    string(REGEX REPLACE "([0-9]+)[.].*" "\\1" osvers "${osvers}")
    set(OS_NAME_VERSION Linux-centos${osvers}-${arch})
  else()
    set(OS_NAME_VERSION Linux-${osid}${osvers}${arch})
  endif()
endif()
#---Build type---------------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE STREQUAL Release)
  string(TOLOWER .${CMAKE_BUILD_TYPE} BUILD_TYPE_FOR_NAME)
endif()

set(CPACK_PACKAGE_RELOCATABLE True)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "root_v${ROOT_VERSION}")
set(CPACK_PACKAGE_FILE_NAME "root_v${ROOT_VERSION}.${OS_NAME_VERSION}${COMPILER_NAME_VERSION}${BUILD_TYPE_FOR_NAME}")
set(CPACK_PACKAGE_EXECUTABLES "root" "ROOT")

if(WIN32)
  set(CPACK_GENERATOR "ZIP;NSIS")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
elseif(APPLE)
  set(CPACK_GENERATOR "TGZ;productbuild")
  set(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
endif()

#----------------------------------------------------------------------------------------------------
# Finally, generate the CPack per-generator options file and include the
# base CPack configuration.
#
configure_file(cmake/modules/CMakeCPackOptions.cmake.in CMakeCPackOptions.cmake @ONLY)
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_BINARY_DIR}/CMakeCPackOptions.cmake)
include(CPack)

#----------------------------------------------------------------------------------------------------
# Define components and installation types (after CPack included!)
#
cpack_add_install_type(full      DISPLAY_NAME "Full Installation")
cpack_add_install_type(minimal   DISPLAY_NAME "Minimal Installation")
cpack_add_install_type(developer DISPLAY_NAME "Developer Installation")

cpack_add_component(applications
    DISPLAY_NAME "ROOT Applications"
    DESCRIPTION "ROOT executables such as root.exe"
     INSTALL_TYPES full minimal developer)

cpack_add_component(libraries
    DISPLAY_NAME "ROOT Libraries"
    DESCRIPTION "All ROOT libraries and dictionaries"
     INSTALL_TYPES full minimal developer)

cpack_add_component(headers
    DISPLAY_NAME "C++ Headers"
    DESCRIPTION "These are needed to do any development"
     INSTALL_TYPES full developer)

cpack_add_component(tests
    DISPLAY_NAME "ROOT Tests and Tutorials"
    DESCRIPTION "These are needed to do any test and tutorial"
     INSTALL_TYPES full developer)
