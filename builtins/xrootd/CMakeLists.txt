# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

include(ExternalProject)

find_package(OpenSSL REQUIRED)

set(XROOTD_VERSION "5.1.0")
set(XROOTD_VERSIONNUM 500010000 CACHE INTERNAL "" FORCE)
set(XROOTD_SRC_URI https://github.com/xrootd/xrootd/archive/v${XROOTD_VERSION}.tar.gz)
set(XROOTD_DESTDIR ${CMAKE_BINARY_DIR}/XROOTD-prefix)
set(XROOTD_ROOTDIR ${XROOTD_DESTDIR})
message(STATUS "Downloading and building XROOTD version ${XROOTD_VERSION}")

  # Guess under which directory XRootD will install its libraires
  set(XROOTD_LIBDIR "lib")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND ${CMAKE_SIZEOF_VOID_P} EQUAL 8
     AND NOT CMAKE_CROSSCOMPILING AND NOT EXISTS "/etc/debian_version")
    set(XROOTD_LIBDIR "lib64")
  endif()

  set(XROOTD_LIBRARIES ${XROOTD_ROOTDIR}/${XROOTD_LIBDIR}/libXrdUtils${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${XROOTD_ROOTDIR}/${XROOTD_LIBDIR}/libXrdClient${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${XROOTD_ROOTDIR}/${XROOTD_LIBDIR}/libXrdCl${CMAKE_SHARED_LIBRARY_SUFFIX})
  ExternalProject_Add(
    XROOTD
    URL ${XROOTD_SRC_URI}
    URL_HASH SHA256=356473ac1652080440c1c2ec80dd8e749ee569003594c3b38869a38999299c5e
    INSTALL_DIR ${XROOTD_ROOTDIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DCMAKE_PREFIX_PATH:STRING=${OPENSSL_PREFIX}
               -DCMAKE_BUILD_TYPE=Release
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${ROOT_EXTERNAL_CXX_FLAGS}
               -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
               -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
               -DENABLE_PYTHON=OFF
               -DCMAKE_INSTALL_RPATH:STRING=${XROOTD_ROOTDIR}/${XROOTD_LIBDIR}
    INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
            COMMAND ${CMAKE_COMMAND} -E copy_directory <INSTALL_DIR>/include/xrootd <INSTALL_DIR>/include
    LOG_DOWNLOAD 1 LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
    BUILD_BYPRODUCTS ${XROOTD_LIBRARIES}
    TIMEOUT 600
  )

if(builtin_openssl)
  add_dependencies(XROOTD OPENSSL)
endif()

list(APPEND XROOTD_LIBRARIES OpenSSL::SSL)

unset(XROOTD_FOUND CACHE)
unset(XROOTD_FOUND PARENT_SCOPE)
set(XROOTD_FOUND TRUE CACHE BOOL "" FORCE)

set(XROOTD_INCLUDE_DIRS ${XROOTD_ROOTDIR}/include/xrootd ${XROOTD_ROOTDIR}/include/xrootd/private)
set(XROOTD_NOMAIN TRUE)
set(XROOTD_CFLAGS "-DROOTXRDVERS=${XROOTD_VERSIONNUM}")

set(XROOTD_INCLUDE_DIR ${XROOTD_ROOTDIR}/include/xrootd ${XROOTD_ROOTDIR}/include/xrootd/private CACHE INTERNAL "" FORCE)
set(XROOTD_INCLUDE_DIRS ${XROOTD_ROOTDIR}/include/xrootd ${XROOTD_ROOTDIR}/include/xrootd/private CACHE INTERNAL "" FORCE)
set(XROOTD_LIBRARY ${XROOTD_PREFIX}/lib/${XROOTD_LIBNAME} CACHE INTERNAL "" FORCE)
set(XROOTD_LIBRARIES ${XROOTD_LIBRARIES} CACHE INTERNAL "" FORCE)

add_library(xrootd INTERFACE)
target_include_directories(xrootd INTERFACE $<BUILD_INTERFACE:${XROOTD_INCLUDE_DIR}>)
target_link_libraries(xrootd INTERFACE $<BUILD_INTERFACE:${XROOTD_LIBRARIES}>)
add_dependencies(xrootd XROOTD)

add_library(Xrootd::Xrootd ALIAS xrootd)

set_property(GLOBAL APPEND PROPERTY ROOT_BUILTIN_TARGETS XROOTD)

install(DIRECTORY ${XROOTD_ROOTDIR}/${XROOTD_LIBDIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries FILES_MATCHING PATTERN "libXrd*")
install(DIRECTORY ${XROOTD_ROOTDIR}/include/xrootd/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers)
if(APPLE)
  # XRootD libraries on mac need the LC_RPATH variable set. The build process already takes care of setting
  #   * BUILD_RPATH = build/XROOTD-prefix/../src
  #   * INSTALL_RPATH = build/lib
  # Since the install directory for the builtin_xrootd target corresponds to the build directory of the main project.
  # Use a post install script to change the LC_RPATH variable of the libraries in the ROOT install folder.
  install(SCRIPT ${CMAKE_CURRENT_LIST_DIR}/XROOTDApplePostInstall.cmake
          CODE "xrootd_libs_change_rpath(${XROOTD_ROOTDIR}/${XROOTD_LIBDIR} ${CMAKE_INSTALL_FULL_LIBDIR})"
  )
endif()
