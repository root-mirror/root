# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT net/rpdutils package
# @author Pere Mato, CERN
############################################################################

include_directories(res ../auth/res)

if(MACOSX_SSL_DEPRECATED)
  set_source_files_properties(src/rpdutils.cxx COMPILE_FLAGS "-Wno-deprecated-declarations")
endif()

if(shadowpw)
  add_definitions(-DR__SHADOWPW)
endif()

if(ssl)
  include_directories(${OPENSSL_INCLUDE_DIR})
endif()

add_library(rpdutil STATIC
  src/daemon.cxx
  src/error.cxx
  src/net.cxx
  src/netpar.cxx
  src/rpdutils.cxx
)

target_include_directories(rpdutil PRIVATE
   ${CMAKE_SOURCE_DIR}/core/clib/inc         # for snprintf.h
   ${CMAKE_SOURCE_DIR}/core/base/inc         # for RConfig.hxx
   ${CMAKE_SOURCE_DIR}/core/foundation/inc   # for TError.hxx
   ${CMAKE_SOURCE_DIR}/core/meta/inc         # for TGenericClassInfo.h
   ${CMAKE_SOURCE_DIR}/net/net/inc           # for NetErrors.h
   ${CMAKE_SOURCE_DIR}/net/auth/inc          # for NetErrors.h
   ${CMAKE_BINARY_DIR}/ginclude              # for RConfigure.h
)

target_link_libraries(rpdutil INTERFACE $<$<PLATFORM_ID:Linux>:crypt>)
set_property(TARGET rpdutil PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(rpdutil move_headers)

ROOT_LINKER_LIBRARY(SrvAuth rpdutils.cxx
  LIBRARIES
    ${OPENSSL_LIBRARIES}
    rpdutil
    rsa
  DEPENDENCIES
    Net
  BUILTINS
    OPENSSL
)

target_include_directories(SrvAuth PRIVATE
   ${CMAKE_SOURCE_DIR}/net/auth/res
   ${CMAKE_SOURCE_DIR}/net/auth/inc
   ${CMAKE_SOURCE_DIR}/rpdutils/res
)
