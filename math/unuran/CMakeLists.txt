# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/unurun package
############################################################################

#---Define package related variables-----------------------------------------------------------------

if(builtin_unuran)

set(UNR_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(UNR_VERSION "1.8.0-root")
set(UNR_TARNAME "unuran-${UNR_VERSION}")
set(UNR_TARGZFILE ${UNR_SRCDIR}/${UNR_TARNAME}.tar.gz)
set(UNR_TARFILE ${UNR_SRCDIR}/${UNR_TARNAME}.tar)
set(UNR_UNTARDIR ${CMAKE_CURRENT_BINARY_DIR}/${UNR_TARNAME})

#---Untar sources at configuration/generation time (needed for listing sources)
if(NOT EXISTS ${UNR_UNTARDIR})
  execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${UNR_TARGZFILE}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
  # This is necessary to replace the config.guess of unuran 1.8.0 as it does not allow to compile the 
  # package on arm64.
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/config.guess_patch1 ${UNR_UNTARDIR}/autoconf/config.guess
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
endif()

if(WIN32)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.win.in  ${UNR_UNTARDIR}/config.h)
else()
  #---Define special compiler settings for unurun-----------------------------------------------------
  set(UNR_CC ${CMAKE_C_COMPILER})
  if(ROOT_ARCHITECTURE MATCHES hpuxia64acc)
    set(UNR_CC "${UNR_CC} +DD64 -Ae")
  elseif(ROOT_ARCHITECTURE MATCHES linuxppc64gcc)
    set(UNR_CC "${UNR_CC} -m64 -fPIC")
  elseif(ROOT_ARCHITECTURE MATCHES linuxx8664gcc)
    set(UNR_CFLAGS "-m64 -fPIC")
  elseif(ROOT_ARCHITECTURE MATCHES linuxicc)
    set(UNR_CFLAGS "-m32")
  elseif(ROOT_ARCHITECTURE MATCHES linuxx8664icc)
    set(UNR_CFLAGS "-m64")
  elseif(ROOT_ARCHITECTURE MATCHES win32)
    set(UNR_CFLAGS "-MD -G5 -GX")
  endif()
  if(CMAKE_OSX_SYSROOT)
    set(UNR_CFLAGS "${UNR_CFLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")
  endif()

  #---configure unuran (required for creating the config.h used by unuran source files)----------------
  add_custom_command(OUTPUT ${UNR_UNTARDIR}/config.h
                     COMMAND GNUMAKE=make ./configure CC=${UNR_CC} CFLAGS=${UNR_CFLAGS} > /dev/null 2>& 1
                     WORKING_DIRECTORY ${UNR_UNTARDIR})
endif()

#---We need to disable some warnings-------------------------------------------------------------------
string(REPLACE -Wall "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  ROOT_ADD_C_FLAG(CMAKE_C_FLAGS -Wno-parentheses-equality)
endif()
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
  ROOT_ADD_C_FLAG(CMAKE_C_FLAGS -Wno-maybe-uninitialized)
endif()

set(unrsources ${UNR_UNTARDIR}/src/utils/*.c
               ${UNR_UNTARDIR}/src/methods/*.c
               ${UNR_UNTARDIR}/src/specfunct/*.c
               ${UNR_UNTARDIR}/src/distr/*.c
               ${UNR_UNTARDIR}/src/distributions/*.c
               ${UNR_UNTARDIR}/src/parser/*.c
               ${UNR_UNTARDIR}/src/tests/*.c
               ${UNR_UNTARDIR}/src/uniform/*.c
               ${UNR_UNTARDIR}/src/urng/*.c )
set(unrconfig ${UNR_UNTARDIR}/config.h)

else()

set(unrsources)
set(unrconfig)

endif(builtin_unuran)


ROOT_STANDARD_LIBRARY_PACKAGE(Unuran
  HEADERS
    TUnuran.h
    TUnuranBaseDist.h
    TUnuranContDist.h
    TUnuranDiscrDist.h
    TUnuranEmpDist.h
    TUnuranMultiContDist.h
    TUnuranSampler.h
  SOURCES
    src/TUnuran.cxx
    src/TUnuranContDist.cxx
    src/TUnuranDiscrDist.cxx
    src/TUnuranEmpDist.cxx
    src/TUnuranMultiContDist.cxx
    src/TUnuranSampler.cxx
    ${unrconfig}
    ${unrsources}
  LIBRARIES
    ${UNURAN_LIBRARIES}
  DEPENDENCIES
    Core
    Hist
    MathCore
)

if(builtin_unuran)
  target_include_directories(Unuran SYSTEM BEFORE PRIVATE ${UNR_UNTARDIR} ${UNR_UNTARDIR}/src ${UNR_UNTARDIR}/src/utils)
  target_compile_definitions(Unuran PRIVATE HAVE_CONFIG_H)
else()
  target_include_directories(Unuran SYSTEM PRIVATE ${UNURAN_INCLUDE_DIRS})
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
