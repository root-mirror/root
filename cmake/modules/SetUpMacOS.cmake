# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(ROOT_ARCHITECTURE macosx)
set(ROOT_PLATFORM macosx)

# https://gitlab.kitware.com/cmake/cmake/issues/19222
if(CMAKE_VERSION VERSION_LESS 3.14.4)
  if(CMAKE_GENERATOR STREQUAL "Ninja")
    find_program(NINJAPROG ninja DOC "looking for Ninja to be used")
    if (NINJAPROG)
      execute_process(COMMAND ${NINJAPROG} --version
        OUTPUT_VARIABLE NINJAVERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE NINJAVERSIONRES)
      if(${NINJAVERSIONRES} EQ 0)
        if(${NINJAVERSION} VERSION_GREATER 1.8.2)
          message(FATAL_ERROR "You have hit https://gitlab.kitware.com/cmake/cmake/issues/19222\n"
            "Your build will be indeterministic, i.e. unreliable for incremental builds."
            "To fix this, please install CMake >= 3.14.4!")
        endif()
      endif()
    endif()
  endif()
endif()

if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  EXECUTE_PROCESS(COMMAND sw_vers "-productVersion"
                  COMMAND cut -d . -f 1-2
                  OUTPUT_VARIABLE MACOSX_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

  MESSAGE(STATUS "Found a macOS system ${MACOSX_VERSION}")

  if(MACOSX_VERSION VERSION_GREATER 10.7 AND ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    set(libcxx ON CACHE BOOL "Build using libc++" FORCE)
  endif()

  if(MACOSX_VERSION VERSION_GREATER 10.4)
    #TODO: check haveconfig and rpath -> set rpath true
    #TODO: check Thread, define link command
    #TODO: more stuff check configure script
    if(CMAKE_SYSTEM_PROCESSOR MATCHES 64)
       if(CMAKE_SYSTEM_PROCESSOR MATCHES arm64)
          MESSAGE(STATUS "Found an AArch64 system")
          set(ROOT_ARCHITECTURE macosxarm64)
       else()
          MESSAGE(STATUS "Found an x86_64 system")
          set(ROOT_ARCHITECTURE macosx64)
       endif()

       SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -m64")
       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -m64")
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
       SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
       SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
    else()
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
       SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
       SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32")
    endif()
  endif()

  if(MACOSX_VERSION VERSION_GREATER 10.6)
    set(MACOSX_SSL_DEPRECATED ON)
  endif()
  if(MACOSX_VERSION VERSION_GREATER 10.7)
    set(MACOSX_ODBC_DEPRECATED ON)
  endif()
  if(MACOSX_VERSION VERSION_GREATER 10.8)
    set(MACOSX_GLU_DEPRECATED ON)
  endif()

  if (CMAKE_COMPILER_IS_GNUCXX)
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -W -Wshadow -Wall -Woverloaded-virtual -fsigned-char -fno-common")
     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -W -Wall -fsigned-char -fno-common")
     SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

     SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -single_module -Wl,-dead_strip_dylibs")
     SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -single_module -Wl,-dead_strip_dylibs")

     set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -bind_at_load -m64")
     set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")

     # Select flags.
     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"          CACHE STRING "Flags for a release build with debug info" FORCE)
     set(CMAKE_CXX_FLAGS_RELEASE        "-O2"             CACHE STRING "Flags for a release build" FORCE)
     set(CMAKE_CXX_FLAGS_DEBUG          "-g"              CACHE STRING "Flags for a debug build" FORCE)
     set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g"          CACHE STRING "Flags for a release build with debug info" FORCE)
     set(CMAKE_C_FLAGS_RELEASE          "-O2"             CACHE STRING "Flags for a release build" FORCE)
     set(CMAKE_C_FLAGS_DEBUG            "-g"              CACHE STRING "Flags for a debug build" FORCE)
  elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -W -Wall -Woverloaded-virtual -fsigned-char -fno-common -Qunused-arguments")
     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -W -Wall -fsigned-char -fno-common -Qunused-arguments")
     if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
     endif()

     SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

     SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -single_module -Wl,-dead_strip_dylibs")
     SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -single_module -Wl,-dead_strip_dylibs")

     set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -bind_at_load -m64")
     set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")

     if(asan)
       # See also core/sanitizer/README.md for what's happening.

       #This should be the right way to do it, but clang 10 seems to have a bug
       #execute_process(COMMAND ${CMAKE_CXX_COMPILER} --print-file-name=libclang_rt.asan_osx_dynamic.dylib OUTPUT_VARIABLE ASAN_RUNTIME_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
       execute_process(COMMAND mdfind -name libclang_rt.asan_osx_dynamic.dylib OUTPUT_VARIABLE ASAN_RUNTIME_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
       set(ASAN_EXTRA_CXX_FLAGS -fsanitize=address -fno-omit-frame-pointer -fsanitize-address-use-after-scope -fsanitize-blacklist=${CMAKE_SOURCE_DIR}/build/ASan_blacklist.txt)
       set(ASAN_EXTRA_SHARED_LINKER_FLAGS "-fsanitize=address -static-libsan")
       set(ASAN_EXTRA_EXE_LINKER_FLAGS "-fsanitize=address -static-libsan -Wl,-u,___asan_default_options -Wl,-u,___lsan_default_options -Wl,-u,___lsan_default_suppressions")
     endif()

     # Select flags.
     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"          CACHE STRING "Flags for a release build with debug info" FORCE)
     set(CMAKE_CXX_FLAGS_RELEASE        "-O2"             CACHE STRING "Flags for a release build" FORCE)
     set(CMAKE_CXX_FLAGS_DEBUG          "-g"              CACHE STRING "Flags for a debug build" FORCE)
     set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g"          CACHE STRING "Flags for a release build with debug info" FORCE)
     set(CMAKE_C_FLAGS_RELEASE          "-O2"             CACHE STRING "Flags for a release build" FORCE)
     set(CMAKE_C_FLAGS_DEBUG            "-g"              CACHE STRING "Flags for a debug build" FORCE)
  else()
    MESSAGE(FATAL_ERROR "There is no setup for this compiler with ID=${CMAKE_CXX_COMPILER_ID} up to now. Don't know what to do. Stop cmake at this point.")
  endif()

  #---Set Linker flags----------------------------------------------------------------------
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=${MACOSX_VERSION}")
else (CMAKE_SYSTEM_NAME MATCHES Darwin)
  MESSAGE(FATAL_ERROR "There is no setup for this this Apple system up to now. Don't know waht to do. Stop cmake at this point.")
endif (CMAKE_SYSTEM_NAME MATCHES Darwin)

#---Avoid puting the libraires and executables in different configuration locations
if(CMAKE_GENERATOR MATCHES Xcode)
  foreach( _conf ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${_conf} _conf )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_conf} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_conf} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_conf} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
  endforeach()
endif()
