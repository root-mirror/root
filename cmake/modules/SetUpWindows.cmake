# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(ROOT_PLATFORM win32)

#----Check the compiler that is used-----------------------------------------------------
if(CMAKE_COMPILER_IS_GNUCXX)

  set(ROOT_ARCHITECTURE win32gcc)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe  -Wall -W -Woverloaded-virtual")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -Wall -W")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"           CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2"              CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG          "-g"               CACHE STRING "Flags for a debug build" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g"           CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_C_FLAGS_RELEASE          "-O2"              CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_C_FLAGS_DEBUG            "-g"               CACHE STRING "Flags for a debug build" FORCE)
elseif(MSVC)
  set(ROOT_ARCHITECTURE win32)

  math(EXPR VC_MAJOR "${MSVC_VERSION} / 100")
  math(EXPR VC_MINOR "${MSVC_VERSION} % 100")

  if(winrtdebug)
    set(BLDCXXFLAGS "-Zc:__cplusplus -MDd -GR")
    set(BLDCFLAGS   "-MDd")
  else()
    set(BLDCXXFLAGS "-Zc:__cplusplus -MD -GR")
    set(BLDCFLAGS   "-MD")
  endif()

  if(CMAKE_PROJECT_NAME STREQUAL ROOT)
    set(CMAKE_CXX_FLAGS "-nologo -I${CMAKE_SOURCE_DIR}/build/win -FIw32pragma.h -FIsehmap.h ${BLDCXXFLAGS} -EHsc- -W3 -wd4141 -wd4291 -wd4244 -wd4049 -D_WIN32 -D_XKEYCHECK_H -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS")
    set(CMAKE_C_FLAGS   "-nologo -I${CMAKE_SOURCE_DIR}/build/win -FIw32pragma.h -FIsehmap.h ${BLDCFLAGS} -EHsc- -W3 -D_WIN32 -DNOMINMAX")
    if(win_broken_tests)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DR__ENABLE_BROKEN_WIN_TESTS")
    endif()
  else()
    set(CMAKE_CXX_FLAGS "-nologo -FIw32pragma.h -FIsehmap.h ${BLDCXXFLAGS} -EHsc- -W3 -wd4244 -D_WIN32")
    set(CMAKE_C_FLAGS   "-nologo -FIw32pragma.h -FIsehmap.h ${BLDCFLAGS} -EHsc- -W3 -D_WIN32")
  endif()

  #---Select compiler flags----------------------------------------------------------------
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -Z7")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
  set(CMAKE_CXX_FLAGS_DEBUG          "-Od -Z7")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -Z7")
  set(CMAKE_C_FLAGS_RELEASE          "-O2")
  set(CMAKE_C_FLAGS_DEBUG            "-Od -Z7")

  #---Set Linker flags----------------------------------------------------------------------
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ignore:4049,4206,4217,4221 -incremental:no")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -ignore:4049,4206,4217,4221 -incremental:no")

  string(TIMESTAMP CURRENT_YEAR "%Y")
  set(ROOT_RC_SCRIPT ${CMAKE_BINARY_DIR}/etc/root.rc)

  foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
  endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
else()
  message(FATAL_ERROR "There is no setup for compiler '${CMAKE_CXX_COMPILER}' on this Windows system up to now. Stop cmake at this point.")
endif()
