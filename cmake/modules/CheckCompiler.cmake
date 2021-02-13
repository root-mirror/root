# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#---------------------------------------------------------------------------------------------------
#  CheckCompiler.cmake
#---------------------------------------------------------------------------------------------------

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "(Apple|)Clang|GNU|Intel|MSVC")
  message(WARNING "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}.")
endif()

if(NOT GENERATOR_IS_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
  if(NOT CMAKE_C_FLAGS AND NOT CMAKE_CXX_FLAGS AND NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Specifies the build type on single-configuration generators" FORCE)
  endif()
endif()

include(CheckLanguage)
#---Enable FORTRAN (unfortunatelly is not not possible in all cases)-------------------------------
if(fortran)
  #--Work-around for CMake issue 0009220
  if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
    set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
  endif()
  if(CMAKE_Fortran_COMPILER)
    # CMAKE_Fortran_COMPILER has already been defined somewhere else, so
    # just check whether it contains a valid compiler
    enable_language(Fortran)
  else()
    # CMAKE_Fortran_COMPILER has not been defined, so first check whether
    # there is a Fortran compiler at all
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
      # Fortran compiler found, however as 'check_language' was executed
      # in a separate process, the result might not be compatible with
      # the C++ compiler, so reset the variable, ...
      unset(CMAKE_Fortran_COMPILER CACHE)
      # ..., and enable Fortran again, this time prefering compilers
      # compatible to the C++ compiler
      enable_language(Fortran)
    endif()
  endif()
else()
  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif()

#----Test if clang setup works----------------------------------------------------------------------
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  exec_program(${CMAKE_CXX_COMPILER} ARGS "--version 2>&1 | grep version" OUTPUT_VARIABLE _clang_version_info)
  string(REGEX REPLACE "^.*[ ]version[ ]([0-9]+)\\.[0-9]+.*" "\\1" CLANG_MAJOR "${_clang_version_info}")
  string(REGEX REPLACE "^.*[ ]version[ ][0-9]+\\.([0-9]+).*" "\\1" CLANG_MINOR "${_clang_version_info}")

  if(CMAKE_GENERATOR STREQUAL "Ninja")
    # LLVM/Clang are automatically checking if we are in interactive terminal mode.
    # We use color output only for Ninja, because Ninja by default is buffering the output,
    # so Clang disables colors as it is sure whether the output goes to a file or to a terminal.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcolor-diagnostics")
  endif()
  if(ccache AND CCACHE_VERSION VERSION_LESS "3.2.0")
    # https://bugzilla.samba.org/show_bug.cgi?id=8118
    # Call to 'ccache clang' is triggering next warning (valid for ccache 3.1.x, fixed in 3.2):
    # "clang: warning: argument unused during compilation: '-c"
    # Adding -Qunused-arguments provides a workaround for the bug.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qunused-arguments")
  endif()
else()
  set(CLANG_MAJOR 0)
  set(CLANG_MINOR 0)
endif()

#---Obtain the major and minor version of the GNU compiler-------------------------------------------
if (CMAKE_COMPILER_IS_GNUCXX)
  string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR ${CMAKE_CXX_COMPILER_VERSION})
  string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR ${CMAKE_CXX_COMPILER_VERSION})
  string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH ${CMAKE_CXX_COMPILER_VERSION})

  if(GCC_PATCH MATCHES "\\.+")
    set(GCC_PATCH "")
  endif()
  if(GCC_MINOR MATCHES "\\.+")
    set(GCC_MINOR "")
  endif()
  if(GCC_MAJOR MATCHES "\\.+")
    set(GCC_MAJOR "")
  endif()
  message(STATUS "Found GCC. Major version ${GCC_MAJOR}, minor version ${GCC_MINOR}")
  if("${GCC_MAJOR}.${GCC_MINOR}" VERSION_GREATER_EQUAL 4.9
      AND CMAKE_GENERATOR STREQUAL "Ninja")
    # GCC checks automatically if we are in interactive terminal mode.
    # We use color output only for Ninja, because Ninja by default is buffering the output,
    # so Clang disables colors as it is sure whether the output goes to a file or to a terminal.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
  endif()
else()
  set(GCC_MAJOR 0)
  set(GCC_MINOR 0)
endif()

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

#---C++ standard----------------------------------------------------------------------

# We want to set the default value of CMAKE_CXX_STANDARD to the compiler default,
# so we check the value of __cplusplus.
# This default value can be overridden by specifying one at the prompt.
if (MSVC)
   set(CXX_STANDARD_STRING 2011)
else()
   execute_process(COMMAND echo __cplusplus
                   COMMAND ${CMAKE_CXX_COMPILER} -E -x c++ -
                   COMMAND tail -n1
                   OUTPUT_VARIABLE CXX_STANDARD_STRING
                   OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
# Lexicographically compare the value of __cplusplus (e.g. "201703L" for C++17) to figure out
# what standard CMAKE_CXX_COMPILER uses by default.
# The standard values that __cplusplus takes are listed e.g. at
# https://en.cppreference.com/w/cpp/preprocessor/replace#Predefined_macros
# but note that compilers might denote partial implementations of new standards (e.g. c++1z)
# with other non-standard values.
if (${CXX_STANDARD_STRING} STRGREATER "201703L")
   set(CXX_STANDARD_STRING 20 CACHE STRING "")
elseif(${CXX_STANDARD_STRING} STRGREATER "201402L")
   set(CXX_STANDARD_STRING 17 CACHE STRING "")
elseif(${CXX_STANDARD_STRING} STRGREATER "201103L")
   set(CXX_STANDARD_STRING 14 CACHE STRING "")
else()
   # We stick to C++11 as a minimum value
   set(CXX_STANDARD_STRING 11 CACHE STRING "")
endif()
set(CMAKE_CXX_STANDARD ${CXX_STANDARD_STRING} CACHE STRING "")
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "")

if(cxx11 OR cxx14 OR cxx17)
  message(DEPRECATION "Options cxx11/14/17 are deprecated. Please use CMAKE_CXX_STANDARD instead.")

  # for backward compatibility
  if(cxx17)
    set(CMAKE_CXX_STANDARD 17 CACHE STRING "" FORCE)
  elseif(cxx14)
    set(CMAKE_CXX_STANDARD 14 CACHE STRING "" FORCE)
  elseif(cxx11)
    set(CMAKE_CXX_STANDARD 11 CACHE STRING "" FORCE)
  endif()

  unset(cxx17 CACHE)
  unset(cxx14 CACHE)
  unset(cxx11 CACHE)
endif()

if(NOT CMAKE_CXX_STANDARD MATCHES "11|14|17")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

# needed by roottest, to be removed once roottest is fixed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION}")

#---Check for libcxx option------------------------------------------------------------
if(libcxx)
  CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" HAS_LIBCXX11)
  if(NOT HAS_LIBCXX11)
    message(STATUS "Current compiler does not suppport -stdlib=libc++ option. Switching OFF libcxx option")
    set(libcxx OFF CACHE BOOL "" FORCE)
  endif()
endif()

#---Need to locate thead libraries and options to set properly some compilation flags----------------
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
  set(CMAKE_THREAD_FLAG -pthread)
else()
  set(CMAKE_THREAD_FLAG)
endif()


#---Setup compiler-specific flags (warning etc)----------------------------------------------
if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  # AppleClang and Clang proper.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wc++11-narrowing -Wsign-compare -Wsometimes-uninitialized -Wconditional-uninitialized -Wheader-guard -Warray-bounds -Wcomment -Wtautological-compare -Wstrncat-size -Wloop-analysis -Wbool-conversion")
elseif(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-implicit-fallthrough -Wno-noexcept-type")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-fallthrough")
  endif()
endif()


#---Setup details depending on the major platform type----------------------------------------------
if(CMAKE_SYSTEM_NAME MATCHES Linux)
  include(SetUpLinux)
elseif(APPLE)
  include(SetUpMacOS)
elseif(WIN32)
  include(SetUpWindows)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_THREAD_FLAG}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_THREAD_FLAG}")

if(libcxx)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

if(gcctoolchain)
  CHECK_CXX_COMPILER_FLAG("--gcc-toolchain=${gcctoolchain}" HAS_GCCTOOLCHAIN)
  if(HAS_GCCTOOLCHAIN)
     set(CMAKE_CXX_FLAGS "--gcc-toolchain=${gcctoolchain} ${CMAKE_CXX_FLAGS}")
  endif()
endif()

if(gnuinstall)
  set(R__HAVE_CONFIG 1)
endif()

#---Check if we use the new libstdc++ CXX11 ABI-----------------------------------------------------
# Necessary to compile check_cxx_source_compiles this early
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
"
#include <string>
#if _GLIBCXX_USE_CXX11_ABI == 0
  #error NOCXX11
#endif
int main() {}
" GLIBCXX_USE_CXX11_ABI)

#---Print the final compiler flags--------------------------------------------------------------------
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
message(STATUS "ROOT Platform: ${ROOT_PLATFORM}")
message(STATUS "ROOT Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "ROOT Processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "ROOT Architecture: ${ROOT_ARCHITECTURE}")
message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}' (flags = '${CMAKE_CXX_FLAGS_${BUILD_TYPE}}')")
message(STATUS "Compiler Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE}}")
