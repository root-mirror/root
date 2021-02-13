# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(ROOT_PLATFORM linux)

if(CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
  if(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
    set(ROOT_ARCHITECTURE linuxx8664icc)
  else()
    set(ROOT_ARCHITECTURE linuxx8664gcc)
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES i686)
  set(FP_MATH_FLAGS "-msse2 -mfpmath=sse")
  if(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
    set(ROOT_ARCHITECTURE linuxicc)
  else()
    set(ROOT_ARCHITECTURE linux)
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES aarch64)
  set(ROOT_ARCHITECTURE linuxarm64)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES arm)
  set(ROOT_ARCHITECTURE linuxarm)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES ppc64)
  set(ROOT_ARCHITECTURE linuxppc64gcc)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES s390x)
  set(ROOT_ARCHITECTURE linuxs390xgcc)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES s390)
  set(ROOT_ARCHITECTURE linuxs390gcc)
else()
  message(FATAL_ERROR "Unknown processor: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

# JIT must be able to resolve symbols from all ROOT binaries.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")

# Set developer flags
if(dev)
  # Warnings are errors.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

  # Do not relink just because a dependent .so has changed.
  # I.e. relink only if a header included by the libs .o-s has changed,
  # whether or not that header "belongs" to a different .so.
  set(CMAKE_LINK_DEPENDS_NO_SHARED On)

  # Split debug info for faster builds.
  if(NOT gnuinstall)
    # We won't install DWARF files.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gsplit-dwarf")
  endif()

  # Try faster linkers, prefer lld then gold.
  execute_process(COMMAND ${CMAKE_C_COMPILER} -fuse-ld=lld -Wl,--version OUTPUT_VARIABLE stdout ERROR_QUIET)
  if("${stdout}" MATCHES "LLD ")
    set(SUPERIOR_LINKER "lld")
  else()
    execute_process(COMMAND ${CMAKE_C_COMPILER} -fuse-ld=gold -Wl,--version OUTPUT_VARIABLE stdout ERROR_QUIET)
    if ("${stdout}" MATCHES "GNU gold")
      set(SUPERIOR_LINKER "gold")
    endif()
  endif()
  # Only lld and gold support --gdb-index
  if(SUPERIOR_LINKER)
    set(LLVM_USE_LINKER "${SUPERIOR_LINKER}")
    if(CMAKE_BUILD_TYPE MATCHES "Deb")
      message(STATUS "Using ${SUPERIOR_LINKER} linker with gdb-index")
      set(GDBINDEX "-Wl,--gdb-index")
    else()
      message(STATUS "Using ${SUPERIOR_LINKER} linker")
    endif()
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=${SUPERIOR_LINKER} ${GDBINDEX}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=${SUPERIOR_LINKER} ${GDBINDEX}")
    set(LLVM_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=${SUPERIOR_LINKER} ${GDBINDEX}")
    set(LLVM_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=${SUPERIOR_LINKER} ${GDBINDEX}")
  endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe ${FP_MATH_FLAGS} -Wshadow -Wall -W -Woverloaded-virtual -fsigned-char")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -Wall -W")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined -Wl,--hash-style=\"both\"")

  if(asan)
    # See also core/sanitizer/README.md for what's happening.
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} --print-file-name=libclang_rt.asan-x86_64.so OUTPUT_VARIABLE ASAN_RUNTIME_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(ASAN_EXTRA_CXX_FLAGS -fsanitize=address -fno-omit-frame-pointer -fsanitize-recover=address)
    set(ASAN_EXTRA_SHARED_LINKER_FLAGS "-fsanitize=address -z undefs")
    set(ASAN_EXTRA_EXE_LINKER_FLAGS "-fsanitize=address -z undefs -Wl,--undefined=__asan_default_options -Wl,--undefined=__lsan_default_options -Wl,--undefined=__lsan_default_suppressions")
  endif()

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g"           CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE        "-O3"              CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG          "-g"               CACHE STRING "Flags for a debug build" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3 -g"           CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_C_FLAGS_RELEASE          "-O3"              CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_C_FLAGS_DEBUG            "-g"               CACHE STRING "Flags for a debug build" FORCE)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe ${FP_MATH_FLAGS} -Wall -W -Woverloaded-virtual -fsigned-char")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -Wall -W")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  endif()

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

  if(asan)
    # See also core/sanitizer/README.md for what's happening.
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} --print-file-name=libclang_rt.asan-x86_64.so OUTPUT_VARIABLE ASAN_RUNTIME_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(ASAN_EXTRA_CXX_FLAGS -fsanitize=address -fno-omit-frame-pointer -fsanitize-address-use-after-scope -fsanitize-blacklist=${CMAKE_SOURCE_DIR}/build/ASan_blacklist.txt)
    set(ASAN_EXTRA_SHARED_LINKER_FLAGS "-fsanitize=address -static-libsan -z undefs")
    set(ASAN_EXTRA_EXE_LINKER_FLAGS "-fsanitize=address -static-libsan -z undefs -Wl,--undefined=__asan_default_options -Wl,--undefined=__lsan_default_options -Wl,--undefined=__lsan_default_suppressions")
  endif()

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"                    CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2"                       CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG          "-g"                        CACHE STRING "Flags for a debug build" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g"                    CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_C_FLAGS_RELEASE          "-O2"                       CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_C_FLAGS_DEBUG            "-g"                        CACHE STRING "Flags for a debug build" FORCE)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1476")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -restrict")

  # Check icc compiler version and set compile flags according to the
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v
                  ERROR_VARIABLE _icc_version_info ERROR_STRIP_TRAILING_WHITESPACE)

  string(REGEX REPLACE "(^V|^icc[ ]v|^icpc[ ]v)ersion[ ]([0-9]+)\\.[0-9]+.*" "\\2" ICC_MAJOR "${_icc_version_info}")
  string(REGEX REPLACE "(^V|^icc[ ]v|^icpc[ ]v)ersion[ ][0-9]+\\.([0-9]+).*" "\\2" ICC_MINOR "${_icc_version_info}")

  if(ICC_MAJOR GREATER 9 OR ICC_MAJOR EQUAL 9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1572")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd1572")
  endif()

  if(ICC_MAJOR GREATER 11 OR ICC_MAJOR EQUAL 11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd279")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd279")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
  endif()

  if(ICC_MAJOR GREATER 14 OR ICC_MAJOR EQUAL 14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd873 -wd2536")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd873 -wd2536")
  endif()

  if(ICC_MAJOR GREATER 15 OR ICC_MAJOR EQUAL 15)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -wd1098 -wd1292 -wd1478 -wd3373")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd597 -wd1098 -wd1292 -wd1478 -wd3373")
  endif()

  # Select flags.
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -fp-model precise -g"          CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -fp-model precise"             CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g"                            CACHE STRING "Flags for a debug build" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -fp-model precise -g"          CACHE STRING "Flags for release build with debug info" FORCE)
  set(CMAKE_C_FLAGS_RELEASE          "-O2 -fp-model precise"             CACHE STRING "Flags for release build" FORCE)
  set(CMAKE_C_FLAGS_DEBUG            "-O0 -g"                            CACHE STRING "Flags for a debug build" FORCE)
endif()
