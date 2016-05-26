# Find the LZ4 includes and library.
#
# This module defines
# LZ4_INCLUDE_DIR, where to locate LZ4 header files
# LZ4_LIBRARIES, the libraries to link against to use LZ4
# LZ4_FOUND.  If false, you cannot build anything that requires LZ4.
#
# inspired by @zzxuanyuan's pull request https://github.com/root-mirror/root/pull/81

if(LZ4_CONFIG_EXECUTABLE)
  set(LZ4_FIND_QUIETLY 1)
endif()
set(LZ4_FOUND 0)

if(NOT LZ4_DIR)
  set(LZ4_DIR $ENV{LZ4_DIR})
endif()

find_path(LZ4_INCLUDE_DIR lz4.h PATHS
  ${LZ4_DIR}/include
  /usr/include
  /usr/local/include
  /opt/lz4/include
  NO_DEFAULT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  DOC "Specify the directory containing lz4.h"
)

find_library(LZ4_LIBRARY NAMES lz4 PATHS
  ${LZ4_DIR}/lib
  /usr/local/lz4/lib
  /usr/local/lib
  /usr/lib/lz4
  /usr/local/lib/lz4
  /usr/lz4/lib /usr/lib
  /usr/lz4 /usr/local/lz4
  /opt/lz4 /opt/lz4/lib
  NO_DEFAULT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  DOC "Specify the lz4 library here."
)

if(LZ4_INCLUDE_DIR)
  message(STATUS "Found LZ4 includes at ${LZ4_INCLUDE_DIR}")
else()
  message(STATUS "LZ4 includes not found")
  message(STATUS "not even in $ENV{LZ4_DIR}")
endif()

if(LZ4_LIBRARY)
  message(STATUS "Found LZ4 library at ${LZ4_LIBRARY}")
else()
  message(STATUS "LZ4 library not found")
endif()


if(LZ4_INCLUDE_DIR AND LZ4_LIBRARY)
  set(LZ4_FOUND 1)
endif()

set(LZ4_LIBRARIES ${LZ4_LIBRARY})
mark_as_advanced(LZ4_FOUND LZ4_LIBRARY LZ4_INCLUDE_DIR)
