# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# Find the native MySQL includes and library
#
#  MYSQL_INCLUDE_DIR - where to find mysql.h, etc.
#  MYSQL_LIBRARIES   - List of libraries when using MySQL.
#  MYSQL_FOUND       - True if MySQL found.

if(MYSQL_INCLUDE_DIR OR MYSQL_)
  # Already in cache, be silent
  SET(MYSQL_FIND_QUIETLY TRUE)
endif()

if(NOT WIN32)
  # Split into two find_program invocations to avoid user-specified
  # mariadb_config being overridden by system mysql_config.
  find_program(MYSQL_CONFIG_EXECUTABLE NAMES mysql_config mariadb_config
    PATH_SUFFIXES bin
    NO_DEFAULT_PATH
    HINTS ${MYSQL_DIR} $ENV{MYSQL_DIR}
    )
  # Will not overwrite if MYSQL_CONFIG_EXECUTABLE is already set.
  find_program(MYSQL_CONFIG_EXECUTABLE NAMES mysql_config mariadb_config
    PATH_SUFFIXES bin
    PATHS /usr/local /usr
    )
endif()

if(MYSQL_CONFIG_EXECUTABLE)
  execute_process(COMMAND ${MYSQL_CONFIG_EXECUTABLE} --cflags OUTPUT_VARIABLE MYSQL_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(MYSQL_CFLAGS)
  string( REGEX MATCHALL "-I[^;]+" MYSQL_INCLUDE_DIRS "${MYSQL_CFLAGS}" )
  string( REPLACE "-I" "" MYSQL_INCLUDE_DIRS "${MYSQL_INCLUDE_DIRS}")
  find_path(MYSQL_INCLUDE_DIR mysql.h PATHS ${MYSQL_INCLUDE_DIRS} NO_DEFAULT_PATH)
  execute_process(COMMAND ${MYSQL_CONFIG_EXECUTABLE} --libs OUTPUT_VARIABLE MYSQL_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  find_path(MYSQL_INCLUDE_DIR mysql.h
    PATH_SUFFIXES include/mysql include/mariadb
    HINTS ${MYSQL_DIR} $ENV{MYSQL_DIR}
    NO_DEFAULT_PATH
    )
  if (NOT MYSQL_INCLUDE_DIR)
    find_path(MYSQL_INCLUDE_DIR mysql.h
      PATH_SUFFIXES include/mysql include/mariadb
      PATHS /usr/local/mysql /usr/local /usr
      )
  endif()
  set(MYSQL_NAMES mysqlclient mysqlclient_r)
  find_library(MYSQL_LIBRARY NAMES ${MYSQL_NAMES}
    PATH_SUFFIXES lib lib/mariadb lib/mysql
    lib/opt lib/opt/mariadb lib/opt/mysql
    HINTS ${MYSQL_DIR} $ENV{MYSQL_DIR}
    PATHS /usr/local/mysql /usr/local /usr
  )
  set(MYSQL_LIBRARIES ${MYSQL_LIBRARY})
endif()

# handle the QUIETLY and REQUIRED arguments and set DCAP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MySQL DEFAULT_MSG MYSQL_INCLUDE_DIR MYSQL_LIBRARIES)

mark_as_advanced(
  MYSQL_CONFIG_EXECUTABLE
  MYSQL_LIBRARY
  MYSQL_INCLUDE_DIR
)
