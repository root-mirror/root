# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

if(pyroot)
  if ((MSVC AND NOT CMAKE_INSTALL_BINDIR STREQUAL CMAKE_INSTALL_PYTHONDIR) OR
      (NOT CMAKE_INSTALL_LIBDIR STREQUAL CMAKE_INSTALL_PYTHONDIR AND (how_many_pythons GREATER 1)))
    message(STATUS "CMAKE_INSTALL_BINDIR = ${CMAKE_INSTALL_BINDIR}; CMAKE_INSTALL_PYTHONDIR = ${CMAKE_INSTALL_PYTHONDIR} ")
    message(FATAL_ERROR "Using a custom CMAKE_INSTALL_PYTHONDIR is not supported for multi-python builds or on Windows. Please check your build configuration.")
  endif()

  if(pyroot_legacy)
    add_subdirectory(pyroot_legacy)
  else()
    # Machinery necessary to loop in parallel over multiple lists
    # This will be used in the CMakeLists of pythonizations,
    # cppyy, cppyy-backend and CPyCppyy to build and install targets
    # for multiple Python versions (only 2 and 3 for now)
    list(LENGTH python_executables len)
    math(EXPR how_many_pythons "${len} - 1")

    add_subdirectory(pyroot)
    add_subdirectory(tpython)
    add_subdirectory(jupyroot)
    add_subdirectory(jsmva)
  endif()
endif()

if(r)
  add_subdirectory(r)
endif()
