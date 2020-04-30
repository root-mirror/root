# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#---Check for Python installation-------------------------------------------------------

message(STATUS "Looking for python")

if(pyroot_experimental)
  unset(PYTHON_INCLUDE_DIR CACHE)
  unset(PYTHON_LIBRARY CACHE)
endif()

# Python is required by header and manpage generation

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14 AND NOT MSVC)

  # Determine whether we should prefer Python 2 or Python 3:
  set(PYTHON_PREFER_VERSION "3")
  # Check whether old `find_package(PythonInterp)` variable was passed.
  # If so, it will be passed to find_package(Python) below. Otherwise,
  # check what `python` points to: Python 2 or 3:
  if(NOT PYTHON_EXECUTABLE)
    find_program(PYTHON_EXECUTABLE "python")
  endif()
  if(PYTHON_EXECUTABLE)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version_info[0])"
                    OUTPUT_VARIABLE PYTHON_PREFER_VERSION
                    ERROR_VARIABLE PYTHON_PREFER_VERSION_ERR)
    if(PYTHON_PREFER_VERSION_ERR)
      message(WARNING "Unable to determine version of ${PYTHON_EXECUTABLE}: ${PYTHON_PREFER_VERSION_ERR}")
    endif()
    string(STRIP "${PYTHON_PREFER_VERSION}" PYTHON_PREFER_VERSION)
  endif()

  message(STATUS "Preferring Python version ${PYTHON_PREFER_VERSION}")

  if("${PYTHON_PREFER_VERSION}" MATCHES "2")
    # Means PYTHON_EXECUTABLE wasn't defined.
    if(PYTHON_INCLUDE_DIRS AND NOT Python2_INCLUDE_DIRS)
      set(Python2_INCLUDE_DIRS "${PYTHON_INCLUDE_DIRS}")
    endif()
    if(PYTHON_LIBRARIES AND NOT Python2_LIBRARIES)
      set(Python2_LIBRARIES "${PYTHON_LIBRARIES}")
    endif()
    find_package(Python2 COMPONENTS Interpreter REQUIRED)
    # Search for NumPy and Development, but not required:
    find_package(Python2 COMPONENTS Development NumPy QUIET)
    # Compat with find_package(PythonInterp), find_package(PythonLibs)
    set(PYTHON_EXECUTABLE "${Python2_EXECUTABLE}" CACHE INTERNAL "" FORCE)
    set(PYTHON_INCLUDE_DIRS "${Python2_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
    set(PYTHON_LIBRARIES "${Python2_LIBRARIES}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_STRING "${Python2_VERSION}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MAJOR "${Python2_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MINOR "${Python2_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
    set(NUMPY_FOUND ${Python2_NumPy_FOUND} CACHE INTERNAL "" FORCE)
    set(NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
  else()
    if(PYTHON_EXECUTABLE AND NOT Python3_EXECUTABLE)
      set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
    endif()
    if(PYTHON_INCLUDE_DIRS AND NOT Python3_INCLUDE_DIRS)
      set(Python3_INCLUDE_DIRS "${PYTHON_INCLUDE_DIRS}")
    endif()
    if(PYTHON_LIBRARIES AND NOT Python3_LIBRARIES)
      set(Python3_LIBRARIES "${PYTHON_LIBRARIES}")
    endif()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    # Search for NumPy and Development, but not required:
    find_package(Python3 COMPONENTS Development NumPy QUIET)
    # Compat with find_package(PythonInterp), find_package(PythonLibs), find_package(NumPy)
    set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE INTERNAL "" FORCE)
    set(PYTHON_INCLUDE_DIRS "${Python3_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
    set(PYTHON_LIBRARIES "${Python3_LIBRARIES}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_STRING "${Python3_VERSION}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MAJOR "${Python3_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MINOR "${Python3_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
    set(NUMPY_FOUND ${Python3_NumPy_FOUND} CACHE INTERNAL "" FORCE)
    set(NUMPY_INCLUDE_DIRS "${Python3_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
  endif()

else()
  find_package(PythonInterp ${python_version} REQUIRED)

  find_package(PythonLibs ${python_version})

  if(NOT "${PYTHONLIBS_VERSION_STRING}" MATCHES "${PYTHON_VERSION_STRING}")
    message(FATAL_ERROR "Version mismatch between Python interpreter (${PYTHON_VERSION_STRING})"
      " and libraries (${PYTHONLIBS_VERSION_STRING}).\nROOT cannot work with this configuration. "
      "Please specify only PYTHON_EXECUTABLE to CMake with an absolute path to ensure matching versions are found.")
  endif()

  find_package(NumPy QUIET)
endif()
