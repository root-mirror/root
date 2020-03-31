# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#---Check for Python installation-------------------------------------------------------

message(STATUS "Looking for Python")

if(pyroot_experimental)
  unset(PYTHON_INCLUDE_DIR CACHE)
  unset(PYTHON_LIBRARY CACHE)
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  # - Check if the PYTHON_EXECUTABLE deprecated variable was passed by the
  # user; if so, check weather it points to Python 2 or 3 and set the
  # appropriate Python{X}_EXECUTABLE variable
  # - Look for Python3 and set the deprecated variables to the ones set
  # automatically by find_package(Python3 ...)
  # - Look for Python2 and set the deprecated variables to the ones set
  # automatically by find_package(Python2 ...) ONLY IF PYTHON3 WASN'T FOUND
  if(PYTHON_EXECUTABLE)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version_info[0])"
                    OUTPUT_VARIABLE PYTHON_PREFER_VERSION
                    ERROR_VARIABLE PYTHON_PREFER_VERSION_ERR)
    if(PYTHON_PREFER_VERSION_ERR)
      message(WARNING "Unable to determine version of ${PYTHON_EXECUTABLE}: ${PYTHON_PREFER_VERSION_ERR}")
    endif()
    string(STRIP "${PYTHON_PREFER_VERSION}" PYTHON_PREFER_VERSION)
    set(Python${PYTHON_PREFER_VERSION}_EXECUTABLE "${PYTHON_EXECUTABLE}")
  endif()

  find_package(Python3 COMPONENTS Interpreter Development NumPy)
  if(Python3_Development_FOUND)
    set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE INTERNAL "" FORCE)
    set(PYTHON_INCLUDE_DIRS "${Python3_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
    set(PYTHON_LIBRARIES "${Python3_LIBRARIES}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_STRING "${Python3_VERSION}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MAJOR "${Python3_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
    set(PYTHON_VERSION_MINOR "${Python3_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
    set(PYTHON_UNDER_VERSION_STRING "${Python3_VERSION_MAJOR}_${Python3_VERSION_MINOR}")
    set(NUMPY_FOUND ${Python3_NumPy_FOUND} CACHE INTERNAL "" FORCE)
    set(NUMPY_INCLUDE_DIRS "${Python3_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
  endif()

  find_package(Python2 COMPONENTS Interpreter Development NumPy)
  if(DEFINED Python2_VERSION AND "${Python2_VERSION}" VERSION_LESS "2.7")
    message(WARNING "Ignoring Python2 installation: unsupported version ${Python2_VERSION} (version>=2.7 required)")
  endif()
  if(Python2_Development_FOUND AND ${Python2_VERSION} VERSION_GREATER_EQUAL "2.7")
    if(NOT Python3_Development_FOUND)
      # Only Python2 was found, set as main
      set(PYTHON_EXECUTABLE "${Python2_EXECUTABLE}" CACHE INTERNAL "" FORCE)
      set(PYTHON_INCLUDE_DIRS "${Python2_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
      set(PYTHON_LIBRARIES "${Python2_LIBRARIES}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_STRING "${Python2_VERSION}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MAJOR "${Python2_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MINOR "${Python2_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
      set(PYTHON_UNDER_VERSION_STRING "${Python2_VERSION_MAJOR}_${Python2_VERSION_MINOR}")
      set(NUMPY_FOUND ${Python2_NumPy_FOUND} CACHE INTERNAL "" FORCE)
      set(NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
    else()
      # Both Python3 and 2 found, set 2 as 'other'
      set(OTHER_PYTHON_EXECUTABLE "${Python2_EXECUTABLE}")
      set(OTHER_PYTHON_INCLUDE_DIRS "${Python2_INCLUDE_DIRS}")
      set(OTHER_PYTHON_LIBRARIES "${Python2_LIBRARIES}")
      set(OTHER_PYTHON_VERSION_STRING "${Python2_VERSION}")
      set(OTHER_PYTHON_VERSION_MAJOR "${Python2_VERSION_MAJOR}")
      set(OTHER_PYTHON_VERSION_MINOR "${Python2_VERSION_MINOR}")
      set(OTHER_PYTHON_UNDER_VERSION_STRING "${Python2_VERSION_MAJOR}_${Python2_VERSION_MINOR}")
      set(OTHER_NUMPY_FOUND ${Python2_NumPy_FOUND})
      set(OTHER_NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}")
    endif()
  endif()

  if(NOT Python3_Development_FOUND AND (NOT Python2_Development_FOUND OR ${Python2_VERSION} VERSION_LESS "2.7"))
    message(FATAL_ERROR "No supported Python 2 or 3 were found")
  endif()

  # Print message saying with which versions of Python are used to build
  if(pyroot_experimental)
    if(NOT Python3_Development_FOUND OR NOT Python2_Development_FOUND)
      message(STATUS "Main Python used to build: ${PYTHON_VERSION_MAJOR}; PyROOT built for version ${PYTHON_VERSION_MAJOR} ")
    elseif(Python3_Development_FOUND AND Python2_Development_FOUND)
      message(STATUS "Main Python used to build: ${PYTHON_VERSION_MAJOR}; PyROOT built for versions ${PYTHON_VERSION_MAJOR} and ${OTHER_PYTHON_VERSION_MAJOR}")
    endif()
  endif()

else()

  find_package(PythonInterp ${python_version} REQUIRED)

  find_package(PythonLibs ${python_version} REQUIRED)

  if(NOT "${PYTHONLIBS_VERSION_STRING}" MATCHES "${PYTHON_VERSION_STRING}")
    message(FATAL_ERROR "Version mismatch between Python interpreter (${PYTHON_VERSION_STRING})"
    " and libraries (${PYTHONLIBS_VERSION_STRING}).\nROOT cannot work with this configuration. "
    "Please specify only PYTHON_EXECUTABLE to CMake with an absolute path to ensure matching versions are found.")
  endif()

  find_package(NumPy)

  set(PYTHON_UNDER_VERSION_STRING "${PYTHON_VERSION_MAJOR}_${PYTHON_VERSION_MINOR}")

endif()

# Create lists of Python 2 and 3 useful variables used to build PyROOT with both versions
# PYTHON_UNDER_VERSION_STRING and OTHER_PYTHON_UNDER_VERSION_STRING in particular are
# introduced because it's not possible to create a library containing '.' in the name
# before the suffix
set(python_executables ${PYTHON_EXECUTABLE} ${OTHER_PYTHON_EXECUTABLE})
set(python_include_dirs ${PYTHON_INCLUDE_DIRS} ${OTHER_PYTHON_INCLUDE_DIRS})
set(python_version_strings ${PYTHON_VERSION_STRING} ${OTHER_PYTHON_VERSION_STRING})
set(python_under_version_strings ${PYTHON_UNDER_VERSION_STRING} ${OTHER_PYTHON_UNDER_VERSION_STRING})
set(python_libraries ${PYTHON_LIBRARIES} ${OTHER_PYTHON_LIBRARIES})
