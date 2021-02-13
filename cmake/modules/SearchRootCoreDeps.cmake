# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#[[
The purpose of this machinery is to search Python and set both the ROOT and PyROOT related
Python variables.

ROOT:
    - PYTHON_EXECUTABLE (set when Interpreter is found)
    - PYTHON_VERSION_STRING        ""
    - PYTHON_VERSION_MAJOR        ""
    - PYTHON_VERSION_MINOR        ""

PyROOT: a set of specific variables is needed for all the Python versions used to build PyROOT.
At the time of writing we only have two, called simply Main and Other.

    - PYTHON_EXECUTABLE_Development_Main (set when Interpreter is found)
    - PYTHON_VERSION_STRING_Development_Main        ""
    - PYTHON_VERSION_MAJOR_Development_Main        ""
    - PYTHON_VERSION_MINOR_Development_Main        ""
    - PYTHON_INCLUDE_DIRS_Development_Main (set when Development is found)
    - PYTHON_LIBRARIES_Development_Main        ""
    - PYTHON_UNDER_VERSION_STRING_Development_Main (used to give a version-dependent name to the libraries, to allow
    multiple builds)
    - Development_Python${version}_FOUND (set if both Interpreter and Development are found for the specified ${version})

    - PYTHON_EXECUTABLE_Development_Other (set when Interpreter is found)
    - PYTHON_VERSION_STRING_Development_Other        ""
    - PYTHON_VERSION_MAJOR_Development_Other        ""
    - PYTHON_VERSION_MINOR_Development_Other        ""
    - PYTHON_INCLUDE_DIRS_Development_Other (set when Development is found)
    - PYTHON_LIBRARIES_Development_Other        ""
    - PYTHON_UNDER_VERSION_STRING_Development_Other (used to give a version-dependent name to the libraries, to allow
    multiple builds)
    - Development_Python${version}_FOUND (set if both Interpreter and Development are found for the specified ${version})

^^^^^^^^

Explanation of the machinery:

The first distinction is based on the CMake version used to build. If it is >= 3.14, than PyROOT can be built
for multiple Python versions. In case CMake >= 3.14, then we check if PYTHON_EXECUTABLE was specified by the user;
if so, PyROOT is built with only that version, otherwise an attempt to build it both with Python2 and 3 is done.

If PYTHON_EXECUTABLE is specified:
    - we check which version it is (from here call X) and call the relative
    find_package(PythonX COMPONENTS Interpreter Development Numpy)
    - if Interpreter is found, we set the ROOT related Python variables
    - if Development is found, we set the PyROOT_Main variables (+ Numpy ones)

If PYTHON_EXECUTABLE is NOT specified:
    - we look for Python3, since we want the highest version to be the preferred one
    - if Python3 Interpreter is found, we set the ROOT related Python variables
    - if Python3 Development is found, we set the PyROOT_Main variables (+ Numpy ones)
    - we then look for Python2, requiring it to be >= 2.7:
        - if Python3 Interpreter wasn't previously found, then it means that Python2 becomes the one used to
        build ROOT and the Main one for PyROOT; we set the ROOT related Python variables and, if Python2 Development
        was found, the PyROOT_Main variables (+ Numpy ones)
        - if Python3 Interpreter was previously found, we then check if both Python2 Interpreter and Development were
        found; if yes, we set the PyROOT_Other variables, otherwise PyROOT will be built with only Python 3
]]

message(STATUS "Looking for Python")

# On macOS, prefer user-provided Pythons.
set(Python_FIND_FRAMEWORK LAST)
set(Python2_FIND_FRAMEWORK LAST)
set(Python3_FIND_FRAMEWORK LAST)

# CMake 3.14 is the minimum because in 3.12 and 3.13, despite find_package(PythonX) being present,
# the Numpy-related variables are not set
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)

  if(PYTHON_EXECUTABLE)
    # - Check if the PYTHON_EXECUTABLE deprecated variable
    # was passed by the user; if so, PyROOT will be built with a single version
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version_info[0])"
                    OUTPUT_VARIABLE PYTHON_PREFER_VERSION
                    ERROR_VARIABLE PYTHON_PREFER_VERSION_ERR)
    if(PYTHON_PREFER_VERSION_ERR)
      message(WARNING "Unable to determine version of ${PYTHON_EXECUTABLE}: ${PYTHON_PREFER_VERSION_ERR}")
    endif()
    string(STRIP "${PYTHON_PREFER_VERSION}" PYTHON_PREFER_VERSION)
    set(Python${PYTHON_PREFER_VERSION}_EXECUTABLE "${PYTHON_EXECUTABLE}")
    find_package(Python${PYTHON_PREFER_VERSION} COMPONENTS Interpreter Development NumPy)
    if(Python${PYTHON_PREFER_VERSION}_Interpreter_FOUND)
      set(PYTHON_EXECUTABLE "${Python${PYTHON_PREFER_VERSION}_EXECUTABLE}")
      set(PYTHON_VERSION_STRING "${Python${PYTHON_PREFER_VERSION}_VERSION}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MAJOR "${Python${PYTHON_PREFER_VERSION}_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MINOR "${Python${PYTHON_PREFER_VERSION}_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
      if(Python${PYTHON_PREFER_VERSION}_Development_FOUND)
        set(PYTHON_INCLUDE_DIRS "${Python${PYTHON_PREFER_VERSION}_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
        set(PYTHON_LIBRARIES "${Python${PYTHON_PREFER_VERSION}_LIBRARIES}" CACHE INTERNAL "" FORCE)
        set(PYTHON_LIBRARY_DIR "${Python${PYTHON_PREFER_VERSION}_LIBRARY_DIRS}" CACHE INTERNAL "" FORCE)
        # Set PyROOT variables
        set(Python${PYTHON_PREFER_VERSION}_Interpreter_Development_FOUND ON) # This means we have both Interpreter and Development, hence we can build PyROOT with this Python
        set(PYTHON_EXECUTABLE_Development_Main "${Python${PYTHON_PREFER_VERSION}_EXECUTABLE}")
        set(PYTHON_VERSION_STRING_Development_Main "${Python${PYTHON_PREFER_VERSION}_VERSION}")
        set(PYTHON_VERSION_MAJOR_Development_Main "${Python${PYTHON_PREFER_VERSION}_VERSION_MAJOR}")
        set(PYTHON_VERSION_MINOR_Development_Main "${Python${PYTHON_PREFER_VERSION}_VERSION_MINOR}")
        set(PYTHON_UNDER_VERSION_STRING_Development_Main "${Python${PYTHON_PREFER_VERSION}_VERSION_MAJOR}_${Python${PYTHON_PREFER_VERSION}_VERSION_MINOR}")
        set(PYTHON_INCLUDE_DIRS_Development_Main "${Python${PYTHON_PREFER_VERSION}_INCLUDE_DIRS}")
        set(PYTHON_LIBRARIES_Development_Main "${Python${PYTHON_PREFER_VERSION}_LIBRARIES}")
      endif()
      if(Python${PYTHON_PREFER_VERSION}_NumPy_FOUND)
        set(NUMPY_FOUND ${Python${PYTHON_PREFER_VERSION}_NumPy_FOUND})
        set(NUMPY_INCLUDE_DIRS "${Python${PYTHON_PREFER_VERSION}_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
      endif()
    endif()
    if(DEFINED Python${PYTHON_PREFER_VERSION}_VERSION AND "${Python${PYTHON_PREFER_VERSION}_VERSION}" VERSION_LESS "2.7")
      message(FATAL_ERROR "Ignoring Python installation: unsupported version ${Python${PYTHON_PREFER_VERSION}_VERSION} (version>=2.7 required)")
    endif()
    if(NOT (Python${PYTHON_PREFER_VERSION}_Interpreter_FOUND AND Python${PYTHON_PREFER_VERSION}_Development_FOUND))
      message(WARNING "No supported Python development package was found for the specified Python executable; PyROOT will not be built")
    endif()

  else()
    # - Look for Python3 and set the deprecated variables to the ones set
    # automatically by find_package(Python3 ...)
    # - Look for Python2 and set the deprecated variables to the ones set
    # automatically by find_package(Python2 ...) ONLY IF PYTHON3 WASN'T FOUND
    find_package(Python3 COMPONENTS Interpreter Development NumPy)
    if(Python3_Interpreter_FOUND)
      set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
      set(PYTHON_VERSION_STRING "${Python3_VERSION}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MAJOR "${Python3_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
      set(PYTHON_VERSION_MINOR "${Python3_VERSION_MINOR}" CACHE INTERNAL "" FORCE )
      if(Python3_Development_FOUND)
        set(PYTHON_INCLUDE_DIRS "${Python3_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
        set(PYTHON_LIBRARIES "${Python3_LIBRARIES}" CACHE INTERNAL "" FORCE)
        set(PYTHON_LIBRARY_DIR "${Python3_LIBRARY_DIRS}" CACHE INTERNAL "" FORCE)
        # Set PyROOT variables
        set(Python3_Interpreter_Development_FOUND ON) # This means we have both Interpreter and Development, hence we can build PyROOT with Python3
        set(PYTHON_EXECUTABLE_Development_Main "${Python3_EXECUTABLE}")
        set(PYTHON_VERSION_STRING_Development_Main "${Python3_VERSION}")
        set(PYTHON_VERSION_MAJOR_Development_Main "${Python3_VERSION_MAJOR}")
        set(PYTHON_VERSION_MINOR_Development_Main "${Python3_VERSION_MINOR}")
        set(PYTHON_UNDER_VERSION_STRING_Development_Main "${Python3_VERSION_MAJOR}_${Python3_VERSION_MINOR}")
        set(PYTHON_INCLUDE_DIRS_Development_Main "${Python3_INCLUDE_DIRS}")
        set(PYTHON_LIBRARIES_Development_Main "${Python3_LIBRARIES}")
      endif()
      if(Python3_NumPy_FOUND)
        set(NUMPY_FOUND ${Python3_NumPy_FOUND})
        set(NUMPY_INCLUDE_DIRS "${Python3_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
      endif()
    endif()

    find_package(Python2 COMPONENTS Interpreter Development NumPy)
    if(DEFINED Python2_VERSION AND "${Python2_VERSION}" VERSION_LESS "2.7")
      message(WARNING "Ignoring Python2 installation: unsupported version ${Python2_VERSION} (version>=2.7 required)")
    endif()
    if("${Python2_VERSION}" VERSION_GREATER_EQUAL "2.7")
      if(NOT Python3_Interpreter_FOUND)
        # Only Python2 was found, set as main
        if(Python2_Interpreter_FOUND)
          set(PYTHON_EXECUTABLE "${Python2_EXECUTABLE}")
          set(PYTHON_VERSION_STRING "${Python2_VERSION}" CACHE INTERNAL "" FORCE)
          set(PYTHON_VERSION_MAJOR "${Python2_VERSION_MAJOR}" CACHE INTERNAL "" FORCE)
          set(PYTHON_VERSION_MINOR "${Python2_VERSION_MINOR}" CACHE INTERNAL "" FORCE)
          if(Python2_Development_FOUND)
            set(PYTHON_INCLUDE_DIRS "${Python2_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
            set(PYTHON_LIBRARIES "${Python2_LIBRARIES}" CACHE INTERNAL "" FORCE)
            set(PYTHON_LIBRARY_DIR "${Python2_LIBRARY_DIRS}" CACHE INTERNAL "" FORCE)
            # Set PyROOT variables
            set(Python2_Interpreter_Development_FOUND ON) # This means we have both Interpreter and Development, hence we can build PyROOT with Python2
            set(PYTHON_EXECUTABLE_Development_Main "${Python2_EXECUTABLE}")
            set(PYTHON_VERSION_STRING_Development_Main "${Python2_VERSION}")
            set(PYTHON_VERSION_MAJOR_Development_Main "${Python2_VERSION_MAJOR}")
            set(PYTHON_VERSION_MINOR_Development_Main "${Python2_VERSION_MINOR}")
            set(PYTHON_UNDER_VERSION_STRING_Development_Main "${Python2_VERSION_MAJOR}_${Python2_VERSION_MINOR}")
            set(PYTHON_INCLUDE_DIRS_Development_Main "${Python2_INCLUDE_DIRS}")
            set(PYTHON_LIBRARIES_Development_Main "${Python2_LIBRARIES}")
          endif()
          if(Python2_NumPy_FOUND)
            set(NUMPY_FOUND ${Python2_NumPy_FOUND})
            set(NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
          endif()
        endif()
      else()
        # Both Python3 and 2 found, set 2 as 'other'
        # In this case, since the 'other' variables are used only for PyROOT (which requires development package),
        # we can simply use the if(Python2_Development_FOUND) condition
        if(Python2_Interpreter_FOUND AND Python2_Development_FOUND)
          set(Python2_Interpreter_Development_FOUND ON) # This means we have both Interpreter and Development, hence we can build PyROOT with Python2
          if(Python3_Interpreter_Development_FOUND)
            set(PYTHON_EXECUTABLE_Development_Other "${Python2_EXECUTABLE}")
            set(PYTHON_VERSION_STRING_Development_Other "${Python2_VERSION}")
            set(PYTHON_VERSION_MAJOR_Development_Other "${Python2_VERSION_MAJOR}")
            set(PYTHON_VERSION_MINOR_Development_Other "${Python2_VERSION_MINOR}")
            set(PYTHON_UNDER_VERSION_STRING_Development_Other "${Python2_VERSION_MAJOR}_${Python2_VERSION_MINOR}")
            set(PYTHON_INCLUDE_DIRS_Development_Other "${Python2_INCLUDE_DIRS}")
            set(PYTHON_LIBRARIES_Development_Other "${Python2_LIBRARIES}")
            if(Python2_NumPy_FOUND)
              set(OTHER_NUMPY_FOUND ${Python2_NumPy_FOUND})
              set(OTHER_NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}")
            endif()
          else()
            set(PYTHON_INCLUDE_DIRS "${Python2_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
            set(PYTHON_LIBRARIES "${Python2_LIBRARIES}" CACHE INTERNAL "" FORCE)
            set(PYTHON_LIBRARY_DIR "${Python2_LIBRARY_DIRS}" CACHE INTERNAL "" FORCE)
            # Set PyROOT variables
            set(PYTHON_EXECUTABLE_Development_Main "${Python2_EXECUTABLE}")
            set(PYTHON_VERSION_STRING_Development_Main "${Python2_VERSION}")
            set(PYTHON_VERSION_MAJOR_Development_Main "${Python2_VERSION_MAJOR}")
            set(PYTHON_VERSION_MINOR_Development_Main "${Python2_VERSION_MINOR}")
            set(PYTHON_UNDER_VERSION_STRING_Development_Main "${Python2_VERSION_MAJOR}_${Python2_VERSION_MINOR}")
            set(PYTHON_INCLUDE_DIRS_Development_Main "${Python2_INCLUDE_DIRS}")
            set(PYTHON_LIBRARIES_Development_Main "${Python2_LIBRARIES}")
            if(Python2_NumPy_FOUND)
              set(NUMPY_FOUND ${Python2_NumPy_FOUND})
              set(NUMPY_INCLUDE_DIRS "${Python2_NumPy_INCLUDE_DIRS}" CACHE INTERNAL "" FORCE)
            endif()
          endif()
        endif()
      endif()
    endif()

    if(NOT Python3_Interpreter_Development_FOUND AND (NOT Python2_Interpreter_Development_FOUND OR "${Python2_VERSION}" VERSION_LESS "2.7"))
      message(WARNING "No supported Python 2 or 3 development packages were found; PyROOT will not be built.")
    endif()

  endif()

else()

  find_package(PythonInterp ${python_version} REQUIRED)

  find_package(PythonLibs ${python_version})

  if(PYTHONLIBS_FOUND)
    if(NOT "${PYTHONLIBS_VERSION_STRING}" MATCHES "${PYTHON_VERSION_STRING}")
      message(FATAL_ERROR "Version mismatch between Python interpreter (${PYTHON_VERSION_STRING})"
      " and libraries (${PYTHONLIBS_VERSION_STRING}).\nROOT cannot work with this configuration. "
      "Please specify only PYTHON_EXECUTABLE to CMake with an absolute path to ensure matching versions are found.")
    endif()
  else()
    message(WARNING "No supported Python development package was found; PyROOT will not be built.")
  endif()

  find_package(NumPy)

  if (PYTHONLIBS_FOUND)
    set(Python${PYTHON_VERSION_MAJOR}_Interpreter_Development_FOUND ON) # This means we have both Interpreter and Development, hence we can build PyROOT with this Python
    set(PYTHON_EXECUTABLE_Development_Main "${PYTHON_EXECUTABLE}")
    set(PYTHON_VERSION_STRING_Development_Main "${PYTHON_VERSION_STRING}")
    set(PYTHON_VERSION_MAJOR_Development_Main "${PYTHON_VERSION_MAJOR}")
    set(PYTHON_VERSION_MINOR_Development_Main "${PYTHON_VERSION_MINOR}")
    set(PYTHON_UNDER_VERSION_STRING_Development_Main "${PYTHON_VERSION_MAJOR}_${PYTHON_VERSION_MINOR}")
    set(PYTHON_INCLUDE_DIRS_Development_Main "${PYTHON_INCLUDE_DIRS}")
    set(PYTHON_LIBRARIES_Development_Main "${PYTHON_LIBRARIES}")
  endif()

endif()

# Create lists of Python 2 and 3 useful variables used to build PyROOT with both versions
# PYTHON_UNDER_VERSION_STRING and OTHER_PYTHON_UNDER_VERSION_STRING in particular are
# introduced because it's not possible to create a library containing '.' in the name
# before the suffix
set(python_executables ${PYTHON_EXECUTABLE_Development_Main} ${PYTHON_EXECUTABLE_Development_Other})
set(python_include_dirs ${PYTHON_INCLUDE_DIRS_Development_Main} ${PYTHON_INCLUDE_DIRS_Development_Other})
set(python_version_strings ${PYTHON_VERSION_STRING_Development_Main} ${PYTHON_VERSION_STRING_Development_Other})
set(python_major_version_strings ${PYTHON_VERSION_MAJOR_Development_Main} ${PYTHON_VERSION_MAJOR_Development_Other})
set(python_under_version_strings ${PYTHON_UNDER_VERSION_STRING_Development_Main} ${PYTHON_UNDER_VERSION_STRING_Development_Other})
set(python_libraries ${PYTHON_LIBRARIES_Development_Main} ${PYTHON_LIBRARIES_Development_Other})
