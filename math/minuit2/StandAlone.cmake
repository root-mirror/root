cmake_minimum_required(VERSION 3.1)

# Tested with and supporting policies up to the following CMake version. 
# Not using ... syntax due to parser bug in MSVC's built-in CMake server mode.
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.12)
endif()

include(FeatureSummary)
include(CMakeDependentOption)

# Check to see if we are inside ROOT and set a smart default
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../build/version_number")
    set(INROOT ON)
else()
    set(INROOT OFF)
endif()

cmake_dependent_option(minuit2_inroot "The source directory is inside the ROOT source" ON "INROOT" OFF)
cmake_dependent_option(minuit2_standalone "Copy in the files from the main ROOT files" OFF "minuit2_inroot" OFF)

if(minuit2_standalone)
    message(STATUS "Copying in files from ROOT sources to make a redistributable source package. You should clean out the new files with make purge or the appropriate git clean command when you are done.")
endif()


# This file adds copy_standalone
include(copy_standalone.cmake)

# Copy these files in if needed
copy_standalone(SOURCE ../../build DESTINATION . OUTPUT VERSION_FILE
                FILES version_number)

copy_standalone(SOURCE ../.. DESTINATION .
                FILES LGPL2_1.txt)

copy_standalone(SOURCE ../.. DESTINATION . OUTPUT LICENSE_FILE
                FILES LICENSE)

file(READ ${VERSION_FILE} versionstr)
string(STRIP ${versionstr} versionstr)
string(REGEX REPLACE "([0-9]+[.][0-9]+)[/]([0-9]+)" "\\1.\\2" versionstr ${versionstr})

project(Minuit2
    VERSION ${versionstr}
    LANGUAGES CXX)


# Inherit default from parent project if not main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS "Minuit2 ${PROJECT_VERSION} standalone")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# Common features to all packages (Math and Minuit2)
# If using this with add_subdirectory, the Minuit2
# namespace does not get automatically prepended,
# so including an alias for that.
add_library(Minuit2Common INTERFACE)
add_library(Minuit2::Common ALIAS Minuit2Common)

# OpenMP support
if(minuit2_omp)
    if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        message(STATUS "Building Minuit2 with OpenMP support")
    endif()
    target_link_libraries(Minuit2Common INTERFACE OpenMP::OpenMP_CXX)
endif()

# MPI support
# Uses the old CXX bindings (deprecated), probably do not activate
if(minuit2_mpi)
    if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        message(STATUS "Building Minuit2 with MPI support")
        message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
    endif()
    target_compile_definitions(Minuit2Common INTERFACE MPIPROC)
    target_link_libraries(Minuit2Common INTERFACE MPI::MPI_CXX)
endif()

# Add the libraries
add_subdirectory(src)

# Exporting targets to allow find_package(Minuit2) to work properly

# Make a config file to make this usable as a CMake Package
# Start by adding the version in a CMake understandable way
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    Minuit2ConfigVersion.cmake
    VERSION ${Minuit2_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# Now, install the Interface targets
install(TARGETS Minuit2Common
        EXPORT Minuit2Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

# Install the export set
install(EXPORT Minuit2Targets
        FILE Minuit2Targets.cmake
        NAMESPACE Minuit2::
        DESTINATION lib/cmake/Minuit2
        )

# Adding the Minuit2Config file
configure_file(Minuit2Config.cmake.in Minuit2Config.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Minuit2Config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/Minuit2ConfigVersion.cmake"
        DESTINATION lib/cmake/Minuit2
        )

# Allow build directory to work for CMake import
export(TARGETS Minuit2Common Minuit2Math Minuit2 NAMESPACE Minuit2:: FILE Minuit2Targets.cmake)
export(PACKAGE Minuit2)

# Only add tests and docs if this is the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    
    # Make adding tests cleaner using this macro
    macro(add_minuit2_test TESTNAME)
        add_executable(${TESTNAME} ${ARGN})
        target_include_directories(${TESTNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
        target_link_libraries(${TESTNAME} PUBLIC Minuit2)
        if(minuit2_mpi)
            add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:${TESTNAME}>" ${MPIEXEC_POSTFLAGS})
        else()
            add_test(NAME ${TESTNAME} COMMAND "$<TARGET_FILE:${TESTNAME}>")
        endif()
    endmacro()

    add_subdirectory(test/MnSim)
    add_subdirectory(test/MnTutorial)

    # Add a build docs target if Doxygen found
    add_subdirectory(doc)
endif()

# Add purge target
if(minuit2_standalone)
    get_property(COPY_STANDALONE_LISTING GLOBAL PROPERTY COPY_STANDALONE_LISTING)
    add_custom_target(purge
        COMMAND ${CMAKE_COMMAND} -E remove ${COPY_STANDALONE_LISTING})
endif()

# Setup package info
add_feature_info(minuit2_omp minuit2_omp "OpenMP (Thread safe FCNs only)")
add_feature_info(minuit2_mpi minuit2_mpi "MPI (Thread safe FCNs only)")
set_package_properties(OpenMP PROPERTIES
    URL "http://www.openmp.org"
    DESCRIPTION "Parallel compiler directives"
    PURPOSE "Parallel FCN calls (Thread safe FCNs only)")
set_package_properties(MPI PROPERTIES
    URL "http://mpi-forum.org"
    DESCRIPTION "Message passing interface"
    PURPOSE "Separate threads (Thread safe FCNs only)")

# Print package info to screen and log if this is the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND)
    feature_summary(FILENAME ${CMAKE_CURRENT_BINARY_DIR}/features.log WHAT ALL)
endif()

# Packaging support
set(CPACK_PACKAGE_VENDOR "root.cern.ch")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Minuit2 standalone fitting tool")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${LICENSE_FILE}")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

# Making a source package should fail if standalone has not been built
# Setting an obvious package name in case a generator is manually specified
if(minuit2_inroot AND NOT minuit2_standalone)
    set(CPACK_SOURCE_GENERATOR "ERROR_MINUIT2_STANDALONE_OFF")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "Minuit2-MISSING_FILES-Source")
else()
    set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
endif()

# CPack collects *everything* except what's listed here.
set(CPACK_SOURCE_IGNORE_FILES
    /test/CMakeLists.txt
    /test/Makefile
    /test/testMinimizer.cxx
    /test/testNdimFit.cxx
    /test/testUnbinGausFit.cxx
    /test/testUserFunc.cxx
    /Module.mk
    /.git
    /dist
    /.*build.*
    /\\\\.DS_Store
    /.*\\\\.egg-info
    /var
    /Pipfile.*$
)
include(CPack)

