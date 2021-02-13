# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

################################################################
# CMakeLists.txt file for building ROOT pythonizations libraries
################################################################

if(dataframe)
    list(APPEND PYROOT_EXTRA_PYSOURCE
        ROOT/pythonization/_rdf_utils.py
        ROOT/pythonization/_rdataframe.py
        ROOT/pythonization/_rtensor.py)
    list(APPEND PYROOT_EXTRA_SOURCE
        src/RDataFramePyz.cxx
        src/RTensorPyz.cxx)
    list(APPEND PYROOT_EXTRA_HEADERS
        inc/RNumpyDS.hxx)
endif()

list(APPEND PYROOT_EXTRA_HEADERS
     inc/TPyDispatcher.h)

set(py_sources
  ROOT/__init__.py
  ROOT/_application.py
  ROOT/_numbadeclare.py
  ROOT/_facade.py
  ROOT/pythonization/__init__.py
  ROOT/pythonization/_cppinstance.py
  ROOT/pythonization/_drawables.py
  ROOT/pythonization/_generic.py
  ROOT/pythonization/_rbdt.py
  ROOT/pythonization/_rooabscollection.py
  ROOT/pythonization/_rooworkspace.py
  ROOT/pythonization/_rvec.py
  ROOT/pythonization/_stl_vector.py
  ROOT/pythonization/_tarray.py
  ROOT/pythonization/_tclonesarray.py
  ROOT/pythonization/_tcollection.py
  ROOT/pythonization/_tcomplex.py
  ROOT/pythonization/_tdirectory.py
  ROOT/pythonization/_tdirectoryfile.py
  ROOT/pythonization/_tfile.py
  ROOT/pythonization/_tgraph.py
  ROOT/pythonization/_th1.py
  ROOT/pythonization/_titer.py
  ROOT/pythonization/_tobject.py
  ROOT/pythonization/_tobjstring.py
  ROOT/pythonization/_tree_inference.py
  ROOT/pythonization/_tseqcollection.py
  ROOT/pythonization/_tstring.py
  ROOT/pythonization/_ttree.py
  ROOT/pythonization/_tvector3.py
  ROOT/pythonization/_tvectort.py
  ${PYROOT_EXTRA_PYSOURCE}
)

set(cpp_sources
    src/PyROOTModule.cxx
    src/PyROOTStrings.cxx
    src/PyROOTWrapper.cxx
    src/RPyROOTApplication.cxx
    src/GenericPyz.cxx
    src/RVecPyz.cxx
    src/TClassPyz.cxx
    src/TClonesArrayPyz.cxx
    src/TDirectoryFilePyz.cxx
    src/TDirectoryPyz.cxx
    src/TFilePyz.cxx
    src/TMemoryRegulator.cxx
    src/TObjectPyz.cxx
    src/TTreePyz.cxx
    src/PyzCppHelpers.cxx
    src/PyzPythonHelpers.cxx
    src/CPPInstancePyz.cxx
    src/FacadeHelpers.cxx
    src/TPyDispatcher.cxx
    inc/TPyDispatcher.h
    ${PYROOT_EXTRA_SOURCE}
)

set(ROOTPySrcDir python/ROOT)
set(ROOT_headers_dir inc)

file(COPY ${ROOTPySrcDir} DESTINATION ${localruntimedir})
# Copy headers inside build_dir/include/ROOT
file(COPY ${ROOT_headers_dir}/ DESTINATION ${CMAKE_BINARY_DIR}/include/ROOT)

foreach(val RANGE ${how_many_pythons})
  list(GET python_under_version_strings ${val} python_under_version_string)
  list(GET python_major_version_strings ${val} python_major_version_string)
  list(GET python_include_dirs ${val} python_include_dir)
  list(GET python_executables ${val} python_executable)

  set(libname ROOTPythonizations${python_under_version_string})

  add_library(${libname} SHARED ${cpp_sources})
  # Set the suffix to '.so' and the prefix to 'lib'
  set_target_properties(${libname} PROPERTIES  ${ROOT_LIBRARY_PROPERTIES})
  if(MSVC)
    set_target_properties(${libname} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set_target_properties(${libname} PROPERTIES SUFFIX ".pyd")
    target_link_libraries(${libname} PUBLIC Core Tree cppyy${python_under_version_string})
  elseif(APPLE)
    target_link_libraries(${libname} PUBLIC -Wl,-bind_at_load -Wl,-undefined -Wl,dynamic_lookup Core Tree cppyy${python_under_version_string})
  else()
    target_link_libraries(${libname} PUBLIC -Wl,--unresolved-symbols=ignore-all Core Tree cppyy${python_under_version_string})
  endif()

  target_include_directories(${libname}
     PRIVATE ${python_include_dir}
     PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>)

  # Disables warnings caused by Py_RETURN_TRUE/Py_RETURN_FALSE
  if(NOT MSVC)
    target_compile_options(${libname} PRIVATE -Wno-strict-aliasing)
  endif()

  # Disables warnings originating from deprecated register keyword in Python
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_STANDARD GREATER_EQUAL 11)
    target_compile_options(${libname} PRIVATE -Wno-register)
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" AND CMAKE_CXX_STANDARD GREATER_EQUAL 11)
    target_compile_options(${libname} PRIVATE -Wno-register)
    target_compile_options(${libname} PRIVATE -Wno-deprecated-register)
  endif()

  # Compile .py files
  foreach(py_source ${py_sources})
    install(CODE "execute_process(COMMAND ${python_executable} -m py_compile ${localruntimedir}/${py_source})")
    install(CODE "execute_process(COMMAND ${python_executable} -O -m py_compile ${localruntimedir}/${py_source})")
  endforeach()

  # Create meta-target PyROOT2 and PyROOT3 (INTERFACE library)
  # Export of targets are not supported for custom targets(add_custom_targets())
  add_library(PyROOT${python_major_version_string} INTERFACE)
  target_link_libraries(PyROOT${python_major_version_string} INTERFACE cppyy_backend${python_under_version_string} cppyy${python_under_version_string} ROOTPythonizations${python_under_version_string})

  # Install library
  install(TARGETS ${libname} EXPORT ${CMAKE_PROJECT_NAME}Exports
                             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries
                             LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHONDIR} COMPONENT libraries
                             ARCHIVE DESTINATION ${CMAKE_INSTALL_PYTHONDIR} COMPONENT libraries)

  # Install meta-target PyROOT2 and PyROOT3 (INTERFACE library)
  # Install library
  install(TARGETS PyROOT${python_major_version_string} EXPORT ${CMAKE_PROJECT_NAME}Exports
                             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries
                             LIBRARY DESTINATION ${CMAKE_INSTALL_PYTHONDIR} COMPONENT libraries
                             ARCHIVE DESTINATION ${CMAKE_INSTALL_PYTHONDIR} COMPONENT libraries)
endforeach()

# Install Python sources and bytecode
install(DIRECTORY ${localruntimedir}/ROOT
        DESTINATION ${CMAKE_INSTALL_PYTHONDIR}
        COMPONENT libraries)

# Install headers required by pythonizations
install(FILES ${PYROOT_EXTRA_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ROOT
        COMPONENT headers)

ROOT_ADD_TEST_SUBDIRECTORY(test)
