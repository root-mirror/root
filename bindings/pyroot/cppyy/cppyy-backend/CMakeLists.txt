# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(py_sources
  cppyy_backend/__init__.py
  cppyy_backend/_cling_config.py
  cppyy_backend/_cppyy_generator.py
  cppyy_backend/_genreflex.py
  cppyy_backend/_rootcling.py
  cppyy_backend/bindings_utils.py
  cppyy_backend/loader.py
  cppyy_backend/_get_cppflags.py
)

set(cppyy_backendPySrcDir cling/python/cppyy_backend)
file(COPY ${cppyy_backendPySrcDir}
     DESTINATION ${localruntimedir}
     PATTERN "cmake" EXCLUDE
     PATTERN "pkg_templates" EXCLUDE)

foreach(val RANGE ${how_many_pythons})
  list(GET python_under_version_strings ${val} python_under_version_string)
  list(GET python_executables ${val} python_executable)

  set(libname cppyy_backend${python_under_version_string})

  add_library(${libname} SHARED clingwrapper/src/clingwrapper.cxx)
  if(MSVC)
    set_target_properties(${libname} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  endif()
  # Set the suffix to '.so' and the prefix to 'lib'
  set_target_properties(${libname} PROPERTIES  ${ROOT_LIBRARY_PROPERTIES})
  target_link_libraries(${libname} Core ${CMAKE_DL_LIBS})

  # cppyy uses ROOT headers from binary directory
  add_dependencies(${libname} move_headers)

  set_property(GLOBAL APPEND PROPERTY ROOT_EXPORTED_TARGETS ${libname})

  # Install library
  install(TARGETS ${libname} EXPORT ${CMAKE_PROJECT_NAME}Exports
                             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries
                             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
                             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)
  if (NOT MSVC AND NOT CMAKE_INSTALL_LIBDIR STREQUAL CMAKE_INSTALL_PYTHONDIR)
    # add a symlink to ${libname} in CMAKE_INSTALL_PYTHONDIR
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink\
            ${CMAKE_INSTALL_FULL_LIBDIR}/lib${libname}.so\
            ${CMAKE_INSTALL_FULL_PYTHONDIR}/lib${libname}.so)")
  endif()

  # Compile .py files
  foreach(py_source ${py_sources})
    install(CODE "execute_process(COMMAND ${python_executable} -m py_compile ${localruntimedir}/${py_source})")
    install(CODE "execute_process(COMMAND ${python_executable} -O -m py_compile ${localruntimedir}/${py_source})")
  endforeach()
endforeach()

# Install Python sources and bytecode
install(DIRECTORY ${localruntimedir}/cppyy_backend
        DESTINATION ${CMAKE_INSTALL_PYTHONDIR}
        COMPONENT libraries)
