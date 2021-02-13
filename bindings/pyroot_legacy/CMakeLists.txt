# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT bindings/pyroot package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(PyROOT
  NO_INSTALL_HEADERS
  HEADERS
    TPyArg.h
    TPyDispatcher.h
    TPyException.h
    TPyFitFunction.h
    TPyReturn.h
    TPyROOTApplication.h
    TPySelector.h
    TPython.h
  SOURCES
    src/Converters.cxx
    src/Cppyy.cxx
    src/Executors.cxx
    src/MethodProxy.cxx
    src/ObjectProxy.cxx
    src/PropertyProxy.cxx
    src/PyRootType.cxx
    src/PyStrings.cxx
    src/Pythonize.cxx
    src/RootModule.cxx
    src/RootWrapper.cxx
    src/TCallContext.cxx
    src/TClassMethodHolder.cxx
    src/TConstructorHolder.cxx
    src/TCustomPyTypes.cxx
    src/TemplateProxy.cxx
    src/TFunctionHolder.cxx
    src/TMemoryRegulator.cxx
    src/TMethodHolder.cxx
    src/TPyArg.cxx
    src/TPyBufferFactory.cxx
    src/TPyClassGenerator.cxx
    src/TPyDispatcher.cxx
    src/TPyException.cxx
    src/TPyFitFunction.cxx
    src/TPyReturn.cxx
    src/TPyROOTApplication.cxx
    src/TPySelector.cxx
    src/TPython.cxx
    src/TSetItemHolder.cxx
    src/TTupleOfInstances.cxx
    src/Utility.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    MathCore
    ROOTVecOps
    Tree
)

ROOT_ADD_CXX_FLAG(_PyROOT_FLAGS -fno-strict-aliasing)
ROOT_ADD_CXX_FLAG(_PyROOT_FLAGS -Wno-parentheses-equality)
ROOT_ADD_CXX_FLAG(_PyROOT_FLAGS -Wno-cast-function-type)
ROOT_ADD_CXX_FLAG(_PyROOT_FLAGS -Wno-deprecated-register)
ROOT_ADD_CXX_FLAG(_PyROOT_FLAGS -Wno-register)

separate_arguments(_PyROOT_FLAGS)

target_compile_options(PyROOT PRIVATE ${_PyROOT_FLAGS})
target_link_libraries(PyROOT PRIVATE ${PYTHON_LIBRARIES_Development_Main})
target_include_directories(PyROOT PRIVATE ${PYTHON_INCLUDE_DIRS_Development_Main})

ROOT_LINKER_LIBRARY(JupyROOT JupyROOT/src/IOHandler.cxx DEPENDENCIES Core CMAKENOEXPORT)

if(MSVC)
  if(${PYTHON_VERSION_MAJOR_Development_Main} GREATER_EQUAL 3)
    set(init_symbol PyInit_libPyROOT)
  else()
    set(init_symbol initlibPyROOT)
  endif()
  add_custom_command(TARGET PyROOT POST_BUILD
                     COMMAND link -dll -nologo -IGNORE:4001 -machine:ix86 -export:${init_symbol}
                     $<TARGET_LINKER_FILE:PyROOT> msvcrt.lib
                     -out:$<TARGET_FILE_DIR:PyROOT>/libPyROOT.pyd)
  install(FILES  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/libPyROOT.pyd
                 DESTINATION ${CMAKE_INSTALL_BINDIR}
                 COMPONENT libraries)
endif()

#---Install python modules--------------------------------------------------

foreach(pyfile cppyy.py _pythonization.py ROOT.py _rdf_utils.py)
  install(FILES ${pyfile} DESTINATION ${runtimedir})
  if(IS_ABSOLUTE ${runtimedir})
    set(absruntimedir ${runtimedir})
  else()
    set(absruntimedir \${CMAKE_INSTALL_PREFIX}/${runtimedir})
  endif()
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE_Development_Main} -m py_compile \$ENV{DESTDIR}${absruntimedir}/${pyfile})")
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE_Development_Main} -O -m py_compile \$ENV{DESTDIR}${absruntimedir}/${pyfile})")
  file(COPY ${pyfile} DESTINATION ${localruntimedir})
endforeach()

if (NOT MSVC AND NOT CMAKE_INSTALL_LIBDIR STREQUAL CMAKE_INSTALL_PYTHONDIR)
  # add a symlink to libPyROOT.so in CMAKE_INSTALL_PYTHONDIR
  install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink\
          ${CMAKE_INSTALL_FULL_LIBDIR}/libPyROOT.so\
          ${CMAKE_INSTALL_FULL_PYTHONDIR}/libPyROOT.so)")
endif()


set(JupyROOTDirName "JupyROOT")
install (DIRECTORY ${JupyROOTDirName} DESTINATION ${runtimedir})
file(COPY ${JupyROOTDirName} DESTINATION ${localruntimedir})

set( JsMVADirName "JsMVA")
install (DIRECTORY ${JsMVADirName} DESTINATION ${runtimedir})
file(COPY ${JsMVADirName} DESTINATION ${localruntimedir})

ROOT_INSTALL_HEADERS()

ROOT_ADD_TEST_SUBDIRECTORY(test)
