# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT main package
# @author Pere Mato, CERN
############################################################################

if(NOT WIN32)
  ROOT_EXECUTABLE(rootn.exe rmain.cxx LIBRARIES New Core MathCore Rint)
  if(CMAKE_SYSTEM_NAME STREQUAL Linux)
    SET_TARGET_PROPERTIES(rootn.exe PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
  endif()
  ROOT_EXECUTABLE(roots.exe roots.cxx LIBRARIES Core MathCore)
  ROOT_EXECUTABLE(xpdtest xpdtest.cxx LIBRARIES Proof Tree Hist RIO Net Thread Matrix MathCore)
endif()
ROOT_EXECUTABLE(root.exe rmain.cxx LIBRARIES Core Rint)
if(MSVC)
  set(root_exports "/EXPORT:_Init_thread_abort /EXPORT:_Init_thread_epoch
      /EXPORT:_Init_thread_footer /EXPORT:_Init_thread_header /EXPORT:_tls_index
      /STACK:4000000")
  set_property(TARGET root.exe APPEND_STRING PROPERTY LINK_FLAGS ${root_exports})
endif()
ROOT_EXECUTABLE(proofserv.exe pmain.cxx LIBRARIES Core MathCore)
if(MSVC)
  ROOT_EXECUTABLE(hadd hadd.cxx LIBRARIES Core RIO Net Hist Graf Graf3d Gpad Tree Matrix MathCore)
else()
  ROOT_EXECUTABLE(hadd hadd.cxx LIBRARIES Core RIO Net Hist Graf Graf3d Gpad Tree Matrix MathCore MultiProc)
endif()
ROOT_EXECUTABLE(rootnb.exe nbmain.cxx LIBRARIES Core)

#---CreateHaddCommandLineOptions------------------------------------------------------------------
generateHeader(hadd
  ${CMAKE_SOURCE_DIR}/main/src/hadd-argparse.py
  ${CMAKE_BINARY_DIR}/ginclude/haddCommandLineOptionsHelp.h
)

if(fortran AND CMAKE_Fortran_COMPILER)
  ROOT_EXECUTABLE(g2root g2root.f LIBRARIES minicern)
  set_target_properties(g2root PROPERTIES COMPILE_FLAGS "-w")
  ROOT_EXECUTABLE(h2root h2root.cxx LIBRARIES Core RIO Net Hist Graf Graf3d Gpad Tree Matrix MathCore Thread minicern)
endif()

file(GLOB utils RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} python/root*)
foreach(rawUtilName ${utils})
  get_filename_component(utilName ${rawUtilName} NAME)
  if(NOT WIN32)
    # We need the .py only on Windows
    string(REPLACE ".py" "" utilName ${utilName})
  endif()
  set(python ${PYTHON_EXECUTABLE_Development_Main})
  configure_file(${rawUtilName} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${utilName} @ONLY)

  install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${utilName}
                               DESTINATION ${CMAKE_INSTALL_BINDIR}
                               RENAME ${utilName}
                               PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                           GROUP_EXECUTE GROUP_READ
                                           WORLD_EXECUTE WORLD_READ
                               COMPONENT applications)
endforeach()

install(FILES python/cmdLineUtils.py DESTINATION ${runtimedir})
if(IS_ABSOLUTE ${runtimedir})
  set(absruntimedir ${runtimedir})
else()
  set(absruntimedir \${CMAKE_INSTALL_PREFIX}/${runtimedir})
endif()
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE_Development_Main} -m py_compile \$ENV{DESTDIR}${absruntimedir}/cmdLineUtils.py)")
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE_Development_Main} -O -m py_compile \$ENV{DESTDIR}${absruntimedir}/cmdLineUtils.py)")
configure_file(python/cmdLineUtils.py ${localruntimedir}/cmdLineUtils.py @ONLY)


if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLING_CXXFLAGS} -Wno-unused-parameter")
else()
  set_source_files_properties(src/rootcling.cxx PROPERTIES COMPILE_FLAGS "-DNOMINMAX -D_XKEYCHECK_H")
endif()

ROOT_EXECUTABLE(rootcling src/rootcling.cxx LIBRARIES RIO Cling Core)

# rootcling includes the ROOT complex header which would build the complex
# dictionary with modules. To make sure that rootcling_stage1 builds this
# dict before we use it, we add a dependency here.
add_dependencies(rootcling complexDict)

target_include_directories(rootcling PRIVATE
        ${CMAKE_SOURCE_DIR}/core/metacling/res
        ${CMAKE_SOURCE_DIR}/core/dictgen/res
        ${CMAKE_SOURCE_DIR}/io/rootpcm/res)
if(WIN32)
  set_target_properties(rootcling PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS 1)
endif()

# Create aliases: rootcint, genreflex.
if(WIN32)
  add_custom_command(TARGET rootcling POST_BUILD
                     COMMAND copy /y rootcling.exe rootcint.exe
                     COMMAND copy /y rootcling.exe genreflex.exe
                     WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
else()
  add_custom_command(TARGET rootcling POST_BUILD
                     COMMAND ln -f rootcling rootcint
                     COMMAND ln -f rootcling genreflex
                     WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                         "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rootcint;${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/genreflex")

if(CMAKE_HOST_UNIX)
  install(CODE "execute_process(COMMAND ln -f rootcling rootcint WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})" COMPONENT applications)
  install(CODE "execute_process(COMMAND ln -f rootcling genreflex WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})" COMPONENT applications)
else()
  if(MSVC)
    install(PROGRAMS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rootcling.exe DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
    install(PROGRAMS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rootcint.exe DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
    install(PROGRAMS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/genreflex.exe DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
  else()
    install(PROGRAMS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rootcint
                      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/genreflex
                      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rlibmap
                      DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
  endif()
endif()
