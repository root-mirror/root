# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# CMakeLists.txt for the ROOT test programs.
# This shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Author: Pere Mato, 25/10/2010
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(test)

# Sergey: make no sence while CMakeLists.txt file cannot be used separately from ROOT
# but variables like ROOT_asimage_FOUND used here and produced in ROOTConfig.cmake
find_package(ROOT REQUIRED)

#---Copy the CTestCustom.cmake file into the build directory---------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
enable_testing()

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  # To avoid to set ld_Library_path to locate the test libraries
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,@loader_path/.")
endif()

if(DEFINED ROOT_SOURCE_DIR)
  # Testing using the binary tree
  set(ROOT_root_CMD ${ROOT_BINDIR}/root.exe)
else()
  # Testing using an installation (assuming access to ROOT CMake modules)
  include_directories(${ROOT_INCLUDE_DIRS})
  include_directories(${ROOT_INCLUDE_DIRS}/../tutorials)
  add_definitions(${ROOT_DEFINITIONS})
  include(RootMacros)
endif()

#---environment-------------------------------------------------------------------------------
ROOT_ADD_TEST(show-environment COMMAND ${CMAKE_COMMAND} -E environment)

#---hworld------------------------------------------------------------------------------------
ROOT_EXECUTABLE(hworld hworld.cxx LIBRARIES Gpad)

#---event-------------------------------------------------------------------------------------
ROOT_STANDARD_LIBRARY_PACKAGE(Event
                              NO_INSTALL_HEADERS
                              HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Event.h
                              SOURCES Event.cxx LINKDEF EventLinkDef.h
                              DEPENDENCIES Hist MathCore)
if(MSVC)
  if(NOT runtime_cxxmodules)
    add_custom_command(TARGET Event POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/libEvent_rdict.pcm
                                       ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/libEvent_rdict.pcm)
  else()
    add_custom_command(TARGET Event POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/Event.pcm
                                     ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/Event.pcm)
  endif()                                                           
endif()
ROOT_EXECUTABLE(eventexe MainEvent.cxx LIBRARIES Event RIO Tree TreePlayer Hist Net)
ROOT_ADD_TEST(test-event COMMAND eventexe)

#---guitest------------------------------------------------------------------------------------
ROOT_EXECUTABLE(guitest guitest.cxx LIBRARIES RIO Gui Gpad Hist)

#---hsimple------------------------------------------------------------------------------------
#ROOT_EXECUTABLE(hsimple hsimple.cxx LIBRARIES RIO Tree Hist)
#ROOT_ADD_TEST(test-hsimple COMMAND hsimple)

#---invalid pointer checks---------------------------------------------------------------------
ROOT_ADD_TEST(test-check-nullptr
  COMMAND ${ROOT_root_CMD} -b -q -l -e "int*ptr=(int*)0" -e "*ptr"
  PASSREGEX "dereference null pointer")

if(NOT MSVC OR win_broken_tests)
  ROOT_ADD_TEST(test-check-invalidptr
    COMMAND ${ROOT_root_CMD} -b -q -l -e "int*ptr=(int*)1" -e "*ptr"
    PASSREGEX "invalid memory address")
endif()

#--minexam-------------------------------------------------------------------------------------
ROOT_EXECUTABLE(minexam minexam.cxx LIBRARIES Minuit)
ROOT_ADD_TEST(test-minexam COMMAND minexam)

#--tstring-------------------------------------------------------------------------------------
ROOT_EXECUTABLE(tstring tstring.cxx LIBRARIES Core)
configure_file(tstring.cxx tstring.cxx COPYONLY)
ROOT_ADD_TEST(test-tstring COMMAND tstring)

#--tcollex-------------------------------------------------------------------------------------
ROOT_EXECUTABLE(tcollex tcollex.cxx LIBRARIES Core)
ROOT_ADD_TEST(test-tcollex COMMAND tcollex)

#--tcollbm------------------------------------------------------------------------------------
ROOT_EXECUTABLE(tcollbm tcollbm.cxx LIBRARIES Core MathCore)
ROOT_ADD_TEST(test-tcollbm COMMAND tcollbm 1000 1000000 LABELS longtest)

#--vvector------------------------------------------------------------------------------------
ROOT_EXECUTABLE(vvector vvector.cxx LIBRARIES Core Matrix RIO)
ROOT_ADD_TEST(test-vvector COMMAND vvector)

#--vmatrix------------------------------------------------------------------------------------
ROOT_EXECUTABLE(vmatrix vmatrix.cxx LIBRARIES Core Matrix RIO)
ROOT_ADD_TEST(test-vmatrix COMMAND vmatrix)

#--vlazy------------------------------------------------------------------------------------
ROOT_EXECUTABLE(vlazy vlazy.cxx LIBRARIES Core Matrix)
ROOT_ADD_TEST(test-vlazy COMMAND vlazy)

#--helloso------------------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(G__Hello ${CMAKE_CURRENT_SOURCE_DIR}/Hello.h MODULE Hello DEPENDENCIES Gpad Graf Hist MathCore Matrix)
ROOT_LINKER_LIBRARY(Hello Hello.cxx G__Hello.cxx LIBRARIES Graf Gpad)

#--Aclockso------------------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(G__Aclock ${CMAKE_CURRENT_SOURCE_DIR}/Aclock.h MODULE Aclock DEPENDENCIES Graf Gpad MathCore)
ROOT_LINKER_LIBRARY(Aclock Aclock.cxx G__Aclock.cxx LIBRARIES Graf Gpad)

#--bench------------------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(G__TBench ${CMAKE_CURRENT_SOURCE_DIR}/TBench.h MODULE TBench LINKDEF benchLinkDef.h DEPENDENCIES MathCore Tree)
ROOT_LINKER_LIBRARY(TBench TBench.cxx G__TBench.cxx LIBRARIES Core MathCore RIO Tree)
ROOT_EXECUTABLE(bench bench.cxx LIBRARIES Core TBench)
ROOT_ADD_TEST(test-bench COMMAND bench -s LABELS longtest)

#--stress------------------------------------------------------------------------------------
  ROOT_EXECUTABLE(stress stress.cxx LIBRARIES Event Core Hist RIO Tree Gpad Postscript)
  ROOT_ADD_TEST(test-stress COMMAND stress -b FAILREGEX "FAILED|Error in"
                            DEPENDS test-event LABELS longtest)

#--stressShapes------------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressShapes stressShapes.cxx LIBRARIES  Geom Tree GenVector Gpad)
ROOT_ADD_TEST(test-stressshapes COMMAND stressShapes -b FAILREGEX "FAILED|Error in")
ROOT_ADD_TEST(test-stressshapes-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressShapes.cxx
              FAILREGEX "FAILED|Error in"
)

#--stressGeometry------------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressGeometry stressGeometry.cxx LIBRARIES Geom Tree GenVector Gpad)
ROOT_ADD_TEST(test-stressgeometry COMMAND stressGeometry -b FAILREGEX "FAILED|Error in" LABELS longtest)
ROOT_ADD_TEST(test-stressgeometry-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressGeometry.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressgeometry LABELS longtest)

#--stressLinear------------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressLinear stressLinear.cxx LIBRARIES Matrix Hist RIO)
ROOT_ADD_TEST(test-stresslinear COMMAND stressLinear FAILREGEX "FAILED|Error in" LABELS longtest)
ROOT_ADD_TEST(test-stresslinear-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressLinear.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stresslinear LABELS longtest)

#--stressGraphics------------------------------------------------------------------------------------
if(ROOT_opengl_FOUND)
  ROOT_EXECUTABLE(stressGraphics stressGraphics.cxx LIBRARIES Graf Gpad Postscript)
  if(MSVC)
    set_property(TARGET stressGraphics APPEND_STRING PROPERTY LINK_FLAGS "/STACK:4000000")
  endif()
  # with ZLIB-CF(builtin_zlib) we have better compression ratio.
  if(ZLIB_CF)
    configure_file(stressGraphics_builtinzlib.ref stressGraphics_builtinzlib.ref COPYONLY)
  else()
    configure_file(stressGraphics.ref stressGraphics.ref COPYONLY)
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../tutorials/graphics/earth.dat earth.dat COPYONLY)
  ROOT_ADD_TEST(test-stressgraphics ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH} COMMAND stressGraphics -b -k FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stressgraphics-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressGraphics.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressgraphics )
endif()

#--stressHistogram------------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressHistogram stressHistogram.cxx LIBRARIES Hist RIO)
ROOT_ADD_TEST(test-stresshistogram COMMAND stressHistogram FAILREGEX "FAILED|Error in" LABELS longtest)
ROOT_ADD_TEST(test-stresshistogram-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressHistogram.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stresshistogram )

#--stressGUI---------------------------------------------------------------------------------------
if(ROOT_asimage_FOUND)
  ROOT_EXECUTABLE(stressGUI stressGUI.cxx LIBRARIES Gui Recorder GuiHtml ASImageGui)
endif()

#--stressSpectrum----------------------------------------------------------------------------------
if (spectrum)
  ROOT_EXECUTABLE(stressSpectrum stressSpectrum.cxx LIBRARIES Hist Spectrum Gpad)
  ROOT_ADD_TEST(test-stressspectrum COMMAND stressSpectrum -b FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stressspectrum-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressSpectrum.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stressspectrum LABELS longtest)
endif()

#--stressVector------------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressVector stressVector.cxx LIBRARIES Physics GenVector)
ROOT_ADD_TEST(test-stressvector COMMAND stressVector FAILREGEX "FAILED|Error in")
ROOT_ADD_TEST(test-stressvector-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressVector.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressvector)

#--stressTMVA--------------------------------------------------------------------------------------
if(CUDA_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDNNCUDA")
endif()

if(BLAS_FOUND AND imt)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDNNCPU")
endif()

if(ROOT_tmva_FOUND)
  ROOT_EXECUTABLE(stressTMVA stressTMVA.cxx LIBRARIES TMVA)
  ROOT_ADD_TEST(test-stresstmva COMMAND stressTMVA -b LABELS longtest TIMEOUT 1800)
  ROOT_ADD_TEST(test-stresstmva-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressTMVA.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stresstmva TIMEOUT 1800)
endif()

#--stressMathMore----------------------------------------------------------------------------------
if(ROOT_mathmore_FOUND)
  ROOT_EXECUTABLE(stressMathMore stressMathMore.cxx LIBRARIES MathMore Smatrix)
  ROOT_ADD_TEST(test-stressmathmore COMMAND stressMathMore FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stressmathmore-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressMathMore.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stressmathmore TIMEOUT 1800)
endif()

#--stressMathCore----------------------------------------------------------------------------------
ROOT_STANDARD_LIBRARY_PACKAGE(TrackMathCoreDict
                              NO_SOURCES
                              NO_INSTALL_HEADERS
                              HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/TrackMathCore.h
                              LINKDEF TrackMathCoreLinkDef.h
                              DEPENDENCIES Core MathCore RIO GenVector Smatrix)
if(MSVC)
  add_custom_command(TARGET TrackMathCoreDict POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/libTrackMathCoreDict.dll
                                     ${CMAKE_CURRENT_BINARY_DIR}/libTrackMathCoreDict.dll)
endif()
ROOT_EXECUTABLE(stressMathCore stressMathCore.cxx LIBRARIES MathCore Hist RIO Tree GenVector)
ROOT_ADD_TEST(test-stressmathcore COMMAND stressMathCore FAILREGEX "FAILED|Error in" LABELS longtest)
ROOT_ADD_TEST(test-stressmathcore-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressMathCore.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressmathcore)

#--stressRooFit----------------------------------------------------------------------------------
if(ROOT_roofit_FOUND)
  if(ROOT_mathmore_FOUND)
    ROOT_EXECUTABLE(stressRooFit stressRooFit.cxx LIBRARIES RooFit RooFitMore)
  else()
    ROOT_EXECUTABLE(stressRooFit stressRooFit.cxx LIBRARIES RooFit )
  endif()  
  configure_file(stressRooFit_ref.root stressRooFit_ref.root COPYONLY)
  ROOT_ADD_TEST(test-stressroofit COMMAND stressRooFit FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stressroofit-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressRooFit.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stressroofit )

  add_subdirectory(fit)
endif()

#--stressRooStats----------------------------------------------------------------------------------
if(ROOT_roofit_FOUND)
  ROOT_EXECUTABLE(stressRooStats stressRooStats.cxx LIBRARIES RooStats)
  configure_file(stressRooStats_ref.root stressRooStats_ref.root COPYONLY)
  ROOT_ADD_TEST(test-stressroostats COMMAND stressRooStats FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stressroostats-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressRooStats.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stressroostats )
  if (ROOT_minuit2_FOUND)
    ROOT_ADD_TEST(test-stressroostats-minuit2 COMMAND stressRooStats -minim Minuit2  FAILREGEX "FAILED|Error in" LABELS longtest)
  endif()
endif()

#--stressHistFactory--------------------------------------------------------------------------------
if(ROOT_roofit_FOUND AND GSL_FOUND AND ROOT_xml_FOUND)
  ROOT_EXECUTABLE(stressHistFactory stressHistFactory.cxx LIBRARIES RooStats HistFactory XMLParser)
  if(NOT MSVC)
    configure_file(HistFactoryTest.tar HistFactoryTest.tar COPYONLY)
    configure_file(stressHistFactory_ref.root stressHistFactory_ref.root COPYONLY)
    ROOT_ADD_TEST(test-stressHistFactory ENVIRONMENT ROOTSYS=${CMAKE_BINARY_DIR} COMMAND stressHistFactory FAILREGEX "FAILED|Error in" LABELS longtest)
    ROOT_ADD_TEST(test-stressHistFactory-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressHistFactory.cxx FAILREGEX "FAILED|Error in" DEPENDS test-stressHistFactory)
  endif()
endif()

#--stressFit---------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressFit stressFit.cxx LIBRARIES MathCore Matrix)
ROOT_ADD_TEST(test-stressfit COMMAND stressFit FAILREGEX "FAILED|Error in" LABELS longtest)
ROOT_ADD_TEST(test-stressfit-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressFit.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressfit)

#--stressHistoFit-----------------------------------------------------------------------------
if(ROOT_unuran_FOUND)
  ROOT_EXECUTABLE(stressHistoFit stressHistoFit.cxx LIBRARIES MathCore Matrix Unuran Tree Gpad)
  ROOT_ADD_TEST(test-stresshistofit COMMAND stressHistoFit FAILREGEX "FAILED|Error in" LABELS longtest)
  ROOT_ADD_TEST(test-stresshistofit-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressHistoFit.cxx
                FAILREGEX "FAILED|Error in" DEPENDS test-stresshistofit )
endif()

#--stressEntryList---------------------------------------------------------------------------
ROOT_EXECUTABLE(stressEntryList stressEntryList.cxx LIBRARIES MathCore Tree Hist)
ROOT_ADD_TEST(test-stressentrylist COMMAND stressEntryList -b FAILREGEX "FAILED|Error in")
ROOT_ADD_TEST(test-stressentrylist-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressEntryList.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressentrylist)

#--stressIterators---------------------------------------------------------------------------
ROOT_EXECUTABLE(stressIterators stressIterators.cxx LIBRARIES Core)
ROOT_ADD_TEST(test-stressiterators COMMAND stressIterators FAILREGEX "FAILED|Error in")
ROOT_ADD_TEST(test-stressiterators-interpreted COMMAND ${ROOT_root_CMD} -b -q -l ${CMAKE_CURRENT_SOURCE_DIR}/stressIterators.cxx
              FAILREGEX "FAILED|Error in" DEPENDS test-stressiterators)

#--stressInterpreter-------------------------------------------------------------------------
ROOT_EXECUTABLE(stressInterpreter stressInterpreter.cxx LIBRARIES Core)
if(WIN32)
  set(cleantmpfiles POSTCMD cmd /c del AutoDict_*)
endif()
ROOT_ADD_TEST(test-stressinterpreter COMMAND stressInterpreter
              PRECMD ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/stressInterpreter.cxx
                                              $<TARGET_FILE_DIR:stressInterpreter>/stressInterpreter.cxx
              ${cleantmpfiles}
              FAILREGEX Error|FAILED)

#--stressHepix-------------------------------------------------------------------------------
ROOT_EXECUTABLE(stressHepix stressHepix.cxx LIBRARIES Core)
#ROOT_ADD_TEST(test-stressHepix COMMAND stressHepix FAILREGEX "FAILED|Error in")

#--stressProof-------------------------------------------------------------------------------
if(proof AND NOT WIN32)
  add_custom_target(TestData COMMAND ${CMAKE_COMMAND} -DDST=${CMAKE_SOURCE_DIR}/files -P ${CMAKE_CURRENT_SOURCE_DIR}/rootDownloadData.cmake)
  ROOT_EXECUTABLE(stressProof stressProof.cxx LIBRARIES Proof ProofPlayer Hist)
  ROOT_ADD_TEST(test-stressproof PRECMD ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target TestData
                                 COMMAND stressProof lite://
                                         -h1 ${CMAKE_SOURCE_DIR}/files/h1 -event ${CMAKE_SOURCE_DIR}/files/event
                                         -l /tmp/stressProof-%d.log -cleanlog -catlog -noprogress
                                         ENVIRONMENT PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
                                         FAILREGEX "FAILED|Error in")
endif()

#--testbits----------------------------------------------------------------------------------
ROOT_EXECUTABLE(testbits testbits.cxx LIBRARIES Core)
ROOT_ADD_TEST(test-testbits COMMAND testbits)

#--ctorture----------------------------------------------------------------------------------
ROOT_EXECUTABLE(ctorture ctorture.cxx LIBRARIES MathCore)
ROOT_ADD_TEST(test-ctorture COMMAND ctorture)

#--qprandom----------------------------------------------------------------------------------
ROOT_EXECUTABLE(QpRandomDriver QpRandomDriver.cxx LIBRARIES Quadp)
ROOT_ADD_TEST(test-qprandom COMMAND QpRandomDriver 10 2 2)

#--threads----------------------------------------------------------------------------------
ROOT_EXECUTABLE(threads threads.cxx LIBRARIES Thread Hist Gpad)
#ROOT_ADD_TEST(test-threads COMMAND threads)

#--stressIOPlugins--------------------------------------------------------------------------
ROOT_EXECUTABLE(stressIOPlugins stressIOPlugins.cxx LIBRARIES Event Core Hist RIO Tree Gpad Postscript)
if(ROOT_xrootd_FOUND)
  ROOT_ADD_TEST(test-stressIOPlugins-xroot COMMAND stressIOPlugins root FAILREGEX "FAILED|Error in")
endif()
if(ROOT_davix_FOUND)
  ROOT_ADD_TEST(test-stressIOPlugins-http COMMAND stressIOPlugins http FAILREGEX "FAILED|Error in")
endif()

#--delaunay----------------------------------------------------------------------------------
ROOT_EXECUTABLE(delaunayTriangulation delaunayTriangulation.cxx LIBRARIES Hist)
ROOT_ADD_TEST(test-delaunay COMMAND delaunayTriangulation)

#--TestTformula------------------------------------------------------------------------------------
ROOT_EXECUTABLE(TFormulaTests TFormulaTests.cxx LIBRARIES Hist RIO)
if(MSVC)
  set_property(TARGET TFormulaTests APPEND_STRING PROPERTY LINK_FLAGS "/STACK:4000000")
endif()
ROOT_ADD_TEST(test-TFormulaTests COMMAND TFormulaTests FAILREGEX "FAILED|Error in")

#--TBB basic test----------------------------------------------------------------------------------
if(ROOT_imt_FOUND)
  include_directories(SYSTEM ${TBB_INCLUDE_DIRS})
  ROOT_EXECUTABLE(testTBB testTBB.cxx LIBRARIES ${TBB_LIBRARIES} BUILTINS TBB)
  set_target_properties(testTBB PROPERTIES COMPILE_FLAGS "${TBB_CXXFLAGS}")
  ROOT_ADD_TEST(test-TBB COMMAND testTBB FAILREGEX "FAILED|Error in")
endif()

#--Vc basic test-----------------------------------------------------------------------------------
if(ROOT_vc_FOUND)
  ROOT_EXECUTABLE(testVc testVc.cxx LIBRARIES ${Vc_LIBRARIES} BUILTINS Vc)
  target_include_directories(testVc SYSTEM BEFORE PRIVATE ${Vc_INCLUDE_DIRS})
  ROOT_ADD_TEST(test-Vc COMMAND testVc FAILREGEX "FAILED|Error in")
endif()

#--VecCore basic test------------------------------------------------------------------------------
if(ROOT_veccore_FOUND)
  ROOT_EXECUTABLE(test-veccore test-veccore.cxx LIBRARIES ${VecCore_LIBRARIES} BUILTINS VECCORE)
  target_include_directories(test-veccore BEFORE PRIVATE ${VecCore_INCLUDE_DIRS})
  target_compile_definitions(test-veccore PRIVATE ${VecCore_DEFINITIONS})
  if(VecCore_Vc_FOUND)
    ROOT_ADD_TEST(VecCore COMMAND test-veccore REGEX "Vc")
  else()
    ROOT_ADD_TEST(VecCore COMMAND test-veccore REGEX "Scalar")
  endif()
endif()

#--Vc GenVector test-----------------------------------------------------------------------------------
if(ROOT_vc_FOUND)
  ROOT_EXECUTABLE(testGenVectorVc testGenVectorVc.cxx LIBRARIES Physics GenVector ${Vc_LIBRARIES} BUILTINS Vc)
  target_include_directories(testGenVectorVc SYSTEM BEFORE PRIVATE ${Vc_INCLUDE_DIRS})
  ROOT_ADD_TEST(test-GenVector-Vc COMMAND testGenVectorVc FAILREGEX "FAILED|Error in")
endif()

#--g2root------------------------------------------------------------------------------------------
if(TARGET g2root)
  ROOT_ADD_TEST(test-g2root COMMAND g2root)
endif()

#--periodic----------------------------------------------------------------------------------------
if (NOT cxxmodules)
  # There are two problems with this test and cxxmodules. First, it picks up a
  # header files from $ROOTSYS/include and builds a library. Naturally, we should
  # build pcm files, however, building module MathCore requires specific -I to
  # build the needed VecCore and Vc components. This test does not have access
  # to the include relevant include paths. Secondly, if we want to reuse the
  # modules from ROOT which makes most sense we can't because we get a hard
  # error such as: fatal error: malformed or corrupted AST file:
  # 'SourceLocation remap refers to unknown module, cannot find include/pcms/1WYSNQV9VBZK7/stl-2OZGQN92C38MI.pcm
  #
  # FIXME: We can fix the first point by moving out all VecCore-related headers
  # such as Math/Types.h and all of its includers in a separate module. Thus,
  # the current test will not require the VecCore (as in the textual case).
  # Alternatively, we can trace the origin of the fatal error and try to remap
  # the source locations.
  if(MSVC)
    set(build_generator_args --build-generator ${CMAKE_GENERATOR})
    if(CMAKE_GENERATOR_PLATFORM)
      list(APPEND build_generator_args --build-generator-platform ${CMAKE_GENERATOR_PLATFORM})
    endif()
    if(CMAKE_GENERATOR_TOOLSET)
      list(APPEND build_generator_args --build-generator-toolset ${CMAKE_GENERATOR_TOOLSET})
    endif()
    ROOT_ADD_TEST(test-periodic-build
      COMMAND ${CMAKE_CTEST_COMMAND} ${build_generator_args}
        --build-and-test ${CMAKE_CURRENT_SOURCE_DIR}/periodic periodic-build
    )
  else()
    ROOT_ADD_TEST(test-periodic-build
      COMMAND
        env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
        ${CMAKE_CTEST_COMMAND} --build-generator ${CMAKE_GENERATOR}
        --build-and-test ${CMAKE_CURRENT_SOURCE_DIR}/periodic periodic-build
    )
  endif()
endif()

#--canary tests------------------------------------------------------------------------------------
if(asserts AND NOT MSVC)
  ROOT_EXECUTABLE(checkAssertsNDEBUG checkAssertsNDEBUG.cxx LIBRARIES Core)
  ROOT_ADD_TEST(checkAssertsNDEBUG COMMAND checkAssertsNDEBUG)
endif()
