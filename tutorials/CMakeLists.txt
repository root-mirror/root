# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# CMakeLists.txt for the ROOT tutorials programs.
# Author: Pere Mato, 25/10/2010
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(tutorials)

# Sergey: make no sence while CMakeLists.txt file cannot be used separately from ROOT
# but variables like ROOT_asimage_FOUND used here and produced in ROOTConfig.cmake
find_package(ROOT REQUIRED)

if(DEFINED ROOT_SOURCE_DIR)  # Testing using the binary tree
  set(ROOT_root_CMD root.exe)
  if(NOT MSVC)  # Ignore environment on Windows
    set(ROOT_environ PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
                     ${ld_library_path}=${CMAKE_BINARY_DIR}/lib:$ENV{${ld_library_path}}
                     ROOTSYS=${CMAKE_BINARY_DIR}
                     PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
  else()
    set(ROOT_environ ROOTSYS=${CMAKE_BINARY_DIR}
                     PYTHONPATH=${CMAKE_BINARY_DIR}/bin;$ENV{PYTHONPATH})
  endif()
else()                       # testing using an installation
  include(${ROOT_USE_FILE})
  if(DEFINED ROOT_CONFIG_EXECUTABLE) #---If ROOT was built with the classic configure/make---
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules)
    include(RootMacros)
    set(ROOT_root_CMD root.exe)
  endif()
  enable_testing()
endif()

#---Copy the CTestCustom.cmake file into the build directory--------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

#---Provide a rootlogon.C file in the current build directory that
#   will affect the way we run all tutorials.
#   This overwrites the existing rootlogon.C and rootalias.C in the
#   tutorials directory which is copied to the build area.
#-------------------------------------------------------------------
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/rootlogon.C "{
  // Needed by ACLiC to use the current dicrectory for scratch area
  gSystem->SetBuildDir(\".\", kTRUE);
  gROOT->SetWebDisplay(\"batch\");
}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/rootalias.C "")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/rootlogoff.C "{}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/.rootrc "
Proof.Sandbox: /tmp/proof
Rint.History:  .root_hist
")

#---Definition of the helper function--------------------------------
function(ROOT_ADD_TUTORIAL macrofile rc)
  string(REPLACE ".C" "" name ${macrofile})
  string(REPLACE "/" "-" name ${name})
  ROOT_ADD_TEST(tutorial-${name} COMMAND ${ROOT_root_CMD} -b -l -n -q ${CMAKE_CURRENT_SOURCE_DIR}/${macrofile}
                PASSRC ${rc} FAILREGEX "Error in" "error:" LABELS tutorial)
endfunction()

#---Tutorials disabled depending on the build components-------------
if(NOT ROOT_minuit2_FOUND)
  set(minuit2_veto fit/fit2dHist.C fit/fit2dHist.C
                   fit/fitCircle.C fit/minuit2FitBench2D.C
                   fit/minuit2FitBench2D.C fit/minuit2FitBench.C
                   fit/minuit2FitBench.C fit/minuit2GausFit.C
                   fit/minuit2GausFit.C fit/NumericalMinimization.C
                   fit/combinedFit.C fit/TestBinomial.C
                   fit/fitNormSum.C fit/vectorizedFit.C
                   tutorials/roostats/rs_bernsteinCorrection.C)
endif()

if (NOT dataframe)
    # RDataFrame
    list(APPEND dataframe_veto dataframe/*.C dataframe/*.py)
    # RDataFrame tutorial in graphs
    list(APPEND dataframe_veto graphs/timeSeriesFromCSV_TDF.C)
    # TMVA tutorials dependent on RDataFrame
    list(APPEND dataframe_veto tmva/tmva*.C)
elseif(MSVC AND NOT win_broken_tests)
    list(APPEND dataframe_veto dataframe/*.py)
endif()

if(NOT sqlite)
    # RDF+SQlite tutorials
    list(APPEND dataframe_veto dataframe/*SQlite*)
endif()

if(MACOSX_VERSION VERSION_EQUAL 10.13)
   list(APPEND dataframe_veto dataframe/df103_NanoAODHiggsAnalysis.*)
endif()

if(NOT ROOT_mlp_FOUND)
   set(mlp_veto legacy/mlp/*.C)
endif()

if(NOT ROOT_spectrum_FOUND)
   set(spectrum_veto spectrum/*.C)
endif()

if(NOT ROOT_roofit_FOUND)
  set(roofit_veto  fit/RoofitDemo.C
                   roofit/*.C roostats/*.C histfactory/*.C
                   roofit/*.py)
else()
  if(MSVC AND NOT win_broken_tests)
    set(roofit_veto roofit/rf104_classfactory.C
                    roofit/rf104_classfactory.py
                    roofit/rf401_importttreethx.C
                    roofit/rf512_wsfactory_oper.C
    )
  endif()
endif()

if(NOT ROOT_unuran_FOUND)
  set(unuran_veto  math/testrandom.C unuran/unuranDemo.C unuran/unuranFoamTest.C
                   math/multidimSampling.C)
endif()

if(NOT ROOT_xml_FOUND)
  set(xml_veto  xml/*.C
                histfactory/*.C    # histfactory requires xml
                unfold/*.C)        # unfold requires xml
endif()

if(NOT ROOT_mpi_FOUND)
  set(mpi_veto io/testTMPIFile.C)
endif()

if(NOT xrootd)
  set(xrootd_veto dataframe/df101_h1Analysis.C
                  dataframe/df102_NanoAODDimuonAnalysis.C
                  dataframe/df103_NanoAODHiggsAnalysis.C
                  tmva/tmva103_Application.C
                  dataframe/df102_NanoAODDimuonAnalysis.py
                  dataframe/df103_NanoAODHiggsAnalysis.py
                  dataframe/df104_HiggsToTwoPhotons.py
                  dataframe/df105_WBosonAnalysis.py
                  dataframe/df106_HiggsToFourLeptons.py
                  dataframe/df107_SingleTopAnalysis.py
                  )
endif()

# variables identifying the package must have the package name  in lower case (it corresponds to the CMake option name)
if(NOT ROOT_r_FOUND)
  set(r_veto  r/*.C)
endif()

if(WIN32)
  set(histfactory_veto histfactory/*.C
                       roostats/StandardFrequentistDiscovery.C) # histfactory doesn't work on Windows
else()
  set(histfactory_veto histfactory/makeExample.C)
endif()


if(NOT ROOT_fitsio_FOUND)
  set(fitsio_veto  fitsio/*.C)
endif()

if(NOT ROOT_mathmore_FOUND)
  set(mathmore_veto
      math/quasirandom.C
      math/exampleMultiRoot.C
      math/Bessel.C
      math/LegendreAssoc.C
      math/Legendre.C
      math/mathmoreIntegration.C
      math/tStudent.C
      math/normalDist.C
      roostats/TestNonCentral.C
      math/Legendre.py
      math/Bessel.py
      math/tStudent.py)
endif()

if(NOT ROOT_fftw3_FOUND)
  set(fftw3_veto roofit/rf208_convolution.C
                 roofit/rf210_angularconv.C
                 roofit/rf211_paramconv.C
                 roofit/rf512_wsfactory_oper.C
                 fft/FFT.C
                 fit/fitConvolution.C)
endif()

if(NOT ROOT_opengl_FOUND)
  set(opengl_veto tree/staff.C
                  gl/*.C)
endif()

if(NOT GRAPHVIZ_FOUND)
  set(gviz_veto graphs/graphstruct.C)
endif()

if(NOT TBB_FOUND AND NOT builtin_tbb)
  set(tbb_veto  multicore/mtbb*.C)
endif()

if(NOT ROOT_imt_FOUND)
  set(imt_veto multicore/imt*.C multicore/mt*.C)
else()
  if(MSVC)
    #---Multiproc is not supported on Windows
    set(imt_veto ${imt_veto} multicore/mp*.C multicore/mtbb201_parallelHistoFill.C)
  endif()
endif()

if(ROOT_CLASSIC_BUILD)
  set(classic_veto multicore/mp104_*.C multicore/mp105_*.C)
endif()

#---These ones requires a display to run-----------------------------
set(gui_veto fit/fitpanel_playback.C
             cocoa/*.C
             geom/building.C geom/cheongwadae.C geom/geom*.C geom/lego.C geom/robot.C geom/south_gate.C geom/station*.C geom/tank.C geom/webdemo.C
             gl/glViewerExercise.C gl/glViewerLOD.C gl/gviz3d.C gl/nucleus.C gl/viewer3DLocal.C gl/viewer3DMaster.C
             gui/*.C
             hist/exec1.C
             hist/exec2.C
             hist/tprofile2polyRealistic.C
             hist/tprofile2polyRealisticModuleError.C
             image/*.C
             graphics/psview.C graphics/gtime.C
             graphics/graph_edit_playback.C
             roostats/ModelInspector.C
             tree/tvdemo.C
             eve/*.C
             webgui/panel/server.cxx webgui/webwindow/server.cxx)

if (NOT ROOT_tmva_FOUND)
  list(APPEND tmva_veto tmva/*.C tmva/*.py tmva/envelope/*.C tmva/keras/*.C tmva/keras/*.py tmva/pytorch/*.py )
else()
  #---These do not need to run for TMVA
  list(APPEND tmva_veto tmva/createData.C)
  if(MSVC AND NOT win_broken_tests)
    list(APPEND tmva_veto tmva/tmva103_Application.C tmva/envelope/classification.C)
  endif()
endif()

if (NOT ROOT_pythia6_FOUND)
  set(pythia_veto pythia/pythiaExample.C)
endif()
if (NOT ROOT_pythia8_FOUND)
  set(pythia_veto ${pythia_veto} pythia/pythia8.C)
else()
  if("$ENV{PYTHIA8}" STREQUAL "")
    get_filename_component(pythia8dir "${PYTHIA8_INCLUDE_DIR}" DIRECTORY)
    list(APPEND ROOT_environ PYTHIA8=${pythia8dir})
  endif()
  if("$ENV{PYTHIA8DATA}" STREQUAL "" AND PYTHIA8_DATA)
    list(APPEND ROOT_environ PYTHIA8DATA=${PYTHIA8_DATA})
  endif()
endif()

if(root7)
  set(root7_veto dataframe/df013_InspectAnalysis.C
                 v7/browser.cxx
                 v7/filedialog.cxx
                 v7/fitpanel.cxx
                 v7/fitpanel6.cxx
      )
  if(NOT davix)
    list(APPEND root7_veto v7/ntuple/ntpl003_lhcbOpenData.C)
    list(APPEND root7_veto v7/ntuple/ntpl004_dimuon.C)
  endif()
  if(NOT dataframe)
    list(APPEND root7_veto v7/ntuple/ntpl004_dimuon.C)
  endif()
else()
  file(GLOB v7_veto_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/ v7/*.cxx v7/*/*.cxx v7/*.C v7/*/*.C)
  list(APPEND root7_veto ${v7_veto_files})
endif()

if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
  set(bits32_veto dataframe/*.C graphs/timeSeriesFrom*.C v7/ntuple/ntpl004_dimuon.C)
endif()

#---These ones are disabled !!! ------------------------------------
set(extra_veto
  legacy/benchmarks.C
  legacy/htmlex.C
  legacy/rootalias.C          # Helper macro
  rootlogon.C          # Helper macro
  rootlogoff.C         # Helper macro
  legacy/rootmarks.C          # Instrumentation. Not a standalone tutorial
  legacy/memstat/memstatExample.C   # Instrumentation. Not a standalone tutorial
  multicore/mp_H1_lambdas.C # not a tutorial; used by mp104_processH1.C et al.
  html/*.C
  net/*.C
  proof/*.C
  sql/*.C
  tree/hsimpleProxy.C # A driver uses this macro which cannot be executed directly
  tree/tree0.C
  tree/tree2a.C
  tree/tree4.C
  roostats/rs401d_FeldmanCousins.C  # Takes too much time
  histfactory/ModifyInterpolation.C
  tree/copytree2.C
  tree/copytree3.C
  tree/copytree.C
  tree/h1analysis*.C # these are not a tutorial but classes used in run_h1analysis.C
  tree/h1chain.C
  http/*.C
  eve7/*.C
  r/rootlogon.C)

set(all_veto hsimple.C
             geom/geometry.C
             ${extra_veto}
             ${gui_veto}
             ${minuit2_veto}
             ${roofit_veto}
             ${unuran_veto}
             ${xml_veto}
             ${mpi_veto}
             ${fitsio_veto}
             ${tmva_veto}
             ${mathmore_veto}
             ${fftw3_veto}
             ${opengl_veto}
             ${gviz_veto}
             ${r_veto}
             ${runtime_cxxmodules_veto}
             ${histfactory_veto}
             ${tbb_veto}
             ${imt_veto}
             ${classic_veto}
             ${pythia_veto}
             ${root7_veto}
             ${bits32_veto}
             ${xrootd_veto}
             ${mlp_veto}
             ${spectrum_veto}
             ${dataframe_veto}
             )

file(GLOB_RECURSE tutorials RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.C)
if(root7 AND webgui)
  file(GLOB_RECURSE tutorials_v7 RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} v7/*.cxx)
  list(APPEND tutorials ${tutorials_v7})
endif()
file(GLOB tutorials_veto RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${all_veto})

list(LENGTH tutorials nTotal)
list(REMOVE_ITEM tutorials ${tutorials_veto})
list(LENGTH tutorials nAfterVeto)
message(STATUS "${nAfterVeto}/${nTotal} C++ tutorials have been activated.")


if(mpi)
  set (temp_list ${tutorials})
  list(FILTER tutorials INCLUDE REGEX "MPI")
  set(mpi_tutorials ${tutorials})
  set(tutorials ${temp_list})
  list(REMOVE_ITEM tutorials ${mpi_tutorials})
endif()

#---Special return code------------------------------------------------
set(returncode_1 fit/fit2a.C fit/graph2dfit.C
                 graphics/arrow.C
                 graphics/crown.C graphics/diamond.C
                 graphics/earth.C graphics/ellipse.C
                 graphics/pavetext.C
                 graphics/tmathtext.C graphics/tmathtext2.C
                 graphs/timeonaxis.C
                 graphs/timeonaxis2.C
                 graphs/timeonaxis3.C
                 graphs/exclusiongraph.C
                 graphs/graphstruct.C
                 hist/ContourList.C
                 hist/hstack.C
                 hist/hbars.C
                 hist/th2polyBoxes.C
                 hist/statsEditing.C
                 hist/cumulative.C
                 hist/hlabels1.C
                 hist/hlabels2.C
                 tree/h1analysis.C
                 math/chi2test.C
                 r/SimpleFitting.C)
#---Dependencies------------------------------------------------------
set(unfold-testUnfold5d-depends tutorial-unfold-testUnfold5c)
set(unfold-testUnfold5c-depends tutorial-unfold-testUnfold5b)
set(unfold-testUnfold5b-depends tutorial-unfold-testUnfold5a)
set(unfold-testUnfold7d-depends tutorial-unfold-testUnfold7c)
set(unfold-testUnfold7c-depends tutorial-unfold-testUnfold7b)
set(unfold-testUnfold7b-depends tutorial-unfold-testUnfold7a)
set(xml-xmlmodifyfile-depends tutorial-xml-xmlnewfile)
set(xml-xmlreadfile-depends tutorial-xml-xmlnewfile)
set(roofit-rf503_wspaceread-depends tutorial-roofit-rf502_wspacewrite)
set(histfactory-example-depends tutorial-roostats-CreateExampleFile)
set(io-readCode-depends tutorial-io-importCode)
set(fit-fit1-depends tutorial-hist-fillrandom)
set(fit-myfit-depends tutorial-fit-fitslicesy)
set(foam-foam_demopers-depends tutorial-foam-foam_demo)
set(tree-staff-depends  tutorial-tree-cernbuild)
set(tree-cernstaff-depends  tutorial-tree-cernbuild)
set(hist-hbars-depends  tutorial-tree-cernbuild)
set(benchmarks-depends tutorial-hsimple
                       tutorial-fit-fit1
                       tutorial-fit-myfit
                       tutorial-hist-h1draw
                       tutorial-hist-FirstContour
                       tutorial-geom-na49view
                       tutorial-tree-ntuple1
                       tutorial-tree-spider
                       tutorial-io-hadd
                       tutorial-io-loopdir
                       tutorial-io-copyFiles)
set(geom-na49view-depends tutorial-geom-geometry)
set(multicore-mt102_readNtuplesFillHistosAndFit-depends tutorial-multicore-mt101_fillNtuples)
set(multicore-mp102_readNtuplesFillHistosAndFit-depends tutorial-multicore-mp101_fillNtuples)
set(multicore-mp105_processEntryList-depends tutorial-multicore-mp104_processH1)

#--many roostats tutorials depending on having creating the file first with histfactory and example_combined_GaussExample_model.root
foreach(tname  ModelInspector OneSidedFrequentistUpperLimitWithBands StandardBayesianMCMCDemo StandardBayesianNumericalDemo
               StandardFeldmanCousinsDemo  StandardFrequentistDiscovery StandardHistFactoryPlotsWithCategories StandardHypoTestDemo
               StandardHypoTestInvDemo StandardProfileInspectorDemo StandardTestStatDistributionDemo
               OneSidedFrequentistUpperLimitWithBands TwoSidedFrequentistUpperLimitWithBands StandardProfileLikelihoodDemo)
  set(roostats-${tname}-depends tutorial-roostats-CreateExampleFile)
endforeach()

#--dependency for TMVA tutorials
set (tmva-TMVAClassificationApplication-depends tutorial-tmva-TMVAClassification)
set (tmva-TMVAClassificationCategory-depends tutorial-tmva-TMVAClassification)
set (tmva-TMVAClassificationCategoryApplication-depends tutorial-tmva-TMVAClassificationCategory)
set (tmva-TMVAMulticlass-depends tutorial-tmva-TMVAMultipleBackgroundExample)
set (tmva-TMVAMulticlassApplication-depends tutorial-tmva-TMVAMulticlass)
set (tmva-TMVARegressionApplication-depends tutorial-tmva-TMVARegression)
set (tmva-tmva101_Training-depends tutorial-tmva-tmva100_DataPreparation-py)
set (tmva-tmva102_Testing-depends tutorial-tmva-tmva101_Training-py)
set (tmva-pytorch-ApplicationClassificationPyTorch-depends tutorial-tmva-pytorch-ClassificationPyTorch-py)
set (tmva-pytorch-ApplicationRegressionPyTorch-depends tutorial-tmva-pytorch-RegressionPyTorch-py)

#--List long-running tutorials to label them as "longtest"
set (long_running
     dataframe/df10[2-7]*
     multicore/mp103*)
file(GLOB long_running RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${long_running})
#--List multithreaded tutorials to run them serially
set (multithreaded
     dataframe/df10[2-7]*
     multicore/mp103*
     tmva/TMVA_CNN_Classification.C
     tmva/TMVA_RNN_Classification.C)
file(GLOB multithreaded RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${multithreaded})

#---Loop over all tutorials and define the corresponding test---------

#---Define the primordial tutorials-----------------------------------
ROOT_ADD_TEST(tutorial-hsimple COMMAND ${ROOT_root_CMD} -b -l -n -q ${CMAKE_CURRENT_SOURCE_DIR}/hsimple.C
                PASSRC 255 FAILREGEX "Error in" "error:" LABELS tutorial)
ROOT_ADD_TEST(tutorial-geom-geometry COMMAND ${ROOT_root_CMD} -b -l -n -q ${CMAKE_CURRENT_SOURCE_DIR}/geom/geometry.C
                FAILREGEX "Error in" "error:" LABELS tutorial)

#---Loop over all tutorials and define the corresponding test---------
foreach(t ${tutorials})
  list(FIND returncode_1 ${t} index)
  if(index EQUAL -1)
    set(rc 0)
  else()
    set(rc 255)
  endif()
  string(REPLACE ".C" "" tname ${t})
  string(REPLACE "/" "-" tname ${tname})

  set(labels tutorial)
  if(${t} IN_LIST long_running)
    list(APPEND labels longtest)
  endif()
  if(${t} IN_LIST multithreaded)
    list(APPEND labels multithreaded)
  endif()

   # These tests on ARM64 need much more than 20 minutes - increase the timeout
   if(ROOT_ARCHITECTURE MATCHES arm64 OR ROOT_ARCHITECTURE MATCHES ppc64)
     set(thisTestTimeout 2400) # 40m
   else()
     set(thisTestTimeout 1200) # 20m
   endif()

  ROOT_ADD_TEST(tutorial-${tname}
                COMMAND ${ROOT_root_CMD} -b -l -q ${CMAKE_CURRENT_SOURCE_DIR}/${t}${${tname}-aclic}
                PASSRC ${rc} FAILREGEX "Error in <" ": error:" "segmentation violation" "FROM HESSE     STATUS=FAILED"
                LABELS ${labels}
                DEPENDS tutorial-hsimple ${${tname}-depends}
                ENVIRONMENT ${ROOT_environ}
                TIMEOUT ${thisTestTimeout})

  if(${t} IN_LIST multithreaded)
    # Makes sure that this doesn't run in parallel with other multithreaded tutorials, and that cmake doesn't start too
    # many other tests. That we use 4 processors is actually a lie, because IMT takes whatever it finds.
    # However, even this poor indication of MT behaviour is a good hint for cmake to reduce congestion.
    set_tests_properties(tutorial-${tname} PROPERTIES RESOURCE_LOCK multithreaded PROCESSORS 4)
  endif()
endforeach()

#---Loop over all MPI tutorials and define the corresponding test---------
foreach(t ${mpi_tutorials})
  list(FIND returncode_1 ${t} index)
  if(index EQUAL -1)
    set(rc 0)
  else()
    set(rc 255)
  endif()
  string(REPLACE ".C" "" tname ${t})
  string(REPLACE "/" "-" tname ${tname})

   # These tests on ARM64 need much more than 20 minutes - increase the timeout
   if(ROOT_ARCHITECTURE MATCHES arm64 OR ROOT_ARCHITECTURE MATCHES ppc64)
     set(thisTestTimeout 2400) # 40m
   else()
     set(thisTestTimeout 1200) # 20m
   endif()

  ROOT_ADD_TEST(tutorial-${tname}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${ROOT_root_CMD} -b -l -q ${CMAKE_CURRENT_SOURCE_DIR}/${t}${${tname}-aclic}
                PASSRC ${rc} FAILREGEX "Error in <" ": error:" "segmentation violation" "FROM HESSE     STATUS=FAILED"
                LABELS tutorial
                DEPENDS tutorial-hsimple ${${tname}-depends}
                ENVIRONMENT ${ROOT_environ}
                TIMEOUT ${thisTestTimeout})
endforeach()

#---Python tutorials-----------------------------------------------------
if(ROOT_pyroot_FOUND)

  file(GLOB_RECURSE pytutorials RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

  # Now python-specific vetos:
  set(pyveto pyroot/demo.py         # requires GUI
             pyroot/fit1_py.py      # not a tutorial
             pyroot/gui_ex.py       # requires GUI
             pyroot/mrt.py          # not really a tutorial
             pyroot/na49geomfile.py # ????
             pyroot/na49visible.py  # ????
             pyroot/parse_CSV_file_with_TTree_ReadStream.py # not a tutorial
             pyroot/numberEntry.py  # requires GUI
             legacy/pyroot/*py      # legacy ...
             histfactory/example.py # not a tutorial
             histfactory/makeQuickModel.py # not a tutorial
             eve/lineset.py         # requires GUI
             sql/sqlcreatedb.py     # same as the C++ case
             sql/sqlfilldb.py       # same as the C++ case
             sql/sqlselect.py       # same as the C++ case
             launcher.py            # Not a tutorial
             )

  if(pyroot_legacy)
    # Disable in old PyROOT, which lacks Import pythonization
    list(APPEND pyveto
        roofit/rf502_wspacewrite.py
        roofit/rf504_simwstool.py
        roofit/rf509_wsinteractive.py
        roofit/rf511_wsfactory_basic.py)
  endif()

  if(NOT dataframe)
    list(APPEND pyveto pyroot/pyroot002_TTreeAsMatrix.py)
  endif()

  if(NOT dataframe
     OR ${PYTHON_VERSION_STRING_Development_Main} VERSION_LESS_EQUAL 2.7.5
     OR (DEFINED ENV{ROOTTEST_IGNORE_NUMBA_PY2} AND PYTHON_VERSION_MAJOR_Development_Main EQUAL 2)
     OR (DEFINED ENV{ROOTTEST_IGNORE_NUMBA_PY3} AND PYTHON_VERSION_MAJOR_Development_Main EQUAL 3)
     OR (MSVC AND NOT win_broken_tests))
    list(APPEND pyveto pyroot/pyroot004_NumbaDeclare.py)
  endif()

  find_python_module(xgboost QUIET)
  if(NOT PY_XGBOOST_FOUND OR NOT dataframe)
    file(GLOB tmva_veto_py RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/tmva10*.py)
    list(APPEND pyveto ${tmva_veto_py})
  endif()

  find_python_module(keras QUIET)
  if(NOT PY_KERAS_FOUND)
    file(GLOB tmva_veto_py RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/keras/*.py)
    list(APPEND pyveto ${tmva_veto_py})
  endif()

  find_python_module(torch QUIET)
  if(NOT PY_TORCH_FOUND)
    file(GLOB tmva_veto_py RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/pytorch/*.py)
    list(APPEND pyveto ${tmva_veto_py})
  endif()
  # disable PyTorch model file used by TMVA_CNN_Classification.C
  list(APPEND pyveto tmva/PyTorch_Generate_CNN_Model.py)

  # Now glob all vetos for pyroot
  file(GLOB pyveto RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${pyveto})

  list(LENGTH pytutorials nTotal)
  # Apply global .C/.py veto from above:
  list(REMOVE_ITEM pytutorials ${tutorials_veto})
  list(REMOVE_ITEM pytutorials ${pyveto})
  list(LENGTH pytutorials nAfterVeto)

  message(STATUS "${nAfterVeto}/${nTotal} python tutorials have been activated.")

  #---Python tutorials dependencies--------------------------------------
  set(pyroot-ntuple1-depends tutorial-pyroot-hsimple-py)
  set(pyroot-h1draw-depends tutorial-pyroot-hsimple-py)
  set(pyroot-benchmarks-depends tutorial-pyroot-hsimple-py
                                tutorial-pyroot-fit1-py
                                tutorial-pyroot-na49view-py
                                tutorial-pyroot-h1draw-py
                                tutorial-pyroot-ntuple1-py)
  set(pyroot-fit1-depends tutorial-pyroot-fillrandom-py)
  set(pyroot-na49view-depends tutorial-pyroot-geometry-py)
  set(roofit-rf104_classfactory-depends tutorial-roofit-rf104_classfactory) #Race condition

  #----------------------------------------------------------------------
  # List requirements for python tutorials.
  # To add a new requirement, add a glob expression that's named requires_<packageName>,
  # and add it to the list "fixtureLists" below.
  file(GLOB requires_numpy   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} dataframe/df026_AsNumpyArrays.py)
  file(GLOB requires_numba   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} pyroot/pyroot004_NumbaDeclare.py)
  file(GLOB requires_pandas  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} dataframe/df026_AsNumpyArrays.py)
  file(GLOB requires_keras   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/keras/*.py)
  file(GLOB requires_torch   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/pytorch/*.py)
  file(GLOB requires_xgboost RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tmva/tmva10*.py)
  set(fixtureLists requires_numpy requires_numba requires_pandas requires_keras requires_xgboost requires_torch)

  # Now set up all the tests
  foreach(t ${pytutorials})
    if (${t} IN_LIST returncode_1)
      set(rc 255)
    else()
      set(rc 0)
    endif()

    set(labels tutorial)
    if(${t} IN_LIST long_running)
      list(APPEND labels longtest)
    endif()
    if(${t} IN_LIST multithreaded)
      list(APPEND labels multithreaded)
    endif()

    string(REPLACE ".py" "" tname ${t})
    string(REPLACE "/" "-" tname ${tname})

    set(tutorial_name tutorial-${tname}-py)

    list(FIND pyexp_fail ${tutorial_name} index)
    if(index EQUAL -1)
      set(py_will_fail "")
    else()
      set(py_will_fail ${PYTESTS_WILLFAIL})
    endif()

    # Test if this tutorial is requiring any fixture
    unset(python_deps)
    foreach(fixtureList ${fixtureLists})
      if(${t} IN_LIST ${fixtureList})
        string(REPLACE "requires_" "" fixture ${fixtureList})
        list(APPEND python_deps ${fixture})
        list(APPEND labels python_runtime_deps)
      endif()
    endforeach()

    ROOT_ADD_TEST(${tutorial_name}
                COMMAND ${PYTHON_EXECUTABLE_Development_Main} ${CMAKE_CURRENT_SOURCE_DIR}/launcher.py ${CMAKE_CURRENT_SOURCE_DIR}/${t}
                PASSRC ${rc} FAILREGEX "Error in" ": error:" "segmentation violation"
                LABELS ${labels}
                DEPENDS ${${tname}-depends}
                ENVIRONMENT ${ROOT_environ}
                PYTHON_DEPS ${python_deps}
                ${py_will_fail})

    if(${t} IN_LIST multithreaded)
      # Makes sure that this doesn't run in parallel with other multithreaded tutorials, and that cmake doesn't start too
      # many other tests. That we use 4 processors is actually a lie, because IMT takes whatever it finds.
      # However, even this poor indication of MT behaviour is a good hint for cmake to reduce congestion.
      set_tests_properties(${tutorial_name} PROPERTIES RESOURCE_LOCK multithreaded PROCESSORS 4)
    endif()
  endforeach()
endif()
