# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT rootfit/histfactory package
# @author Pere Mato, CERN
############################################################################

if(mathmore OR ROOT_mathmore_FOUND)
  set(ROOT_MATHMORE_LIBRARY MathMore)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(HistFactory
  HEADERS
    RooStats/HistFactory/Asimov.h
    RooStats/HistFactory/Channel.h
    RooStats/HistFactory/ConfigParser.h
    RooStats/HistFactory/Data.h
    RooStats/HistFactory/EstimateSummary.h
    RooStats/HistFactory/FlexibleInterpVar.h
    RooStats/HistFactory/HistFactoryException.h
    RooStats/HistFactory/HistFactoryModelUtils.h
    RooStats/HistFactory/HistFactoryNavigation.h
    RooStats/HistFactory/HistFactorySimultaneous.h
    RooStats/HistFactory/HistoToWorkspaceFactoryFast.h
    RooStats/HistFactory/HistoToWorkspaceFactory.h
    RooStats/HistFactory/HistRef.h
    RooStats/HistFactory/LinInterpVar.h
    RooStats/HistFactory/MakeModelAndMeasurementsFast.h
    RooStats/HistFactory/Measurement.h
    RooStats/HistFactory/ParamHistFunc.h
    RooStats/HistFactory/PiecewiseInterpolation.h
    RooStats/HistFactory/PreprocessFunction.h
    RooStats/HistFactory/RooBarlowBeestonLL.h
    RooStats/HistFactory/Sample.h
    RooStats/HistFactory/Systematics.h
  SOURCES
    src/Asimov.cxx
    src/Channel.cxx
    src/ConfigParser.cxx
    src/Data.cxx
    src/EstimateSummary.cxx
    src/FlexibleInterpVar.cxx
    src/Helper.cxx
    src/Helper.h
    src/HistFactoryModelUtils.cxx
    src/HistFactoryNavigation.cxx
    src/HistFactorySimultaneous.cxx
    src/HistoToWorkspaceFactory.cxx
    src/HistoToWorkspaceFactoryFast.cxx
    src/HistRef.cxx
    src/LinInterpVar.cxx
    src/MakeModelAndMeasurements.cxx
    src/MakeModelAndMeasurementsFast.cxx
    src/Measurement.cxx
    src/ParamHistFunc.cxx
    src/PiecewiseInterpolation.cxx
    src/PreprocessFunction.cxx
    src/RooBarlowBeestonLL.cxx
    src/Sample.cxx
    src/Systematics.cxx
  DICTIONARY_OPTIONS
    "-writeEmptyRootPCM"
  DEPENDENCIES
    RooFit
    RooFitCore
    Tree
    RIO
    Hist
    Matrix
    MathCore
    Graf
    Gpad
    RooStats
    XMLParser
)

# For recent clang, this can facilitate auto-vectorisation.
# In RooFit, the errno side effect is not needed, anyway:
target_compile_options(HistFactory PUBLIC -fno-math-errno)

ROOT_EXECUTABLE(hist2workspace
  MakeModelAndMeasurements.cxx
  hist2workspace.cxx
  LIBRARIES
    HistFactory
    RooStats
    RooFit
    RooFitCore
    RIO
)

#---Createhist2workspaceCommandLineOptions------------------------------------------------------------------
generateHeader(hist2workspace
  ${CMAKE_CURRENT_SOURCE_DIR}/src/hist2workspace-argparse.py
  ${CMAKE_BINARY_DIR}/ginclude/hist2workspaceCommandLineOptionsHelp.h
)

file(COPY config/prepareHistFactory DESTINATION  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/prepareHistFactory
                PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                            GROUP_EXECUTE GROUP_READ
                            WORLD_EXECUTE WORLD_READ
                DESTINATION ${CMAKE_INSTALL_BINDIR})

ROOT_ADD_TEST_SUBDIRECTORY(test)
