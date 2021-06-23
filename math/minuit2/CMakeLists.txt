# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

cmake_minimum_required(VERSION 3.1)

if(NOT CMAKE_PROJECT_NAME STREQUAL ROOT)
  project(Minuit2 LANGUAGES CXX)
endif(NOT CMAKE_PROJECT_NAME STREQUAL ROOT)

option(minuit2_mpi "Enable support for MPI in Minuit2")
option(minuit2_omp "Enable support for OpenMP in Minuit2")

# This package can be built separately
# or as part of ROOT.
if(CMAKE_PROJECT_NAME STREQUAL ROOT)
  ROOT_STANDARD_LIBRARY_PACKAGE(Minuit2
    HEADERS
      Minuit2/ABObj.h
      Minuit2/ABProd.h
      Minuit2/ABSum.h
      Minuit2/ABTypes.h
      Minuit2/AnalyticalGradientCalculator.h
      Minuit2/BFGSErrorUpdator.h
      Minuit2/CombinedMinimizer.h
      Minuit2/CombinedMinimumBuilder.h
      Minuit2/ContoursError.h
      Minuit2/DavidonErrorUpdator.h
      Minuit2/ExternalInternalGradientCalculator.h
      Minuit2/FCNAdapter.h
      Minuit2/FCNBase.h
      Minuit2/FCNGradAdapter.h
      Minuit2/FCNGradientBase.h
      Minuit2/FumiliBuilder.h
      Minuit2/FumiliChi2FCN.h
      Minuit2/FumiliErrorUpdator.h
      Minuit2/FumiliFCNAdapter.h
      Minuit2/FumiliFCNBase.h
      Minuit2/FumiliGradientCalculator.h
      Minuit2/FumiliMaximumLikelihoodFCN.h
      Minuit2/FumiliMinimizer.h
      Minuit2/FumiliStandardChi2FCN.h
      Minuit2/FumiliStandardMaximumLikelihoodFCN.h
      Minuit2/FunctionGradient.h
      Minuit2/FunctionMinimizer.h
      Minuit2/FunctionMinimum.h
      Minuit2/GenericFunction.h
      Minuit2/GradientCalculator.h
      Minuit2/HessianGradientCalculator.h
      Minuit2/InitialGradientCalculator.h
      Minuit2/LASymMatrix.h
      Minuit2/LAVector.h
      Minuit2/LaInverse.h
      Minuit2/LaOuterProduct.h
      Minuit2/LaProd.h
      Minuit2/LaSum.h
      Minuit2/MPIProcess.h
      Minuit2/MatrixInverse.h
      Minuit2/MinimumBuilder.h
      Minuit2/MinimumError.h
      Minuit2/MinimumErrorUpdator.h
      Minuit2/MinimumParameters.h
      Minuit2/MinimumSeed.h
      Minuit2/MinimumSeedGenerator.h
      Minuit2/MinimumState.h
      Minuit2/MinosError.h
      Minuit2/Minuit2Minimizer.h
      Minuit2/MinuitParameter.h
      Minuit2/MnApplication.h
      Minuit2/MnConfig.h
      Minuit2/MnContours.h
      Minuit2/MnCovarianceSqueeze.h
      Minuit2/MnCross.h
      Minuit2/MnEigen.h
      Minuit2/MnFcn.h
      Minuit2/MnFumiliMinimize.h
      Minuit2/MnFunctionCross.h
      Minuit2/MnGlobalCorrelationCoeff.h
      Minuit2/MnHesse.h
      Minuit2/MnLineSearch.h
      Minuit2/MnMachinePrecision.h
      Minuit2/MnMatrix.h
      Minuit2/MnMigrad.h
      Minuit2/MnMinimize.h
      Minuit2/MnMinos.h
      Minuit2/MnParabola.h
      Minuit2/MnParabolaFactory.h
      Minuit2/MnParabolaPoint.h
      Minuit2/MnParameterScan.h
      Minuit2/MnPlot.h
      Minuit2/MnPosDef.h
      Minuit2/MnPrint.h
      Minuit2/MnScan.h
      Minuit2/MnSeedGenerator.h
      Minuit2/MnSimplex.h
      Minuit2/MnStrategy.h
      Minuit2/MnTiny.h
      Minuit2/MnTraceObject.h
      Minuit2/MnUserCovariance.h
      Minuit2/MnUserFcn.h
      Minuit2/MnUserParameterState.h
      Minuit2/MnUserParameters.h
      Minuit2/MnUserTransformation.h
      Minuit2/MnVectorTransform.h
      Minuit2/ModularFunctionMinimizer.h
      Minuit2/NegativeG2LineSearch.h
      Minuit2/Numerical2PGradientCalculator.h
      Minuit2/ParametricFunction.h
      Minuit2/ScanBuilder.h
      Minuit2/ScanMinimizer.h
      Minuit2/SimplexBuilder.h
      Minuit2/SimplexMinimizer.h
      Minuit2/SimplexParameters.h
      Minuit2/SimplexSeedGenerator.h
      Minuit2/SinParameterTransformation.h
      Minuit2/SqrtLowParameterTransformation.h
      Minuit2/SqrtUpParameterTransformation.h
      Minuit2/StackAllocator.h
      Minuit2/VariableMetricBuilder.h
      Minuit2/VariableMetricEDMEstimator.h
      Minuit2/VariableMetricMinimizer.h
      Minuit2/VectorOuterProduct.h
      TMinuit2TraceObject.h
    SOURCES
      src/AnalyticalGradientCalculator.cxx
      src/BFGSErrorUpdator.cxx
      src/CMakeLists.txt
      src/CombinedMinimumBuilder.cxx
      src/DavidonErrorUpdator.cxx
      src/ExternalInternalGradientCalculator.cxx
      src/FitterUtil.h
      src/FumiliBuilder.cxx
      src/FumiliErrorUpdator.cxx
      src/FumiliGradientCalculator.cxx
      src/FumiliMinimizer.cxx
      src/FumiliStandardChi2FCN.cxx
      src/FumiliStandardMaximumLikelihoodFCN.cxx
      src/HessianGradientCalculator.cxx
      src/InitialGradientCalculator.cxx
      src/LaEigenValues.cxx
      src/LaInnerProduct.cxx
      src/LaInverse.cxx
      src/LaOuterProduct.cxx
      src/LaSumOfElements.cxx
      src/LaVtMVSimilarity.cxx
      src/MPIProcess.cxx
      src/MinimumBuilder.cxx
      src/Minuit2Minimizer.cxx
      src/MnApplication.cxx
      src/MnContours.cxx
      src/MnCovarianceSqueeze.cxx
      src/MnEigen.cxx
      src/MnFcn.cxx
      src/MnFumiliMinimize.cxx
      src/MnFunctionCross.cxx
      src/MnGlobalCorrelationCoeff.cxx
      src/MnHesse.cxx
      src/MnLineSearch.cxx
      src/MnMachinePrecision.cxx
      src/MnMinos.cxx
      src/MnParabolaFactory.cxx
      src/MnParameterScan.cxx
      src/MnPlot.cxx
      src/MnPosDef.cxx
      src/MnPrint.cxx
      src/MnPrintImpl.cxx
      src/MnScan.cxx
      src/MnSeedGenerator.cxx
      src/MnStrategy.cxx
      src/MnTiny.cxx
      src/MnTraceObject.cxx
      src/MnUserFcn.cxx
      src/MnUserParameterState.cxx
      src/MnUserParameters.cxx
      src/MnUserTransformation.cxx
      src/ModularFunctionMinimizer.cxx
      src/NegativeG2LineSearch.cxx
      src/Numerical2PGradientCalculator.cxx
      src/ParametricFunction.cxx
      src/ScanBuilder.cxx
      src/SimplexBuilder.cxx
      src/SimplexParameters.cxx
      src/SimplexSeedGenerator.cxx
      src/SinParameterTransformation.cxx
      src/SqrtLowParameterTransformation.cxx
      src/SqrtUpParameterTransformation.cxx
      src/TMinuit2TraceObject.cxx
      src/VariableMetricBuilder.cxx
      src/VariableMetricEDMEstimator.cxx
      src/mnbins.cxx
      src/mndasum.cxx
      src/mndaxpy.cxx
      src/mnddot.cxx
      src/mndscal.cxx
      src/mndspmv.cxx
      src/mndspr.cxx
      src/mnlsame.cxx
      src/mnteigen.cxx
      src/mntplot.cxx
      src/mnvert.cxx
      src/mnxerbla.cxx
    DICTIONARY_OPTIONS
      -writeEmptyRootPCM
    DEPENDENCIES
      MathCore
      Hist
)
endif()

if(minuit2_omp)
  find_package(OpenMP REQUIRED)
  find_package(Threads REQUIRED)

  if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_DEFINITIONS -D_GLIBCXX_PARALLEL)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_LIBRARIES})
  endif()

  if(CMAKE_PROJECT_NAME STREQUAL ROOT)
    target_link_libraries(Minuit2 PUBLIC Threads::Threads OpenMP::OpenMP_CXX)
  endif()
endif()

if(minuit2_mpi)
  find_package(MPI REQUIRED)

  if(NOT TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_COMPILE_DEFINITIONS ${MPI_CXX_COMPILE_DEFINITIONS})
    set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_FLAGS})
    set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MPI_CXX_INCLUDE_DIRS})
    set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
  endif()

  if(CMAKE_PROJECT_NAME STREQUAL ROOT)
    target_compile_definitions(Minuit2 PUBLIC MPIPROC)
    target_link_libraries(Minuit2 PUBLIC MPI::MPI_CXX)
  endif()
endif()

if(CMAKE_PROJECT_NAME STREQUAL ROOT)
  add_definitions(-DUSE_ROOT_ERROR)
  ROOT_ADD_TEST_SUBDIRECTORY(test)
else()
  include(StandAlone.cmake)
endif()
