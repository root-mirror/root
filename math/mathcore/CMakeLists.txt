# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/mathcore package
############################################################################

set(HEADERS
  Fit/BasicFCN.h
  Fit/BinData.h
  Fit/Chi2FCN.h
  Fit/DataOptions.h
  Fit/DataRange.h
  Fit/FcnAdapter.h
  Fit/FitConfig.h
  Fit/FitData.h
  Fit/FitExecutionPolicy.h
  Fit/FitResult.h
  Fit/FitUtil.h
  Fit/Fitter.h
  Fit/LogLikelihoodFCN.h
  Fit/ParameterSettings.h
  Fit/PoissonLikelihoodFCN.h
  Fit/SparseData.h
  Fit/UnBinData.h
  Math/AdaptiveIntegratorMultiDim.h
  Math/AllIntegrationTypes.h
  Math/BasicMinimizer.h
  Math/BrentMethods.h
  Math/BrentMinimizer1D.h
  Math/BrentRootFinder.h
  Math/ChebyshevPol.h
  Math/Delaunay2D.h
  Math/DistFuncMathCore.h
  Math/DistSampler.h
  Math/DistSamplerOptions.h
  Math/Error.h
  Math/Factory.h
  Math/FitMethodFunction.h
  Math/Functor.h
  Math/GaussIntegrator.h
  Math/GaussLegendreIntegrator.h
  Math/GenAlgoOptions.h
  Math/GoFTest.h
  Math/IFunction.h
  Math/IFunctionfwd.h
  Math/IMinimizer1D.h
  Math/IOptions.h
  Math/IParamFunction.h
  Math/IParamFunctionfwd.h
  Math/IRootFinderMethod.h
  Math/Integrator.h
  Math/IntegratorMultiDim.h
  Math/IntegratorOptions.h
  Math/KDTree.h
  Math/LCGEngine.h
  Math/Math.h
  Math/MersenneTwisterEngine.h
  Math/MinimTransformFunction.h
  Math/MinimTransformVariable.h
  Math/Minimizer.h
  Math/MinimizerOptions.h
  Math/MinimizerVariableTransformation.h
  Math/MixMaxEngine.h
  Math/MultiDimParamFunctionAdapter.h
  Math/OneDimFunctionAdapter.h
  Math/ParamFunctor.h
  Math/PdfFuncMathCore.h
  Math/ProbFuncMathCore.h
  Math/QuantFuncMathCore.h
  Math/Random.h
  Math/RandomFunctions.h
  Math/RanluxppEngine.h
  Math/RichardsonDerivator.h
  Math/RootFinder.h
  Math/SpecFuncMathCore.h
  Math/StdEngine.h
  Math/TDataPoint.h
  Math/TDataPointN.h
  Math/TRandomEngine.h
  Math/Types.h
  Math/Util.h
  Math/VirtualIntegrator.h
  Math/WrappedFunction.h
  Math/WrappedParamFunction.h
  TComplex.h
  TKDTree.h
  TKDTreeBinning.h
  TMath.h
  TRandom.h
  TRandom1.h
  TRandom2.h
  TRandom3.h
  TRandomGen.h
  TStatistic.h
  VectorizedTMath.h
)

if(runtime_cxxmodules)
  # Adding Math/Math.h to the MathCore module makes the generation of GenVector
  # dictionary hand due to a bug in rootcling/clang. The redeclaration chain of
  # 'abs' shadow declaration (from stl math.h) is broken.
  # FIXME: Revise after a llvm upgrade or reproduce it outside rootcling.
  list(REMOVE_ITEM HEADERS "Math/Math.h")
endif()

ROOT_ADD_C_FLAG(_flags -Wno-strict-overflow)  # Avoid what it seems a compiler false positive warning
ROOT_ADD_C_FLAG(_flags -Wno-maybe-uninitialized)  # Avoid what it seems a compiler false positive warning
ROOT_ADD_C_FLAG(_flags -Wno-parentheses-equality)

set_source_files_properties(src/triangle.c COMPILE_FLAGS "${_flags}")

if(imt)
  set(MATHCORE_DEPENDENCIES Imt)
endif()

if(veccore)
  set(MATHCORE_BUILTINS VECCORE)
  set(MATHCORE_LIBRARIES ${VecCore_LIBRARIES})
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(MathCore
  HEADERS
    ${HEADERS}
  SOURCES
    src/AdaptiveIntegratorMultiDim.cxx
    src/BasicMinimizer.cxx
    src/BinData.cxx
    src/BrentMethods.cxx
    src/BrentMinimizer1D.cxx
    src/BrentRootFinder.cxx
    src/ChebyshevPol.cxx
    src/DataRange.cxx
    src/Delaunay2D.cxx
    src/DistSampler.cxx
    src/DistSamplerOptions.cxx
    src/Factory.cxx
    src/FitConfig.cxx
    src/FitData.cxx
    src/FitResult.cxx
    src/FitUtil.cxx
    src/Fitter.cxx
    src/GaussIntegrator.cxx
    src/GaussLegendreIntegrator.cxx
    src/GenAlgoOptions.cxx
    src/GoFTest.cxx
    src/Integrator.cxx
    src/IntegratorOptions.cxx
    src/MersenneTwisterEngine.cxx
    src/MinimTransformFunction.cxx
    src/MinimizerOptions.cxx
    src/MinimizerVariableTransformation.cxx
    src/MixMaxEngineImpl17.cxx
    src/MixMaxEngineImpl240.cxx
    src/MixMaxEngineImpl256.cxx
    src/ParameterSettings.cxx
    src/PdfFuncMathCore.cxx
    src/ProbFuncMathCore.cxx
    src/QuantFuncMathCore.cxx
    src/RandomFunctions.cxx
    src/RanluxppEngineImpl.cxx
    src/RichardsonDerivator.cxx
    src/RootFinder.cxx
    src/SparseData.cxx
    src/SpecFuncCephes.cxx
    src/SpecFuncCephes.h
    src/SpecFuncCephesInv.cxx
    src/SpecFuncMathCore.cxx
    src/StdEngine.cxx
    src/TComplex.cxx
    src/TKDTree.cxx
    src/TKDTreeBinning.cxx
    src/TMath.cxx
    src/TRandom.cxx
    src/TRandom1.cxx
    src/TRandom2.cxx
    src/TRandom3.cxx
    src/TRandomGen.cxx
    src/TStatistic.cxx
    src/UnBinData.cxx
    src/triangle.c
    src/VectorizedTMath.cxx
  LIBRARIES
    ${MATHCORE_LIBRARIES}
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    ${MATHCORE_DEPENDENCIES}
  BUILTINS
    ${MATHCORE_BUILTINS}
)

target_include_directories(MathCore PRIVATE ${Vc_INCLUDE_DIR})
target_include_directories(MathCore PRIVATE ${VecCore_INCLUDE_DIRS})

list(APPEND math_incl ${Vc_INCLUDE_DIR})
list(APPEND math_incl ${VecCore_INCLUDE_DIRS})

foreach(incl ${math_incl})
   target_include_directories(MathCore PUBLIC $<BUILD_INTERFACE:${incl}>)
endforeach()

target_compile_definitions(MathCore INTERFACE ${VecCore_DEFINITIONS})
target_link_libraries(MathCore PRIVATE ${CMAKE_THREAD_LIBS_INIT})

ROOT_ADD_TEST_SUBDIRECTORY(test)
