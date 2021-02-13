# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/mathmore package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(MathMore
  HEADERS
    Math/ChebyshevApprox.h
    Math/Derivator.h
    Math/DistFuncMathMore.h
    Math/GSLIntegrator.h
    Math/GSLMCIntegrator.h
    Math/GSLMinimizer.h
    Math/GSLMinimizer1D.h
    Math/GSLMultiRootFinder.h
    Math/GSLNLSMinimizer.h
    Math/GSLQuasiRandom.h
    Math/GSLRandom.h
    Math/GSLRandomFunctions.h
    Math/GSLRndmEngines.h
    Math/GSLRootFinder.h
    Math/GSLRootFinderDeriv.h
    Math/GSLSimAnMinimizer.h
    Math/InterpolationTypes.h
    Math/Interpolator.h
    Math/KelvinFunctions.h
    Math/MCParameters.h
    Math/PdfFuncMathMore.h
    Math/Polynomial.h
    Math/ParamFunction.h
    Math/QuasiRandom.h
    Math/RootFinderAlgorithms.h
    Math/SpecFuncMathMore.h
    Math/Vavilov.h
    Math/VavilovAccurate.h
    Math/VavilovAccurateCdf.h
    Math/VavilovAccuratePdf.h
    Math/VavilovAccurateQuantile.h
    Math/VavilovFast.h
SOURCES
    src/ChebyshevApprox.cxx
    src/Derivator.cxx
    src/GSLDerivator.cxx
    src/GSLError.cxx
    src/GSLIntegrator.cxx
    src/GSLInterpolator.cxx
    src/GSLMCIntegrator.cxx
    src/GSLMinimizer.cxx
    src/GSLMinimizer1D.cxx
    src/GSLMultiRootFinder.cxx
    src/GSLNLSMinimizer.cxx
    src/GSLQuasiRandom.cxx
    src/GSLRndmEngines.cxx
    src/GSLRootFinder.cxx
    src/GSLRootFinderDeriv.cxx
    src/GSLRootHelper.cxx
    src/GSLSimAnMinimizer.cxx
    src/GSLSimAnnealing.cxx
    src/Interpolator.cxx
    src/KelvinFunctions.cxx
    src/MCParameters.cxx
    src/MultiNumGradFunction.cxx
    src/PdfFuncMathMore.cxx
    src/Polynomial.cxx
    src/QuantFuncMathMore.cxx
    src/RootFinderAlgorithms.cxx
    src/SpecFuncMathMore.cxx
    src/Vavilov.cxx
    src/VavilovAccurate.cxx
    src/VavilovAccurateCdf.cxx
    src/VavilovAccuratePdf.cxx
    src/VavilovAccurateQuantile.cxx
    src/VavilovFast.cxx
    src/cblas.cxx
  LINKDEF
    Math/LinkDef.h
  DEPENDENCIES
    MathCore
  BUILTINS
    GSL
)

target_include_directories(MathMore SYSTEM PRIVATE ${GSL_INCLUDE_DIR})
target_link_libraries(MathMore PRIVATE ${GSL_LIBRARIES})

ROOT_ADD_TEST_SUBDIRECTORY(test)
