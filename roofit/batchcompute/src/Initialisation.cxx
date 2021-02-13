#include "RooBatchCompute.h"

#include "TEnv.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <exception>


/**
 * The dispatch pointer points to the instance of the compute library in use, provided it has been loaded. 
 * The pointer is of type RooBatchComputeInterface*, so that calling functions through it are always virtual calls.
 * \see RooBatchComputeInterface, RooBatchComputeClass, RF_ARCH
 */
RooBatchCompute::RooBatchComputeInterface* RooBatchCompute::dispatch=nullptr;

namespace {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inspect cpu capabilities, and load the optimal library for RooFit computations.
void loadComputeLibrary() {

  std::string libName="libRooBatchCompute_GENERIC";

#ifdef R__RF_ARCHITECTURE_SPECIFIC_LIBS
  
  __builtin_cpu_init();
  if (gEnv->GetValue("RooFit.LoadOptimisedComputationLibrary", 1) == 0) {
    if (gDebug>0) {
      std::cout << "In roofitcore/InitUtils.cxx:loadComputeLibrary(): RooFit.LoadOptimisedComputationLibrary is set to 0, using generic RooBatchCompute library." << std::endl;
    }
  }
  
  #if __GNUC__ > 5 || defined(__clang__)
  //skylake-avx512 support
  else if (__builtin_cpu_supports("avx512cd") && __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512dq"))  {
    libName = "libRooBatchCompute_AVX512";
  }
  #endif 
  
  else if (__builtin_cpu_supports("avx2")) {
    libName = "libRooBatchCompute_AVX2";
  } else if (__builtin_cpu_supports("avx")) {
    libName = "libRooBatchCompute_AVX";
  } else if (__builtin_cpu_supports("sse4.1")) {
    libName = "libRooBatchCompute_SSE4.1";
  }

#else //R__RF_ARCHITECTURE_SPECIFIC_LIBS not defined

  if (gDebug>0) {
    std::cout << "In roofitcore/InitUtils.cxx:loadComputeLibrary(): Architecture specifics libraries not supported." << std::endl;
  }

#endif //R__RF_ARCHITECTURE_SPECIFIC_LIBS

  const auto returnValue = gSystem->Load(libName.c_str());
  if (returnValue == -1 || returnValue == -2) {
    throw std::runtime_error("RooFit was unable to load its computation library " + libName);
  } else if (returnValue == 1) {
    // Library should not have been loaded before we tried to do it.
    throw std::logic_error("RooFit computation library " + libName + " was loaded before RooFit initialisation began.");
  } else if (gDebug>0) {
    std::cout << "In roofitcore/InitUtils.cxx:loadComputeLibrary(): Library " + libName + " was loaded successfully" << std::endl;
  }
}

} //end anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// A RAII that performs RooFit's static initialisation.
static struct RooBatchComputeInitialiser {
  RooBatchComputeInitialiser() {
    loadComputeLibrary();
  }
} __RooBatchComputeInitialiser;

