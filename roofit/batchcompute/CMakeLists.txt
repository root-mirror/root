# Library which powers fast batch computations in Roofit.
ROOT_LINKER_LIBRARY(RooBatchCompute
    src/Initialisation.cxx
    src/RunContext.cxx
  DEPENDENCIES
    Core
    MathCore
)

if(vdt AND NOT builtin_vdt)
  target_include_directories(RooBatchCompute PUBLIC ${VDT_INCLUDE_DIR})
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # gcc has a bug in that it fails to put its internal -lgcc into the right place when linking.
  # We need it to check cpu flags in src/Initialisation.cxx
  # Here, we add an explicit link instruction according to the workaround posted here:
  # https://bugs.launchpad.net/ubuntu/+source/gcc-5/+bug/1568899
  target_link_libraries(RooBatchCompute PRIVATE -lgcc_s -lgcc)
endif()


############################################################################################################################################
# Instantiations of the shared objects which provide the actual computation functions.

# Generic implementation for CPUs that don't support vector instruction sets.
ROOT_LINKER_LIBRARY(RooBatchCompute_GENERIC src/RooBatchCompute.cxx TYPE SHARED DEPENDENCIES RooFitCore RooBatchCompute)

# Windows platform and ICC compiler need special code and testing, thus the feature has not been implemented yet for these.
if (ROOT_PLATFORM MATCHES "linux|macosx" AND CMAKE_SYSTEM_PROCESSOR MATCHES x86_64 AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

  target_compile_options(RooBatchCompute PRIVATE -DR__RF_ARCHITECTURE_SPECIFIC_LIBS)

  ROOT_LINKER_LIBRARY(RooBatchCompute_SSE4.1  src/RooBatchCompute.cxx TYPE SHARED DEPENDENCIES RooFitCore RooBatchCompute)
  ROOT_LINKER_LIBRARY(RooBatchCompute_AVX     src/RooBatchCompute.cxx TYPE SHARED DEPENDENCIES RooFitCore RooBatchCompute)
  ROOT_LINKER_LIBRARY(RooBatchCompute_AVX2    src/RooBatchCompute.cxx TYPE SHARED DEPENDENCIES RooFitCore RooBatchCompute)

  # Flags -fno-signaling-nans, -fno-trapping-math and -O3 are necessary to enable autovectorization (especially for GCC).
  set(common-flags $<$<CXX_COMPILER_ID:GNU>:-fno-signaling-nans>)
  list(APPEND common-flags $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>: -fno-trapping-math -O3>)

  target_compile_options(RooBatchCompute_SSE4.1  PRIVATE ${common-flags} -msse4    -DRF_ARCH=SSE4)
  target_compile_options(RooBatchCompute_AVX     PRIVATE ${common-flags} -mavx     -DRF_ARCH=AVX )
  target_compile_options(RooBatchCompute_AVX2    PRIVATE ${common-flags} -mavx2    -DRF_ARCH=AVX2)

  # AVX512 is only supported in gcc 6+
  # We focus on AVX512 capable processors that support at least the skylake-avx512 instruction sets.
  if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 6)
    ROOT_LINKER_LIBRARY(RooBatchCompute_AVX512  src/RooBatchCompute.cxx TYPE SHARED DEPENDENCIES RooFitCore RooBatchCompute)
    target_compile_options(RooBatchCompute_AVX512  PRIVATE ${common-flags} -march=skylake-avx512 -DRF_ARCH=AVX512)
  endif()

endif()

ROOT_INSTALL_HEADERS()
