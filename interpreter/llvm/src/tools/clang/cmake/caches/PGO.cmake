set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "")
set(CLANG_ENABLE_BOOTSTRAP ON CACHE BOOL "")
set(LLVM_BUILD_EXTERNAL_COMPILER_RT ON CACHE BOOL "")

set(LLVM_TARGETS_TO_BUILD X86 CACHE STRING "")
set(BOOTSTRAP_LLVM_BUILD_INSTRUMENTED ON CACHE BOOL "")
set(CLANG_BOOTSTRAP_TARGETS
  generate-profdata
  stage2
  stage2-check-all
  stage2-check-llvm
  stage2-check-clang
  stage2-test-suite CACHE STRING "")

set(CLANG_BOOTSTRAP_CMAKE_ARGS
  -C ${CMAKE_CURRENT_LIST_DIR}/PGO-stage2-instrumented.cmake
  CACHE STRING "")
