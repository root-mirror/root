# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#---Custom CTest settings---------------------------------------------------

set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "100000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "64000")

if(CTEST_BUILD_NAME MATCHES arm64)
  set(CTEST_TEST_TIMEOUT 2400)
else()
  set(CTEST_TEST_TIMEOUT 1200)
endif()

set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
        "Warning: Rank mismatch in argument"
        "Warning: Actual argument contains too few elements"
        "has no symbols"                                         # library.a(object.c.o) has no symbols
        "note: variable tracking size limit exceeded"            # vc/tests/sse_blend.cpp
        "warning is a GCC extension"
        "bug in GCC 4.8.1"
        "warning: please use fgets or getline instead"           # deprecated use of std functions cint/ROOT
        "is dangerous, better use"                               # deprecated use of std functions cint/ROOT
        "function is dangerous and should not be used"           # deprecated use of std functions cint/ROOT
    )
set(CTEST_CUSTOM_ERROR_EXCEPTION ${CTEST_CUSTOM_ERROR_EXCEPTION}
        "fatal error: cannot open file"
        "remark: ")

#---Include other CTest Custom files----------------------------------------
if(DEFINED CTEST_BINARY_DIRECTORY)
  set(dir ${CTEST_BINARY_DIRECTORY})
else()
  set(dir .)
endif()
include(${dir}/test/CTestCustom.cmake OPTIONAL)
include(${dir}/roottest/CTestCustom.cmake OPTIONAL)
include(${dir}/rootbench/CTestCustom.cmake OPTIONAL)
include(${dir}/tutorials/CTestCustom.cmake OPTIONAL)

#---Load custom environment variables for test configurations at runtime------
include(${dir}/CTestEnvVars.cmake OPTIONAL)
