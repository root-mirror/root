# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/multiproc package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(MultiProc STAGE1
  HEADERS
    MPCode.h
    MPSendRecv.h
    PoolUtils.h
    TMPClient.h
    TMPWorker.h
    TMPWorkerExecutor.h
    TProcPool.h
    ROOT/TProcessExecutor.hxx
  SOURCES
    src/MPSendRecv.cxx
    src/TMPClient.cxx
    src/TMPWorker.cxx
    src/TProcessExecutor.cxx
  LIBRARIES
    ${CMAKE_DL_LIBS}
  DEPENDENCIES
    Core
    Net
)
