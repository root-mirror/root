# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT proof/proofd package
# @author Pere Mato, CERN
############################################################################

set(PROOFD_SOURCES
  src/XpdObject.cxx
  src/XProofProtUtils.cxx
  src/XrdProofConn.cxx
  src/XrdProofdAdmin.cxx
  src/XrdProofdAux.cxx
  src/XrdProofdClient.cxx
  src/XrdProofdClientMgr.cxx
  src/XrdProofdConfig.cxx
  src/XrdProofdManager.cxx
  src/XrdProofdNetMgr.cxx
  src/XrdProofdPriorityMgr.cxx
  src/XrdProofdProofServ.cxx
  src/XrdProofdProofServMgr.cxx
  src/XrdProofdProtocol.cxx
  src/XrdProofdResponse.cxx
  src/XrdProofdSandbox.cxx
  src/XrdProofGroup.cxx
  src/XrdProofPhyConn.cxx
  src/XrdProofSched.cxx
  src/XrdProofWorker.cxx
  src/XrdROOT.cxx
  ../../net/rpdutils/src/rpdconn.cxx
)

include_directories(${XROOTD_INCLUDE_DIRS})
include_directories(AFTER ${CMAKE_CURRENT_SOURCE_DIR}/inc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../net/rpdutils/res)
add_definitions(${XROOTD_CFLAGS})
#
# This should not be needed any longer after fix 039136791216af4620b9041e89d60c3ac0bab2f7
# Uncomment if not the case (G Ganis, 18/5/2016)
# ROOT_ADD_CXX_FLAG(CMAKE_CXX_FLAGS -Wno-undefined-bool-conversion)

ROOT_ADD_CXX_FLAG(CMAKE_CXX_FLAGS -Wno-nonnull-compare)

if(WIN32)
  ROOT_LINKER_LIBRARY(XrdProofd
    XProofProtUtils.cxx
    LIBRARIES
      ${XROOTD_LIBRARIES}
  )
else()
  if (NOT XROOTD_NOMAIN)
    ROOT_EXECUTABLE(xproofd
      ${PROOFD_SOURCES}
      LIBRARIES
        ${XROOTD_LIBRARIES}
      BUILTINS
        XROOTD
    )
  endif()

  ROOT_LINKER_LIBRARY(XrdProofd
    ${PROOFD_SOURCES}
    LIBRARIES
      ${XROOTD_LIBRARIES}
      rpdutil
      ${CMAKE_THREAD_LIBS_INIT}
      ${CMAKE_DL_LIBS}
    BUILTINS
      XROOTD
  )
  set_target_properties(XrdProofd PROPERTIES INTERFACE_LINK_LIBRARIES "")

  ROOT_EXECUTABLE(proofexecv
    proofexecv.cxx
    ${CMAKE_SOURCE_DIR}/net/rpdutils/src/rpdconn.cxx
    ${CMAKE_SOURCE_DIR}/net/rpdutils/src/rpdpriv.cxx
  )

  target_include_directories(proofexecv PRIVATE
    ${CMAKE_SOURCE_DIR}/core/base/inc    # for RConfig.hxx
    ${CMAKE_BINARY_DIR}/ginclude         # for RConfigure.h
  )

endif()
