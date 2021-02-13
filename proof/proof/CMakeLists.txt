# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT proof/proof package
# @author Pere Mato, CERN
############################################################################

if(alien)
  add_definitions(-DALIENDSMGR)
  set(PROOF_EXTRA_HEADERS TDataSetManagerAliEn.h)
  set(PROOF_EXTRA_SOURCES src/TDataSetManagerAliEn.cxx)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(Proof
  HEADERS
    TCondor.h
    TDataSetManagerAliEn.h
    TDataSetManagerFile.h
    TDataSetManager.h
    TDSet.h
    TDSetProxy.h
    TLockPath.h
    TPackMgr.h
    TProofChain.h
    TProofCondor.h
    TProofDebug.h
    TProof.h
    TProofLite.h
    TProofLog.h
    TProofMgr.h
    TProofMgrLite.h
    TProofNodeInfo.h
    TProofOutputFile.h
    TProofOutputList.h
    TProofProgressStatus.h
    TProofQueryResult.h
    TProofResources.h
    TProofResourcesStatic.h
    TProofServ.h
    TProofServLite.h
    TProofSuperMaster.h
    TQueryResultManager.h
    TSelVerifyDataSet.h
    TSlave.h
    TSlaveLite.h
    TVirtualProofPlayer.h
    TVirtualPacketizer.h
    ${PROOF_EXTRA_HEADERS}
  SOURCES
    src/TCondor.cxx
    src/TDataSetManager.cxx
    src/TDataSetManagerFile.cxx
    src/TDSet.cxx
    src/TDSetProxy.cxx
    src/TLockPath.cxx
    src/TPackMgr.cxx
    src/TProofChain.cxx
    src/TProofCondor.cxx
    src/TProof.cxx
    src/TProofDebug.cxx
    src/TProofLite.cxx
    src/TProofLog.cxx
    src/TProofMgr.cxx
    src/TProofMgrLite.cxx
    src/TProofNodeInfo.cxx
    src/TProofOutputFile.cxx
    src/TProofOutputList.cxx
    src/TProofProgressStatus.cxx
    src/TProofQueryResult.cxx
    src/TProofResourcesStatic.cxx
    src/TProofServ.cxx
    src/TProofServLite.cxx
    src/TProofSuperMaster.cxx
    src/TQueryResultManager.cxx
    src/TSelVerifyDataSet.cxx
    src/TSlave.cxx
    src/TSlaveLite.cxx
    src/TVirtualProofPlayer.cxx
    src/TVirtualPacketizer.cxx
    ${PROOF_EXTRA_SOURCES}
  DEPENDENCIES
    MathCore
    Net
    RIO
    Thread
    Tree
  LIBRARIES
    ${ROOT_ATOMIC_LIBS}
)
