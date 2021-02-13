# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT proof/proofplayer package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ProofPlayer
  NO_INSTALL_HEADERS
  HEADERS
    TDrawFeedback.h
    TEventIter.h
    TOutputListSelectorDataMap.h
    TPacketizerAdaptive.h
    TPacketizerFile.h
    TPacketizer.h
    TPacketizerMulti.h
    TPacketizerUnit.h
    TPerfStats.h
    TProofLimitsFinder.h
    TProofMonSender.h
    TProofMonSenderML.h
    TProofMonSenderSQL.h
    TProofPlayer.h
    TProofPlayerLite.h
    TStatsFeedback.h
    TStatus.h
  SOURCES
    TDrawFeedback.cxx
    TEventIter.cxx
    TOutputListSelectorDataMap.cxx
    TPacketizerAdaptive.cxx
    TPacketizer.cxx
    TPacketizerFile.cxx
    TPacketizerMulti.cxx
    TPacketizerUnit.cxx
    TPerfStats.cxx
    TProofLimitsFinder.cxx
    TProofMonSender.cxx
    TProofMonSenderML.cxx
    TProofMonSenderSQL.cxx
    TProofPlayer.cxx
    TProofPlayerLite.cxx
    TStatsFeedback.cxx
    TStatus.cxx
  DEPENDENCIES
    Hist
    Net
    MathCore
    Proof
    RIO
    Thread
    Tree
    TreePlayer
)

ROOT_STANDARD_LIBRARY_PACKAGE(ProofDraw
  NO_INSTALL_HEADERS
  HEADERS
    TProofDraw.h
  SOURCES
    TProofDraw.cxx
  LINKDEF
    LinkDefDraw.h
  DEPENDENCIES
    ProofPlayer
    TreePlayer
)

ROOT_INSTALL_HEADERS()
