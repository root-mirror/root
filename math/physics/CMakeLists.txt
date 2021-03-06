# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT math/physics package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(Physics
  HEADERS
    TFeldmanCousins.h
    TDecay.h
    TGenDecay.h
    TGenFoamDecay.h
    TGenInterface.h
    TGenPhaseSpace.h
    TLorentzRotation.h
    TLorentzVector.h
    TQuaternion.h
    TRobustEstimator.h
    TRolke.h
    TRotation.h
    TVector2.h
    TVector3.h
  SOURCES
    src/TFeldmanCousins.cxx
    src/TDecay.cxx
    src/TGenDecay.cxx
    src/TGenFoamDecay.cxx
    src/TGenPhaseSpace.cxx
    src/TLorentzRotation.cxx
    src/TLorentzVector.cxx
    src/TQuaternion.cxx
    src/TRobustEstimator.cxx
    src/TRolke.cxx
    src/TRotation.cxx
    src/TVector2.cxx
    src/TVector3.cxx
  DEPENDENCIES
    Foam
    Matrix
    MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
