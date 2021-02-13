/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : BDTEventWrapper                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *                                                                                *
 *                                                                                *
 * Author: Doug Schouten (dschoute@sfu.ca)                                        *
 *                                                                                *
 *                                                                                *
 * Copyright (c) 2007:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Texas at Austin, USA                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TMVA/BDTEventWrapper.h"

/*! \class TMVA::BDTEventWrapper
\ingroup TMVA
*/

#include "RtypesCore.h"

using namespace TMVA;

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

BDTEventWrapper::BDTEventWrapper(const Event* e) : fEvent(e) {

   fBkgWeight = 0.0;
   fSigWeight = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

BDTEventWrapper::~BDTEventWrapper() {
}

////////////////////////////////////////////////////////////////////////////////
/// Set the accumulated weight, for sorted signal/background events
///
/// @param type - true for signal, false for background
/// @param weight - the total weight

void BDTEventWrapper::SetCumulativeWeight(Bool_t type, Double_t weight) {


   if(type) fSigWeight = weight;
   else fBkgWeight = weight;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the accumulated weight

Double_t BDTEventWrapper::GetCumulativeWeight(Bool_t type) const {
   if(type) return fSigWeight;
   return fBkgWeight;
}
