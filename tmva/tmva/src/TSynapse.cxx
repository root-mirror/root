// @(#)root/tmva $Id$
// Author: Matt Jachowski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TSynapse                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Matt Jachowski  <jachowski@stanford.edu> - Stanford University, USA       *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::TSynapse
\ingroup TMVA
Synapse class used by TMVA artificial neural network methods
*/

#include "TMVA/TSynapse.h"

#include "TMVA/TNeuron.h"

#include "TMVA/MsgLogger.h"

#include "TMVA/Types.h"

#include "ThreadLocalStorage.h"

static const Int_t fgUNINITIALIZED = -1;

ClassImp(TMVA::TSynapse);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::TSynapse::TSynapse()
   : fWeight( 0 ),
     fLearnRate( 0 ),
     fDelta( 0 ),
     fDEDw( 0 ),
     fCount( 0 ),
     fPreNeuron( NULL ),
     fPostNeuron( NULL )
{
   fWeight     = fgUNINITIALIZED;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::TSynapse::~TSynapse()
{
}

////////////////////////////////////////////////////////////////////////////////
/// set synapse weight

void TMVA::TSynapse::SetWeight(Double_t weight)
{
   fWeight = weight;
}

////////////////////////////////////////////////////////////////////////////////
/// get output of pre-neuron weighted by synapse weight

Double_t TMVA::TSynapse::GetWeightedValue()
{
   if (fPreNeuron == NULL)
      Log() << kFATAL << "<GetWeightedValue> synapse not connected to neuron" << Endl;

   return (fWeight * fPreNeuron->GetActivationValue());
}

////////////////////////////////////////////////////////////////////////////////
/// get error field of post-neuron weighted by synapse weight

Double_t TMVA::TSynapse::GetWeightedDelta()
{
   if (fPostNeuron == NULL)
      Log() << kFATAL << "<GetWeightedDelta> synapse not connected to neuron" << Endl;

   return fWeight * fPostNeuron->GetDelta();
}

////////////////////////////////////////////////////////////////////////////////
/// adjust the weight based on the error field all ready calculated by CalculateDelta

void TMVA::TSynapse::AdjustWeight()
{
   Double_t wDelta = fDelta / fCount;
   fWeight += -fLearnRate * wDelta;
   InitDelta();
}

////////////////////////////////////////////////////////////////////////////////
/// calculate/adjust the error field for this synapse

void TMVA::TSynapse::CalculateDelta()
{
   fDelta += fPostNeuron->GetDelta() * fPreNeuron->GetActivationValue();
   fCount++;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::MsgLogger& TMVA::TSynapse::Log() const
{
   TTHREAD_TLS_DECL_ARG(MsgLogger,logger,"TSynapse");  //! message logger, static to save resources
   return logger;
}
