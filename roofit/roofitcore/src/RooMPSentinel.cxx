/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooMPSentinel.cxx
\class RooMPSentinel
\ingroup Roofitcore

RooMPSentinel is a singleton class that keeps track of all
parellel execution processes for goodness-of-fit calculations.
The primary task of RooMPSentinel is to terminate all server processes
when the main ROOT process is exiting.
**/


#include "RooFit.h"

#include "RooMPSentinel.h"
#include "RooRealMPFE.h"

using namespace std;

ClassImp(RooMPSentinel);
  ;


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooMPSentinel::RooMPSentinel() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor. Terminate all parallel processes still registered with
/// the sentinel

RooMPSentinel::~RooMPSentinel() 
{
  TIterator *iter = _mpfeSet.createIterator() ;
  RooRealMPFE* mpfe ;
  while((mpfe=(RooRealMPFE*)iter->Next())) {
    mpfe->standby() ;
  }
  delete iter ;
}
 


////////////////////////////////////////////////////////////////////////////////
/// Register given multi-processor front-end object with the sentinel

void RooMPSentinel::add(RooRealMPFE& mpfe) 
{
  _mpfeSet.add(mpfe,kTRUE) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove given multi-processor front-end object from the sentinel

void RooMPSentinel::remove(RooRealMPFE& mpfe) 
{
  _mpfeSet.remove(mpfe,kTRUE) ;
}
