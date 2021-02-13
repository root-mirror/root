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
\file RooAbsBinning.cxx
\class RooAbsBinning
\ingroup Roofitcore

RooAbsBinning is the abstract base class for RooRealVar binning definitions.
This class defines the interface to retrieve bin boundaries, ranges etc.
**/

#include "RooAbsBinning.h"

#include "RooAbsReal.h"
#include "RooFit.h"
#include "TBuffer.h"
#include "TClass.h"

#include "Riostream.h"

using namespace std;

ClassImp(RooAbsBinning);
;


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooAbsBinning::RooAbsBinning(const char* name) : TNamed(name,name)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsBinning::~RooAbsBinning() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Print binning name

void RooAbsBinning::printName(ostream& os) const 
{
  os << GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print binning title

void RooAbsBinning::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print binning class name

void RooAbsBinning::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print binning arguments (the RooAbsReal objects represening
/// the variable bin boundaries for parameterized binning implementations

void RooAbsBinning::printArgs(ostream& os) const 
{
  os << "[ " ;    
  if (lowBoundFunc()) {
    os << "lowerBound=" << lowBoundFunc()->GetName() ;
  }
  if (highBoundFunc()) {
    if (lowBoundFunc()) {
      os << " " ;
    }
    os << "upperBound=" << highBoundFunc()->GetName() ;
  }
  os << " ]" ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Print binning value, i.e the bin boundary positions

void RooAbsBinning::printValue(ostream &os) const
{
  Int_t n = numBins() ;
  os << "B(" ;
  
  Int_t i ;
  for (i=0 ; i<n ; i++) {
    if (i>0) {
      os << " : " ;
    }
    os << binLow(i) ;
  }
  os << " : " << binHigh(n-1) ;
  os << ")" ;

}



////////////////////////////////////////////////////////////////////////////////
/// Custom streamer implementing schema evolution between V1 and V2 persistent binnings

void RooAbsBinning::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v==1) {
	TObject::Streamer(R__b);
      } else {
	TNamed::Streamer(R__b);
      }
      RooPrintable::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooAbsBinning::IsA());
   } else {
      R__c = R__b.WriteVersion(RooAbsBinning::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      RooPrintable::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

