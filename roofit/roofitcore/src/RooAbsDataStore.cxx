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
\file RooAbsDataStore.cxx
\class RooAbsDataStore
\ingroup Roofitcore

RooAbsDataStore is the abstract base class for data collection that
use a TTree as internal storage mechanism
**/

#include "RooFit.h"
#include "RooMsgService.h"
#include "RooAbsDataStore.h"

#include "Riostream.h"
#include <iomanip>
#include "TClass.h"

using namespace std ;

ClassImp(RooAbsDataStore);
;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooAbsDataStore::RooAbsDataStore() 
{
  _doDirtyProp = kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooAbsDataStore::RooAbsDataStore(const char* name, const char* title, const RooArgSet& vars) : 
  TNamed(name,title)
{
  // clone the fundamentals of the given data set into internal buffer
  _vars.add(vars) ;

  _doDirtyProp = kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////

RooAbsDataStore::RooAbsDataStore(const RooAbsDataStore& other, const char* newname) : TNamed(other), RooPrintable(other)
{
  if (newname) {
    SetName(newname) ;
  }
  _vars.add(other._vars) ;
  _doDirtyProp = other._doDirtyProp ;
}



////////////////////////////////////////////////////////////////////////////////

RooAbsDataStore::RooAbsDataStore(const RooAbsDataStore& other, const RooArgSet& vars, const char* newname) : TNamed(other), RooPrintable(other)
{
  if (newname) {
    SetName(newname) ;
  }
  _vars.add(vars) ;
  _doDirtyProp = other._doDirtyProp ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsDataStore::~RooAbsDataStore()
{

}



////////////////////////////////////////////////////////////////////////////////
/// Return true if currently loaded coordinate is considered valid within
/// the current range definitions of all observables

Bool_t RooAbsDataStore::valid() const 
{
  return kTRUE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print name of dataset

void RooAbsDataStore::printName(ostream& os) const 
{
  os << GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print title of dataset

void RooAbsDataStore::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print class name of dataset

void RooAbsDataStore::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print value of the dataset, i.e. the sum of weights contained in the dataset

void RooAbsDataStore::printValue(ostream& os) const 
{
  os << numEntries() << " entries" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print argument of dataset, i.e. the observable names

void RooAbsDataStore::printArgs(ostream& os) const 
{
  os << "[" ;    
  Bool_t first(kTRUE) ;
  for (const auto arg : _vars) {
    if (first) {
      first=kFALSE ;
    } else {
      os << "," ;
    }
    os << arg->GetName() ;
  }
  os << "]" ;
}






////////////////////////////////////////////////////////////////////////////////
/// Define default print options, for a given print style

Int_t RooAbsDataStore::defaultPrintContents(Option_t* /*opt*/) const 
{
  return kName|kClassName|kArgs|kValue ;
}





////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooAbsDataStore::printMultiline(ostream& os, Int_t /*content*/, Bool_t verbose, TString indent) const 
{
  os << indent << "DataStore " << GetName() << " (" << GetTitle() << ")" << endl ;
  os << indent << "  Contains " << numEntries() << " entries" << endl;

  if (!verbose) {
    os << indent << "  Observables " << _vars << endl ;
  } else {
    os << indent << "  Observables: " << endl ;
    _vars.printStream(os,kName|kValue|kExtras|kTitle,kVerbose,indent+"  ") ;
  }

  if(verbose) {
    if (_cachedVars.getSize()>0) {
      os << indent << "  Caches " << _cachedVars << endl ;
    }
//     if(_truth.getSize() > 0) {
//       os << indent << "  Generated with ";
//       TString deeper(indent) ;
//       deeper += "   " ;
//       _truth.printStream(os,kName|kValue,kStandard,deeper) ;
//     }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Get the weights of the events in the range [first, last[.
/// This is a slow default implementation that will fill a vector with every
/// event retrieved one by one.
/// Derived classes may just return the range of the original data.

RooSpan<const double> RooAbsDataStore::getWeightBatch(std::size_t first, std::size_t last) const {

    std::vector<double> ret;//TODO try to align storage
    ret.reserve(last-first);

    for (auto i = first; i < last; ++i) {
      ret.push_back(weight(i));
    }

    return RooSpan<const double>(std::move(ret));
}

