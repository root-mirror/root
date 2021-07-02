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

#include "RooAbsDataStore.h"

#include "TClass.h"


////////////////////////////////////////////////////////////////////////////////
/// Print class name of dataset
void RooAbsDataStore::printClassName(std::ostream& os) const { os << IsA()->GetName() ; }


////////////////////////////////////////////////////////////////////////////////
/// Print argument of dataset, i.e. the observable names
void RooAbsDataStore::printArgs(std::ostream& os) const  { _vars.printValue(os); }


////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface
void RooAbsDataStore::printMultiline(std::ostream& os, Int_t /*content*/, Bool_t verbose, TString indent) const
{
  os << indent << "DataStore " << GetName() << " (" << GetTitle() << ")" << std::endl ;
  os << indent << "  Contains " << numEntries() << " entries" << std::endl;

  if (!verbose) {
    os << indent << "  Observables " << _vars << std::endl ;
  } else {
    os << indent << "  Observables: " << std::endl ;
    _vars.printStream(os,kName|kValue|kExtras|kTitle,kVerbose,indent+"  ") ;
  }

  if(verbose && !_cachedVars.empty()) {
    os << indent << "  Caches " << _cachedVars << std::endl ;
  }
}
