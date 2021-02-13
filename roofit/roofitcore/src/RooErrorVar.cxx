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
\file RooErrorVar.cxx
\class RooErrorVar
\ingroup Roofitcore

RooErrorVar is an auxilary class that represents the error
of a RooRealVar as a seperate object. The main reason of
existence of this class is to facilitate the reuse of existing
techniques to perform calculations that involve a RooRealVars
error, such as calculating the pull value.
**/

#include "RooErrorVar.h"
#include "RooAbsBinning.h"
#include "RooStreamParser.h"
#include "RooRangeBinning.h"
#include "RooMsgService.h"
#include "RooUniformBinning.h"

using namespace std;

ClassImp(RooErrorVar);
;



////////////////////////////////////////////////////////////////////////////////
/// Construct an lvalue variable representing the error of RooRealVar input

RooErrorVar::RooErrorVar(const char *name, const char *title, const RooRealVar& input) :
  RooAbsRealLValue(name,title),
  _realVar("realVar","RooRealVar with error",this,(RooAbsReal&)input)
{
  _binning = new RooUniformBinning(-1,1,100) ;
}



////////////////////////////////////////////////////////////////////////////////

RooErrorVar::RooErrorVar(const RooErrorVar& other, const char* name) :
  RooAbsRealLValue(other,name),
  _realVar("realVar",this,other._realVar)
{
  _binning = other._binning->clone() ;

  // Copy constructor

  TIterator* iter = other._altBinning.MakeIterator() ;
  RooAbsBinning* binning ;
  while((binning=(RooAbsBinning*)iter->Next())) {
    _altBinning.Add(binning->clone()) ;
  }
  delete iter ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor 

RooErrorVar::~RooErrorVar()
{
  delete _binning ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return value, i.e. error on input variable

Double_t RooErrorVar::getValV(const RooArgSet*) const 
{ 
  return evaluate();
}



////////////////////////////////////////////////////////////////////////////////
/// Return true if we have binning with given name

Bool_t RooErrorVar::hasBinning(const char* name) const
{
  return _altBinning.FindObject(name) ? kTRUE : kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return binning with given name. If no binning exists with such a name, clone the default
/// binning on the fly if so requested

const RooAbsBinning& RooErrorVar::getBinning(const char* name, Bool_t verbose, Bool_t createOnTheFly) const 
{
  return const_cast<RooErrorVar*>(this)->getBinning(name,verbose,createOnTheFly) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return binning with given name. If no binning exists with such a name, clone the default
/// binning on the fly if so requested

RooAbsBinning& RooErrorVar::getBinning(const char* name, Bool_t /*verbose*/, Bool_t createOnTheFly) 
{
  // Return default (normalization) binning and range if no name is specified
  if (name==0) {
    return *_binning ;
  }
  
  // Check if binning with this name has been created already
  RooAbsBinning* binning = (RooAbsBinning*) _altBinning.FindObject(name) ;
  if (binning) {
    return *binning ;
  }

  // Return default binning if binning is not found and no creation is requested
  if (!createOnTheFly) {
    return *_binning ;
  }

  // Create a new RooRangeBinning with this name with default range
  binning = new RooRangeBinning(getMin(),getMax(),name) ;
  coutI(Contents) << "RooErrorVar::getBinning(" << GetName() << ") new range named '" 
		  << name << "' created with default bounds" << endl ;

  _altBinning.Add(binning) ;

  return *binning ;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a list of all binning names. An empty name implies the default binning.
/// A 0 pointer should be passed to getBinning in this case.

std::list<std::string> RooErrorVar::getBinningNames() const
{
  std::list<std::string> binningNames(1, "");

  RooFIter iter = _altBinning.fwdIterator();
  const RooAbsArg* binning = 0;
  while((binning = iter.next())) {
    const char* name = binning->GetName();
    binningNames.push_back(name);
  }
  return binningNames;
}


/// Remove lower bound from named binning, or default binning if name is null
void RooErrorVar::removeMin(const char* name) {
  getBinning(name).setMin(-RooNumber::infinity()) ;
}


/// Remove upper bound from named binning, or default binning if name is null
void RooErrorVar::removeMax(const char* name) {
  getBinning(name).setMax(RooNumber::infinity()) ;
}


/// Remove both upper and lower bounds from named binning, or
/// default binning if name is null
void RooErrorVar::removeRange(const char* name) {
  getBinning(name).setRange(-RooNumber::infinity(),RooNumber::infinity()) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Store given binning with this variable under the given name

void RooErrorVar::setBinning(const RooAbsBinning& binning, const char* name) 
{
  if (!name) {
    if (_binning) delete _binning ;
    _binning = binning.clone() ;
  } else {

    // Remove any old binning with this name
    RooAbsBinning* oldBinning = (RooAbsBinning*) _altBinning.FindObject(name) ;
    if (oldBinning) {
      _altBinning.Remove(oldBinning) ;
      delete oldBinning ;
    }

    // Insert new binning in list of alternative binnings
    RooAbsBinning* newBinning = binning.clone() ;
    newBinning->SetName(name) ;
    newBinning->SetTitle(name) ;
    _altBinning.Add(newBinning) ;

  }
  

}



////////////////////////////////////////////////////////////////////////////////
/// Set the lower bound of the range with the given name to the given value
/// If name is a null pointer, set the lower bound of the default range

void RooErrorVar::setMin(const char* name, Double_t value) 
{
  // Set new minimum of fit range 
  RooAbsBinning& binning = getBinning(name) ;

  // Check if new limit is consistent
  if (value >= getMax()) {
    coutW(InputArguments) << "RooErrorVar::setMin(" << GetName() 
			  << "): Proposed new fit min. larger than max., setting min. to max." << endl ;
    binning.setMin(getMax()) ;
  } else {
    binning.setMin(value) ;
  }

  // Clip current value in window if it fell out
  if (!name) {
    Double_t clipValue ;
    if (!inRange(_value,0,&clipValue)) {
      setVal(clipValue) ;
    }
  }
    
  setShapeDirty() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set the upper bound of the range with the given name to the given value
/// If name is a null pointer, set the upper bound of the default range

void RooErrorVar::setMax(const char* name, Double_t value)
{
  // Set new maximum of fit range 
  RooAbsBinning& binning = getBinning(name) ;

  // Check if new limit is consistent
  if (value < getMin()) {
    coutW(InputArguments) << "RooErrorVar::setMax(" << GetName() 
			  << "): Proposed new fit max. smaller than min., setting max. to min." << endl ;
    binning.setMax(getMin()) ;
  } else {
    binning.setMax(value) ;
  }

  // Clip current value in window if it fell out
  if (!name) {
    Double_t clipValue ;
    if (!inRange(_value,0,&clipValue)) {
      setVal(clipValue) ;
    }
  }

  setShapeDirty() ;
}

/// Set default binning to nBins uniform bins
void RooErrorVar::setBins(Int_t nBins) {
  setBinning(RooUniformBinning(getMin(),getMax(),nBins)) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the upper and lower lower bound of the range with the given name to the given values
/// If name is a null pointer, set the upper and lower bounds of the default range

void RooErrorVar::setRange( const char* name, Double_t min, Double_t max) 
{
  Bool_t exists = name ? (_altBinning.FindObject(name)?kTRUE:kFALSE) : kTRUE ;

  // Set new fit range 
  RooAbsBinning& binning = getBinning(name,kFALSE) ;

  // Check if new limit is consistent
  if (min>max) {
    coutW(InputArguments) << "RooErrorVar::setRange(" << GetName() 
			  << "): Proposed new fit max. smaller than min., setting max. to min." << endl ;
    binning.setRange(min,min) ;
  } else {
    binning.setRange(min,max) ;
  }

  if (!exists) {
    coutI(InputArguments) << "RooErrorVar::setRange(" << GetName() 
			  << ") new range named '" << name << "' created with bounds [" 
			  << min << "," << max << "]" << endl ;
  }

  setShapeDirty() ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream

Bool_t RooErrorVar::readFromStream(istream& is, Bool_t /*compact*/, Bool_t verbose) 
{
  TString token,errorPrefix("RooErrorVar::readFromStream(") ;
  errorPrefix.Append(GetName()) ;
  errorPrefix.Append(")") ;
  RooStreamParser parser(is,errorPrefix) ;
  Double_t value(0) ;

    // Compact mode: Read single token
  if (parser.readDouble(value,verbose)) return kTRUE ;
  if (isValidReal(value,verbose)) {
    setVal(value) ;
    return kFALSE ;
  } else {
    return kTRUE ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Write value to stream

void RooErrorVar::writeToStream(ostream& os, Bool_t /*compact*/) const
{
  os << getVal() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Force the internal value cache to be up to date

void RooErrorVar::syncCache(const RooArgSet*) 
{ 
  _value = evaluate() ; 
}



