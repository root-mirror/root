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
\file RooAbsHiddenReal.cxx
\class RooAbsHiddenReal
\ingroup Roofitcore

RooAbsHiddenReal is a base class for objects that want to hide
their return value from interactive use, e.g. for implementations
of parameter unblinding functions. This class overrides all
printing methods with versions that do not reveal the objects value
and it has a protected version of getVal()
**/

#include "RooFit.h"

#include "Riostream.h"

#include "RooArgSet.h"
#include "RooAbsHiddenReal.h"
#include "RooCategory.h"
#include "RooMsgService.h"

using namespace std;

ClassImp(RooAbsHiddenReal);
;

RooCategory* RooAbsHiddenReal::_dummyBlindState = 0;


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooAbsHiddenReal::RooAbsHiddenReal(const char *name, const char *title, const char* unit)
  : RooAbsReal(name,title,unit),
    _state("state","Blinding state",this,dummyBlindState())
{  
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooAbsHiddenReal::RooAbsHiddenReal(const char *name, const char *title, RooAbsCategory& blindState, const char* unit)
  : RooAbsReal(name,title,unit),
  _state("state","Blinding state",this,blindState)
{  
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsHiddenReal::RooAbsHiddenReal(const RooAbsHiddenReal& other, const char* name) : 
  RooAbsReal(other, name),
  _state("state",this,other._state)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor 

RooAbsHiddenReal::~RooAbsHiddenReal() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Special version of printValue that doesn't reveal the objects value

void RooAbsHiddenReal::printValue(ostream& os) const
{
  os << "(hidden)" ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Special version of readFromStream that disallows reading from stream

Bool_t RooAbsHiddenReal::readFromStream(istream& is, Bool_t compact, Bool_t verbose)
{
  if (isHidden()) {
    // No-op version of readFromStream 
    coutE(InputArguments) << "RooAbsHiddenReal::readFromStream(" << GetName() << "): not allowed" << endl ;
    return kTRUE ;
  } else {
    return readFromStream(is,compact,verbose) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Special version of writeToStream that disallows reading from stream

void RooAbsHiddenReal::writeToStream(ostream& os, Bool_t compact) const
{
  if (isHidden()) {
    // No-op version of writeToStream 
    coutE(InputArguments) << "RooAbsHiddenReal::writeToStream(" << GetName() << "): not allowed" << endl ;
  } else {
    RooAbsReal::writeToStream(os,compact) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Return reference to internal dummy RooCategory implementation
/// blinding state switch

RooAbsCategory& RooAbsHiddenReal::dummyBlindState() const 
{
  if (!_dummyBlindState) {
    _dummyBlindState = new RooCategory("dummyBlindState","dummy blinding state") ;
    _dummyBlindState->defineType("Normal",0) ;
    _dummyBlindState->defineType("Blind",1) ;
    _dummyBlindState->setIndex(1) ;
  }
  return *_dummyBlindState ;
}


