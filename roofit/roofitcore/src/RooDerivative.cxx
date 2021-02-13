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
\file RooDerivative.cxx
\class RooDerivative
\ingroup Roofitcore

RooDerivative represents the first, second, or third order derivative
of any RooAbsReal as calculated (numerically) by the MathCore Richardson
derivator class.
**/


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooDerivative.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooFunctor.h"

#include "Math/WrappedFunction.h"
#include "Math/RichardsonDerivator.h"

using namespace std;

ClassImp(RooDerivative);



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooDerivative::RooDerivative() : _order(1), _eps(1e-7), _ftor(0), _rd(0)
{
}



////////////////////////////////////////////////////////////////////////////////

RooDerivative::RooDerivative(const char* name, const char* title, RooAbsReal& func, RooRealVar& x, Int_t orderIn, Double_t epsIn) :
  RooAbsReal(name, title),
  _order(orderIn),
  _eps(epsIn),
  _nset("nset","nset",this,kFALSE,kFALSE),
  _func("function","function",this,func),
  _x("x","x",this,x),
  _ftor(0),
  _rd(0)
{
  if (_order<0 || _order>3 ) {
    throw std::string(Form("RooDerivative::ctor(%s) ERROR, derivation order must be 1,2 or 3",name)) ;
  }
}

////////////////////////////////////////////////////////////////////////////////

RooDerivative::RooDerivative(const char* name, const char* title, RooAbsReal& func, RooRealVar& x, const RooArgSet& nset, Int_t orderIn, Double_t epsIn) :
  RooAbsReal(name, title),
  _order(orderIn),
  _eps(epsIn),
  _nset("nset","nset",this,kFALSE,kFALSE),
  _func("function","function",this,func),
  _x("x","x",this,x),
  _ftor(0),
  _rd(0)
{
  if (_order<0 || _order>3) { 
    throw std::string(Form("RooDerivative::ctor(%s) ERROR, derivation order must be 1,2 or 3",name)) ;
  }
  _nset.add(nset) ;
}



////////////////////////////////////////////////////////////////////////////////

RooDerivative::RooDerivative(const RooDerivative& other, const char* name) :
  RooAbsReal(other, name), 
  _order(other._order),  
  _eps(other._eps),
  _nset("nset",this,other._nset),
  _func("function",this,other._func),
  _x("x",this,other._x),
  _ftor(0),
  _rd(0)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooDerivative::~RooDerivative() 
{
  if (_rd) delete _rd ;
  if (_ftor) delete _ftor ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate value

Double_t RooDerivative::evaluate() const 
{
  if (!_ftor) {
    _ftor = _func.arg().functor(_x.arg(),RooArgSet(),_nset)  ;    
    ROOT::Math::WrappedFunction<RooFunctor&> wf(*_ftor);
    _rd = new ROOT::Math::RichardsonDerivator(wf,_eps*(_x.max()-_x.min()),kTRUE) ;
  }
  
  switch (_order) {    
  case 1: return _rd->Derivative1(_x);
  case 2: return _rd->Derivative2(_x);
  case 3: return _rd->Derivative3(_x);
  }
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Zap functor and derivator ;

Bool_t RooDerivative::redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) 
{
  delete _ftor ;
  delete _rd ; 
  _ftor = 0 ;
  _rd = 0 ;
  return kFALSE ;
}
