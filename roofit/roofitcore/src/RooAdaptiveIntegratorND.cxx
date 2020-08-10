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
\file RooAdaptiveIntegratorND.cxx
\class RooAdaptiveIntegratorND
\ingroup Roofitcore

RooAdaptiveIntegratorND implements an adaptive one-dimensional 
numerical integration algorithm.
**/


#include "RooFit.h"
#include "Riostream.h"

#include "TClass.h"
#include "RooAdaptiveIntegratorND.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNumber.h"
#include "RooMsgService.h"
#include "RooNumIntFactory.h"
#include "RooMultiGenFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

#include <assert.h>



using namespace std;

ClassImp(RooAdaptiveIntegratorND);
;

// Register this class with RooNumIntConfig

////////////////////////////////////////////////////////////////////////////////
/// Register RooAdaptiveIntegratorND, its parameters, dependencies and capabilities with RooNumIntFactory

void RooAdaptiveIntegratorND::registerIntegrator(RooNumIntFactory& fact)
{
  RooRealVar maxEval2D("maxEval2D","Max number of function evaluations for 2-dim integrals",100000) ;
  RooRealVar maxEval3D("maxEval3D","Max number of function evaluations for 3-dim integrals",1000000) ;
  RooRealVar maxEvalND("maxEvalND","Max number of function evaluations for >3-dim integrals",10000000) ;
  RooRealVar maxWarn("maxWarn","Max number of warnings on precision not reached that is printed",5) ;

  fact.storeProtoIntegrator(new RooAdaptiveIntegratorND(),RooArgSet(maxEval2D,maxEval3D,maxEvalND,maxWarn)) ;
}
 


////////////////////////////////////////////////////////////////////////////////
/// Default ctor

RooAdaptiveIntegratorND::RooAdaptiveIntegratorND()
{
  _xmin = 0 ;
  _xmax = 0 ;
  _epsRel = 1e-7 ;
  _epsAbs = 1e-7 ;
  _nmax = 10000 ;
  _func = 0 ;
  _integrator = 0 ;
  _nError = 0 ;
  _nWarn = 0 ;
  _useIntegrandLimits = kTRUE ;
  _intName = "(none)" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of integral on given function binding and with given configuration. The
/// integration limits are taken from the definition in the function binding
///_func = function.

RooAdaptiveIntegratorND::RooAdaptiveIntegratorND(const RooAbsFunc& function, const RooNumIntConfig& config) :
  RooAbsIntegrator(function)
{

  _func = new RooMultiGenFunction(function) ;  
  _nWarn = static_cast<Int_t>(config.getConfigSection("RooAdaptiveIntegratorND").getRealValue("maxWarn")) ;
  switch (_func->NDim()) {
  case 1: throw string(Form("RooAdaptiveIntegratorND::ctor ERROR dimension of function must be at least 2")) ;
  case 2: _nmax = static_cast<Int_t>(config.getConfigSection("RooAdaptiveIntegratorND").getRealValue("maxEval2D")) ; break ; 
  case 3: _nmax = static_cast<Int_t>(config.getConfigSection("RooAdaptiveIntegratorND").getRealValue("maxEval3D")) ; break ;
  default: _nmax = static_cast<Int_t>(config.getConfigSection("RooAdaptiveIntegratorND").getRealValue("maxEvalND")) ; break ;
  }
  // by default do not use absolute tolerance (see https://root.cern/phpBB3/viewtopic.php?f=15&t=20071 )
  _epsAbs = 0.0;
  _epsRel = config.epsRel();      
  _integrator = new ROOT::Math::AdaptiveIntegratorMultiDim(_epsAbs,_epsRel,_nmax) ;
  _integrator->SetFunction(*_func) ;
  _useIntegrandLimits=kTRUE ;

  _xmin = 0 ;
  _xmax = 0 ;
  _nError = 0 ;
  _nWarn = 0 ;
  checkLimits() ;
  _intName = function.getName() ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Virtual constructor with given function and configuration. Needed by RooNumIntFactory

RooAbsIntegrator* RooAdaptiveIntegratorND::clone(const RooAbsFunc& function, const RooNumIntConfig& config) const
{
  RooAbsIntegrator* ret = new RooAdaptiveIntegratorND(function,config) ;
  
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAdaptiveIntegratorND::~RooAdaptiveIntegratorND()
{
  delete[] _xmin ;
  delete[] _xmax ;
  delete _integrator ;
  delete _func ;
  if (_nError>_nWarn) {
    coutW(NumIntegration) << "RooAdaptiveIntegratorND::dtor(" << _intName 
			  << ") WARNING: Number of suppressed warningings about integral evaluations where target precision was not reached is " << _nError-_nWarn << endl ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Check that our integration range is finite and otherwise return kFALSE.
/// Update the limits from the integrand if requested.

Bool_t RooAdaptiveIntegratorND::checkLimits() const 
{
  if (!_xmin) {
    _xmin = new Double_t[_func->NDim()] ;
    _xmax = new Double_t[_func->NDim()] ;
  }
  
  if (_useIntegrandLimits) {
    for (UInt_t i=0 ; i<_func->NDim() ; i++) {
      _xmin[i]= integrand()->getMinLimit(i);
      _xmax[i]= integrand()->getMaxLimit(i);
    }
  }

  return kTRUE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Change our integration limits. Return kTRUE if the new limits are
/// ok, or otherwise kFALSE. Always returns kFALSE and does nothing
/// if this object was constructed to always use our integrand's limits.

Bool_t RooAdaptiveIntegratorND::setLimits(Double_t *xmin, Double_t *xmax) 
{
  if(_useIntegrandLimits) {
    oocoutE((TObject*)0,Integration) << "RooAdaptiveIntegratorND::setLimits: cannot override integrand's limits" << endl;
    return kFALSE;
  }
  for (UInt_t i=0 ; i<_func->NDim() ; i++) {
    _xmin[i]= xmin[i];
    _xmax[i]= xmax[i];
  }

  return checkLimits();
}




////////////////////////////////////////////////////////////////////////////////
/// Evaluate integral at given function binding parameter values

Double_t RooAdaptiveIntegratorND::integral(const Double_t* /*yvec*/) 
{
  Double_t ret = _integrator->Integral(_xmin,_xmax) ;  
  if (_integrator->Status()==1) {
    _nError++ ;
    if (_nError<=_nWarn) {
      coutW(NumIntegration) << "RooAdaptiveIntegratorND::integral(" << integrand()->getName() << ") WARNING: target rel. precision not reached due to nEval limit of "
			    << _nmax << ", estimated rel. precision is " << Form("%3.1e",_integrator->RelError()) << endl ;
    } 
    if (_nError==_nWarn) {
      coutW(NumIntegration) << "RooAdaptiveIntegratorND::integral(" << integrand()->getName() 
			    << ") Further warnings on target precision are suppressed conform specification in integrator specification" << endl ;
    }    
  }  
  return ret ;
}

