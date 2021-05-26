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
#ifndef ROO_FUNCTOR
#define ROO_FUNCTOR

#include "RooArgSet.h"
#include "RooAbsReal.h"

class RooAbsFunc ;
class RooAbsPdf ;

class RooFunctor {

public:
  RooFunctor(const RooAbsFunc& func) ;
  RooFunctor(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters) ;
  RooFunctor(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters, const RooArgSet& nset) ;
  RooFunctor(const RooFunctor&) ;
  virtual ~RooFunctor() ;

  Int_t nObs() const {
    // Number of observables
    return _nobs ;
  }
  Int_t nPar() const {
    // Number of parameters;
    return _npar ;
  }

  //Int_t numCall() const { return _binding->numCall() ; }
  //void resetNumCall() const { _binding->resetNumCall() ; }


  Double_t operator()(Double_t x) const { return eval(x) ; }
  Double_t operator()(const Double_t* x, const Double_t* p) const { return eval(x,p) ; }
  Double_t operator()(const Double_t* x) const { return eval(x) ; }

  Double_t eval(const Double_t* /*x*/, const Double_t* /*p*/) const ;
  Double_t eval(const Double_t* /*x*/) const ;
  Double_t eval(Double_t  /*x*/) const ;

  RooAbsFunc& binding() { return *_binding ; }

protected:

  Bool_t          _ownBinding ; // Do we own the binding function
  RooArgSet       _nset ; // Normalization observables
  RooAbsFunc*     _binding ; // Function binding
  Double_t*       _x ; //! Transfer array ;
  Int_t           _npar ; //! Number of parameters ;
  Int_t           _nobs ; //! Number of observables ;

  ClassDef(RooFunctor,0) // Export RooAbsReal as functor
};

#endif

