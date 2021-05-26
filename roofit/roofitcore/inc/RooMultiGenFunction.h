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
#ifndef ROO_MULTI_GEN_FUNCTION
#define ROO_MULTI_GEN_FUNCTION

#include "RooFunctor.h"
#include "Math/IFunction.h"

class RooAbsReal ;
class RooArgList ;
class RooArgSet ;
class RooAbsFunc ;

class RooMultiGenFunction : public ROOT::Math::IMultiGenFunction {

public:
  RooMultiGenFunction(const RooAbsFunc& func) ;
  RooMultiGenFunction(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters) ;
  RooMultiGenFunction(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters, const RooArgSet& nset) ;
  RooMultiGenFunction(const RooMultiGenFunction& other) ;
  virtual ~RooMultiGenFunction() ;

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const {
    return new RooMultiGenFunction(*this) ;
  }

/*   Int_t numCall() const { return _ftor.numCall() ; } */
/*   void resetNumCall() const { _ftor.resetNumCall() ; } */

  unsigned int NDim() const { return _ftor.nObs() ; }

protected:

  double DoEval(const double*) const ;

  RooFunctor _ftor ;

  ClassDef(RooMultiGenFunction,0) // Export RooAbsReal as functor
};

#endif

