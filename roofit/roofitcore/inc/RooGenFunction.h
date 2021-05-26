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
#ifndef ROO_GEN_FUNCTION
#define ROO_GEN_FUNCTION

#include "RooFunctor.h"
#include "Math/IFunction.h"

class RooGenFunction : public ROOT::Math::IGenFunction {

public:
  RooGenFunction(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters) ;
  RooGenFunction(const RooAbsReal& func, const RooArgList& observables, const RooArgList& parameters, const RooArgSet& nset) ;
  RooGenFunction(const RooGenFunction& other) ;
  virtual ~RooGenFunction() ;

  virtual ROOT::Math::IBaseFunctionOneDim* Clone() const {
    return new RooGenFunction(*this) ;
  }

/*   Int_t numCall() const { return _ftor.numCall() ; } */
/*   void resetNumCall() const { _ftor.resetNumCall() ; } */

protected:

  double DoEval(double) const ;

  RooFunctor _ftor ;

  ClassDef(RooGenFunction,0) // Export RooAbsReal as functor
};

#endif

