/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooConvIntegrandBinding.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_CONV_INTEGRAND_BINDING
#define ROO_CONV_INTEGRAND_BINDING

#include "RooAbsFunc.h"

class RooAbsRealLValue;
class RooAbsReal;
class RooArgSet;

class RooConvIntegrandBinding : public RooAbsFunc {
public:
  RooConvIntegrandBinding(const RooAbsReal& func, const RooAbsReal& model,
                RooAbsReal& x, RooAbsReal& xprime,
                     const RooArgSet* nset=0, Bool_t clipInvalid=kFALSE);
  virtual ~RooConvIntegrandBinding();

  virtual Double_t operator()(const Double_t xvector[]) const;
  virtual Double_t getMinLimit(UInt_t dimension) const;
  virtual Double_t getMaxLimit(UInt_t dimension) const;
  inline void setNormalizationSet(const RooArgSet* nset) {
    // Use the supplied nset as normalization set for calls to func and model
    _nset = nset ;
  }

protected:
  void loadValues(const Double_t xvector[], Bool_t clipInvalid=kFALSE) const;

  const RooAbsReal *_func;   // Pointer to input function
  const RooAbsReal *_model ; // Pointer to input resolution model

  RooAbsRealLValue **_vars;  // Array of pointers to variables
  const RooArgSet *_nset;    // Normalization set to be used for function evaluations
  mutable Bool_t _xvecValid; // If true _xvec defines a valid point
  Bool_t _clipInvalid ;      // If true, invalid x values are clipped into their valid range

  ClassDef(RooConvIntegrandBinding,0) // RooAbsFunc representation of convolution integrands
};

#endif

