/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_DERIVATIVE
#define ROO_DERIVATIVE

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"

namespace ROOT{ namespace Math {
class RichardsonDerivator;
}}

class RooRealVar;
class RooArgList ;

class RooDerivative : public RooAbsReal {
public:

  RooDerivative() ;
  RooDerivative(const char *name, const char *title, RooAbsReal& func, RooRealVar& x, Int_t order=1, Double_t eps=0.001) ;
  RooDerivative(const char *name, const char *title, RooAbsReal& func, RooRealVar& x, const RooArgSet& nset, Int_t order=1, Double_t eps=0.001) ;
  virtual ~RooDerivative() ;

  RooDerivative(const RooDerivative& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooDerivative(*this, newname); }

  Int_t order() const { return _order ; }
  Double_t eps() const { return _eps ; }
  void setEps(Double_t e) { _eps = e ; }

  Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) ;

protected:

  Int_t _order ;                         // Derivation order
  Double_t _eps ;                        // Precision
  RooSetProxy  _nset ;                   // Normalization set (optional)
  RooRealProxy _func ;                   // Input function
  RooRealProxy _x     ;                  // Observable
  mutable RooFunctor*  _ftor ;                   //! Functor binding of RooAbsReal
  mutable ROOT::Math::RichardsonDerivator *_rd ; //! Derivator

  Double_t evaluate() const;

  ClassDef(RooDerivative,1) // Representation of derivative of any RooAbsReal
};

#endif
