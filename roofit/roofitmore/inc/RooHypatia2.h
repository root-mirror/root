// Author: Stephan Hageboeck, CERN, Oct 2019
// Based on RooIpatia2 by Diego Martinez Santos, Nikhef, Diego.Martinez.Santos@cern.ch
/*****************************************************************************
 * Project: RooFit                                                           *
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOHYPATIA2
#define ROOHYPATIA2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooAbsReal;
 
class RooHypatia2 : public RooAbsPdf {
public:
  RooHypatia2() {} ;
  RooHypatia2(const char *name, const char *title,
	      RooAbsReal& x, RooAbsReal& lambda, RooAbsReal& zeta, RooAbsReal& beta,
	      RooAbsReal& sigma, RooAbsReal& mu, RooAbsReal& a, RooAbsReal& n, RooAbsReal& a2, RooAbsReal& n2);
  RooHypatia2(const RooHypatia2& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override { return new RooHypatia2(*this,newname); }
  inline virtual ~RooHypatia2() { }

  /* Analytical integrals need testing.

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  double analyticalIntegral(Int_t code, const char* rangeName=0) const override;

  */

  
private:
  RooRealProxy _x;
  RooRealProxy _lambda;
  RooRealProxy _zeta;
  RooRealProxy _beta;
  RooRealProxy _sigma;
  RooRealProxy _mu;
  RooRealProxy _a;
  RooRealProxy _n;
  RooRealProxy _a2;
  RooRealProxy _n2;

  Double_t evaluate() const override;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const override;

  /// \cond CLASS_DEF_DOXY
  ClassDefOverride(RooHypatia2, 1);
  /// \endcond
};
 
#endif
