// Author: Stephan Hageboeck, CERN, May 2019
/*****************************************************************************
 * Project: RooFit                                                           *
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
#ifndef ROO_JOHNSON
#define ROO_JOHNSON

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooConstVar.h"

class RooRealVar;

class RooJohnson final : public RooAbsPdf {
public:
  RooJohnson() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300

  RooJohnson(const char *name, const char *title,
            RooAbsReal& mass, RooAbsReal& mu, RooAbsReal& lambda,
            RooAbsReal& gamma, RooAbsReal& delta,
            double massThreshold = -std::numeric_limits<double>::max());

  RooJohnson(const RooJohnson& other, const char* newName = nullptr);

  virtual ~RooJohnson() = default;

  TObject* clone(const char* newname) const override {
    return new RooJohnson(*this,newname);
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const override;
  void generateEvent(Int_t code) override;

private:
  enum AnaInt_t {kMass = 1, kMean, kLambda, kGamma, kDelta};

  RooRealProxy _mass;
  RooRealProxy _mu;
  RooRealProxy _lambda;

  RooRealProxy _gamma;
  RooRealProxy _delta;

  double _massThreshold{-1.E300};

  Double_t evaluate() const override;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const override;

  ClassDefOverride(RooJohnson,1)
};

#endif
