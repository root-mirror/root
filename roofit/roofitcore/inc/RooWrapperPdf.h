// Author: Stephan Hageboeck, CERN
/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2018, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_WRAPPER_PDF
#define ROO_WRAPPER_PDF

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooAbsPdf.h"
#include <list>

class RooWrapperPdf final : public RooAbsPdf {
public:

  RooWrapperPdf() { };
  /// Construct a new RooWrapperPdf.
  /// \param[in] name A name to identify this object.
  /// \param[in] title Title (for e.g. plotting)
  /// \param[in] inputFunction Any RooAbsReal that should be converted into a PDF. Although it's possible
  /// to pass a PDF, it only makes sense for non-PDF functions.
  RooWrapperPdf(const char *name, const char *title, RooAbsReal& inputFunction) :
    RooAbsPdf(name, title),
    _func("inputFunction", "Function to be converted into a PDF", this, inputFunction) { }
  virtual ~RooWrapperPdf() {};

  RooWrapperPdf(const RooWrapperPdf& other, const char* name = 0) :
    RooAbsPdf(other, name),
    _func("inputFunction", this, other._func) { }

  virtual TObject* clone(const char* newname) const override {
    return new RooWrapperPdf(*this, newname);
  }

  // Analytical Integration handling
  Bool_t forceAnalyticalInt(const RooAbsArg& dep) const override {
    return _func.arg().forceAnalyticalInt(dep);
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet,
      const char* rangeName=0) const override {
    return _func.arg().getAnalyticalIntegralWN(allVars, analVars, normSet, rangeName);
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& numVars,
      const char* rangeName=0) const override {
    return _func.arg().getAnalyticalIntegral(allVars, numVars, rangeName);
  }
  double analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const override {
    return _func.arg().analyticalIntegralWN(code, normSet, rangeName);
  }
  double analyticalIntegral(Int_t code, const char* rangeName=0) const override {
    return _func.arg().analyticalIntegral(code, rangeName);
  }


  // Internal toy generation. Since our _func is not a PDF (if it is, it doesn't make sense to use this wrapper),
  // we cannot do anything.
  /// Get specialised generator. Since the underlying function is not a PDF, this will always return zero.
//  Int_t getGenerator(const RooArgSet& /*directVars*/, RooArgSet& /*generateVars*/,
//      bool /*staticInitOK = true*/) const override { return 0; }
//  void initGenerator(Int_t /*code*/) override { }
//  void generateEvent(Int_t /*code*/) override { }
//  Bool_t isDirectGenSafe(const RooAbsArg& /*arg*/) const override { return false; }


  // Hints for optimized brute-force sampling
  Int_t getMaxVal(const RooArgSet& vars) const override {
    return _func.arg().getMaxVal(vars);
  }
  Double_t maxVal(Int_t code) const override {
    return _func.arg().maxVal(code);
  }
  Int_t minTrialSamples(const RooArgSet& arGenObs) const override {
    return _func.arg().minTrialSamples(arGenObs);
  }

  // Plotting and binning hints
  Bool_t isBinnedDistribution(const RooArgSet& obs) const override {
    return _func.arg().isBinnedDistribution(obs);
  }
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override {
    return _func.arg().binBoundaries(obs, xlo, xhi);
  }
  std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override {
    return _func.arg().plotSamplingHint(obs, xlo, xhi);
  }



private:
  RooRealProxy _func;

  double evaluate() const override {
    return _func;
  }

  ClassDefOverride(RooWrapperPdf,1)
};

#endif
