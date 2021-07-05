/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooExponential.h,v 1.10 2007/07/12 20:30:49 wouter Exp $
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
#ifndef ROO_EXPONENTIAL
#define ROO_EXPONENTIAL

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooExponential : public RooAbsPdf {
public:
  RooExponential() {} ;
  RooExponential(const char *name, const char *title,
       RooAbsReal& _x, RooAbsReal& _c);
  RooExponential(const RooExponential& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override { return new RooExponential(*this,newname); }
  inline virtual ~RooExponential() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

protected:
  RooRealProxy x;
  RooRealProxy c;

  Double_t evaluate() const override;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const override;
  
private:
  ClassDefOverride(RooExponential,1) // Exponential PDF
};

#endif
