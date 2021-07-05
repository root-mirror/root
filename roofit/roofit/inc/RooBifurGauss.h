/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBifurGauss.h,v 1.12 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   Abi Soffer, Colorado State University, abi@slac.stanford.edu            *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California,         *
 *                          Colorado State University                        *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BIFUR_GAUSS
#define ROO_BIFUR_GAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooBifurGauss : public RooAbsPdf {
public:
  RooBifurGauss() {} ;
  RooBifurGauss(const char *name, const char *title, RooAbsReal& _x,
      RooAbsReal& _mean, RooAbsReal& _sigmaL, RooAbsReal& _sigmaR);

  RooBifurGauss(const RooBifurGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBifurGauss(*this,newname); }
  inline virtual ~RooBifurGauss() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;


protected:

  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;

  Double_t evaluate() const;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;

private:

  ClassDef(RooBifurGauss,1) // Bifurcated Gaussian PDF
};

#endif
