/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBreitWigner.h,v 1.9 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   AS, Abi Soffer, Colorado State University, abi@slac.stanford.edu        *
 *   TS, Thomas Schietinger, SLAC, schieti@slac.stanford.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          Colorado State University                        *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BREITWIGNER
#define ROO_BREITWIGNER

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooBreitWigner : public RooAbsPdf {
public:
  RooBreitWigner() {} ;
  RooBreitWigner(const char *name, const char *title,
         RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _width);
  RooBreitWigner(const RooBreitWigner& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBreitWigner(*this,newname); }
  inline virtual ~RooBreitWigner() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;

  Double_t evaluate() const ;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;
  inline bool canComputeBatchWithCuda() const { return true; }

//   void initGenerator();
//   Int_t generateDependents();

private:

  ClassDef(RooBreitWigner,1) // Breit Wigner PDF
};

#endif
