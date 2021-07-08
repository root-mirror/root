/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooVoigtian.h,v 1.7 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   TS, Thomas Schietinger, SLAC,          schieti@slac.stanford.edu        *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_VOIGTIAN
#define ROO_VOIGTIAN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooVoigtian : public RooAbsPdf {
public:
  RooVoigtian() {} ;
  RooVoigtian(const char *name, const char *title,
         RooAbsReal& _x, RooAbsReal& _mean,
              RooAbsReal& _width, RooAbsReal& _sigma,
              Bool_t doFast = kFALSE);
  RooVoigtian(const RooVoigtian& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooVoigtian(*this,newname); }
  inline virtual ~RooVoigtian() { }

// These methods allow the user to select the fast evaluation
// of the complex error function using look-up tables
// (default is the "slow" CERNlib algorithm)

  inline void selectFastAlgorithm()    { _doFast = kTRUE;  }
  inline void selectDefaultAlgorithm() { _doFast = kFALSE; }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy sigma ;

  Double_t evaluate() const ;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;
  inline bool canComputeBatchWithCuda() const { return true; }

private:

  Bool_t _doFast;
  ClassDef(RooVoigtian,2) // Voigtian PDF (Gauss (x) BreitWigner)
};

#endif

