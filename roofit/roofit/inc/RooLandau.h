/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooLandau.h,v 1.5 2007/07/12 20:30:49 wouter Exp $
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
#ifndef ROO_LANDAU
#define ROO_LANDAU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooLandau : public RooAbsPdf {
public:
  RooLandau() {} ;
  RooLandau(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma);
  RooLandau(const RooLandau& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooLandau(*this,newname); }
  inline virtual ~RooLandau() { }

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);
  
protected:
  
  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy sigma ;
  
  Double_t evaluate() const ;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;
  inline bool canComputeBatchWithCuda() const { return true; }
  
private:
  
  ClassDef(RooLandau,1) // Landau Distribution PDF
};

#endif
