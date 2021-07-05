/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooCBShape.h,v 1.11 2007/07/12 20:30:49 wouter Exp $
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
#ifndef ROO_CB_SHAPE
#define ROO_CB_SHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooCBShape : public RooAbsPdf {
public:
  RooCBShape() {} ;
  RooCBShape(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0, RooAbsReal& _sigma,
        RooAbsReal& _alpha, RooAbsReal& _n);

  RooCBShape(const RooCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooCBShape(*this,newname); }

  inline virtual ~RooCBShape() { }

  virtual Int_t getAnalyticalIntegral( RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0 ) const;
  virtual Double_t analyticalIntegral( Int_t code, const char* rangeName=0 ) const;

  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;

protected:

  Double_t ApproxErf(Double_t arg) const ;

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alpha;
  RooRealProxy n;

  Double_t evaluate() const;
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const;


private:

  ClassDef(RooCBShape,1) // Crystal Ball lineshape PDF
};

#endif
