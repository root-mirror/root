/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooDstD0BG.h,v 1.12 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   UE, Ulrik Egede,     RAL,               U.Egede@rl.ac.uk                *
 *   MT, Max Turri,       UC Santa Cruz      turri@slac.stanford.edu         *
 *   CC, Chih-hsiang Cheng, Stanford         chcheng@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          RAL and Stanford University. All rights reserved.*
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_DstD0_BG
#define ROO_DstD0_BG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooDstD0BG : public RooAbsPdf {
public:
  RooDstD0BG() {} ;
  RooDstD0BG(const char *name, const char *title,
        RooAbsReal& _dm, RooAbsReal& _dm0, RooAbsReal& _c,
        RooAbsReal& _a, RooAbsReal& _b);

  RooDstD0BG(const RooDstD0BG& other, const char *name=0) ;
  virtual TObject *clone(const char *newname) const {
    return new RooDstD0BG(*this,newname); }
  inline virtual ~RooDstD0BG() { };

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy dm ;
  RooRealProxy dm0 ;
  RooRealProxy C,A,B ;

  Double_t evaluate() const;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const;

private:

  ClassDef(RooDstD0BG,1) // D*-D0 mass difference background PDF
};

#endif
