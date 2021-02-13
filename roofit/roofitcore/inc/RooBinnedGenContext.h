/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_BINNED_GEN_CONTEXT
#define ROO_BINNED_GEN_CONTEXT

#include "RooAbsGenContext.h"
#include "RooArgSet.h"


class RooDataSet;
class RooRealIntegral;
class RooAcceptReject;
class RooDataHist;
class TRandom;
class TIterator;

class RooBinnedGenContext : public RooAbsGenContext {
public:
  RooBinnedGenContext(const RooAbsPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0,
                   const RooArgSet* auxProto=0, Bool_t _verbose= kFALSE);
  virtual ~RooBinnedGenContext();
  
  RooDataSet* generate(Double_t nEvents=0, Bool_t skipInit=kFALSE, Bool_t extendedMode=kFALSE) ;

  virtual void setProtoDataOrder(Int_t*)  {}

  virtual void attach(const RooArgSet& params) ;

  virtual void printMultiline(std::ostream &os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const ;

  virtual void setExpectedData(Bool_t) ;

protected:

  virtual void initGenerator(const RooArgSet &theEvent);
  virtual void generateEvent(RooArgSet &theEvent, Int_t remaining);

  RooBinnedGenContext(const RooBinnedGenContext& other) ;

  const RooArgSet* _vars ;       
  RooArgSet* _pdfSet ;           //  Set owned all nodes of internal clone of p.d.f
  RooAbsPdf *_pdf ;              //  Pointer to cloned p.d.f
  RooDataHist* _hist ;          // Histogram 
  Bool_t _expectedData ;        // Asimov?

  ClassDef(RooBinnedGenContext,0) // Specialized context for generating a dataset from a binned pdf
};

#endif
