/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooProdGenContext.h,v 1.15 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_PROD_GEN_CONTEXT
#define ROO_PROD_GEN_CONTEXT

#include "RooAbsGenContext.h"
#include "RooArgSet.h"
#include <list>

class RooProdPdf;
class RooDataSet;
class RooRealIntegral;
class RooAcceptReject;
class TRandom;
class TIterator;
class RooSuperCategory ;

class RooProdGenContext : public RooAbsGenContext {
public:
  RooProdGenContext(const RooProdPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0,
          const RooArgSet* auxProto=0, Bool_t _verbose= kFALSE);
  virtual ~RooProdGenContext();

  virtual void setProtoDataOrder(Int_t* lut) ;
  virtual void printMultiline(std::ostream &os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const ;

  virtual void attach(const RooArgSet& params) ;

protected:

  virtual void initGenerator(const RooArgSet &theEvent);
  virtual void generateEvent(RooArgSet &theEvent, Int_t remaining);

  void updateCCDTable() ;


  RooProdGenContext(const RooProdGenContext& other) ;

  RooArgSet _commonCats ;        // Common category dependents
  RooArgSet* _ccdCloneSet ;
  RooSuperCategory* _ccdSuper ;  // SuperCategory of Common category dependents
  RooArgSet* _pdfCloneSet ;
  RooAbsPdf* _pdfClone ;
  RooRealIntegral* _pdfCcdInt ;
  RooArgSet _uniObs ;            // Observable to be generated with flat distribution
  TIterator* _uniIter ;          // Iterator over uniform observables
  Bool_t _ccdRefresh ;
  Double_t * _ccdTable ;
  const RooProdPdf *_pdf ;       //  Original PDF
  std::list<RooAbsGenContext*>  _gcList ; //  List of component generator contexts
  RooArgSet _ownedMultiProds ;   //  Owned auxiliary multi-term product PDFs

  ClassDef(RooProdGenContext,0) // Context for efficient generation of a a dataset from a RooProdPdf
};

#endif
