/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooSimultaneous.h,v 1.42 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_SIMULTANEOUS
#define ROO_SIMULTANEOUS

//#include "THashList.h"
#include "TList.h"
#include "RooAbsPdf.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooAICRegistry.h"
#include "RooObjCacheManager.h"
#include "RooAbsCacheElement.h"
#include "RooArgList.h"
#include <map>
#include <string>
class RooAbsCategoryLValue ;
class RooFitResult ;
class RooPlot ;
class RooAbsData ;
class RooLinkedList ;

class RooSimultaneous : public RooAbsPdf {
public:

  // Constructors, assignment etc
  inline RooSimultaneous() : _plotCoefNormRange(0) { }
  RooSimultaneous(const char *name, const char *title, RooAbsCategoryLValue& indexCat) ;
  RooSimultaneous(const char *name, const char *title, std::map<std::string,RooAbsPdf*> pdfMap, RooAbsCategoryLValue& inIndexCat) ;
  RooSimultaneous(const char *name, const char *title, const RooArgList& pdfList, RooAbsCategoryLValue& indexCat) ;
  RooSimultaneous(const RooSimultaneous& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSimultaneous(*this,newname) ; }
  virtual ~RooSimultaneous() ;

  virtual Double_t evaluate() const ;
  virtual Bool_t selfNormalized() const { return kTRUE ; }
  Bool_t addPdf(const RooAbsPdf& pdf, const char* catLabel) ;

  virtual ExtendMode extendMode() const ;

  virtual Double_t expectedEvents(const RooArgSet* nset) const ;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset) ; }

  virtual Bool_t forceAnalyticalInt(const RooAbsArg&) const { return kTRUE ; }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;

  using RooAbsPdf::plotOn ;
  virtual RooPlot* plotOn(RooPlot* frame,
           const RooCmdArg& arg1            , const RooCmdArg& arg2=RooCmdArg(),
           const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
           const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
           const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg(),
           const RooCmdArg& arg9=RooCmdArg(), const RooCmdArg& arg10=RooCmdArg()) const {
    return RooAbsReal::plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10) ;
  }
  virtual RooPlot* plotOn(RooPlot* frame, RooLinkedList& cmdList) const ;

  // Backward compatibility function
  virtual RooPlot *plotOn(RooPlot *frame, Option_t* drawOptions, Double_t scaleFactor=1.0,
           ScaleType stype=Relative, const RooAbsData* projData=0, const RooArgSet* projSet=0,
           Double_t precision=1e-3, Bool_t shiftToZero=kFALSE, const RooArgSet* projDataSet=0,
           Double_t rangeLo=0, Double_t rangeHi=0, RooCurve::WingMode wmode=RooCurve::Extended) const;

  RooAbsPdf* getPdf(const char* catName) const ;
  const RooAbsCategoryLValue& indexCat() const { return (RooAbsCategoryLValue&) _indexCat.arg() ; }


  virtual RooDataSet* generateSimGlobal(const RooArgSet& whatVars, Int_t nEvents) ;

  virtual RooDataHist* fillDataHist(RooDataHist *hist, const RooArgSet* nset, Double_t scaleFactor,
                Bool_t correctForBinVolume=kFALSE, Bool_t showProgress=kFALSE) const ;

protected:

  void initialize(RooAbsCategoryLValue& inIndexCat, std::map<std::string,RooAbsPdf*> pdfMap) ;

  virtual void selectNormalization(const RooArgSet* depSet=0, Bool_t force=kFALSE) ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;
  mutable RooSetProxy _plotCoefNormSet ;
  const TNamed* _plotCoefNormRange ;

  class CacheElem : public RooAbsCacheElement {
  public:
    virtual ~CacheElem() {} ;
    RooArgList containedArgs(Action) { return RooArgList(_partIntList) ; }
    RooArgList _partIntList ;
  } ;
  mutable RooObjCacheManager _partIntMgr ; // Component normalization manager


  friend class RooSimGenContext ;
  friend class RooSimSplitGenContext ;
  virtual RooAbsGenContext* autoGenContext(const RooArgSet &vars, const RooDataSet* prototype=0, const RooArgSet* auxProto=0,
                  Bool_t verbose=kFALSE, Bool_t autoBinned=kTRUE, const char* binnedTag="") const ;
  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                  const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;

  RooCategoryProxy _indexCat ; // Index category
  TList    _pdfProxyList ;     // List of PDF proxies (named after applicable category state)
  Int_t    _numPdf ;           // Number of registered PDFs

  ClassDef(RooSimultaneous,2)  // Simultaneous operator p.d.f, functions like C++  'switch()' on input p.d.fs operating on index category5A
};

#endif
