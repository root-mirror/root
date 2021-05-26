/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooProdPdf.h,v 1.44 2007/07/16 21:04:28 wouter Exp $
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
#ifndef ROO_PROD_PDF
#define ROO_PROD_PDF

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooLinkedList.h"
#include "RooAICRegistry.h"
#include "RooObjCacheManager.h"
#include "RooCmdArg.h"

#include <vector>
#include <list>
#include <string>

namespace RooBatchCompute {
struct RunContext;
}

typedef RooArgList* pRooArgList ;
typedef RooLinkedList* pRooLinkedList ;

class RooProdPdf : public RooAbsPdf {
public:
  RooProdPdf() ;
  RooProdPdf(const char *name, const char *title, Double_t cutOff=0);
  RooProdPdf(const char *name, const char *title,
       RooAbsPdf& pdf1, RooAbsPdf& pdf2, Double_t cutOff=0) ;
  RooProdPdf(const char* name, const char* title, const RooArgList& pdfList, Double_t cutOff=0) ;
  RooProdPdf(const char* name, const char* title, const RooArgSet& fullPdfSet, const RooLinkedList& cmdArgList) ;

  RooProdPdf(const char* name, const char* title, const RooArgSet& fullPdfSet,
           const RooCmdArg& arg1            , const RooCmdArg& arg2=RooCmdArg(),
             const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
             const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
             const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg()) ;

  RooProdPdf(const char* name, const char* title,
             const RooCmdArg& arg1,             const RooCmdArg& arg2=RooCmdArg(),
             const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
             const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
             const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg()) ;

  RooProdPdf(const RooProdPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooProdPdf(*this,newname) ; }
  virtual ~RooProdPdf() ;

  virtual Bool_t checkObservables(const RooArgSet* nset) const ;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& dep) const ;
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Bool_t selfNormalized() const { return _selfNorm ; }

  virtual ExtendMode extendMode() const ;
  virtual Double_t expectedEvents(const RooArgSet* nset) const ;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset) ; }

  const RooArgList& pdfList() const { return _pdfList ; }

  virtual Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  virtual void initGenerator(Int_t code) ;
  virtual void generateEvent(Int_t code);
  virtual Bool_t isDirectGenSafe(const RooAbsArg& arg) const ;

  // Constraint management
  virtual RooArgSet* getConstraints(const RooArgSet& observables, RooArgSet& constrainedParams, Bool_t stripDisconnected) const ;

  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ;
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  Bool_t isBinnedDistribution(const RooArgSet& obs) const  ;

  void printMetaArgs(std::ostream& os) const ;

  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;
  void fixRefRange(const char* rangeName) ;

  void setSelfNormalized(Bool_t flag) { _selfNorm = flag ; }
  void setDefNormSet(const RooArgSet& nset) { _defNormSet.removeAll() ; _defNormSet.addClone(nset) ; }


  Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) ;

  RooArgSet* getConnectedParameters(const RooArgSet& observables) const ;

  RooArgSet* findPdfNSet(RooAbsPdf& pdf) const ;


private:

  Double_t evaluate() const ;
  virtual RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const;

  RooAbsReal* makeCondPdfRatioCorr(RooAbsReal& term, const RooArgSet& termNset, const RooArgSet& termImpSet, const char* normRange, const char* refRange) const ;

  virtual void getParametersHook(const RooArgSet* /*nset*/, RooArgSet* /*list*/, Bool_t stripDisconnected) const ;

  void initializeFromCmdArgList(const RooArgSet& fullPdfSet, const RooLinkedList& l) ;

  void factorizeProduct(const RooArgSet& normSet, const RooArgSet& intSet,
                        RooLinkedList& termList,   RooLinkedList& normList,
                        RooLinkedList& impDepList, RooLinkedList& crossDepList,
                        RooLinkedList& intList) const;
  std::string makeRGPPName(const char* pfx, const RooArgSet& term, const RooArgSet& iset, const RooArgSet& nset, const char* isetRangeName) const ;
  void groupProductTerms(RooLinkedList& groupedTerms, RooArgSet& outerIntDeps,
                         const RooLinkedList& terms, const RooLinkedList& norms,
                         const RooLinkedList& imps, const RooLinkedList& ints, const RooLinkedList& cross) const ;



  Int_t getPartIntList(const RooArgSet* nset, const RooArgSet* iset, const char* isetRangeName=0) const ;

  std::vector<RooAbsReal*> processProductTerm(const RooArgSet* nset, const RooArgSet* iset, const char* isetRangeName,
                     const RooArgSet* term,const RooArgSet& termNSet, const RooArgSet& termISet,
                     Bool_t& isOwned, Bool_t forceWrap=kFALSE) const ;


  virtual CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised ; } ;
  virtual void setCacheAndTrackHints(RooArgSet&) ;

  // The cache object
  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem() : _isRearranged(kFALSE) { }
    virtual ~CacheElem() = default;
    // Payload
    RooArgList _partList ;
    RooArgList _numList ;
    RooArgList _denList ;
    RooArgList _ownedList ;
    std::vector<std::unique_ptr<RooArgSet>> _normList;
    Bool_t _isRearranged ;
    std::unique_ptr<RooAbsReal> _rearrangedNum{};
    std::unique_ptr<RooAbsReal> _rearrangedDen{};
    // Cache management functions
    virtual RooArgList containedArgs(Action) ;
    virtual void printCompactTreeHook(std::ostream&, const char *, Int_t, Int_t) ;
  private:
    CacheElem(const CacheElem&) ;
  } ;
  mutable RooObjCacheManager _cacheMgr ; // The cache manager

  void rearrangeProduct(CacheElem&) const;
  RooAbsReal* specializeIntegral(RooAbsReal& orig, const char* targetRangeName) const ;
  RooAbsReal* specializeRatio(RooFormulaVar& input, const char* targetRangeName) const ;
  Double_t calculate(const RooProdPdf::CacheElem& cache, Bool_t verbose=kFALSE) const ;


  friend class RooProdGenContext ;
  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                  const RooArgSet *auxProto=0, Bool_t verbose= kFALSE) const ;


  mutable RooAICRegistry _genCode ; //! Registry of composite direct generator codes

  Double_t _cutOff ;       //  Cutoff parameter for running product
  RooListProxy _pdfList ;  //  List of PDF components
  RooLinkedList _pdfNSetList ; // List of PDF component normalization sets
  Int_t _extendedIndex ;   //  Index of extended PDF (if any)

  void useDefaultGen(Bool_t flag=kTRUE) { _useDefaultGen = flag ; }
  Bool_t _useDefaultGen ; // Use default or distributed event generator

  mutable TNamed* _refRangeName ; // Reference range name for interpretation of conditional products

  Bool_t _selfNorm ; // Is self-normalized
  RooArgSet _defNormSet ; // Default normalization set

private:

  ClassDef(RooProdPdf,4) // PDF representing a product of PDFs
};


#endif
