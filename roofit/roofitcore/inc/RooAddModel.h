/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAddPdf.h,v 1.46 2007/07/12 20:30:28 wouter Exp $
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
#ifndef ROO_ADD_MODEL
#define ROO_ADD_MODEL

#include "RooResolutionModel.h"
#include "RooListProxy.h"
#include "RooSetProxy.h"
#include "RooAICRegistry.h"
#include "RooNormSetCache.h"
#include "RooNameSet.h"
#include "RooObjCacheManager.h"

class RooAddModel : public RooResolutionModel {
public:

  RooAddModel() ;
  RooAddModel(const char *name, const char *title, const RooArgList& pdfList, const RooArgList& coefList, Bool_t ownPdfList=kFALSE) ;
  RooAddModel(const RooAddModel& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooAddModel(*this,newname) ; }
  virtual RooResolutionModel* convolution(RooFormulaVar* basis, RooAbsArg* owner) const ;
  virtual ~RooAddModel() ;

  Double_t evaluate() const ;
  virtual Bool_t checkObservables(const RooArgSet* nset) const ;

  virtual Int_t basisCode(const char* name) const ;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const {
    // Force RooRealIntegral to offer all observables for internal integration
    return kTRUE ;
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Bool_t selfNormalized() const {
    // Model is self normalized when used as p.d.f
    return _basisCode==0 ? kTRUE : kFALSE ;
  }

  virtual ExtendMode extendMode() const {
    // Return extended mode capabilities
    return (_haveLastCoef || _allExtendable) ? MustBeExtended : CanNotBeExtended;
  }
  virtual Double_t expectedEvents(const RooArgSet* nset) const ;
  virtual Double_t expectedEvents(const RooArgSet& nset) const {
    // Return expected number of events for extended likelihood calculation
    // which is the sum of all coefficients
    return expectedEvents(&nset) ;
  }

  const RooArgList& pdfList() const {
    // Return list of component p.d.fs
    return _pdfList ;
  }
  const RooArgList& coefList() const {
    // Return list of coefficients of component p.d.f.s
    return _coefList ;
  }

  Bool_t isDirectGenSafe(const RooAbsArg& arg) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);


  void fixCoefNormalization(const RooArgSet& refCoefNorm) ;
  void fixCoefRange(const char* rangeName) ;
  virtual void resetErrorCounters(Int_t resetValue=10) ;

  void printMetaArgs(std::ostream& os) const ;

protected:

  friend class RooAddGenContext ;
  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                       const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;

  virtual void selectNormalization(const RooArgSet* depSet=0, Bool_t force=kFALSE) ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;

  mutable RooSetProxy _refCoefNorm ;   //! Reference observable set for coefficient interpretation
  mutable TNamed* _refCoefRangeName ;  //! Reference range name for coefficient interpretation

  Bool_t _projectCoefs ;         // If true coefficients need to be projected for use in evaluate()
  mutable Double_t* _coefCache ; //! Transiet cache with transformed values of coefficients


  class CacheElem : public RooAbsCacheElement {
  public:
    virtual ~CacheElem() {} ;

    RooArgList _suppNormList ; // Supplemental normalization list

    RooArgList _projList ; // Projection integrals to be multiplied with coefficients
    RooArgList _suppProjList ; // Projection integrals to be multiplied with coefficients for supplemental normalization terms
    RooArgList _refRangeProjList ; // Range integrals to be multiplied with coefficients (reference range)
    RooArgList _rangeProjList ; // Range integrals to be multiplied with coefficients (target range)

    virtual RooArgList containedArgs(Action) ;

  } ;
  mutable RooObjCacheManager _projCacheMgr ;  // Manager of cache with coefficient projections and transformations
  CacheElem* getProjCache(const RooArgSet* nset, const RooArgSet* iset=0, const char* rangeName=0) const ;
  void updateCoefficients(CacheElem& cache, const RooArgSet* nset) const ;

  typedef RooArgList* pRooArgList ;
  void getCompIntList(const RooArgSet* nset, const RooArgSet* iset, pRooArgList& compIntList, Int_t& code, const char* isetRangeName) const ;
  class IntCacheElem : public RooAbsCacheElement {
  public:
    virtual ~IntCacheElem() {} ;
    RooArgList _intList ; // List of component integrals
    virtual RooArgList containedArgs(Action) ;
  } ;

  mutable RooObjCacheManager _intCacheMgr ; // Manager of cache with integrals

  mutable RooAICRegistry _codeReg ;  //! Registry of component analytical integration codes

  RooListProxy _pdfList ;   //  List of component PDFs
  RooListProxy _coefList ;  //  List of coefficients
  mutable RooArgList* _snormList{nullptr};  //!  List of supplemental normalization factors

  Bool_t _haveLastCoef ;    //  Flag indicating if last PDFs coefficient was supplied in the ctor
  Bool_t _allExtendable ;   //  Flag indicating if all PDF components are extendable

  mutable Int_t _coefErrCount ; //! Coefficient error counter

  mutable RooArgSet _ownedComps ; //! Owned components

private:

  ClassDef(RooAddModel,1) // Resolution model representing a sum of resolution models
};

#endif
