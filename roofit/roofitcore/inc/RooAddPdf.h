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
#ifndef ROO_ADD_PDF
#define ROO_ADD_PDF

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooSetProxy.h"
#include "RooAICRegistry.h"
#include "RooNormSetCache.h"
#include "RooNameSet.h"
#include "RooObjCacheManager.h"
#include "RooNameReg.h"

#include <vector>
#include <list>
#include <utility>

class RooAddPdf : public RooAbsPdf {
public:

  RooAddPdf() ;
  RooAddPdf(const char *name, const char *title=0);
  RooAddPdf(const char *name, const char *title,
       RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooAbsReal& coef1) ;
  RooAddPdf(const char *name, const char *title, const RooArgList& pdfList) ;
  RooAddPdf(const char *name, const char *title, const RooArgList& pdfList, const RooArgList& coefList, Bool_t recursiveFraction=kFALSE) ;

  RooAddPdf(const RooAddPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooAddPdf(*this,newname) ; }
  virtual ~RooAddPdf() ;

  virtual Bool_t checkObservables(const RooArgSet* nset) const ;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const {
    // Force RooRealIntegral to offer all observables for internal integration
    return kTRUE ;
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Bool_t selfNormalized() const {
    // P.d.f is self normalized
    return kTRUE ;
  }

  virtual ExtendMode extendMode() const {
    // Return extended mode capabilities
    return ((_haveLastCoef&&!_recursive) || _allExtendable) ? MustBeExtended : CanNotBeExtended;
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

  void fixCoefNormalization(const RooArgSet& refCoefNorm) ;
  void fixCoefRange(const char* rangeName) ;

  const RooArgSet& getCoefNormalization() const { return _refCoefNorm ; }
  const char* getCoefRange() const { return _refCoefRangeName?RooNameReg::str(_refCoefRangeName):"" ; }

  virtual Double_t getValV(const RooArgSet *set = 0) const;

  virtual void resetErrorCounters(Int_t resetValue=10) ;

  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ;
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  Bool_t isBinnedDistribution(const RooArgSet& obs) const  ;

  void printMetaArgs(std::ostream& os) const ;

  virtual CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised ; } ;
  virtual void setCacheAndTrackHints(RooArgSet&) ;

protected:

  virtual void selectNormalization(const RooArgSet* depSet=0, Bool_t force=kFALSE) ;
  virtual void selectNormalizationRange(const char* rangeName=0, Bool_t force=kFALSE) ;

  mutable RooSetProxy _refCoefNorm ;   // Reference observable set for coefficient interpretation
  mutable TNamed* _refCoefRangeName ;  // Reference range name for coefficient interpreation

  Bool_t _projectCoefs ;         // If true coefficients need to be projected for use in evaluate()
  std::vector<double> _coefCache; //! Transient cache with transformed values of coefficients


  class CacheElem : public RooAbsCacheElement {
  public:
    virtual ~CacheElem() {} ;

    RooArgList _suppNormList ; // Supplemental normalization list
    Bool_t    _needSupNorm ; // Does the above list contain any non-unit entries?

    RooArgList _projList ; // Projection integrals to be multiplied with coefficients
    RooArgList _suppProjList ; // Projection integrals to be multiplied with coefficients for supplemental normalization terms
    RooArgList _refRangeProjList ; // Range integrals to be multiplied with coefficients (reference range)
    RooArgList _rangeProjList ; // Range integrals to be multiplied with coefficients (target range)

    virtual RooArgList containedArgs(Action) ;

  } ;
  mutable RooObjCacheManager _projCacheMgr ;  // Manager of cache with coefficient projections and transformations
  CacheElem* getProjCache(const RooArgSet* nset, const RooArgSet* iset=0, const char* rangeName=0) const ;
  void updateCoefficients(CacheElem& cache, const RooArgSet* nset) const ;


  friend class RooAddGenContext ;
  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype=0,
                                       const RooArgSet* auxProto=0, Bool_t verbose= kFALSE) const ;


  Double_t evaluate() const;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const;


  mutable RooAICRegistry _codeReg ;  //! Registry of component analytical integration codes

  RooListProxy _pdfList ;   //  List of component PDFs
  RooListProxy _coefList ;  //  List of coefficients
  mutable RooArgList* _snormList{nullptr};  //!  List of supplemental normalization factors

  Bool_t _haveLastCoef ;    //  Flag indicating if last PDFs coefficient was supplied in the ctor
  Bool_t _allExtendable ;   //  Flag indicating if all PDF components are extendable
  Bool_t _recursive ;       //  Flag indicating is fractions are treated recursively

  mutable Int_t _coefErrCount ; //! Coefficient error counter

private:
  std::pair<const RooArgSet*, CacheElem*> getNormAndCache(const RooArgSet* defaultNorm = nullptr) const;

  ClassDef(RooAddPdf,3) // PDF representing a sum of PDFs
};

#endif
