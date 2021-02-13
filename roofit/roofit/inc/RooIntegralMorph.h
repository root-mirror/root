/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2007, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOLINEARMORPH
#define ROOLINEARMORPH

#include "RooAbsCachedPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <vector>

class RooBrentRootFinder;

class TH1D;

class RooIntegralMorph : public RooAbsCachedPdf {
public:
   RooIntegralMorph() :  _cache(nullptr)  {
    // coverity[UNINIT_CTOR]
  } ;
  RooIntegralMorph(const char *name, const char *title,
         RooAbsReal& _pdf1,
         RooAbsReal& _pdf2,
           RooAbsReal& _x,
         RooAbsReal& _alpha, Bool_t cacheAlpha=kFALSE);
  RooIntegralMorph(const RooIntegralMorph& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooIntegralMorph(*this,newname); }
  inline virtual ~RooIntegralMorph() { }

  Bool_t selfNormalized() const {
    // P.d.f is self normalized
    return kTRUE ;
  }
  void setCacheAlpha(Bool_t flag) {
    // Activate caching of p.d.f. shape for all values of alpha as well
    _cacheMgr.sterilize() ; _cacheAlpha = flag ;
  }
  Bool_t cacheAlpha() const {
    // If true caching of p.d.f for all alpha values is active
    return _cacheAlpha ;
  }

  virtual void preferredObservableScanOrder(const RooArgSet& obs, RooArgSet& orderedObs) const ;

  class MorphCacheElem : public RooAbsCachedPdf::PdfCacheElem {
  public:
    MorphCacheElem(RooIntegralMorph& self, const RooArgSet* nset) ;
    ~MorphCacheElem() ;
    void calculate(TIterator* iter) ;
    virtual RooArgList containedArgs(Action) ;

  protected:

    void findRange() ;
    Double_t calcX(Double_t y, Bool_t& ok) ;
    Int_t binX(Double_t x) ;
    void fillGap(Int_t ixlo, Int_t ixhi,Double_t splitPoint=0.5) ;
    void interpolateGap(Int_t ixlo, Int_t ixhi) ;

    RooIntegralMorph* _self ; //
    RooArgSet* _nset ;
    RooAbsPdf* _pdf1 ; // PDF1
    RooAbsPdf* _pdf2 ; // PDF2
    RooRealVar* _x   ; // X
    RooAbsReal* _alpha ; // ALPHA
    RooAbsReal* _c1 ; // CDF of PDF 1
    RooAbsReal* _c2 ; // CDF of PDF 2
    RooAbsFunc* _cb1 ; // Binding of CDF1
    RooAbsFunc* _cb2 ; // Binding of CDF2
    RooBrentRootFinder* _rf1 ; // ROOT finder on CDF1
    RooBrentRootFinder* _rf2 ; // ROOT finder of CDF2 ;

    std::vector<Double_t> _yatX ; //
    std::vector<Double_t> _calcX; //
    Int_t _yatXmin, _yatXmax ;
    Int_t _ccounter ;

    Double_t _ycutoff ;

  } ;

protected:

  friend class MorphCacheElem ;
  virtual PdfCacheElem* createCache(const RooArgSet* nset) const ;
  virtual const char* inputBaseName() const ;
  virtual RooArgSet* actualObservables(const RooArgSet& nset) const ;
  virtual RooArgSet* actualParameters(const RooArgSet& nset) const ;
  virtual void fillCacheObject(PdfCacheElem& cache) const ;

  RooRealProxy pdf1 ; // First input shape
  RooRealProxy pdf2 ; // Second input shape
  RooRealProxy x ;    // Observable
  RooRealProxy alpha ; // Interpolation parameter
  Bool_t _cacheAlpha ; // If true, both (x,alpha) are cached
  mutable MorphCacheElem* _cache ; // Current morph cache element in use


  Double_t evaluate() const ;

private:

  ClassDef(RooIntegralMorph,1) // Linear shape interpolation operator p.d.f
};

#endif
