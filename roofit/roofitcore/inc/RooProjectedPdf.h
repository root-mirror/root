/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOPROJECTEDPDF
#define ROOPROJECTEDPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooObjCacheManager.h"
#include "RooSetProxy.h" 

class RooProjectedPdf : public RooAbsPdf {
public:

  RooProjectedPdf() ;
  RooProjectedPdf(const char *name, const char *title,  RooAbsReal& _intpdf, const RooArgSet& intObs);
  RooProjectedPdf(const RooProjectedPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooProjectedPdf(*this,newname); }
  inline virtual ~RooProjectedPdf() { }

  // Analytical integration support
  virtual Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;
  virtual Bool_t forceAnalyticalInt(const RooAbsArg& dep) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void initGenerator(Int_t /*code*/) {} ; // optional pre-generation initialization
  void generateEvent(Int_t code);

  virtual Bool_t selfNormalized() const { return kTRUE ; }

  // Handle projection of projection explicitly
  virtual RooAbsPdf* createProjection(const RooArgSet& iset) ;  

  void printMetaArgs(std::ostream& os) const ;


protected:

  RooRealProxy intpdf ; // p.d.f that is integrated
  RooSetProxy intobs ;  // observables that p.d.f is integrated over
  RooSetProxy deps ;    // dependents of this p.d.f

  class CacheElem : public RooAbsCacheElement {
  public:
    virtual ~CacheElem() { delete _projection ; } ;
    // Payload
    RooAbsReal* _projection ;
    // Cache management functions
    virtual RooArgList containedArgs(Action) ; 
    virtual void printCompactTreeHook(std::ostream&, const char *, Int_t, Int_t) ;
  } ;
  mutable RooObjCacheManager _cacheMgr ; //! The cache manager

  Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) ;

  const RooAbsReal* getProjection(const RooArgSet* iset, const RooArgSet* nset, const char* rangeName, int& code) const ;
  Double_t evaluate() const ;

private:

  ClassDef(RooProjectedPdf,1) // Operator p.d.f calculating projection of another p.d.f
};
 
#endif
