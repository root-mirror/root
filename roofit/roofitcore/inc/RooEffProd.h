/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooEffProd.h,v 1.2 2007/05/11 10:14:56 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven, NIKHEF/VU                                            *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_EFF_PROD
#define ROO_EFF_PROD

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooObjCacheManager.h"

class RooEffProd: public RooAbsPdf {
public:
  // Constructors, assignment etc
  inline RooEffProd() : _nset(0), _fixedNset(0) { };
  virtual ~RooEffProd();
  RooEffProd(const char *name, const char *title, RooAbsPdf& pdf, RooAbsReal& efficiency);
  RooEffProd(const RooEffProd& other, const char* name=0);

  virtual TObject* clone(const char* newname) const { return new RooEffProd(*this,newname); }

  virtual RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                       const RooArgSet* auxProto, Bool_t verbose) const;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg& /*dep*/) const {
    // Return kTRUE to force RooRealIntegral to offer all observables for internal integration
    return kTRUE ;
  }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;

protected:

  const RooAbsPdf* pdf() const {
    // Return pointer to pdf in product
    return (RooAbsPdf*) _pdf.absArg() ;
  }
  const RooAbsReal* eff() const {
    // Return pointer to efficiency function in product
    return (RooAbsReal*) _eff.absArg() ;
  }

  // Function evaluation
  virtual Double_t evaluate() const ;

  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem() : _clone(0), _int(0) {}
    virtual ~CacheElem() { delete _int ; delete _clone ; }
    // Payload
    RooArgSet   _intObs ;
    RooEffProd* _clone ;
    RooAbsReal* _int ;
    // Cache management functions
    virtual RooArgList containedArgs(Action) ;
  } ;
  mutable RooObjCacheManager _cacheMgr ; //! The cache manager


  // the real stuff...
  RooRealProxy _pdf ;               // Probability Density function
  RooRealProxy _eff;                // Efficiency function
  mutable const RooArgSet* _nset  ; //! Normalization set to be used in evaluation

  RooArgSet* _fixedNset ; //! Fixed normalization set overriding default normalization set (if provided)

  ClassDef(RooEffProd,2) // Product operator p.d.f of (PDF x efficiency) implementing optimized generator context
};

#endif
