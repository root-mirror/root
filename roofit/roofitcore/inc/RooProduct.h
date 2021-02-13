/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooProduct.h,v 1.5 2007/05/11 09:11:30 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   GR, Gerhard Raven,   VU Amsterdan,     graven@nikhef.nl                 *
 *                                                                           *
 * Copyright (c) 2000-2007, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_PRODUCT
#define ROO_PRODUCT

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooCacheManager.h"
#include "RooObjCacheManager.h"

#include <list>

class RooRealVar;
class RooArgList;


class RooProduct : public RooAbsReal {
public:

  RooProduct() ;
  RooProduct(const char *name, const char *title, const RooArgList& _prodSet) ;

  RooProduct(const RooProduct& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooProduct(*this, newname); }
  virtual Bool_t forceAnalyticalInt(const RooAbsArg& dep) const ;
  virtual Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars,
                                                   const RooArgSet* normSet,
                                                   const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;


  RooArgList components() { RooArgList tmp(_compRSet) ; tmp.add(_compCSet) ; return tmp ; }

  virtual ~RooProduct() ;

  class ProdMap ;

  void printMetaArgs(std::ostream& os) const ;

  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const ;

  virtual CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised ; } ;
  virtual void setCacheAndTrackHints(RooArgSet&) ;

protected:

  RooListProxy _compRSet ;
  RooListProxy _compCSet ;

  class CacheElem : public RooAbsCacheElement {
  public:
      virtual ~CacheElem();
      // Payload
      RooArgList _prodList ;
      RooArgList _ownedList ;
      virtual RooArgList containedArgs(Action) ;
  };
  mutable RooObjCacheManager _cacheMgr ; // The cache manager
                                                                                                                                                             

  Double_t calculate(const RooArgList& partIntList) const;
  Double_t evaluate() const;
  const char* makeFPName(const char *pfx,const RooArgSet& terms) const ;
  ProdMap* groupProductTerms(const RooArgSet&) const;
  Int_t getPartIntList(const RooArgSet* iset, const char *rangeName=0) const;
    



  ClassDef(RooProduct,2) // Product of RooAbsReal and/or RooAbsCategory terms
};

#endif
