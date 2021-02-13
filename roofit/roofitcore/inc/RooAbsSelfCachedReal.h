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

#ifndef ROOABSSELFCACHEDREAL
#define ROOABSSELFCACHEDREAL

#include "RooAbsCachedReal.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooHistPdf.h"


class RooAbsSelfCachedReal : public RooAbsCachedReal {
public:

  RooAbsSelfCachedReal() {} ;
  RooAbsSelfCachedReal(const char *name, const char *title, Int_t ipOrder=0);
  RooAbsSelfCachedReal(const RooAbsSelfCachedReal& other, const char* name=0) ;
  virtual ~RooAbsSelfCachedReal() ;

protected:

  virtual const char* inputBaseName() const { 
    // Use own name as base name for caches
    return GetName() ; 
  }
  virtual RooArgSet* actualObservables(const RooArgSet& nset) const ;
  virtual RooArgSet* actualParameters(const RooArgSet& nset) const ;
  virtual void fillCacheObject(FuncCacheElem& cache) const ;  

private:

  ClassDef(RooAbsSelfCachedReal,0) // Abstract base class for self-caching functions
};
 
#endif
