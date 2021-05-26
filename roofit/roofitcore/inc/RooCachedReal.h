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

#ifndef ROOCACHEDREAL
#define ROOCACHEDREAL

#include "RooAbsCachedReal.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"


class RooCachedReal : public RooAbsCachedReal {
public:
  RooCachedReal() : _cacheSource(kFALSE) {
    // coverity[UNINIT_CTOR]
  }
  RooCachedReal(const char *name, const char *title, RooAbsReal& _func, const RooArgSet& cacheObs);
  RooCachedReal(const char *name, const char *title, RooAbsReal& _func);
  RooCachedReal(const RooCachedReal& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCachedReal(*this,newname); }
  virtual ~RooCachedReal() ;

  void setCdfBoundaries(Bool_t flag) {
    // If flag is true the RooHistFunc that represent the cache histogram
    // will use special boundary conditions for use with cumulative distribution
    // functions: at the lower bound the function is forced to converge at zero and the upper
    // bound is the function is forced to converge at 1.0
    _useCdfBoundaries = flag ;
  }
  Bool_t getCdfBoundaries() const {
    // If true the c.d.f boundary mode is active
    return _useCdfBoundaries ;
  }

  Bool_t cacheSource() const { return _cacheSource ; }
  void setCacheSource(Bool_t flag) { _cacheSource = flag ; }

protected:

  virtual const char* inputBaseName() const {
    // Return base name for caches, i.e. the name of the cached function
    return func.arg().GetName() ;
  } ;
  virtual RooArgSet* actualObservables(const RooArgSet& nset) const ;
  virtual RooArgSet* actualParameters(const RooArgSet& nset) const ;
  virtual void fillCacheObject(FuncCacheElem& cacheFunc) const ;
  virtual Double_t evaluate() const {
    // Dummy evaluate, it is never called
    return func ;
  }

  void operModeHook() ;

  virtual FuncCacheElem* createCache(const RooArgSet* nset) const ;

  virtual const char* payloadUniqueSuffix() const { return func.arg().aggregateCacheUniqueSuffix() ; }

  RooRealProxy func ;           // Proxy to function being cached
  RooSetProxy  _cacheObs ;      // Variables to be cached
  Bool_t _useCdfBoundaries ;    // Are c.d.f boundary conditions used by the RooHistFuncs?
  Bool_t _cacheSource ;         // Keep an attached clone of the source in the cache for fast operation

private:

  ClassDef(RooCachedReal,2) // P.d.f class that wraps another p.d.f and caches its output

};

#endif
