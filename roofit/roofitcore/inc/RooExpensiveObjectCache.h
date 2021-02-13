/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_EXPENSIVE_OBJECT_CACHE
#define ROO_EXPENSIVE_OBJECT_CACHE

#include "TObject.h"
#include "RooArgSet.h"
#include "TString.h"
#include <map>

class RooExpensiveObjectCache : public TObject {
public:

  RooExpensiveObjectCache() ;
  RooExpensiveObjectCache(const RooExpensiveObjectCache&) ;
  virtual ~RooExpensiveObjectCache() ;

  Bool_t registerObject(const char* ownerName, const char* objectName, TObject& cacheObject, TIterator* paramIter) ;
  Bool_t registerObject(const char* ownerName, const char* objectName, TObject& cacheObject, const RooArgSet& params) ;
  const TObject* retrieveObject(const char* name, TClass* tclass, const RooArgSet& params) ;

  const TObject* getObj(Int_t uniqueID) ;
  Bool_t clearObj(Int_t uniqueID) ;
  Bool_t setObj(Int_t uniqueID, TObject* obj) ;
  void clearAll() ;

  void importCacheObjects(RooExpensiveObjectCache& other, const char* ownerName, Bool_t verbose=kFALSE) ;

  static RooExpensiveObjectCache& instance() ;

  Int_t size() const { return _map.size() ; }

  void print() const ;

  class ExpensiveObject {
  public:
    ExpensiveObject() { _uid = 0 ; _payload = 0 ; } ;
    ExpensiveObject(Int_t uid, const char* ownerName, TObject& payload, TIterator* paramIter) ;
    ExpensiveObject(Int_t uid, const ExpensiveObject& other) ;
    virtual ~ExpensiveObject() ;
    Bool_t matches(TClass* tc, const RooArgSet& params) ;

    Int_t uid() const { return _uid ; }
    const TObject* payload() const { return _payload ; }
    TObject* payload() { return _payload ; }
    void setPayload(TObject* obj) { _payload = obj ; }
    const char* ownerName() const { return _ownerName.Data() ; }

    void print() const;

  protected:
    
    Int_t _uid ; // Unique element ID ;
    TObject* _payload ; // Payload
    std::map<TString,Double_t> _realRefParams ; // Names and values of real-valued reference parameters
    std::map<TString,Int_t> _catRefParams ; // Names and values of discrete-valued reference parameters 
    TString _ownerName ; // Name of RooAbsArg object that is associated to cache contents
  
    ClassDef(ExpensiveObject,2) ; // Cache element containing expensive object and parameter values for which object is valid
} ;

 
protected:

  Int_t _nextUID ; 

  std::map<TString,ExpensiveObject*> _map ;
 
  
  ClassDef(RooExpensiveObjectCache,2) // Singleton class that serves as session repository for expensive objects
};

#endif
