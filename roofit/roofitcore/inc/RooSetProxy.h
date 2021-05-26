/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooSetProxy.h,v 1.21 2007/07/13 21:24:36 wouter Exp $
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
#ifndef ROO_SET_PROXY
#define ROO_SET_PROXY

#include "RooAbsProxy.h"
#include "RooAbsArg.h"
#include "RooArgSet.h"

class RooSetProxy final : public RooArgSet, public RooAbsProxy  {
public:

#ifdef USEMEMPOOLFORARGSET
  void* operator new (size_t bytes);
  void operator delete (void *ptr);
#endif

  // Constructors, assignment etc.
  RooSetProxy() :
  _owner{nullptr},
  _defValueServer{false},
  _defShapeServer{false} { }

  RooSetProxy(const char* name, const char* desc, RooAbsArg* owner,
         Bool_t defValueServer=kTRUE, Bool_t defShapeServer=kFALSE) ;
  RooSetProxy(const char* name, RooAbsArg* owner, const RooSetProxy& other) ;
  virtual ~RooSetProxy() ;

  virtual const char* name() const override { return GetName() ; }

  // List content management (modified for server hooks)
  using RooAbsCollection::add;
  virtual Bool_t add(const RooAbsArg& var, Bool_t silent=kFALSE) override;
  virtual Bool_t add(const RooAbsArg& var, Bool_t valueServer, Bool_t shapeServer, Bool_t silent) ;

  using RooAbsCollection::addOwned;
  virtual Bool_t addOwned(RooAbsArg& var, Bool_t silent=kFALSE) override;

  using RooAbsCollection::addClone;
  virtual RooAbsArg *addClone(const RooAbsArg& var, Bool_t silent=kFALSE) override;

  virtual Bool_t replace(const RooAbsArg& var1, const RooAbsArg& var2) override;
  virtual Bool_t remove(const RooAbsArg& var, Bool_t silent=kFALSE, Bool_t matchByNameOnly=kFALSE) override;
  Bool_t remove(const RooAbsCollection& list, Bool_t silent=kFALSE, Bool_t matchByNameOnly=kFALSE) ;
  virtual void removeAll() override;

  virtual void print(std::ostream& os, Bool_t addContents=kFALSE) const override;

  RooSetProxy& operator=(const RooArgSet& other) ;

protected:

  RooAbsArg* _owner ;
  Bool_t _defValueServer ;
  Bool_t _defShapeServer ;

  virtual Bool_t changePointer(const RooAbsCollection& newServerSet, Bool_t nameChange=kFALSE, Bool_t factoryInitMode=kFALSE) override;

  ClassDefOverride(RooSetProxy,1) // Proxy class for a RooArgSet
};

#endif

