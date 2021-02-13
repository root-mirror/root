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
#ifndef ROO_OBJ_CACHE_MANAGER
#define ROO_OBJ_CACHE_MANAGER

#include "Rtypes.h"

#include "RooNormSetCache.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsCache.h"
#include "RooAbsCacheElement.h"
#include "RooCacheManager.h"

class RooNameSet;


class RooObjCacheManager : public RooCacheManager<RooAbsCacheElement> {

public:

  RooObjCacheManager(RooAbsArg* owner=0, Int_t maxSize=2, Bool_t clearCacheOnServerRedirect=kTRUE, Bool_t allowOptimize=kFALSE) ;
  RooObjCacheManager(const RooObjCacheManager& other, RooAbsArg* owner=0) ;
  virtual ~RooObjCacheManager() ;
  
  virtual Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) ;
  virtual void operModeHook() ;
  virtual void optimizeCacheMode(const RooArgSet& /*obs*/, RooArgSet& /*optSet*/, RooLinkedList& /*processedNodes*/) ;
  virtual void printCompactTreeHook(std::ostream&, const char *) ;
  virtual void findConstantNodes(const RooArgSet& /*obs*/, RooArgSet& /*cacheList*/, RooLinkedList& /*processedNodes*/) ;

  virtual void insertObjectHook(RooAbsCacheElement&) ;

  void sterilize() ;

  static void doClearObsList(Bool_t flag) { _clearObsList = flag ; }
  static Bool_t clearObsList() { return _clearObsList ; }

  void setClearOnRedirect(Bool_t flag) { _clearOnRedirect = flag ; }
 
protected:

  Bool_t _clearOnRedirect ;
  Bool_t _allowOptimize ; 
  Bool_t _optCacheModeSeen  ;              //! 

  RooArgSet* _optCacheObservables ;        //! current optCacheObservables 

  static Bool_t _clearObsList ; // Clear obslist on sterilize?
  
  ClassDef(RooObjCacheManager,3) // Cache manager for generic caches that contain RooAbsArg objects
} ;



#endif 
