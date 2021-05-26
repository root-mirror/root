/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooArgProxy.h,v 1.21 2007/07/12 20:30:28 wouter Exp $
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
#ifndef ROO_ARG_PROXY
#define ROO_ARG_PROXY

#include "TNamed.h"
#include "RooAbsProxy.h"
#include "RooAbsArg.h"


class RooArgProxy : public TNamed, public RooAbsProxy  {
public:

  // Constructors, assignment etc.
  RooArgProxy() : _owner(0), _arg(0), _valueServer(kFALSE), _shapeServer(kFALSE), _isFund(kTRUE), _ownArg(kFALSE) {
    // Default constructor
  } ;
  RooArgProxy(const char* name, const char* desc, RooAbsArg* owner,
         Bool_t valueServer, Bool_t shapeServer, Bool_t proxyOwnsArg=kFALSE) ;
  RooArgProxy(const char* name, const char* desc, RooAbsArg* owner, RooAbsArg& arg,
         Bool_t valueServer, Bool_t shapeServer, Bool_t proxyOwnsArg=kFALSE) ;
  RooArgProxy(const char* name, RooAbsArg* owner, const RooArgProxy& other) ;
  virtual ~RooArgProxy() ;
  inline RooAbsArg* absArg() const {
    // Return pointer to contained argument
    return _arg ;
  }

  virtual const char* name() const {
    // Return name of proxy
    return GetName() ;
  }
  virtual void print(std::ostream& os, Bool_t addContents=kFALSE) const ;

protected:

  friend class RooSimultaneous ;
  RooAbsArg* _owner ;       // Pointer to owner of proxy
  RooAbsArg* _arg ;         // Pointer to content of proxy

  Bool_t _valueServer ;     // If true contents is value server of owner
  Bool_t _shapeServer ;     // If true contents is shape server of owner
  Bool_t _isFund ;          // If true proxy contains an lvalue
  Bool_t _ownArg ;          // If true proxy owns contents

  friend class RooAbsArg ;

  inline Bool_t isValueServer() const {
    // Returns true of contents is value server of owner
    return _valueServer ;
  }
  inline Bool_t isShapeServer() const {
    // Returns true if contents is shape server of owner
    return _shapeServer ;
  }
  virtual Bool_t changePointer(const RooAbsCollection& newServerSet, Bool_t nameChange=kFALSE, Bool_t factoryInitMode=kFALSE) ;

  virtual void changeDataSet(const RooArgSet* newNormSet) ;

  ClassDef(RooArgProxy,1) // Abstract proxy for RooAbsArg objects
};

#endif

