/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooSharedProperties.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_SHARED_PROPERTY
#define ROO_ABS_SHARED_PROPERTY

#include "TObject.h"
#include "TUUID.h"
#include "TString.h"

class RooSharedProperties : public TObject {
public:

  RooSharedProperties() ;
  RooSharedProperties(const char* uuidstr) ;
  virtual ~RooSharedProperties() ;
  Bool_t operator==(const RooSharedProperties& other) const ;

  virtual RooSharedProperties* clone() = 0 ;

  virtual void Print(Option_t* opts=0) const ;

  void increaseRefCount() { _refCount++ ; }
  void decreaseRefCount() { if (_refCount>0) _refCount-- ; }
  Int_t refCount() const { return _refCount ; }

  void setInSharedList() { _inSharedList = kTRUE ; }
  Bool_t inSharedList() const { return _inSharedList ; }

   TString asString() const { return TString(_uuid.AsString()); }

protected:

  TUUID _uuid ; // Unique object ID
  Int_t _refCount ; //! Use count 
  Int_t _inSharedList ; //! Is in shared list

  ClassDef(RooSharedProperties,1) // Abstract interface for shared property implementations
};


#endif
