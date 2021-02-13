/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsHiddenReal.h,v 1.10 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_HIDDEN_REAL
#define ROO_ABS_HIDDEN_REAL

class RooArgSet ;
#include "RooAbsReal.h"
#include "RooCategoryProxy.h"

class RooCategory ;

class RooAbsHiddenReal : public RooAbsReal {
public:
  // Constructors, assignment etc.
  inline RooAbsHiddenReal() { 
    // Default constructor
  }
  RooAbsHiddenReal(const char *name, const char *title, const char *unit= "") ;
  RooAbsHiddenReal(const char *name, const char *title, RooAbsCategory& blindState, const char *unit= "") ;
  RooAbsHiddenReal(const RooAbsHiddenReal& other, const char* name=0) ;
  virtual ~RooAbsHiddenReal();
  
  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // Printing interface (human readable)
  virtual void printValue(std::ostream& stream) const ;
  
  inline Bool_t isHidden() const { 
    // If true, hiding mode is active
    return _state.arg().getCurrentIndex()!=0 ; 
  }

  Double_t getHiddenVal(const RooArgSet* nset=0) const { 
    // Bypass accessor to function value that also works in hidden mode
    return RooAbsReal::getVal(nset) ; 
  }

protected:

  // This is dubious from a C++ point of view, but it blocks the interactive user
  // from accidentally calling getVal() without explicit cast, which is the whole
  // point of this class
  virtual Double_t getValV(const RooArgSet* nset=0) const { 
    // Forward call to RooAbsReal
    return RooAbsReal::getValV(nset) ; 
  }

  static RooCategory* _dummyBlindState ; 
  RooAbsCategory& dummyBlindState() const ;

  RooCategoryProxy _state ; // Proxy to hiding state category

  ClassDef(RooAbsHiddenReal,1) // Abstract hidden real-valued variable
};

#endif
