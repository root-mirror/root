/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooErrorVar.h,v 1.16 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ERROR_VAR
#define ROO_ERROR_VAR

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

#include <list>
#include <string>

class RooVectorDataStore;

class RooErrorVar : public RooAbsRealLValue {
public:
  // Constructors, assignment etc.
  inline RooErrorVar() {
    // Default constructor
  }
  RooErrorVar(const char *name, const char *title, const RooRealVar& input) ;
  RooErrorVar(const RooErrorVar& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooErrorVar(*this,newname); }
  virtual ~RooErrorVar() ;

  virtual Double_t getValV(const RooArgSet* set=0) const ;

  virtual Double_t evaluate() const {
    // return error of input RooRealVar
    return ((RooRealVar&)_realVar.arg()).getError() ;
  }

  virtual void setVal(Double_t value) {
    // Set error of input RooRealVar to value
    ((RooRealVar&)_realVar.arg()).setVal(value) ;
  }

  inline virtual Bool_t isFundamental() const {
    // Return kTRUE as we implement a fundamental type of AbsArg that can be stored in a dataset
    return kTRUE ;
  }

  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // Set/get finite fit range limits
  inline void setMin(Double_t value) {
    // Set lower bound of default range to value
    setMin(0,value) ;
  }
  inline void setMax(Double_t value) {
    // Set upper bound of default range to value
    setMax(0,value) ;
  }
  inline void setRange(Double_t min, Double_t max) {
    // Set default ranges to [min,max]
    setRange(0,min,max) ;
  }
  void setMin(const char* name, Double_t value) ;
  void setMax(const char* name, Double_t value) ;
  void setRange(const char* name, Double_t min, Double_t max) ;

  void setBins(Int_t nBins);
  void setBinning(const RooAbsBinning& binning, const char* name=0) ;
  const RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE) const ;
  RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE) ;
  Bool_t hasBinning(const char* name) const ;
  std::list<std::string> getBinningNames() const ;

  // Set infinite fit range limits
  void removeMin(const char* name=0);
  void removeMax(const char* name=0);
  void removeRange(const char* name=0);

  using RooAbsRealLValue::operator= ;
  using RooAbsRealLValue::setVal ;

protected:

  RooLinkedList _altBinning ;  //! Optional alternative ranges and binnings

  void syncCache(const RooArgSet* set=0) ;

  RooRealProxy _realVar ; // RealVar with the original error
  RooAbsBinning* _binning ; //! Pointer to default binning definition

  ClassDef(RooErrorVar,1) // RooAbsRealLValue representation of an error of a RooRealVar
};

#endif
