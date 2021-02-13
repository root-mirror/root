/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooLinearVar.h,v 1.19 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_LINEAR_VAR
#define ROO_LINEAR_VAR

#include <cmath>
#include <cfloat>
#include <string>
#include <list>
#include "RooAbsRealLValue.h"
#include "RooRealProxy.h"
#include "RooFormula.h"
#include "RooLinTransBinning.h"

class RooArgSet ;

class RooLinearVar : public RooAbsRealLValue {
public:
  // Constructors, assignment etc.
  RooLinearVar() {} ;
  RooLinearVar(const char *name, const char *title, RooAbsRealLValue& variable, const RooAbsReal& slope, const RooAbsReal& offset, const char *unit= "") ;
  RooLinearVar(const RooLinearVar& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooLinearVar(*this,newname); }
  virtual ~RooLinearVar() ;
  
  // Parameter value and error accessors
  virtual void setVal(Double_t value) ;

  // Jacobian and limits
  virtual Bool_t hasBinning(const char* name) const ;
  virtual const RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE) const ;
  virtual RooAbsBinning& getBinning(const char* name=0, Bool_t verbose=kTRUE, Bool_t createOnTheFly=kFALSE)  ;
  virtual std::list<std::string> getBinningNames() const;

  virtual Double_t jacobian() const ;
  virtual Bool_t isJacobianOK(const RooArgSet& depList) const ;

  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // Printing interface (human readable)

  using RooAbsRealLValue::operator= ;
  using RooAbsRealLValue::setVal ;

protected:

  virtual Double_t evaluate() const ;

  mutable RooLinTransBinning _binning ;
  RooLinkedList _altBinning ; //!
  RooRealProxy _var ;         // Input observable
  RooRealProxy _slope ;       // Slope of transformation
  RooRealProxy _offset ;      // Offset of transformation

  ClassDef(RooLinearVar,1) // Lvalue linear transformation function 
};

#endif
