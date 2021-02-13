/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsFunc.h,v 1.9 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_FUNC
#define ROO_ABS_FUNC

#include "Rtypes.h"
#include "RooSpan.h"

#include <list>
#include <vector>

class RooAbsRealLValue ;

class RooAbsFunc {
public:
  inline RooAbsFunc(UInt_t dimension) : _ncall(0), _dimension(dimension), _valid(kTRUE) { }
  inline RooAbsFunc(const RooAbsFunc& other) : _ncall(0), _dimension(other._dimension), _valid(kTRUE) { }

  inline virtual ~RooAbsFunc() { }
  inline UInt_t getDimension() const { 
    // Dimension of function
    return _dimension; 
  }
  inline Bool_t isValid() const { 
    // Is function in valid state
    return _valid; 
  }

  virtual Double_t operator()(const Double_t xvector[]) const = 0;
  virtual RooSpan<const double> getValues(std::vector<RooSpan<const double>> /*coordinates*/) const {
    throw std::logic_error("Not implemented.");
  }
  virtual Double_t getMinLimit(UInt_t dimension) const = 0;
  virtual Double_t getMaxLimit(UInt_t dimension) const = 0;

  Int_t numCall() const { 
    // Return number of function calls since last reset
    return _ncall ; 
  }
  void resetNumCall() const { 
    // Reset function call counter
    _ncall = 0 ; 
  }

  virtual void saveXVec() const {
    // Interface to save current values of observables (if supported by binding implementation)
  } ;
  virtual void restoreXVec() const {
    // Interface to restore observables to saved values (if supported
    // by binding implementation)
  } ;

  virtual const char* getName() const { 
    // Name of function binding
    return "(unnamed)" ; 
  }  

  virtual std::list<Double_t>* binBoundaries(Int_t) const { return 0 ; }

  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const {
    // Interface for returning an optional hint for initial sampling points when constructing a curve 
    // projected on observable.  
    return 0 ; 
  }

protected:
  mutable Int_t _ncall ; // Function call counter
  UInt_t _dimension;     // Number of observables
  Bool_t _valid;         // Is binding in valid state?
   ClassDef(RooAbsFunc,0) // Abstract real-valued function interface
};

#endif

