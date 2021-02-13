/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsBinning.h,v 1.13 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_BINNING
#define ROO_ABS_BINNING

#include "Rtypes.h"
#include "RooPrintable.h"
#include "TNamed.h" 
class TIterator ;
class RooAbsRealLValue ;
class RooAbsReal ;

class RooAbsBinning : public TNamed, public RooPrintable {
public:

  RooAbsBinning(const char* name=0) ;
  RooAbsBinning(const RooAbsBinning& other, const char* name=0) : TNamed(name,name), RooPrintable(other) {
    // Copy constructor
  }
  virtual TObject* Clone(const char* newname=0) const { return clone(newname) ; }
  virtual RooAbsBinning* clone(const char* name=0) const = 0 ;
  virtual ~RooAbsBinning() ;

  /// Return number of bins.
  Int_t numBins() const { 
    return numBoundaries()-1 ; 
  }
  virtual Int_t numBoundaries() const = 0 ;
  virtual Int_t binNumber(Double_t x) const = 0 ;
  virtual Int_t rawBinNumber(Double_t x) const { return binNumber(x) ; }
  virtual Double_t binCenter(Int_t bin) const = 0 ;
  virtual Double_t binWidth(Int_t bin) const = 0 ;
  virtual Double_t binLow(Int_t bin) const = 0 ;
  virtual Double_t binHigh(Int_t bin) const = 0 ;
  virtual Bool_t isUniform() const { return kFALSE ; }

  virtual void setRange(Double_t xlo, Double_t xhi) = 0 ;
  /// Change lower bound to xlo.
  virtual void setMin(Double_t xlo) { 
    setRange(xlo,highBound()) ; 
  }
  /// Change upper bound to xhi.
  virtual void setMax(Double_t xhi) { 
    setRange(lowBound(),xhi) ; 
  }

  virtual Double_t lowBound() const = 0 ;
  virtual Double_t highBound() const = 0 ;
  virtual Double_t averageBinWidth() const = 0 ;


  virtual Double_t* array() const = 0 ;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printArgs(std::ostream& os) const ;
  virtual void printValue(std::ostream& os) const ;
  
  /// Interface function. If true, min/max of binning is parameterized by external RooAbsReals.
  /// Default to `false`, unless overridden by a sub class.
  virtual Bool_t isParameterized() const { 
    return kFALSE ; 
  }
  /// Return pointer to RooAbsReal parameterized lower bound, if any.
  virtual RooAbsReal* lowBoundFunc() const { 
    return 0 ; 
  }
  /// Return pointer to RooAbsReal parameterized upper bound, if any.
  virtual RooAbsReal* highBoundFunc() const { 
    return 0 ; 
  }
  /// If true (default), the range definition can be shared across clones of a RooRealVar.
  virtual Bool_t isShareable() const { 
    return kTRUE ; 
  }
  /// Hook interface function to execute code upon insertion into a RooAbsRealLValue.
  virtual void insertHook(RooAbsRealLValue&) const {  }
  /// Hook interface function to execute code upon removal from a RooAbsRealLValue.
  virtual void removeHook(RooAbsRealLValue&) const {  }


protected:  

  ClassDef(RooAbsBinning,2) // Abstract base class for binning specification
};

#endif
