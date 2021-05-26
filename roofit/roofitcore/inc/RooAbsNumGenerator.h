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
#ifndef ROO_ABS_NUM_GENERATOR
#define ROO_ABS_NUM_GENERATOR

#include "TNamed.h"
#include "RooPrintable.h"
#include "RooArgSet.h"
#include "RooArgList.h"

class RooAbsReal;
class RooRealVar;
class RooDataSet;
class RooRealBinding;
class RooNumGenConfig ;

class RooAbsNumGenerator : public TNamed, public RooPrintable {
public:
  RooAbsNumGenerator() : _cloneSet(0), _funcClone(0), _funcMaxVal(0), _verbose(kFALSE), _isValid(kFALSE), _funcValStore(0), _funcValPtr(0), _cache(0) {} ;
  RooAbsNumGenerator(const RooAbsReal &func, const RooArgSet &genVars, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0);
  virtual RooAbsNumGenerator* clone(const RooAbsReal&, const RooArgSet& genVars, const RooArgSet& condVars,
                const RooNumGenConfig& config, Bool_t verbose=kFALSE, const RooAbsReal* maxFuncVal=0) const = 0 ;

  Bool_t isValid() const {
    // If true, generator is in a valid state
    return _isValid;
  }
  virtual ~RooAbsNumGenerator();

  inline void setVerbose(Bool_t verbose= kTRUE) {
    // If flag is true, verbose messaging will be active during generation
    _verbose= verbose;
  }
  inline Bool_t isVerbose() const {
    // Return status of verbose messaging flag
    return _verbose;
  }

  virtual const RooArgSet *generateEvent(UInt_t remaining, Double_t& resampleRatio) = 0;
  virtual Double_t getFuncMax() { return 0 ; }

   inline virtual void Print(Option_t *options= 0) const {
     // ascii printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printArgs(std::ostream& os) const ;

  void attachParameters(const RooArgSet& vars) ;

  // Advertisement of capabilities
  virtual Bool_t canSampleCategories() const { return kFALSE ; }
  virtual Bool_t canSampleConditional() const { return kFALSE ; } // Must implement getFuncMax()

protected:

  RooArgSet *_cloneSet;                // Set owning clone of input function
  RooAbsReal *_funcClone;              // Pointer to top level node of cloned function
  const RooAbsReal *_funcMaxVal ;      // Container for maximum function value
  RooArgSet _catVars,_realVars ;       // Sets of discrete and real valued observables
  Bool_t _verbose, _isValid;           // Verbose and valid flag
  RooRealVar *_funcValStore,*_funcValPtr; // RRVs storing function value in context and in output dataset

  RooDataSet *_cache;                  // Dataset holding generated values of observables

  ClassDef(RooAbsNumGenerator,0) // Abstract base class for numeric event generator algorithms
};

#endif
