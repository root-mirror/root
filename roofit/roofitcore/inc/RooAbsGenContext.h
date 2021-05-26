/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsGenContext.h,v 1.15 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ABS_GEN_CONTEXT
#define ROO_ABS_GEN_CONTEXT

#include "TNamed.h"
#include "RooPrintable.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"

class RooDataSet;

class RooAbsGenContext : public TNamed, public RooPrintable {
public:
  RooAbsGenContext(const RooAbsPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0, const RooArgSet* auxProto=0,
         Bool_t _verbose= kFALSE) ;
  virtual ~RooAbsGenContext();

  virtual RooDataSet *generate(Double_t nEvents= 0, Bool_t skipInit=kFALSE, Bool_t extendedMode=kFALSE);

  Bool_t isValid() const {
    // If true generator context is in a valid state
    return _isValid;
  }

  inline void setVerbose(Bool_t verbose= kTRUE) {
    // Set/clear verbose messaging
    _verbose= verbose;
  }
  inline Bool_t isVerbose() const {
    // If true verbose messaging is active
    return _verbose;
  }

  virtual void setProtoDataOrder(Int_t* lut) ;

   inline virtual void Print(Option_t *options= 0) const {
     // Print context information on stdout
     printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  virtual void attach(const RooArgSet& params) ;

  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printArgs(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const;

  virtual Int_t defaultPrintContents(Option_t* opt) const ;
  virtual StyleOption defaultPrintStyle(Option_t* opt) const ;

  virtual void setExpectedData(Bool_t) {} ;

  virtual void generateEvent(RooArgSet &theEvent, Int_t remaining) = 0;
  virtual void initGenerator(const RooArgSet &theEvent);

protected:

  virtual RooDataSet* createDataSet(const char* name, const char* title, const RooArgSet& obs) ;

  void resampleData(Double_t& ratio) ;

  const RooDataSet *_prototype; // Pointer to prototype dataset
  RooArgSet *_theEvent;         // Pointer to observable event being generated
  Bool_t _isValid;              // Is context in valid state?
  Bool_t _verbose;              // Verbose messaging?
  UInt_t _expectedEvents;       // Number of expected events from extended p.d.f
  RooArgSet _protoVars;         // Prototype observables
  Int_t _nextProtoIndex;        // Next prototype event to load according to LUT
  RooAbsPdf::ExtendMode _extendMode ;  // Extended mode capabilities of p.d.f.
  Int_t* _protoOrder ;          // LUT with traversal order of prototype data
  TString _normRange ;          // Normalization range of pdf

  RooDataSet* _genData ;        //! Data being generated

  ClassDef(RooAbsGenContext,0) // Abstract context for generating a dataset from a PDF
};

#endif
