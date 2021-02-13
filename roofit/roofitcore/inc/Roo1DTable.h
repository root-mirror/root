/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: Roo1DTable.h,v 1.19 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_1D_TABLE
#define ROO_1D_TABLE

#include "TObjArray.h"
#include "RooTable.h"
#include <vector>

class Roo1DTable : public RooTable {
public:

  // Constructors, cloning and assignment
  Roo1DTable() {
    // Default constructor
    // coverity[UNINIT_CTOR]
  } ;
  virtual ~Roo1DTable();
  Roo1DTable(const char *name, const char *title, const RooAbsCategory &cat);
  Roo1DTable(const Roo1DTable& other) ;

  virtual void fill(RooAbsCategory& cat, Double_t weight=1.0) ;
  Double_t get(const char* label, Bool_t silent=kFALSE) const ;
  Double_t getFrac(const char* label, Bool_t silent=kFALSE) const ;
  Double_t get(const int index, Bool_t silent=kFALSE) const ;
  Double_t getFrac(const int index, Bool_t silent=kFALSE) const ;
  Double_t getOverflow() const ;

  // Printing interface (human readable)
  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printValue(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;
  virtual Int_t defaultPrintContents(Option_t* opt) const ;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface (human readable)
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  virtual Bool_t isIdentical(const RooTable& other) ;

protected:

  
  TObjArray _types ;             // Array of defined category states
  std::vector<Double_t> _count ; // Array of counters for each state
  Double_t  _total ;             // Total number of entries
  Double_t  _nOverflow ;         // Number of overflow entries

  ClassDef(Roo1DTable,1) // 1-dimensional table
};

#endif
