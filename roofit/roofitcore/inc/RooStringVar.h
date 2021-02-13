/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooStringVar.h,v 1.23 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_STRING_VAR
#define ROO_STRING_VAR

#include "RooAbsArg.h"

#include <string>

class RooStringVar final : public RooAbsArg {
public:
  // Constructors, assignment etc.
  RooStringVar() { }
  RooStringVar(const char *name, const char *title, const char* value, Int_t size=1024) ;
  RooStringVar(const RooStringVar& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override { return new RooStringVar(*this,newname); }
  virtual ~RooStringVar() = default;

  // Parameter value and error accessors
  virtual operator TString() {return TString(_string.c_str()); }
  const char* getVal() const { clearValueDirty(); return _string.c_str(); }
  void setVal(const char* newVal) { _string = newVal ? newVal : ""; setValueDirty(); }
  virtual RooAbsArg& operator=(const char* newVal) { setVal(newVal); return *this; }

  // We implement a fundamental type of AbsArg that can be stored in a dataset
  bool isFundamental() const override { return true; }

  // I/O streaming interface (machine readable)
  bool readFromStream(std::istream& is, Bool_t compact, Bool_t verbose) override;
  void writeToStream(std::ostream& os, Bool_t /*compact*/) const override { os << _string; }

  // Return value and unit accessors
  bool operator==(const char* val) const { return _string == val; }
  bool operator==(const RooAbsArg& other) const override {
    auto otherStr = dynamic_cast<const RooStringVar*>(&other);
    return otherStr && _string == otherStr->_string;
  }
  bool isIdentical(const RooAbsArg& other, Bool_t /*assumeSameType*/) const override { return *this == other; }

  // Printing interface (human readable)
  virtual void printValue(std::ostream& os) const override { os << _string; }


  RooAbsArg *createFundamental(const char* newname=0) const override {
    return new RooStringVar(newname ? newname : GetName(), GetTitle(), "", 1);
  }

protected:
  // Internal consistency checking (needed by RooDataSet)
  virtual Bool_t isValid() const override { return true; }
  virtual Bool_t isValidString(const char*, Bool_t /*printError=kFALSE*/) const { return true; }

  virtual void syncCache(const RooArgSet* /*nset*/ = nullptr) override { }
  void copyCache(const RooAbsArg* source, Bool_t valueOnly=kFALSE, Bool_t setValDiry=kTRUE) override;
  virtual void attachToTree(TTree& t, Int_t bufSize=32000) override;
  virtual void attachToVStore(RooVectorDataStore&) override { }
  virtual void fillTreeBranch(TTree& t) override;
  virtual void setTreeBranchStatus(TTree& t, Bool_t active) override;

private:
  std::string _string;
  ClassDefOverride(RooStringVar,2) // String-valued variable
};

#endif
