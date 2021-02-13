/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooFormula.h,v 1.34 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_FORMULA
#define ROO_FORMULA

#include "RooPrintable.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "TFormula.h"
#include "RooSpan.h"

#include <memory>
#include <vector>
#include <string>

namespace RooBatchCompute {
  struct RunContext;
}
class RooAbsReal;

class RooFormula : public TNamed, public RooPrintable {
public:
  // Constructors etc.
  RooFormula() ;
  RooFormula(const char* name, const char* formula, const RooArgList& varList, bool checkVariables = true);
  RooFormula(const RooFormula& other, const char* name=0);
  virtual TObject* Clone(const char* newName = nullptr) const {return new RooFormula(*this, newName);}
	
  ////////////////////////////////////////////////////////////////////////////////
  /// Return list of arguments which are used in the formula.
  RooArgSet actualDependents() const {return usedVariables();}
  Bool_t changeDependents(const RooAbsCollection& newDeps, Bool_t mustReplaceAll, Bool_t nameChange) ;

  /// Return pointer to the parameter with given name.
  /// \return Parameter if in use, nullptr if not in use.
  RooAbsArg* getParameter(const char* name) const {
    return usedVariables().find(name);
  }

  /// Return pointer to parameter at given index. This returns
  /// irrespective of whether the parameter is in use.
  RooAbsArg* getParameter(Int_t index) const {
    return _origList.at(index);
  }

  Bool_t ok() const { return _tFormula != nullptr; }
  /// Evalute all parameters/observables, and then evaluate formula.
  Double_t eval(const RooArgSet* nset=0) const;
  RooSpan<double> evaluateSpan(const RooAbsReal* dataOwner, RooBatchCompute::RunContext& inputData, const RooArgSet* nset = nullptr) const;

  /// DEBUG: Dump state information
  void dump() const;
  Bool_t reCompile(const char* newFormula) ;


  virtual void printValue(std::ostream& os) const ;
  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printArgs(std::ostream& os) const ;
  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;

  virtual void Print(Option_t *options= 0) const {
    // Printing interface (human readable)
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  std::string formulaString() const {
    return _tFormula ? _tFormula->GetTitle() : "";
  }

private:
  RooFormula& operator=(const RooFormula& other);
  std::string processFormula(std::string origFormula) const;
  RooArgList  usedVariables() const;
  std::string reconstructFormula(std::string internalRepr) const;
  void installFormulaOrThrow(const std::string& formulaa);

  RooArgList _origList; //! Original list of dependents
  std::vector<bool> _isCategory; //! Whether an element of the _origList is a category.
  std::unique_ptr<TFormula> _tFormula; //! The formula used to compute values

  ClassDef(RooFormula,0)
};

#endif
