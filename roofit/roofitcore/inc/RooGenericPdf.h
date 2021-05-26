/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooGenericPdf.h,v 1.20 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_GENERIC_PDF
#define ROO_GENERIC_PDF

#include "RooAbsPdf.h"
#include "RooFormula.h"
#include "RooListProxy.h"

class RooArgList ;

class RooGenericPdf : public RooAbsPdf {
public:
  // Constructors, assignment etc
  inline RooGenericPdf(){}
  RooGenericPdf(const char *name, const char *title, const char* formula, const RooArgList& dependents);
  RooGenericPdf(const char *name, const char *title, const RooArgList& dependents);
  RooGenericPdf(const RooGenericPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooGenericPdf(*this,newname); }

  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

  // Printing interface (human readable)
  void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const ;
  void printMetaArgs(std::ostream& os) const ;

  // Debugging
  void dumpFormula() { formula().dump() ; }

protected:

  RooFormula& formula() const ;

  // Function evaluation
  RooListProxy _actualVars ;
  virtual Double_t evaluate() const ;
  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext& inputData, const RooArgSet* normSet) const;

  Bool_t setFormula(const char* formula) ;

  // Post-processing of server redirection
  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

  virtual Bool_t isValidReal(Double_t value, Bool_t printError) const ;

  std::unique_ptr<RooFormula> _formula{nullptr}; //! Formula engine
  TString _formExpr ;            // Formula expression string

  ClassDef(RooGenericPdf,1) // Generic PDF defined by string expression and list of variables
};

#endif
