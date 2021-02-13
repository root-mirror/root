/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Author:                                                                   *
 *   Tristan du Pree, Nikhef, Amsterdam, tdupree@nikhef.nl                   *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MULTIBINOMIAL
#define ROO_MULTIBINOMIAL

#include "RooAbsReal.h"
#include "RooListProxy.h"

class RooArgList;


class RooMultiBinomial : public RooAbsReal {
 public:
  // Constructors, assignment etc
  inline RooMultiBinomial() {
  }

  RooMultiBinomial(const char *name, const char *title, const RooArgList& effFuncList, const RooArgList& catList, Bool_t ignoreNonVisible);
  RooMultiBinomial(const RooMultiBinomial& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooMultiBinomial(*this,newname); }
  virtual ~RooMultiBinomial();

 protected:

  // Function evaluation
  virtual Double_t evaluate() const ;

 private:

  RooListProxy _catList ; // Accept/reject categories
  RooListProxy _effFuncList ; // Efficiency functions per category
  Bool_t _ignoreNonVisible ; // Ignore combination of only rejects (since invisible)

  ClassDef(RooMultiBinomial,1) // Simultaneous pdf of N Binomial distributions with associated efficiency functions
  };

#endif
