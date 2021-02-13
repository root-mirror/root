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
#ifndef ROO_RECURSIVE_FRACTION
#define ROO_RECURSIVE_FRACTION

#include "RooAbsReal.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooRecursiveFraction : public RooAbsReal {
public:

  RooRecursiveFraction() ;
  RooRecursiveFraction(const char *name, const char *title, const RooArgList& fracSet) ;
  virtual ~RooRecursiveFraction() ;

  RooRecursiveFraction(const RooRecursiveFraction& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooRecursiveFraction(*this, newname); }

protected:

  RooListProxy _list ;

  Double_t evaluate() const;

  ClassDef(RooRecursiveFraction,1) // Recursive fraction formula f1*(1-f2)*(1-f3) etc...
} ;

#endif
