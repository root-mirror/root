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
#ifndef ROO_RANGE_BOOLEAN
#define ROO_RANGE_BOOLEAN

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "TString.h"

class RooRealVar;
class RooArgList ;
#include <list>

class RooRangeBoolean : public RooAbsReal {
public:

  RooRangeBoolean() ;
  RooRangeBoolean(const char* name, const char* title, RooAbsRealLValue& x, const char* rangeName) ;
  RooRangeBoolean(const RooRangeBoolean& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooRangeBoolean(*this, newname); }
  virtual ~RooRangeBoolean() ;


  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ;

protected:

  RooRealProxy _x;
  TString _rangeName ;

  Double_t evaluate() const;

  ClassDef(RooRangeBoolean,1) // Polynomial function
};

#endif
