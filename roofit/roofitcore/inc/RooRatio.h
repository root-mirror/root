// Author: Rahul Balasubramanian, Nikhef 01 Apr 2021

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_RATIO
#define ROO_RATIO

#include "RooAbsReal.h"
#include "RooProduct.h"
#include "RooRealProxy.h"

#include <list>

class RooRealVar;
class RooArgList;
class RooProduct;

class RooRatio : public RooAbsReal {
public:
  RooRatio();
  RooRatio(const char *name, const char *title, Double_t numerator,
           Double_t denominator);
  RooRatio(const char *name, const char *title, Double_t numerator,
           RooAbsReal &denominator);
  RooRatio(const char *name, const char *title, RooAbsReal &numerator,
           Double_t denominator);
  RooRatio(const char *name, const char *title, RooAbsReal &numerator,
           RooAbsReal &denominator);
  RooRatio(const char *name, const char *title,
           const RooArgList &num, const RooArgList &denom);

  RooRatio(const RooRatio &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const {
    return new RooRatio(*this, newname);
  }
  virtual ~RooRatio();

protected:
  RooRealProxy _numerator;
  RooRealProxy _denominator;

  Double_t evaluate() const;

  ClassDef(RooRatio, 2) // Ratio of two RooAbsReal and/or numbers
};

#endif
