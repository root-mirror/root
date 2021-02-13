/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooPullVar.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_PULL_VAR
#define ROO_PULL_VAR

#include "RooTemplateProxy.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"


class RooPullVar : public RooAbsReal {
public:

  RooPullVar() ;
  RooPullVar(const char *name, const char *title, RooRealVar& measurement, RooAbsReal& truth) ;
  virtual ~RooPullVar() ;

  RooPullVar(const RooPullVar& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooPullVar(*this, newname); }


protected:

  RooTemplateProxy<RooRealVar> _meas ;
  RooRealProxy _true ;

  Double_t evaluate() const;

  ClassDef(RooPullVar,1) // Calculation of pull of measurement w.r.t a truth value
};

#endif
