/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooPolyVar.h,v 1.7 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_POLY_VAR
#define ROO_POLY_VAR

#include <vector>

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooPolyVar : public RooAbsReal {
public:

  RooPolyVar() ;
  RooPolyVar(const char* name, const char* title, RooAbsReal& x) ;
  RooPolyVar(const char *name, const char *title,
      RooAbsReal& _x, const RooArgList& _coefList, Int_t lowestOrder=0) ;

  RooPolyVar(const RooPolyVar& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooPolyVar(*this, newname); }
  virtual ~RooPolyVar() ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy _x;
  RooListProxy _coefList ;
  Int_t _lowestOrder ;

  mutable std::vector<Double_t> _wksp; //! do not persist

  Double_t evaluate() const;

  ClassDef(RooPolyVar,1) // Polynomial function
};

#endif
