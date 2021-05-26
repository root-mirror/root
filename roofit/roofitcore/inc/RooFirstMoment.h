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
#ifndef ROO_FIRST_MOMENT
#define ROO_FIRST_MOMENT

#include "RooAbsMoment.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"


class RooRealVar;
class RooArgList ;

class RooFirstMoment : public RooAbsMoment {
public:

  RooFirstMoment() ;
  RooFirstMoment(const char *name, const char *title, RooAbsReal& func, RooRealVar& x) ;
  RooFirstMoment(const char *name, const char *title, RooAbsReal& func, RooRealVar& x, const RooArgSet& nset, Bool_t intNSet=kFALSE) ;
  virtual ~RooFirstMoment() ;

  RooFirstMoment(const RooFirstMoment& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooFirstMoment(*this, newname); }

  const RooAbsReal& xF() { return _xf.arg() ; }
  const RooAbsReal& ixF() { return _ixf.arg() ; }
  const RooAbsReal& iF() { return _if.arg() ; }

protected:

  RooRealProxy _xf ;                     // X*F
  RooRealProxy _ixf ;                    // Int(X*F(X))dx ;
  RooRealProxy _if ;                     // Int(F(x))dx ;
  Double_t evaluate() const;

  ClassDef(RooFirstMoment,1) // Representation of moment in a RooAbsReal in a given RooRealVar
};

#endif
