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
#ifndef ROO_ABS_MOMENT
#define ROO_ABS_MOMENT

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"


class RooRealVar;
class RooArgList ;

class RooAbsMoment : public RooAbsReal {
public:

  RooAbsMoment() ;
  RooAbsMoment(const char *name, const char *title, RooAbsReal& func, RooRealVar& x, Int_t order=1, Bool_t takeRoot=kFALSE) ;
  RooAbsMoment(const RooAbsMoment& other, const char* name = 0);
  virtual ~RooAbsMoment() ;

  Int_t order() const { return _order ; }
  Bool_t central() const { return _mean.absArg() ? kTRUE : kFALSE ; }
  RooAbsReal* mean() { return (RooAbsReal*) _mean.absArg() ; }


protected:

  Int_t _order ;                         ///< Moment order
  Int_t _takeRoot ;                      ///< Return n-order root of moment
  RooSetProxy  _nset ;                   ///< Normalization set (optional)
  RooRealProxy _func ;                   ///< Input function
  RooRealProxy _x     ;                  ///< Observable
  RooRealProxy _mean ;                   ///< Mean (if calculated for central moment)

  ClassDef(RooAbsMoment,1) // Abstract representation of moment in a RooAbsReal in a given RooRealVar
};

#endif
