/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
\file RooRealAnalytic.cxx
\class RooRealAnalytic
\ingroup Roofitcore

Lightweight RooAbsFunc interface adaptor that binds an analytic integral of a
RooAbsReal object (specified by a code) to a set of dependent variables.
**/


#include "RooFit.h"

#include "RooRealAnalytic.h"
#include "RooAbsReal.h"
#include "RooAbsRealLValue.h"

#include <assert.h>

using namespace std;

ClassImp(RooRealAnalytic);
;


////////////////////////////////////////////////////////////////////////////////
/// Evaluate our analytic integral at the specified values of the dependents.

Double_t RooRealAnalytic::operator()(const Double_t xvector[]) const 
{
  assert(isValid());
  loadValues(xvector);  
  _ncall++ ;
  return _code ? _func->analyticalIntegralWN(_code,_nset,_rangeName?_rangeName->GetName():0):_func->getVal(_nset) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate the analytic integral of the function at the specified values of the dependents.
RooSpan<const double> RooRealAnalytic::getValues(std::vector<RooSpan<const double>> coordinates) const {
  assert(isValid());
  _ncall += coordinates.front().size();

  if (!_batchBuffer)
    _batchBuffer.reset(new std::vector<double>());
  _batchBuffer->resize(coordinates.front().size());
  RooSpan<double> results(*_batchBuffer);

  for (std::size_t i=0; i < coordinates.front().size(); ++i) {
    for (unsigned int dim=0; dim < coordinates.size(); ++dim) {
      _vars[dim]->setVal(coordinates[dim][i]);
    }

    if (_code == 0) {
      results[i] = _func->getVal(_nset);
    } else {
      results[i] = _func->analyticalIntegralWN(_code,_nset,_rangeName?_rangeName->GetName():0);
    }
  }

  return results;
}

