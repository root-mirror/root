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
\file RooPullVar.cxx
\class RooPullVar
\ingroup Roofitcore

RooPullVar represents the pull of a measurement w.r.t. the true value
using the measurement and its error. Both the true value and
the measured value (with error) are taken from two user-supplied
RooRealVars. If the measured parameter has an asymmetric error, the proper
side of that error will be used:
\f[
 \mathrm{Pull}_x = \frac{x_\mathrm{meas}-x_\mathrm{true}}{\Delta_x}
\f]
**/

#include "RooPullVar.h"

#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

using namespace std;

ClassImp(RooPullVar);
;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooPullVar::RooPullVar()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Construct the pull of the RooRealVar 'meas'.
///
/// \param[in] name  Name of the pull variable.
/// \param[in] title The title (for plotting).
/// \param[in] meas  The measurement. This variable needs to have an error, and therefore is a RooRealVar.
/// \param[in] truth The true value.
RooPullVar::RooPullVar(const char* name, const char* title, RooRealVar& meas, RooAbsReal& truth) :
  RooAbsReal(name, title),
  _meas("meas","Measurement",this,meas),
  _true("true","Truth",this,truth)
{
}





////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooPullVar::RooPullVar(const RooPullVar& other, const char* name) :
  RooAbsReal(other, name), 
  _meas("meas",this,other._meas),
  _true("true",this,other._true)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooPullVar::~RooPullVar() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate pull. Use asymmetric error if defined in measurement,
/// otherwise use symmetric error. If measurement has no error
/// return zero.

Double_t RooPullVar::evaluate() const 
{
  const auto& meas = _meas.arg();
  if (meas.hasAsymError()) {
    Double_t delta = _meas-_true ;
    if (delta<0) {
      return delta/meas.getAsymErrorHi() ;
    } else {
      return -delta/meas.getAsymErrorLo() ;
    }
  } else if (meas.hasError()) {
    return (_meas-_true)/meas.getError() ;    
  } else {
    return 0 ;
  }
}


