/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
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

/** \class RooLandau
    \ingroup Roofit

Landau distribution p.d.f
\image html RF_Landau.png "PDF of the Landau distribution."
**/

#include "RooLandau.h"
#include "RooHelpers.h"
#include "RooFit.h"
#include "RooRandom.h"
#include "RooBatchCompute.h"

#include "TMath.h"

ClassImp(RooLandau);

////////////////////////////////////////////////////////////////////////////////

RooLandau::RooLandau(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma)
{
  RooHelpers::checkRangeOfParameters(this, {&_sigma}, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

RooLandau::RooLandau(const RooLandau& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooLandau::evaluate() const
{
  return TMath::Landau(x, mean, sigma);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of Landau distribution.  
RooSpan<double> RooLandau::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {
  return RooBatchCompute::dispatch->computeLandau(this, evalData, x->getValues(evalData, normSet), mean->getValues(evalData, normSet), sigma->getValues(evalData, normSet));
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooLandau::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooLandau::generateEvent(Int_t code)
{
  assert(1 == code); (void)code;
  Double_t xgen ;
  while(1) {
    xgen = RooRandom::randomGenerator()->Landau(mean,sigma);
    if (xgen<x.max() && xgen>x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}
