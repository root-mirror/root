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

/** \class RooCBShape
    \ingroup Roofit

PDF implementing the Crystal Ball line shape.
**/

#include "RooCBShape.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "BatchHelpers.h"
#include "RooVDTHeaders.h"

#include "TMath.h"

#include <cmath>

using namespace std;

ClassImp(RooCBShape);

////////////////////////////////////////////////////////////////////////////////

Double_t RooCBShape::ApproxErf(Double_t arg) const
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;

  return RooMath::erf(arg);
}

////////////////////////////////////////////////////////////////////////////////

RooCBShape::RooCBShape(const char *name, const char *title,
             RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
             RooAbsReal& _alpha, RooAbsReal& _n) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alpha("alpha", "Alpha", this, _alpha),
  n("n", "Order", this, _n)
{
}

////////////////////////////////////////////////////////////////////////////////

RooCBShape::RooCBShape(const RooCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alpha("alpha", this, other.alpha),
  n("n", this, other.n)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooCBShape::evaluate() const {
  Double_t t = (m-m0)/sigma;
  if (alpha < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)alpha);

  if (t >= -absAlpha) {
    return exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t b= n/absAlpha - absAlpha;

    return a/TMath::Power(b - t, n);
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace {
//Author: Emmanouil Michalainas, CERN 21 August 2019

template<class Tm, class Tm0, class Tsigma, class Talpha, class Tn>
void compute(	size_t batchSize,
	double * __restrict output,
	Tm M, Tm0 M0, Tsigma S, Talpha A, Tn N)
{
  for (size_t i=0; i<batchSize; i++) {
    const double t = (M[i]-M0[i]) / S[i];
    if ((A[i]>0 && t>=-A[i])   ||   (A[i]<0 && -t>=A[i])) {
      output[i] = -0.5*t*t;
    } else {
      output[i] = N[i] / (N[i] -A[i]*A[i] -A[i]*t);
      output[i] = _rf_fast_log(output[i]);
      output[i] *= N[i];
      output[i] -= 0.5*A[i]*A[i];
    }
  }
  
  for (size_t i=0; i<batchSize; i++) {
    output[i] = _rf_fast_exp(output[i]);
  }
}
};

RooSpan<double> RooCBShape::evaluateBatch(std::size_t begin, std::size_t batchSize) const {
  using namespace BatchHelpers;

  EvaluateInfo info = getInfo( {&m, &m0, &sigma, &alpha, &n}, begin, batchSize );
  if (info.nBatches == 0) {
    return {};
  }
  auto output = _batchData.makeWritableBatchUnInit(begin, batchSize);
  auto mData = m.getValBatch(begin, info.size);

  if (info.nBatches==1 && !mData.empty()) {
    compute(info.size, output.data(), mData.data(),
    BracketAdapter<double> (m0),
    BracketAdapter<double> (sigma),
    BracketAdapter<double> (alpha),
    BracketAdapter<double> (n));
  }
  else {
    compute(info.size, output.data(),
    BracketAdapterWithMask (m,m.getValBatch(begin,info.size)),
    BracketAdapterWithMask (m0,m0.getValBatch(begin,info.size)),
    BracketAdapterWithMask (sigma,sigma.getValBatch(begin,info.size)),
    BracketAdapterWithMask (alpha,alpha.getValBatch(begin,info.size)),
    BracketAdapterWithMask (n,n.getValBatch(begin,info.size)));
  }
  return output;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static const double sqrtPiOver2 = 1.2533141373;
  static const double sqrt2 = 1.4142135624;

  R__ASSERT(code==1);
  double result = 0.0;
  bool useLog = false;

  if( fabs(n-1.0) < 1.0e-05 )
    useLog = true;

  double sig = fabs((Double_t)sigma);

  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;

  if(alpha < 0) {
    double tmp = tmin;
    tmin = -tmax;
    tmax = -tmp;
  }

  double absAlpha = fabs((Double_t)alpha);

  if( tmin >= -absAlpha ) {
    result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                - ApproxErf(tmin/sqrt2) );
  }
  else if( tmax <= -absAlpha ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;

    if(useLog) {
      result += a*sig*( log(b-tmin) - log(b-tmax) );
    }
    else {
      result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                                - 1.0/(TMath::Power(b-tmax,n-1.0)) );
    }
  }
  else {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;

    double term1 = 0.0;
    if(useLog) {
      term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
    }
    else {
      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                              - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }

    double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                     - ApproxErf(-absAlpha/sqrt2) );


    result += term1 + term2;
  }

  return result != 0. ? result : 1.E-300;
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise that we know the maximum of self for given (m0,alpha,n,sigma)

Int_t RooCBShape::getMaxVal(const RooArgSet& vars) const
{
  RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
    return 1 ;
  }
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooCBShape::maxVal(Int_t code) const
{
  R__ASSERT(code==1) ;

  // The maximum value for given (m0,alpha,n,sigma)
  return 1.0 ;
}
