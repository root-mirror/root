/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   TS, Thomas Schietinger, SLAC,           schieti@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooVoigtian
    \ingroup Roofit

RooVoigtian is an efficient implementation of the convolution of a
Breit-Wigner with a Gaussian, making use of the complex error function.
RooFitCore provides two algorithms for the evaluation of the complex error
function (the default CERNlib C335 algorithm, and a faster, look-up-table
based method). By default, RooVoigtian employs the default (CERNlib)
algorithm. Select the faster algorithm either in the constructor, or with
the selectFastAlgorithm() method.
**/

#include "RooVoigtian.h"
#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "RooBatchCompute.h"

#include <cmath>
#include <complex>
using namespace std;

ClassImp(RooVoigtian);

////////////////////////////////////////////////////////////////////////////////

RooVoigtian::RooVoigtian(const char *name, const char *title,
          RooAbsReal& _x, RooAbsReal& _mean,
          RooAbsReal& _width, RooAbsReal& _sigma,
              Bool_t doFast) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  width("width","Breit-Wigner Width",this,_width),
  sigma("sigma","Gauss Width",this,_sigma),
  _doFast(doFast)
{

}

////////////////////////////////////////////////////////////////////////////////

RooVoigtian::RooVoigtian(const RooVoigtian& other, const char* name) :
  RooAbsPdf(other,name), x("x",this,other.x), mean("mean",this,other.mean),
  width("width",this,other.width),sigma("sigma",this,other.sigma),
  _doFast(other._doFast)
{

}

////////////////////////////////////////////////////////////////////////////////

Double_t RooVoigtian::evaluate() const
{
  Double_t s = (sigma>0) ? sigma : -sigma ;
  Double_t w = (width>0) ? width : -width ;

  Double_t coef= -0.5/(s*s);
  Double_t arg = x - mean;

  // return constant for zero width and sigma
  if (s==0. && w==0.) return 1.;

  // Breit-Wigner for zero sigma
  if (s==0.) return (1./(arg*arg+0.25*w*w));

  // Gauss for zero width
  if (w==0.) return exp(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  Double_t c = 1./(sqrt(2.)*s);
  Double_t a = 0.5*c*w;
  Double_t u = c*arg;
  std::complex<Double_t> z(u,a) ;
  std::complex<Double_t> v(0.) ;

  if (_doFast) {
    v = RooMath::faddeeva_fast(z);
  } else {
    v = RooMath::faddeeva(z);
  }
  return c * v.real();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of Voigtian distribution. 
void RooVoigtian::computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const
{
  rbc::dispatch->compute(rbc::Voigtian, output, nEvents, dataMap, {&*x,&*mean,&*width,&*sigma,&*_norm});
}
