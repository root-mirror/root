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

/** \class RooGaussian
    \ingroup Roofit

Plain Gaussian p.d.f
**/

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooGaussian.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"

using namespace std;

ClassImp(RooGaussian);

////////////////////////////////////////////////////////////////////////////////

RooGaussian::RooGaussian(const char *name, const char *title,
          RooAbsReal& _x, RooAbsReal& _mean,
          RooAbsReal& _sigma) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma)
{
}

////////////////////////////////////////////////////////////////////////////////

RooGaussian::RooGaussian(const RooGaussian& other, const char* name) :
  RooAbsPdf(other,name), x("x",this,other.x), mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooGaussian::evaluate() const
{
  const double arg = x - mean;
  const double sig = sigma;
  return exp(-0.5*arg*arg/(sig*sig)) ;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate and return the negative log-likelihood of the Poisson

Double_t RooGaussian::getLogVal(const RooArgSet* set) const
{
  return RooAbsPdf::getLogVal(set) ;
//   Double_t prob = getVal(set) ;
//   return log(prob) ;

  Double_t arg= x - mean;
  Double_t sig = sigma ;

  //static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  //Double_t extra = -0.5*arg*arg/(sig*sig) - log(2*rootPiBy2*sig) ;
  Double_t extra = -0.5*arg*arg/(sig*sig) - log(analyticalIntegral(1,0)) ;

  return extra ;

}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGaussian::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  if (matchArgs(allVars,analVars,mean)) return 2 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooGaussian::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1 || code==2);

  //The normalisation constant 1./sqrt(2*pi*sigma^2) is left out in evaluate().
  //Therefore, the integral is scaled up by that amount to make RooFit normalise
  //correctly.
  const double resultScale = sqrt(TMath::TwoPi()) * sigma;

  //Here everything is scaled and shifted into a standard normal distribution:
  const double xscale = TMath::Sqrt2() * sigma;
  double max = 0.;
  double min = 0.;
  if (code == 1){
    max = (x.max(rangeName)-mean)/xscale;
    min = (x.min(rangeName)-mean)/xscale;
  } else { //No == 2 test because of assert
    max = (mean.max(rangeName)-x)/xscale;
    min = (mean.min(rangeName)-x)/xscale;
  }


  //Here we go for maximum precision: We compute all integrals in the UPPER
  //tail of the Gaussian, because erfc has the highest precision there.
  //Therefore, the different cases for range limits in the negative hemisphere are mapped onto
  //the equivalent points in the upper hemisphere using erfc(-x) = 2. - erfc(x)
  const double ecmin = std::erfc(std::abs(min));
  const double ecmax = std::erfc(std::abs(max));


  return resultScale * 0.5 * (
      min*max < 0.0 ? 2.0 - (ecmin + ecmax)
                    : max <= 0. ? ecmax - ecmin : ecmin - ecmax
  );
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGaussian::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  if (matchArgs(directVars,generateVars,mean)) return 2 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooGaussian::generateEvent(Int_t code)
{
  assert(code==1 || code==2) ;
  Double_t xgen ;
  if(code==1){
    while(1) {
      xgen = RooRandom::randomGenerator()->Gaus(mean,sigma);
      if (xgen<x.max() && xgen>x.min()) {
   x = xgen ;
   break;
      }
    }
  } else if(code==2){
    while(1) {
      xgen = RooRandom::randomGenerator()->Gaus(x,sigma);
      if (xgen<mean.max() && xgen>mean.min()) {
   mean = xgen ;
   break;
      }
    }
  } else {
    cout << "error in RooGaussian generateEvent"<< endl;
  }

  return;
}
