 /*****************************************************************************
  * Project: RooFit                                                           *
  *                                                                           *
  * Simple Poisson PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           *
  *****************************************************************************/

/** \class RooPoisson
    \ingroup Roofit

Poisson pdf
**/

#include "RooPoisson.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooNaNPacker.h"
#include "RooBatchCompute.h"

#include "Math/ProbFuncMathCore.h"

ClassImp(RooPoisson);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooPoisson::RooPoisson(const char *name, const char *title,
             RooAbsReal& _x,
             RooAbsReal& _mean,
             Bool_t noRounding) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  _noRounding(noRounding)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

 RooPoisson::RooPoisson(const RooPoisson& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   mean("mean",this,other.mean),
   _noRounding(other._noRounding),
   _protectNegative(other._protectNegative)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation in terms of the TMath::Poisson() function.

Double_t RooPoisson::evaluate() const
{
  Double_t k = _noRounding ? x : floor(x);
  if(_protectNegative && mean<0) {
    RooNaNPacker np;
    np.setPayload(-mean);
    return np._payload;
  }
  return TMath::Poisson(k,mean) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute multiple values of the Poisson distribution.  
void RooPoisson::computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const
{
  rbc::dispatch->compute(rbc::Poisson, output, nEvents, dataMap, {&*x,&*mean,&*_norm}, 
    {static_cast<double>(_protectNegative), static_cast<double>(_noRounding)});
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooPoisson::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  if (matchArgs(allVars, analVars, mean)) return 2;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooPoisson::analyticalIntegral(Int_t code, const char* rangeName) const
{
  R__ASSERT(code == 1 || code == 2) ;

  if(_protectNegative && mean<0)
    return exp(-2*mean); // make it fall quickly

  if (code == 1) {
    // Implement integral over x as summation. Add special handling in case
    // range boundaries are not on integer values of x
    const double xmin = std::max(0., x.min(rangeName));
    const double xmax = x.max(rangeName);

    if (xmax < 0. || xmax < xmin) {
      return 0.;
    }
    if (!x.hasMax() || RooNumber::isInfinite(xmax)) {
      //Integrating the full Poisson distribution here 
      return 1.;
    }

    // The range as integers. ixmin is included, ixmax outside.
    const unsigned int ixmin = xmin;
    const unsigned int ixmax = std::min(xmax + 1.,
                                        (double)std::numeric_limits<unsigned int>::max());
    
    // Sum from 0 to just before the bin outside of the range.
    if (ixmin == 0) {
      return ROOT::Math::poisson_cdf(ixmax - 1, mean);
    }
    else {
      // If necessary, subtract from 0 to the beginning of the range
      if (ixmin <= mean) {
        return ROOT::Math::poisson_cdf(ixmax - 1, mean) - ROOT::Math::poisson_cdf(ixmin - 1, mean);
      }
      else {
        //Avoid catastrophic cancellation in the high tails:
        return ROOT::Math::poisson_cdf_c(ixmin - 1, mean) - ROOT::Math::poisson_cdf_c(ixmax - 1, mean);
      }

    }

  } else if(code == 2) {

    // the integral with respect to the mean is the integral of a gamma distribution
    Double_t mean_min = mean.min(rangeName);
    Double_t mean_max = mean.max(rangeName);

    Double_t ix;
    if(_noRounding) ix = x + 1;
    else ix = Int_t(TMath::Floor(x)) + 1.0; // negative ix does not need protection (gamma returns 0.0)

    return ROOT::Math::gamma_cdf(mean_max, ix, 1.0) - ROOT::Math::gamma_cdf(mean_min, ix, 1.0);
  }

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Advertise internal generator in x

Int_t RooPoisson::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////
/// Implement internal generator using TRandom::Poisson

void RooPoisson::generateEvent(Int_t code)
{
  R__ASSERT(code==1) ;
  Double_t xgen ;
  while(1) {
    xgen = RooRandom::randomGenerator()->Poisson(mean);
    if (xgen<=x.max() && xgen>=x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}
