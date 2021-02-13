/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *    T. Skwarnicki modify RooCBShape to Symmetrical Double-Sided CB         *
 *    Michael Wilkinson add to RooFit source                                 *
 *****************************************************************************/

/** \class RooCBShape
    \ingroup Roofit

PDF implementing the Symmetrical Double-Sided Crystall Ball line shape.
\f[
  f(m;m_0,\sigma,\alpha,n) =
  \begin{cases}
    \exp \left( - \frac{1}{2} \cdot \left[ \frac{m - m_0}{\sigma} \right]^2 \right), & \mbox{for }\left| \frac{m - m_0}{\sigma} \right| < |\alpha| \\
    A \cdot (B + \left| \frac{m - m_0}{\sigma} \right|)^{-n}, & \mbox{for }\left| \frac{m - m_0}{\sigma} \right| \geq |\alpha|, \\
  \end{cases}
\f]
times some normalization factor,
where
\f[
  \begin{align}
    A &= \left(\frac{n}{\left| \alpha \right|}\right)^n \cdot \exp\left(- \frac {\left| \alpha \right|^2}{2}\right) \\
    B &= \frac{n}{\left| \alpha \right|}  - \left| \alpha \right| \\
  \end{align}
\f]
**/

#include "RooSDSCBShape.h"
#include "RooFit.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"

#include "TMath.h"

#include <exception>
#include <math.h>

using namespace std;

ClassImp(RooSDSCBShape);

////////////////////////////////////////////////////////////////////////////////

Double_t RooSDSCBShape::ApproxErf(Double_t arg) const
{
  constexpr double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;

  return RooMath::erf(arg);
}


////////////////////////////////////////////////////////////////////////////////

RooSDSCBShape::RooSDSCBShape(const char *name, const char *title,
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

RooSDSCBShape::RooSDSCBShape(const RooSDSCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alpha("alpha", this, other.alpha),
  n("n", this, other.n)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooSDSCBShape::evaluate() const {

  Double_t t = fabs( (Double_t)( (m-m0)/sigma ) );
  Double_t absAlpha = fabs( (Double_t)alpha);

  if (t < absAlpha) {
    return exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t b= n/absAlpha - absAlpha;

    return a/TMath::Power(b + t, n);
  }

}

////////////////////////////////////////////////////////////////////////////////

Int_t RooSDSCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooSDSCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{
  constexpr double sqrtPiOver2 = 1.2533141373;
  constexpr double sqrt2 = 1.4142135624;

  R__ASSERT(code==1);
  (void)code;  // suppress unused warning
  double result = 0.0;
  bool useLog = false;

  if( fabs(n-1.0) < 1.0e-05 )
    useLog = true;

  double sig = fabs((Double_t)sigma);

  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;

  double absAlpha = fabs((Double_t)alpha);

  if( ( tmin >= -absAlpha ) && ( tmax <= absAlpha ) ) {
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
  else if( tmin >= absAlpha ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;

    if(useLog) {
      result += a*sig*( log(b+tmax) - log(b+tmin) );
    }
    else {
      result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b+tmax,n-1.0))
                                - 1.0/(TMath::Power(b+tmin,n-1.0)) );
    }
  }
  else if( ( tmin < -absAlpha ) && ( tmax <= absAlpha ) ) {
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
  else if( ( tmin < -absAlpha ) && ( tmax > absAlpha ) ) {
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

    double term2 = sig*sqrtPiOver2*(   ApproxErf( absAlpha/sqrt2)
                                     - ApproxErf(-absAlpha/sqrt2) );


    double term3 = 0.0;
    if(useLog) {
      term3 = a*sig*(  log(b+tmax) - log(n/absAlpha));
    }
    else {
      term3 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b+tmax,n-1.0))
                              - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }

    result += term1 + term2 + term3;
  }
  else if( ( tmin >= -absAlpha ) && ( tmax > absAlpha ) ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;

    double term1 = 0.0;
    if(useLog) {
      term1 = a*sig*(  log(b+tmax) - log(n/absAlpha));
    }
    else {
      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b+tmax,n-1.0))
                              - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }

    double term2 = sig*sqrtPiOver2*(   ApproxErf(-tmin/sqrt2)
                                     - ApproxErf(-absAlpha/sqrt2) );


    result += term1 + term2;
  }
  else {
    throw std::logic_error("This should never happen! Error in RooSDSCBShape logic. Please file a bug report.");
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise that we know the maximum of self for given (m0,alpha,n,sigma)

Int_t RooSDSCBShape::getMaxVal(const RooArgSet& vars) const
{
  RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
    return 1 ;
  }
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooSDSCBShape::maxVal(Int_t code) const
{
  R__ASSERT(code==1) ;

  // The maximum value for given (m0,alpha,n,sigma)
  return 1.0 ;
}
