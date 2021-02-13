/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *    T. Skwarnicki modify RooCBShape to Asymmetrical Double-Sided CB                     *
 *    Michael Wilkinson add to RooFit source                                 *
 *****************************************************************************/

/** \class RooCBShape
    \ingroup Roofit

PDF implementing the Asymmetrical Double-Sided Crystall Ball line shape.
\f[
  f(m;m_0,\sigma,\alpha_L,n_L,\alpha_R,n_R) =
  \begin{cases}
    A_L \cdot (B_L - \frac{m - m_0}{\sigma})^{-n_L}, & \mbox{for }\frac{m - m_0}{\sigma} < -\alpha_L \\
    \exp \left( - \frac{1}{2} \cdot \left[ \frac{m - m_0}{\sigma} \right]^2 \right), & \mbox{for }\frac{m - m_0}{\sigma} \leq \alpha_R \\
    A_R \cdot (B_R + \frac{m - m_0}{\sigma})^{-n_R}, & \mbox{otherwise}, \\
  \end{cases}
\f]
times some normalization factor,
where
\f[
  \begin{align}
    A_i &= \left(\frac{n_i}{\left| \alpha_i \right|}\right)^{n_i} \cdot \exp\left(- \frac {\left| \alpha_i \right|^2}{2}\right) \\
    B_i &= \frac{n_i}{\left| \alpha_i \right|}  - \left| \alpha_i \right| \\
  \end{align}
\f]
**/

#include "RooDSCBShape.h"
#include "RooFit.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"

#include "TMath.h"

#include <exception>
#include <math.h>

using namespace std;

ClassImp(RooDSCBShape);

////////////////////////////////////////////////////////////////////////////////

Double_t RooDSCBShape::ApproxErf(Double_t arg) const
{
  constexpr double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;

  return RooMath::erf(arg);
}

////////////////////////////////////////////////////////////////////////////////

RooDSCBShape::RooDSCBShape(const char *name, const char *title,
		       RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
		       RooAbsReal& _alphaL, RooAbsReal& _nL,
                       RooAbsReal& _alphaR, RooAbsReal& _nR) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaL("alphaL", "Left Alpha", this, _alphaL),
  nL("nL", "Left Order", this, _nL),
  alphaR("alphaR", "Right Alpha", this, _alphaR),
  nR("nR", "Right Order", this, _nR)
{
}

////////////////////////////////////////////////////////////////////////////////

RooDSCBShape::RooDSCBShape(const RooDSCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), alphaL("alphaL", this, other.alphaL),
  nL("nL", this, other.nL), alphaR("alphaR", this, other.alphaR),
  nR("nR", this, other.nR)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooDSCBShape::evaluate() const {

  Double_t t = (Double_t)( (m-m0)/sigma );
  Double_t absAlphaL = fabs( (Double_t)alphaL);
  Double_t absAlphaR = fabs( (Double_t)alphaR);


  if (t < -alphaL) {
    Double_t a =  TMath::Power(nL/absAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t b= nL/absAlphaL - absAlphaL;

    return a/TMath::Power(b - t, nL);
  }
  else if(t<=alphaR){
    return exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(nR/absAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t b= nR/absAlphaR - absAlphaR;

    return a/TMath::Power(b + t, nR);
  }

}

////////////////////////////////////////////////////////////////////////////////

Int_t RooDSCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooDSCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{
  constexpr double sqrtPiOver2 = 1.2533141373;
  constexpr double sqrt2 = 1.4142135624;

  R__ASSERT(code==1);
  double result = 0.0;
  bool useLogL = false;
  bool useLogR = false;

  if( fabs(nL-1.0) < 1.0e-05 )
    useLogL = true;
  if( fabs(nR-1.0) < 1.0e-05 )
    useLogR = true;

  double sig = fabs((Double_t)sigma);

  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;

  double absAlphaL = fabs((Double_t)alphaL);
  double absAlphaR = fabs((Double_t)alphaR);

  if( ( tmin >= -absAlphaL ) && ( tmax <= absAlphaR ) ) {
    result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                - ApproxErf(tmin/sqrt2) );
  }
  else if( tmax <= -absAlphaL ) {
    double a = TMath::Power(nL/absAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    double b = nL/absAlphaL - absAlphaL;

    if(useLogL) {
      result += a*sig*( log(b-tmin) - log(b-tmax) );
    }
    else {
      result += a*sig/(1.0-nL)*(   1.0/(TMath::Power(b-tmin,nL-1.0))
                                - 1.0/(TMath::Power(b-tmax,nL-1.0)) );
    }
  }
  else if( tmin >= absAlphaR ) {
    double a = TMath::Power(nR/absAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    double b = nR/absAlphaR - absAlphaR;

    if(useLogR) {
      result += a*sig*( log(b+tmax) - log(b+tmin) );
    }
    else {
      result += a*sig/(1.0-nR)*(   1.0/(TMath::Power(b+tmax,nR-1.0))
                                - 1.0/(TMath::Power(b+tmin,nR-1.0)) );
    }
  }
  else if( ( tmin < -absAlphaL ) && ( tmax <= absAlphaR ) ) {
    double a = TMath::Power(nL/absAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    double b = nL/absAlphaL - absAlphaL;

    double term1 = 0.0;
    if(useLogL) {
      term1 = a*sig*(  log(b-tmin) - log(nL/absAlphaL));
    }
    else {
      term1 = a*sig/(1.0-nL)*(   1.0/(TMath::Power(b-tmin,nL-1.0))
                              - 1.0/(TMath::Power(nL/absAlphaL,nL-1.0)) );
    }

    double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                     - ApproxErf(-absAlphaL/sqrt2) );

    result += term1 + term2;
  }
  else if( ( tmin < -absAlphaL ) && ( tmax > absAlphaR ) ) {
    double aL = TMath::Power(nL/absAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    double aR = TMath::Power(nR/absAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    double bL = nL/absAlphaL - absAlphaL;
    double bR = nR/absAlphaR - absAlphaR;

    double term1 = 0.0;
    if(useLogL) {
      term1 = aL*sig*(  log(bL-tmin) - log(nL/absAlphaL));
    }
    else {
      term1 = aL*sig/(1.0-nL)*(   1.0/(TMath::Power(bL-tmin,nL-1.0))
                              - 1.0/(TMath::Power(nL/absAlphaL,nL-1.0)) );
    }
    double term2 = sig*sqrtPiOver2*(   ApproxErf( absAlphaR/sqrt2)
                                     - ApproxErf(-absAlphaL/sqrt2) );

    double term3 = 0.0;
    if(useLogR) {
      term3 = aR*sig*(  log(bR+tmax) - log(nR/absAlphaR));
    }
    else {
      term3 = aR*sig/(1.0-nR)*(   1.0/(TMath::Power(bR+tmax,nR-1.0))
                              - 1.0/(TMath::Power(nR/absAlphaR,nR-1.0)) );
    }
    result += term1 + term2 + term3;
  }
  else if( ( tmin >= -absAlphaL ) && ( tmax > absAlphaR ) ) {
    double a = TMath::Power(nR/absAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    double b = nR/absAlphaR - absAlphaR;

    double term1 = 0.0;
    if(useLogR) {
      term1 = a*sig*(  log(b+tmax) - log(nR/absAlphaR));
    }
    else {
      term1 = a*sig/(1.0-nR)*(   1.0/(TMath::Power(b+tmax,nR-1.0))
                              - 1.0/(TMath::Power(nR/absAlphaR,nR-1.0)) );
    }

    double term2 = sig*sqrtPiOver2*(   ApproxErf(absAlphaR/sqrt2)
                                     - ApproxErf(tmin/sqrt2) );


    result += term1 + term2;
  }
  else {
    throw std::logic_error("This should never happen! Error in RooDSCBShape logic. Please file a bug report.");
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Advertise that we know the maximum of self for given (m0,alpha,n,sigma)

Int_t RooDSCBShape::getMaxVal(const RooArgSet& vars) const
{
  RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
    return 1 ;
  }
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooDSCBShape::maxVal(Int_t code) const
{
  R__ASSERT(code==1) ;

  // The maximum value for given (m0,alpha,n,sigma)
  return 1.0 ;
}
