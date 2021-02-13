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

/** \class RooGaussModel
    \ingroup Roofit

Class RooGaussModel implements a RooResolutionModel that models a Gaussian
distribution. Object of class RooGaussModel can be used
for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
**/

#include "RooFit.h"

#include "TMath.h"
#include "Riostream.h"
#include "RooGaussModel.h"
#include "RooRealConstant.h"
#include "RooRandom.h"

#include "TError.h"

using namespace std;

ClassImp(RooGaussModel);

////////////////////////////////////////////////////////////////////////////////

RooGaussModel::RooGaussModel(const char *name, const char *title, RooAbsRealLValue& xIn,
              RooAbsReal& _mean, RooAbsReal& _sigma) :
  RooResolutionModel(name,title,xIn),
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,(RooRealVar&)RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,(RooRealVar&)RooRealConstant::value(1))
{
}

////////////////////////////////////////////////////////////////////////////////

RooGaussModel::RooGaussModel(const char *name, const char *title, RooAbsRealLValue& xIn,
              RooAbsReal& _mean, RooAbsReal& _sigma,
              RooAbsReal& _msSF) :
  RooResolutionModel(name,title,xIn),
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_msSF),
  ssf("ssf","Sigma Scale Factor",this,_msSF)
{
}

////////////////////////////////////////////////////////////////////////////////

RooGaussModel::RooGaussModel(const char *name, const char *title, RooAbsRealLValue& xIn,
              RooAbsReal& _mean, RooAbsReal& _sigma,
              RooAbsReal& _meanSF, RooAbsReal& _sigmaSF) :
  RooResolutionModel(name,title,xIn),
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{
}

////////////////////////////////////////////////////////////////////////////////

RooGaussModel::RooGaussModel(const RooGaussModel& other, const char* name) :
  RooResolutionModel(other,name),
  _flatSFInt(other._flatSFInt),
  _asympInt(other._asympInt),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  msf("msf",this,other.msf),
  ssf("ssf",this,other.ssf)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooGaussModel::~RooGaussModel()
{
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGaussModel::basisCode(const char* name) const
{
  if (!TString("exp(-@0/@1)").CompareTo(name)) return expBasisPlus ;
  if (!TString("exp(@0/@1)").CompareTo(name)) return expBasisMinus ;
  if (!TString("exp(-abs(@0)/@1)").CompareTo(name)) return expBasisSum ;
  if (!TString("exp(-@0/@1)*sin(@0*@2)").CompareTo(name)) return sinBasisPlus ;
  if (!TString("exp(@0/@1)*sin(@0*@2)").CompareTo(name)) return sinBasisMinus ;
  if (!TString("exp(-abs(@0)/@1)*sin(@0*@2)").CompareTo(name)) return sinBasisSum ;
  if (!TString("exp(-@0/@1)*cos(@0*@2)").CompareTo(name)) return cosBasisPlus ;
  if (!TString("exp(@0/@1)*cos(@0*@2)").CompareTo(name)) return cosBasisMinus ;
  if (!TString("exp(-abs(@0)/@1)*cos(@0*@2)").CompareTo(name)) return cosBasisSum ;
  if (!TString("(@0/@1)*exp(-@0/@1)").CompareTo(name)) return linBasisPlus ;
  if (!TString("(@0/@1)*(@0/@1)*exp(-@0/@1)").CompareTo(name)) return quadBasisPlus ;
  if (!TString("exp(-@0/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisPlus;
  if (!TString("exp(@0/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisMinus;
  if (!TString("exp(-abs(@0)/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisSum;
  if (!TString("exp(-@0/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisPlus;
  if (!TString("exp(@0/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisMinus;
  if (!TString("exp(-abs(@0)/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisSum;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooGaussModel::evaluate() const
{
  // *** 1st form: Straight Gaussian, used for unconvoluted PDF or expBasis with 0 lifetime ***
  static Double_t root2(std::sqrt(2.)) ;
  static Double_t root2pi(std::sqrt(2.*std::atan2(0.,-1.))) ;
  static Double_t rootpi(std::sqrt(std::atan2(0.,-1.))) ;

  BasisType basisType = (BasisType)( (_basisCode == 0) ? 0 : (_basisCode/10) + 1 );
  BasisSign basisSign = (BasisSign)( _basisCode - 10*(basisType-1) - 2 ) ;

  Double_t tau = (_basisCode!=noBasis)?((RooAbsReal*)basis().getParameter(1))->getVal():0 ;
  if (basisType == coshBasis && _basisCode!=noBasis ) {
     Double_t dGamma = ((RooAbsReal*)basis().getParameter(2))->getVal();
     if (dGamma==0) basisType = expBasis;
  }

  if (basisType==none || ((basisType==expBasis || basisType==cosBasis) && tau==0.)) {
    Double_t xprime = (x-(mean*msf))/(sigma*ssf) ;
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 1st form" << endl ;

    Double_t result = std::exp(-0.5*xprime*xprime)/(sigma*ssf*root2pi) ;
    if (_basisCode!=0 && basisSign==Both) result *= 2 ;
    return result ;
  }

  // *** 2nd form: 0, used for sinBasis, linBasis, and quadBasis with tau=0 ***
  if (tau==0) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  // *** 3nd form: Convolution with exp(-t/tau), used for expBasis and cosBasis(omega=0) ***
  Double_t omega =  (basisType==sinBasis  || basisType==cosBasis)  ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t dgamma = (basisType==sinhBasis || basisType==coshBasis) ? ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t _x = omega *tau ;
  Double_t _y = tau*dgamma/2;
  Double_t xprime = (x-(mean*msf))/tau ;
  Double_t c = (sigma*ssf)/(root2*tau) ;
  Double_t u = xprime/(2*c) ;

  if (basisType==expBasis || (basisType==cosBasis && _x==0.)) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 3d form tau=" << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += evalCerf(0,-u,c).real();
    if (basisSign!=Plus)  result += evalCerf(0, u,c).real();
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::getVal(" << GetName() << ") got nan during case 1 " << endl; }
    return result ;
  }

  // *** 4th form: Convolution with exp(-t/tau)*sin(omega*t), used for sinBasis(omega<>0,tau<>0) ***
  if (basisType==sinBasis) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 4th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (_x==0.) return result ;
    if (basisSign!=Minus) result += -evalCerf(-_x,-u,c).imag();
    if (basisSign!=Plus)  result += -evalCerf( _x, u,c).imag();
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooGaussModel::getVal(" << GetName() << ") got nan during case 3 " << endl;
    return result ;
  }

  // *** 5th form: Convolution with exp(-t/tau)*cos(omega*t), used for cosBasis(omega<>0) ***
  if (basisType==cosBasis) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 5th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += evalCerf(-_x,-u,c).real();
    if (basisSign!=Plus)  result += evalCerf( _x, u,c).real();
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooGaussModel::getVal(" << GetName() << ") got nan during case 4 " << endl;
    return result ;
  }

  // ***8th form: Convolution with exp(-|t|/tau)*cosh(dgamma*t/2), used for         coshBasisSum ***
  if (basisType==coshBasis || basisType ==sinhBasis) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 8th form tau = " << tau << endl ;
    Double_t result(0);
    int sgn = ( basisType == coshBasis ? +1 : -1 );
    if (basisSign!=Minus) result += 0.5*(    evalCerf(0,-u,c*(1-_y)).real()+sgn*evalCerf(0,-u,c*(1+_y)).real()) ;
    if (basisSign!=Plus)  result += 0.5*(sgn*evalCerf(0, u,c*(1-_y)).real()+    evalCerf(0, u,c*(1+_y)).real()) ;
    if (TMath::IsNaN(result)) cxcoutE(Tracing) << "RooGaussModel::getVal(" << GetName() << ") got nan during case 8 " << endl;
    return result ;
  }

  // *** 6th form: Convolution with (t/tau)*exp(-t/tau), used for linBasis ***
  if (basisType==linBasis) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 6th form tau = " << tau << endl ;
    R__ASSERT(basisSign==Plus);  // This should only be for positive times

    Double_t f0 = std::exp(-xprime+c*c) * RooMath::erfc(-u+c);
    Double_t f1 = std::exp(-u*u);
    return (xprime - 2*c*c)*f0 + (2*c/rootpi)*f1 ;
  }

  // *** 7th form: Convolution with (t/tau)^2*exp(-t/tau), used for quadBasis ***
  if (basisType==quadBasis) {
    if (verboseEval()>2) cout << "RooGaussModel::evaluate(" << GetName() << ") 7th form tau = " << tau << endl ;
    R__ASSERT(basisSign==Plus);  // This should only be for positive times

    Double_t f0 = std::exp(-xprime+c*c) * RooMath::erfc(-u+c);
    Double_t f1 = std::exp(-u*u);
    Double_t x2c2 = xprime - 2*c*c;
    return ( x2c2*x2c2*f0 + (2*c/rootpi)*x2c2*f1 + 2*c*c*f0 );
  }

  R__ASSERT(0) ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGaussModel::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  switch(_basisCode) {

  // Analytical integration capability of raw PDF
  case noBasis:

    // Optionally advertise flat integral over sigma scale factor
    if (_flatSFInt) {
      if (matchArgs(allVars,analVars,RooArgSet(convVar(),ssf.arg()))) return 2 ;
    }

    if (matchArgs(allVars,analVars,convVar())) return 1 ;
    break ;

  // Analytical integration capability of convoluted PDF
  case expBasisPlus:
  case expBasisMinus:
  case expBasisSum:
  case sinBasisPlus:
  case sinBasisMinus:
  case sinBasisSum:
  case cosBasisPlus:
  case cosBasisMinus:
  case cosBasisSum:
  case linBasisPlus:
  case quadBasisPlus:
  case coshBasisMinus:
  case coshBasisPlus:
  case coshBasisSum:
  case sinhBasisMinus:
  case sinhBasisPlus:
  case sinhBasisSum:

    // Optionally advertise flat integral over sigma scale factor
    if (_flatSFInt) {

      if (matchArgs(allVars,analVars,RooArgSet(convVar(),ssf.arg()))) {
   return 2 ;
      }
    }

    if (matchArgs(allVars,analVars,convVar())) return 1 ;
    break ;
  }

  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooGaussModel::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static const Double_t root2 = std::sqrt(2.) ;
  //static Double_t rootPiBy2 = std::sqrt(std::atan2(0.0,-1.0)/2.0);
  static const Double_t rootpi = std::sqrt(std::atan2(0.0,-1.0));
  Double_t ssfInt(1.0) ;

  // Code must be 1 or 2
  R__ASSERT(code==1||code==2) ;
  if (code==2) ssfInt = (ssf.max(rangeName)-ssf.min(rangeName)) ;

  BasisType basisType = (BasisType)( (_basisCode == 0) ? 0 : (_basisCode/10) + 1 );
  BasisSign basisSign = (BasisSign)( _basisCode - 10*(basisType-1) - 2 ) ;

  // *** 1st form: Straight Gaussian, used for unconvoluted PDF or expBasis with 0 lifetime ***
  Double_t tau = (_basisCode!=noBasis)?((RooAbsReal*)basis().getParameter(1))->getVal():0 ;
  if (basisType == coshBasis && _basisCode!=noBasis ) {
     Double_t dGamma = ((RooAbsReal*)basis().getParameter(2))->getVal();
     if (dGamma==0) basisType = expBasis;
  }
  if (basisType==none || ((basisType==expBasis || basisType==cosBasis) && tau==0.)) {
    Double_t xscale = root2*(sigma*ssf);
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 1st form" << endl ;

    Double_t xpmin = (x.min(rangeName)-(mean*msf))/xscale ;
    Double_t xpmax = (x.max(rangeName)-(mean*msf))/xscale ;

    Double_t result ;
    if (_asympInt) { // modified FMV, 07/24/03
        result = 1.0 ;
    } else {
       result = 0.5*(RooMath::erf(xpmax)-RooMath::erf(xpmin)) ;
    }

    if (_basisCode!=0 && basisSign==Both) result *= 2 ;
    //cout << "Integral 1st form " << " result= " << result*ssfInt << endl;
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 1 " << endl; }
    return result*ssfInt ;
  }


  Double_t omega = ((basisType==sinBasis)||(basisType==cosBasis)) ?  ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;
  Double_t dgamma =((basisType==sinhBasis)||(basisType==coshBasis)) ?  ((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;

  // *** 2nd form: unity, used for sinBasis and linBasis with tau=0 (PDF is zero) ***
  if (tau==0) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  // *** 3rd form: Convolution with exp(-t/tau), used for expBasis and cosBasis(omega=0) ***
  Double_t c = (sigma*ssf)/(root2*tau) ;
  Double_t xpmin = (x.min(rangeName)-(mean*msf))/tau ;
  Double_t xpmax = (x.max(rangeName)-(mean*msf))/tau ;
  Double_t umin = xpmin/(2*c) ;
  Double_t umax = xpmax/(2*c) ;

  if (basisType==expBasis || (basisType==cosBasis && omega==0.)) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 3d form tau=" << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += evalCerfInt(+1,0,tau,-umin,-umax,c).real();
    if (basisSign!=Plus)  result += evalCerfInt(-1,0,tau, umin, umax,c).real();
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 3 " << endl; }
    return result*ssfInt ;
  }

  // *** 4th form: Convolution with exp(-t/tau)*sin(omega*t), used for sinBasis(omega<>0,tau<>0) ***
  Double_t _x = omega * tau ;
  Double_t _y = tau*dgamma/2;

  if (basisType==sinBasis) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 4th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (_x==0) return result*ssfInt ;
    if (basisSign!=Minus) result += -1*evalCerfInt(+1,-_x,tau,-umin,-umax,c).imag();
    if (basisSign!=Plus)  result += -1*evalCerfInt(-1, _x,tau, umin, umax,c).imag();
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 4 " << endl; }
    return result*ssfInt ;
  }

  // *** 5th form: Convolution with exp(-t/tau)*cos(omega*t), used for cosBasis(omega<>0) ***
  if (basisType==cosBasis) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 5th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += evalCerfInt(+1,-_x,tau,-umin,-umax,c).real();
    if (basisSign!=Plus)  result += evalCerfInt(-1, _x,tau, umin, umax,c).real();
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 5 " << endl; }
    return result*ssfInt ;
  }

  // *** 8th form: Convolution with exp(-|t|/tau)*cosh(dgamma*t/2), used for coshBasis ***
  // *** 9th form: Convolution with exp(-|t|/tau)*sinh(dgamma*t/2), used for sinhBasis ***
  if (basisType==coshBasis || basisType == sinhBasis) {
    if (verboseEval()>0) {cout << "RooGaussModel::analyticalIntegral(" << GetName()                             << ") 8th form tau=" << tau << endl ; }
    Double_t result(0) ;
    int sgn = ( basisType == coshBasis ? +1 : -1 );
    if (basisSign!=Minus) result += 0.5*(    evalCerfInt(+1,0,tau/(1-_y),-umin,-umax,c*(1-_y)).real()+ sgn*evalCerfInt(+1,0,tau/(1+_y),-umin,-umax,c*(1+_y)).real());
    if (basisSign!=Plus)  result += 0.5*(sgn*evalCerfInt(-1,0,tau/(1-_y), umin, umax,c*(1-_y)).real()+     evalCerfInt(-1,0,tau/(1+_y), umin, umax,c*(1+_y)).real());
    if (TMath::IsNaN(result)) { cxcoutE(Tracing) << "RooGaussModel::analyticalIntegral(" << GetName() << ") got nan during case 6 " << endl; }
    return result*ssfInt ;
  }

  // *** 6th form: Convolution with (t/tau)*exp(-t/tau), used for linBasis ***
  if (basisType==linBasis) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 6th form tau=" << tau << endl ;

    Double_t f0 = RooMath::erf(-umax) - RooMath::erf(-umin);
    Double_t f1 = std::exp(-umax*umax) - std::exp(-umin*umin);

    Double_t tmp1 = std::exp(-xpmax)*RooMath::erfc(-umax + c);
    Double_t tmp2 = std::exp(-xpmin)*RooMath::erfc(-umin + c);

    Double_t f2 = tmp1 - tmp2;
    Double_t f3 = xpmax*tmp1 - xpmin*tmp2;

    Double_t expc2 = std::exp(c*c);

    return -tau*(              f0 +
        (2*c/rootpi)*f1 +
        (1 - 2*c*c)*expc2*f2 +
          expc2*f3
      )*ssfInt;
  }

  // *** 7th form: Convolution with (t/tau)*(t/tau)*exp(-t/tau), used for quadBasis ***
  if (basisType==quadBasis) {
    if (verboseEval()>0) cout << "RooGaussModel::analyticalIntegral(" << GetName() << ") 7th form tau=" << tau << endl ;

    Double_t f0 = RooMath::erf(-umax) - RooMath::erf(-umin);

    Double_t tmpA1 = std::exp(-umax*umax);
    Double_t tmpA2 = std::exp(-umin*umin);

    Double_t f1 = tmpA1 - tmpA2;
    Double_t f2 = umax*tmpA1 - umin*tmpA2;

    Double_t tmpB1 = std::exp(-xpmax)*RooMath::erfc(-umax + c);
    Double_t tmpB2 = std::exp(-xpmin)*RooMath::erfc(-umin + c);

    Double_t f3 = tmpB1 - tmpB2;
    Double_t f4 = xpmax*tmpB1 - xpmin*tmpB2;
    Double_t f5 = xpmax*xpmax*tmpB1 - xpmin*xpmin*tmpB2;

    Double_t expc2 = std::exp(c*c);

    return -tau*( 2*f0 +
        (4*c/rootpi)*((1-c*c)*f1 + c*f2) +
        (2*c*c*(2*c*c-1) + 2)*expc2*f3 - (4*c*c-2)*expc2*f4 + expc2*f5
                )*ssfInt;
  }
  R__ASSERT(0) ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////
/// use the approximation: erf(z) = exp(-z*z)/(std::sqrt(pi)*z)
/// to explicitly cancel the divergent exp(y*y) behaviour of
/// CWERF for z = x + i y with large negative y

std::complex<Double_t> RooGaussModel::evalCerfApprox(Double_t _x, Double_t u, Double_t c)
{
  static const Double_t rootpi= std::sqrt(std::atan2(0.,-1.));
  const std::complex<Double_t> z(_x * c, u + c);
  const std::complex<Double_t> zc(u + c, - _x * c);
  const std::complex<Double_t> zsq((z.real() + z.imag()) * (z.real() - z.imag()),
     2. * z.real() * z.imag());
  const std::complex<Double_t> v(-zsq.real() - u*u, -zsq.imag());
  const std::complex<Double_t> ev = std::exp(v);
  const std::complex<Double_t> mez2zcrootpi = -std::exp(zsq)/(zc*rootpi);

  return 2. * (ev * (mez2zcrootpi + 1.));
}

////////////////////////////////////////////////////////////////////////////////

std::complex<Double_t> RooGaussModel::evalCerfInt(Double_t sign, Double_t _x, Double_t tau, Double_t umin, Double_t umax, Double_t c) const
{
  std::complex<Double_t> diff(2., 0.);
  if (!_asympInt) {
    diff = evalCerf(_x,umin,c);
    diff -= evalCerf(_x,umax,c);
    diff += RooMath::erf(umin) - RooMath::erf(umax);
    diff *= sign;
  }
  diff *= std::complex<Double_t>(1., _x);
  diff *= tau / (1.+_x*_x);
  return diff;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGaussModel::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooGaussModel::generateEvent(Int_t code)
{
  R__ASSERT(code==1) ;
  Double_t xmin = x.min();
  Double_t xmax = x.max();
  TRandom *generator = RooRandom::randomGenerator();
  while(true) {
    Double_t xgen = generator->Gaus(mean*msf,sigma*ssf);
    if (xgen<xmax && xgen>xmin) {
      x = xgen ;
      return ;
    }
  }
}
