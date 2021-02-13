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

/** \class RooGExpModel
    \ingroup Roofit

The RooGExpModel is a RooResolutionModel implementation that models
a resolution function that is the convolution of a Gaussian with
a one-sided exponential. Such objects can be used
for analytical convolutions with classes inheriting from RooAbsAnaConvPdf.
\f[
  \mathrm{GExp} = \exp \left( -\frac{1}{2} \left(\frac{x-\mu}{\sigma} \right)^2 \right)^2
    \otimes \exp\left( -\frac{x}{\tau} \right)
\f]

**/

#include "RooGExpModel.h"

#include "RooFit.h"
#include "RooMath.h"
#include "RooRealConstant.h"
#include "RooRandom.h"
#include "TMath.h"


using namespace std;

ClassImp(RooGExpModel);

////////////////////////////////////////////////////////////////////////////////
/// Create a Gauss (x) Exp model with mean, sigma and tau parameters and scale factors for each parameter.
///
/// \note If scale factors for the parameters are not needed, `RooConst(1.)` can be passed.
///
/// \param[in] name Name of this instance.
/// \param[in] title Title (e.g. for plotting)
/// \param[in] xIn The convolution observable.
/// \param[in] meanIn The mean of the Gaussian.
/// \param[in] sigmaIn Width of the Gaussian.
/// \param[in] rlifeIn Lifetime constant \f$ \tau \f$.
/// \param[in] meanSF  Scale factor for mean.
/// \param[in] sigmaSF Scale factor for sigma.
/// \param[in] rlifeSF Scale factor for rlife.
/// \param[in] nlo   Include next-to-leading order for higher accuracy of convolution.
/// \param[in] type  Switch between normal and flipped model.
RooGExpModel::RooGExpModel(const char *name, const char *title, RooAbsRealLValue& xIn,
    RooAbsReal& meanIn, RooAbsReal& sigmaIn, RooAbsReal& rlifeIn,
    RooAbsReal& meanSF, RooAbsReal& sigmaSF, RooAbsReal& rlifeSF,
    Bool_t nlo, Type type) :
  RooResolutionModel(name, title, xIn),
  _mean("mean", "Mean of Gaussian component", this, meanIn),
  sigma("sigma", "Width", this, sigmaIn),
  rlife("rlife", "Life time", this, rlifeIn),
  _meanSF("meanSF", "Scale factor for mean", this, meanSF),
  ssf("ssf", "Sigma Scale Factor", this, sigmaSF),
  rsf("rsf", "RLife Scale Factor", this, rlifeSF),
  _flip(type==Flipped),
  _nlo(nlo),
  _flatSFInt(false),
  _asympInt(false)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Create a Gauss (x) Exp model with sigma and tau parameters.
///
/// \param[in] name Name of this instance.
/// \param[in] title Title (e.g. for plotting)
/// \param[in] xIn The convolution observable.
/// \param[in] _sigma Width of the Gaussian.
/// \param[in] _rlife Lifetime constant \f$ \tau \f$.
/// \param[in] nlo   Include next-to-leading order for higher accuracy of convolution.
/// \param[in] type  Switch between normal and flipped model.
RooGExpModel::RooGExpModel(const char *name, const char *title, RooAbsRealLValue& xIn,
            RooAbsReal& _sigma, RooAbsReal& _rlife,
            Bool_t nlo, Type type) :
  RooResolutionModel(name,title,xIn),
  _mean("mean", "Mean of Gaussian component", this, RooRealConstant::value(0.)),
  sigma("sigma","Width",this,_sigma),
  rlife("rlife","Life time",this,_rlife),
  _meanSF("meanSF", "Scale factor for mean", this, RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,(RooRealVar&)RooRealConstant::value(1)),
  rsf("rsf","RLife Scale Factor",this,(RooRealVar&)RooRealConstant::value(1)),
  _flip(type==Flipped),_nlo(nlo), _flatSFInt(kFALSE), _asympInt(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a Gauss (x) Exp model with sigma and tau parameters.
///
/// \param[in] name Name of this instance.
/// \param[in] title Title (e.g. for plotting)
/// \param[in] xIn The convolution observable.
/// \param[in] _sigma Width of the Gaussian.
/// \param[in] _rlife Lifetime constant \f$ \tau \f$.
/// \param[in] _rsSF Scale factor for both sigma and tau.
/// \param[in] nlo   Include next-to-leading order for higher accuracy of convolution.
/// \param[in] type  Switch between normal and flipped model.
RooGExpModel::RooGExpModel(const char *name, const char *title, RooAbsRealLValue& xIn,
            RooAbsReal& _sigma, RooAbsReal& _rlife,
            RooAbsReal& _rsSF,
            Bool_t nlo, Type type) :
  RooResolutionModel(name,title,xIn),
  _mean("mean", "Mean of Gaussian component", this, RooRealConstant::value(0.)),
  sigma("sigma","Width",this,_sigma),
  rlife("rlife","Life time",this,_rlife),
  _meanSF("meanSF", "Scale factor for mean", this, RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,_rsSF),
  rsf("rsf","RLife Scale Factor",this,_rsSF),
  _flip(type==Flipped),
  _nlo(nlo),
  _flatSFInt(kFALSE),
  _asympInt(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a Gauss (x) Exp model with sigma and tau parameters and scale factors.
///
/// \param[in] name Name of this instance.
/// \param[in] title Title (e.g. for plotting)
/// \param[in] xIn The convolution observable.
/// \param[in] _sigma Width of the Gaussian.
/// \param[in] _rlife Lifetime constant \f$ \tau \f$.
/// \param[in] _sigmaSF Scale factor for sigma.
/// \param[in] _rlifeSF Scale factor for rlife.
/// \param[in] nlo   Include next-to-leading order for higher accuracy of convolution.
/// \param[in] type  Switch between normal and flipped model.
RooGExpModel::RooGExpModel(const char *name, const char *title, RooAbsRealLValue& xIn,
            RooAbsReal& _sigma, RooAbsReal& _rlife,
            RooAbsReal& _sigmaSF, RooAbsReal& _rlifeSF,
            Bool_t nlo, Type type) :
  RooResolutionModel(name,title,xIn),
  _mean("mean", "Mean of Gaussian component", this, RooRealConstant::value(0.)),
  sigma("sigma","Width",this,_sigma),
  rlife("rlife","Life time",this,_rlife),
  _meanSF("meanSF", "Scale factor for mean", this, RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF),
  rsf("rsf","RLife Scale Factor",this,_rlifeSF),
  _flip(type==Flipped),
  _nlo(nlo),
  _flatSFInt(kFALSE),
  _asympInt(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////

RooGExpModel::RooGExpModel(const RooGExpModel& other, const char* name) :
  RooResolutionModel(other,name),
  _mean("mean", this, other._mean),
  sigma("sigma",this,other.sigma),
  rlife("rlife",this,other.rlife),
  _meanSF("meanSf", this, other._meanSF),
  ssf("ssf",this,other.ssf),
  rsf("rsf",this,other.rsf),
  _flip(other._flip),
  _nlo(other._nlo),
  _flatSFInt(other._flatSFInt),
  _asympInt(other._asympInt)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooGExpModel::~RooGExpModel()
{
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGExpModel::basisCode(const char* name) const
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
  if (!TString("exp(-@0/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisPlus;
  if (!TString("exp(@0/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisMinus;
  if (!TString("exp(-abs(@0)/@1)*sinh(@0*@2/2)").CompareTo(name)) return sinhBasisSum;
  if (!TString("exp(-@0/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisPlus;
  if (!TString("exp(@0/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisMinus;
  if (!TString("exp(-abs(@0)/@1)*cosh(@0*@2/2)").CompareTo(name)) return coshBasisSum;
  return 0 ;
}


namespace {
////////////////////////////////////////////////////////////////////////////////
/// Approximation of the log of the complex error function
Double_t logErfC(Double_t xx)
{
  Double_t t,z,ans;
  z=fabs(xx);
  t=1.0/(1.0+0.5*z);

  if(xx >= 0.0)
    ans=log(t)+(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+
   t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  else
    ans=log(2.0-t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+
        t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277))))))))));

  return ans;
}

////////////////////////////////////////////////////////////////////////////////
/// use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z)
/// to explicitly cancel the divergent exp(y*y) behaviour of
/// CWERF for z = x + i y with large negative y

std::complex<Double_t> evalCerfApprox(Double_t swt, Double_t u, Double_t c)
{
  static Double_t rootpi= sqrt(atan2(0.,-1.));
  std::complex<Double_t> z(swt*c,u+c);
  std::complex<Double_t> zc(u+c,-swt*c);
  std::complex<Double_t> zsq= z*z;
  std::complex<Double_t> v= -zsq - u*u;

  return std::exp(v)*(-std::exp(zsq)/(zc*rootpi) + 1.)*2.;
}


// Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
std::complex<Double_t> evalCerf(Double_t swt, Double_t u, Double_t c)
{
  std::complex<Double_t> z(swt*c,u+c);
  return (z.imag()>-4.0) ? RooMath::faddeeva_fast(z)*std::exp(-u*u) : evalCerfApprox(swt,u,c) ;
}


// Calculate Re(exp(-u^2) cwerf(i(u+c)))
// added FMV, 08/17/03
inline Double_t evalCerfRe(Double_t u, Double_t c) {
  Double_t expArg = u*2*c+c*c ;
  if (expArg<300) {
     return exp(expArg) * RooMath::erfc(u+c);
  } else {
     return exp(expArg+logErfC(u+c));
  }
}

}



////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::evaluate() const
{
  static Double_t root2(sqrt(2.)) ;

  BasisType basisType = (BasisType)( (_basisCode == 0) ? 0 : (_basisCode/10) + 1 );
  BasisSign basisSign = (BasisSign)( _basisCode - 10*(basisType-1) - 2 ) ;

  Double_t fsign = _flip?-1:1 ;

  Double_t sig = sigma*ssf ;
  Double_t rtau = rlife*rsf ;

  Double_t tau = (_basisCode!=noBasis)?((RooAbsReal*)basis().getParameter(1))->getVal():0. ;
  // added, FMV 07/27/03
  if (basisType == coshBasis && _basisCode!=noBasis ) {
     Double_t dGamma = ((RooAbsReal*)basis().getParameter(2))->getVal();
     if (dGamma==0) basisType = expBasis;
  }

  // *** 1st form: Straight GExp, used for unconvoluted PDF or expBasis with 0 lifetime ***
  if (basisType==none || ((basisType==expBasis || basisType==cosBasis) && tau==0.)) {
    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName() << ") 1st form" << endl ;

    Double_t expArg = sig*sig/(2*rtau*rtau) + fsign*(x - _mean*_meanSF)/rtau ;

    Double_t result ;
    if (expArg<300) {
      result = 1/(2*rtau) * exp(expArg) * RooMath::erfc(sig/(root2*rtau) + fsign*(x - _mean*_meanSF)/(root2*sig));
    } else {
      // If exponent argument is very large, bring canceling RooMath::erfc() term inside exponent
      // to avoid floating point over/underflows of intermediate calculations
      result = 1/(2*rtau) * exp(expArg + logErfC(sig/(root2*rtau) + fsign*(x - _mean*_meanSF)/(root2*sig))) ;
    }

//     Double_t result = 1/(2*rtau)
//                     * exp(sig*sig/(2*rtau*rtau) + fsign*x/rtau)
//                     * RooMath::erfc(sig/(root2*rtau) + fsign*x/(root2*sig));

    // equivalent form, added FMV, 07/24/03
    //Double_t xprime = x/rtau ;
    //Double_t c = sig/(root2*rtau) ;
    //Double_t u = xprime/(2*c) ;
    //Double_t result = 0.5*evalCerf(fsign*u,c).real() ;  // sign=-1 !

    if (_basisCode!=0 && basisSign==Both) result *= 2 ;
    //cout << "1st form " << "x= " << x << " result= " << result << endl;
    return result ;
  }

  // *** 2nd form: 0, used for sinBasis and cosBasis with tau=0 ***
  if (tau==0) {
    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  Double_t omega = (basisType!=expBasis)?((RooAbsReal*)basis().getParameter(2))->getVal():0. ;

  // *** 3nd form: Convolution with exp(-t/tau), used for expBasis and cosBasis(omega=0) ***
  if (basisType==expBasis || (basisType==cosBasis && omega==0.)) {
    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName() << ") 3d form tau=" << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += calcDecayConv(+1,tau,sig,rtau,fsign) ;  // modified FMV,08/13/03
    if (basisSign!=Plus)  result += calcDecayConv(-1,tau,sig,rtau,fsign) ;  // modified FMV,08/13/03
    //cout << "3rd form " << "x= " << x << " result= " << result << endl;
    return result ;
  }

  // *** 4th form: Convolution with exp(-t/tau)*sin(omega*t), used for sinBasis(omega<>0,tau<>0) ***
  Double_t wt = omega *tau ;
  if (basisType==sinBasis) {
    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName() << ") 4th form omega = "
              << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (wt==0.) return result ;
    if (basisSign!=Minus) result += -1*calcSinConv(+1,sig,tau,omega,rtau,fsign).imag() ;
    if (basisSign!=Plus) result += -1*calcSinConv(-1,sig,tau,omega,rtau,fsign).imag() ;
    //cout << "4th form " << "x= " << x << " result= " << result << endl;
    return result ;
  }

  // *** 5th form: Convolution with exp(-t/tau)*cos(omega*t), used for cosBasis(omega<>0) ***
  if (basisType==cosBasis) {
    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName()
              << ") 5th form omega = " << omega << ", tau = " << tau << endl ;
    Double_t result(0) ;
    if (basisSign!=Minus) result += calcSinConv(+1,sig,tau,omega,rtau,fsign).real() ;
    if (basisSign!=Plus) result += calcSinConv(-1,sig,tau,omega,rtau,fsign).real() ;
    //cout << "5th form " << "x= " << x << " result= " << result << endl;
    return result ;
  }


  // *** 6th form: Convolution with exp(-t/tau)*sinh(dgamma*t/2), used for sinhBasis ***
  if (basisType==sinhBasis) {
    Double_t dgamma = ((RooAbsReal*)basis().getParameter(2))->getVal();

    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName()
              << ") 6th form = " << dgamma << ", tau = " << tau << endl;
    Double_t result(0);
    //if (basisSign!=Minus) result += calcSinhConv(+1,+1,-1,tau,dgamma,sig,rtau,fsign);
    //if (basisSign!=Plus) result += calcSinhConv(-1,-1,+1,tau,dgamma,sig,rtau,fsign);
    // better form, since it also accounts for the numerical divergence region, added FMV, 07/24/03
    Double_t tau1 = 1/(1/tau-dgamma/2) ;
    Double_t tau2 = 1/(1/tau+dgamma/2) ;
    if (basisSign!=Minus) result += 0.5*(calcDecayConv(+1,tau1,sig,rtau,fsign)-calcDecayConv(+1,tau2,sig,rtau,fsign));
          // modified FMV,08/13/03
    if (basisSign!=Plus) result += 0.5*(calcDecayConv(-1,tau2,sig,rtau,fsign)-calcDecayConv(-1,tau1,sig,rtau,fsign));
          // modified FMV,08/13/03
    //cout << "6th form " << "x= " << x << " result= " << result << endl;
    return result;
  }

  // *** 7th form: Convolution with exp(-t/tau)*cosh(dgamma*t/2), used for coshBasis ***
  if (basisType==coshBasis) {
    Double_t dgamma = ((RooAbsReal*)basis().getParameter(2))->getVal();

    if (verboseEval()>2) cout << "RooGExpModel::evaluate(" << GetName()
               << ") 7th form = " << dgamma << ", tau = " << tau << endl;
    Double_t result(0);
    //if (basisSign!=Minus) result += calcCoshConv(+1,tau,dgamma,sig,rtau,fsign);
    //if (basisSign!=Plus) result += calcCoshConv(-1,tau,dgamma,sig,rtau,fsign);
    // better form, since it also accounts for the numerical divergence region, added FMV, 07/24/03
    Double_t tau1 = 1/(1/tau-dgamma/2) ;
    Double_t tau2 = 1/(1/tau+dgamma/2) ;
    if (basisSign!=Minus) result += 0.5*(calcDecayConv(+1,tau1,sig,rtau,fsign)+calcDecayConv(+1,tau2,sig,rtau,fsign));
          // modified FMV,08/13/03
    if (basisSign!=Plus) result += 0.5*(calcDecayConv(-1,tau1,sig,rtau,fsign)+calcDecayConv(-1,tau2,sig,rtau,fsign));
          // modified FMV,08/13/03
    //cout << "7th form " << "x= " << x << " result= " << result << endl;
    return result;
  }
  R__ASSERT(0) ;
  return 0 ;
  }


////////////////////////////////////////////////////////////////////////////////

std::complex<Double_t> RooGExpModel::calcSinConv(Double_t sign, Double_t sig, Double_t tau, Double_t omega, Double_t rtau, Double_t fsign) const
{
  static Double_t root2(sqrt(2.)) ;

  Double_t s1= -sign*(x - _mean*_meanSF)/tau;
  //Double_t s1= x/tau;
  Double_t c1= sig/(root2*tau);
  Double_t u1= s1/(2*c1);
  Double_t s2= (x - _mean*_meanSF)/rtau;
  Double_t c2= sig/(root2*rtau);
  Double_t u2= fsign*s2/(2*c2) ;
  //Double_t u2= s2/(2*c2) ;

  std::complex<Double_t> eins(1,0);
  std::complex<Double_t> k(1/tau,sign*omega);
  //return (evalCerf(-sign*omega*tau,u1,c1)+evalCerf(0,u2,c2)*fsign*sign) / (eins + k*fsign*sign*rtau) ;

  return (evalCerf(-sign*omega*tau,u1,c1)+std::complex<Double_t>(evalCerfRe(u2,c2),0)*fsign*sign) / (eins + k*fsign*sign*rtau) ;
  // equivalent form, added FMV, 07/24/03
  //return (evalCerf(-sign*omega*tau,-sign*u1,c1)+evalCerf(0,fsign*u2,c2)*fsign*sign) / (eins + k*fsign*sign*rtau) ;
}

// added FMV,08/18/03

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::calcSinConv(Double_t sign, Double_t sig, Double_t tau, Double_t rtau, Double_t fsign) const
{
  static Double_t root2(sqrt(2.)) ;

  Double_t s1= -sign*(x - _mean*_meanSF)/tau;
  //Double_t s1= x/tau;
  Double_t c1= sig/(root2*tau);
  Double_t u1= s1/(2*c1);
  Double_t s2= (x - _mean*_meanSF)/rtau;
  Double_t c2= sig/(root2*rtau);
  Double_t u2= fsign*s2/(2*c2) ;
  //Double_t u2= s2/(2*c2) ;

  Double_t eins(1);
  Double_t k(1/tau);
  return (evalCerfRe(u1,c1)+evalCerfRe(u2,c2)*fsign*sign) / (eins + k*fsign*sign*rtau) ;
  // equivalent form, added FMV, 07/24/03
  //return (evalCerf(-sign*u1,c1).real()+evalCerf(fsign*u2,c2).real()*fsign*sign) / (eins + k*fsign*sign*rtau) ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::calcDecayConv(Double_t sign, Double_t tau, Double_t sig, Double_t rtau, Double_t fsign) const
// modified FMV,08/13/03
{
  static Double_t root2(sqrt(2.)) ;
  static Double_t root2pi(sqrt(2*atan2(0.,-1.))) ;
  static Double_t rootpi(sqrt(atan2(0.,-1.)));

  // Process flip status
  Double_t xp(x - _mean*_meanSF) ;
  //if (_flip) {
  //  xp   *= -1 ;
  //  sign *= -1 ;
  //}
  xp *= fsign ;    // modified FMV,08/13/03
  sign *= fsign ;  // modified FMV,08/13/03

  Double_t cFly;
  if ((sign<0)&&(fabs(tau-rtau)<tau/260)) {

    Double_t MeanTau=0.5*(tau+rtau);
    if (fabs(xp/MeanTau)>300) {
      return 0 ;
    }

    cFly=1./(MeanTau*MeanTau*root2pi) *
      exp(-(-xp/MeanTau-sig*sig/(2*MeanTau*MeanTau)))
      *(sig*exp(-1/(2*sig*sig)*TMath::Power((sig*sig/MeanTau+xp),2))
   -(sig*sig/MeanTau+xp)*(rootpi/root2)*RooMath::erfc(sig/(root2*MeanTau)+xp/(root2*sig)));

    if(_nlo) {
      Double_t epsilon=0.5*(tau-rtau);
      Double_t a=sig/(root2*MeanTau)+xp/(root2*sig);
      cFly += 1./(MeanTau*MeanTau)
   *exp(-(-xp/MeanTau-sig*sig/(2*MeanTau*MeanTau)))
   *0.5/MeanTau*epsilon*epsilon*
   (exp(-a*a)*(sig/MeanTau*root2/rootpi
          -(4*a*sig*sig)/(2*rootpi*MeanTau*MeanTau)
          +(-4/rootpi+8*a*a/rootpi)/6
          *TMath::Power(sig/(root2*MeanTau),3)
          +2/rootpi*(sig*sig/(MeanTau*MeanTau)+xp/MeanTau)*
          (sig/(root2*MeanTau)-a*(sig*sig)/(2*MeanTau*MeanTau))
          +2/rootpi*((3*sig*sig)/(2*MeanTau*MeanTau)+xp/MeanTau+
                0.5*TMath::Power(sig*sig/(MeanTau*MeanTau)+xp/MeanTau,2))*sig/(root2*MeanTau))
    -(2*sig*sig/(MeanTau*MeanTau)+xp/MeanTau+(sig*sig/(MeanTau*MeanTau)+xp/MeanTau)*
      (3*sig*sig/(2*MeanTau*MeanTau)+xp/MeanTau)
      +TMath::Power(sig*sig/(MeanTau*MeanTau)+xp/MeanTau,3)/6)*RooMath::erfc(a));
    }

  } else {

    Double_t expArg1 = sig*sig/(2*tau*tau)-sign*xp/tau ;
    Double_t expArg2 = sig*sig/(2*rtau*rtau)+xp/rtau ;

    Double_t term1, term2 ;
    if (expArg1<300) {
      term1 = exp(expArg1) *RooMath::erfc(sig/(root2*tau)-sign*xp/(root2*sig)) ;
    } else {
      term1 = exp(expArg1+logErfC(sig/(root2*tau)-sign*xp/(root2*sig))) ; ;
    }
    if (expArg2<300) {
      term2 = exp(expArg2) *RooMath::erfc(sig/(root2*rtau)+xp/(root2*sig)) ;
    } else {
      term2 = exp(expArg2+logErfC(sig/(root2*rtau)+xp/(root2*sig))) ; ;
    }

    cFly=(term1+sign*term2)/(2*(tau+sign*rtau));

    // WVE prevent numeric underflows
    if (cFly<1e-100) {
      cFly = 0 ;
    }

    // equivalent form, added FMV, 07/24/03
    //cFly = calcSinConv(sign, sig, tau, rtau, fsign)/(2*tau) ;
  }

  return cFly*2*tau ;
}

/* commented FMV, 07/24/03

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::calcCoshConv(Double_t sign, Double_t tau, Double_t dgamma, Double_t sig, Double_t rtau, Double_t fsign) const
{


  static Double_t root2(sqrt(2.)) ;
  static Double_t root2pi(sqrt(2*atan2(0.,-1.))) ;
  static Double_t rootpi(sqrt(atan2(0.,-1.)));
  Double_t tau1 = 1/(1/tau-dgamma/2);
  Double_t tau2 = 1/(1/tau+dgamma/2);
  Double_t cFly;
  Double_t xp(x);

  //if (_flip) {
  //  xp   *= -1 ;
  //  sign *= -1 ;
  //}
  xp *= fsign ;    // modified FMV,08/13/03
  sign *= fsign ;  // modified FMV,08/13/03

  cFly=tau1*(exp(sig*sig/(2*tau1*tau1)-sign*xp/tau1)
     *RooMath::erfc(sig/(root2*tau1)-sign*xp/(root2*sig))
     +sign*exp(sig*sig/(2*rtau*rtau)+xp/rtau)
     *RooMath::erfc(sig/(root2*rtau)+xp/(root2*sig)))/(2*(tau1+sign*rtau))
    +tau2*(exp(sig*sig/(2*tau2*tau2)-sign*xp/tau2)
     *RooMath::erfc(sig/(root2*tau2)-sign*xp/(root2*sig))
     +sign*exp(sig*sig/(2*rtau*rtau)+xp/rtau)
     *RooMath::erfc(sig/(root2*rtau)+xp/(root2*sig)))/(2*(tau2+sign*rtau));;
  return cFly;
}
*/

/* commented FMV, 07/24/03

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::calcSinhConv(Double_t sign, Double_t sign1, Double_t sign2, Double_t tau, Double_t dgamma, Double_t sig, Double_t rtau, Double_t fsign) const
{
  static Double_t root2(sqrt(2.)) ;
  static Double_t root2pi(sqrt(2*atan2(0.,-1.))) ;
  static Double_t rootpi(sqrt(atan2(0.,-1.)));
  Double_t tau1 = 1/(1/tau-dgamma/2);
  Double_t tau2 = 1/(1/tau+dgamma/2);
  Double_t cFly;
  Double_t xp(x);

  //if (_flip) {
  //  xp   *= -1 ;
  //  sign1 *= -1 ;
  //  sign2 *= -1 ;
  //}
  xp *= fsign ;    // modified FMV,08/13/03
  sign1 *= fsign ;  // modified FMV,08/13/03
  sign2 *= fsign ;  // modified FMV,08/13/03

  cFly=sign1*tau1*(exp(sig*sig/(2*tau1*tau1)-sign*xp/tau1)
     *RooMath::erfc(sig/(root2*tau1)-sign*xp/(root2*sig))
     +sign*exp(sig*sig/(2*rtau*rtau)+xp/rtau)
     *RooMath::erfc(sig/(root2*rtau)+xp/(root2*sig)))/(2*(tau1+sign*rtau))
    +sign2*tau2*(exp(sig*sig/(2*tau2*tau2)-sign*xp/tau2)
     *RooMath::erfc(sig/(root2*tau2)-sign*xp/(root2*sig))
     +sign*exp(sig*sig/(2*rtau*rtau)+xp/rtau)
     *RooMath::erfc(sig/(root2*rtau)+xp/(root2*sig)))/(2*(tau2+sign*rtau));;
  return cFly;
}
*/

////////////////////////////////////////////////////////////////////////////////

Int_t RooGExpModel::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  switch(_basisCode) {

  // Analytical integration capability of raw PDF
  case noBasis:
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
  case sinhBasisPlus:
  case sinhBasisMinus:
  case sinhBasisSum:
  case coshBasisPlus:
  case coshBasisMinus:
  case coshBasisSum:

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

Double_t RooGExpModel::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static Double_t root2 = sqrt(2.) ;
//   static Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  Double_t ssfInt(1.0) ;

  // Code must be 1 or 2
  R__ASSERT(code==1||code==2) ;
  if (code==2) {
    ssfInt = (ssf.max(rangeName)-ssf.min(rangeName)) ;
  }

  BasisType basisType = (BasisType)( (_basisCode == 0) ? 0 : (_basisCode/10) + 1 );
  BasisSign basisSign = (BasisSign)( _basisCode - 10*(basisType-1) - 2 ) ;

  Double_t tau = (_basisCode!=noBasis)?((RooAbsReal*)basis().getParameter(1))->getVal():0 ;

  // added FMV, 07/24/03
  if (basisType == coshBasis && _basisCode!=noBasis ) {
     Double_t dGamma = ((RooAbsReal*)basis().getParameter(2))->getVal();
     if (dGamma==0) basisType = expBasis;
  }
  Double_t fsign = _flip?-1:1 ;
  Double_t sig = sigma*ssf ;
  Double_t rtau = rlife*rsf ;

  // *** 1st form????
  if (basisType==none || ((basisType==expBasis || basisType==cosBasis) && tau==0.)) {
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName() << ") 1st form" << endl ;

    //Double_t result = 1.0 ; // WVE inferred from limit(tau->0) of cosBasisNorm
    // finite+asymtotic normalization, added FMV, 07/24/03
    Double_t xpmin = (x.min(rangeName) - _mean*_meanSF)/rtau ;
    Double_t xpmax = (x.max(rangeName) - _mean*_meanSF)/rtau ;
    Double_t c = sig/(root2*rtau) ;
    Double_t umin = xpmin/(2*c) ;
    Double_t umax = xpmax/(2*c) ;
    Double_t result ;
    if (_asympInt) {
      result = 1.0 ;
    } else {
      result = 0.5*evalCerfInt(-fsign,rtau,fsign*umin,fsign*umax,c)/rtau ;  //WVEFIX add 1/rtau
    }

    if (_basisCode!=0 && basisSign==Both) result *= 2 ;
    //cout << "Integral 1st form " << " result= " << result*ssfInt << endl;
    return result*ssfInt ;
  }

  Double_t omega = (basisType!=expBasis) ?((RooAbsReal*)basis().getParameter(2))->getVal() : 0 ;

  // *** 2nd form: unity, used for sinBasis and cosBasis with tau=0 (PDF is zero) ***
  //if (tau==0&&omega!=0) {
  if (tau==0) {  // modified, FMV 07/24/03
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName() << ") 2nd form" << endl ;
    return 0. ;
  }

  // *** 3rd form: Convolution with exp(-t/tau), used for expBasis and cosBasis(omega=0) ***
  if (basisType==expBasis || (basisType==cosBasis && omega==0.)) {
    //Double_t result = 2*tau ;
    //if (basisSign==Both) result *= 2 ;
    // finite+asymtotic normalization, added FMV, 07/24/03
    Double_t result(0.);
    if (basisSign!=Minus) result += calcSinConvNorm(+1,tau,sig,rtau,fsign,rangeName);
    if (basisSign!=Plus) result += calcSinConvNorm(-1,tau,sig,rtau,fsign,rangeName);
    //cout << "Integral 3rd form " << " result= " << result*ssfInt << endl;
    return result*ssfInt ;
  }

  // *** 4th form: Convolution with exp(-t/tau)*sin(omega*t), used for sinBasis(omega<>0,tau<>0) ***
  Double_t wt = omega * tau ;
  if (basisType==sinBasis) {
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName() << ") 4th form omega = "
              << omega << ", tau = " << tau << endl ;
    //cout << "sin integral" << endl;
    Double_t result(0) ;
    if (wt==0) return result ;
    //if (basisSign!=Minus) result += calcSinConvNorm(+1,tau,omega).imag() ;
    //if (basisSign!=Plus) result += calcSinConvNorm(-1,tau,omega).imag() ;
    // finite+asymtotic normalization, added FMV, 07/24/03
    if (basisSign!=Minus) result += -1*calcSinConvNorm(+1,tau,omega,sig,rtau,fsign,rangeName).imag();
    if (basisSign!=Plus) result += -1*calcSinConvNorm(-1,tau,omega,sig,rtau,fsign,rangeName).imag();
    //cout << "Integral 4th form " << " result= " << result*ssfInt << endl;
    return result*ssfInt ;
  }

  // *** 5th form: Convolution with exp(-t/tau)*cos(omega*t), used for cosBasis(omega<>0) ***
  if (basisType==cosBasis) {
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName()
              << ") 5th form omega = " << omega << ", tau = " << tau << endl ;
    //cout << "cos integral" << endl;
    Double_t result(0) ;
    //if (basisSign!=Minus) result += calcSinConvNorm(+1,tau,omega).real() ;
    //if (basisSign!=Plus) result += calcSinConvNorm(-1,tau,omega).real() ;
    // finite+asymtotic normalization, added FMV, 07/24/03
    if (basisSign!=Minus) result += calcSinConvNorm(+1,tau,omega,sig,rtau,fsign,rangeName).real();
    if (basisSign!=Plus) result += calcSinConvNorm(-1,tau,omega,sig,rtau,fsign,rangeName).real();
    //cout << "Integral 5th form " << " result= " << result*ssfInt << endl;
    return result*ssfInt ;
  }

  Double_t dgamma = ((basisType==coshBasis)||(basisType==sinhBasis))?((RooAbsReal*)basis().getParameter(2))->getVal():0 ;

  // *** 6th form: Convolution with exp(-t/tau)*sinh(dgamma*t/2), used for sinhBasis ***
  if (basisType==sinhBasis) {
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName()
              << ") 6th form dgamma = " << dgamma << ", tau = " << tau << endl ;
    Double_t tau1 = 1/(1/tau-dgamma/2);
    Double_t tau2 = 1/(1/tau+dgamma/2);
    //cout << "sinh integral" << endl;
    Double_t result(0) ;
    //if (basisSign!=Minus) result += tau1-tau2 ;
    //if (basisSign!=Plus) result += tau2-tau1 ;
    // finite+asymtotic normalization, added FMV, 07/24/03
    if (basisSign!=Minus) result += 0.5*(calcSinConvNorm(+1,tau1,sig,rtau,fsign,rangeName)-
                calcSinConvNorm(+1,tau2,sig,rtau,fsign,rangeName));
    if (basisSign!=Plus) result += 0.5*(calcSinConvNorm(-1,tau2,sig,rtau,fsign,rangeName)-
               calcSinConvNorm(-1,tau1,sig,rtau,fsign,rangeName));
    //cout << "Integral 6th form " << " result= " << result*ssfInt << endl;
    return result;
    }

  // ** 7th form: Convolution with exp(-t/tau)*cosh(dgamma*t/2), used for coshBasis ***
  if (basisType==coshBasis) {
    if (verboseEval()>0) cout << "RooGExpModel::analyticalIntegral(" << GetName()
              << ") 6th form dgamma = " << dgamma << ", tau = " << tau << endl ;
    //cout << "cosh integral" << endl;
    Double_t tau1 = 1/(1/tau-dgamma/2);
    Double_t tau2 = 1/(1/tau+dgamma/2);
    //Double_t result = (tau1+tau2) ;
    //if (basisSign==Both) result *= 2 ;
    // finite+asymtotic normalization, added FMV, 07/24/03
    Double_t result(0);
    if (basisSign!=Minus) result += 0.5*(calcSinConvNorm(+1,tau1,sig,rtau,fsign,rangeName)+
                calcSinConvNorm(+1,tau2,sig,rtau,fsign,rangeName));
    if (basisSign!=Plus) result += 0.5*(calcSinConvNorm(-1,tau1,sig,rtau,fsign,rangeName)+
               calcSinConvNorm(-1,tau2,sig,rtau,fsign,rangeName));
    //cout << "Integral 7th form " << " result= " << result*ssfInt << endl;
    return result;

    }

  R__ASSERT(0) ;
  return 1 ;
}

// modified FMV, 07/24/03. Finite+asymtotic normalization

////////////////////////////////////////////////////////////////////////////////
///  old code (asymptotic normalization only)
///  std::complex<Double_t> z(1/tau,sign*omega);
///  return z*2/(omega*omega+1/(tau*tau));

std::complex<Double_t> RooGExpModel::calcSinConvNorm(Double_t sign, Double_t tau, Double_t omega,
                Double_t sig, Double_t rtau, Double_t fsign, const char* rangeName) const
{
  static Double_t root2(sqrt(2.)) ;

  Double_t smin1= (x.min(rangeName) - _mean*_meanSF)/tau;
  Double_t smax1= (x.max(rangeName) - _mean*_meanSF)/tau;
  Double_t c1= sig/(root2*tau);
  Double_t umin1= smin1/(2*c1);
  Double_t umax1= smax1/(2*c1);
  Double_t smin2= (x.min(rangeName) - _mean*_meanSF)/rtau;
  Double_t smax2= (x.max(rangeName) - _mean*_meanSF)/rtau;
  Double_t c2= sig/(root2*rtau);
  Double_t umin2= smin2/(2*c2) ;
  Double_t umax2= smax2/(2*c2) ;

  std::complex<Double_t> eins(1,0);
  std::complex<Double_t> k(1/tau,sign*omega);
  std::complex<Double_t> term1 = evalCerfInt(sign,-sign*omega*tau, tau, -sign*umin1, -sign*umax1, c1);
  //std::complex<Double_t> term2 = evalCerfInt(-fsign,0., rtau, fsign*umin2, fsign*umax2, c2)*std::complex<Double_t>(fsign*sign,0);
  std::complex<Double_t> term2 = std::complex<Double_t>(evalCerfInt(-fsign, rtau, fsign*umin2, fsign*umax2, c2)*fsign*sign,0);
  return (term1+term2)/(eins + k*fsign*sign*rtau) ;
}

// added FMV, 08/17/03

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::calcSinConvNorm(Double_t sign, Double_t tau, Double_t sig, Double_t rtau, Double_t fsign, const char* rangeName) const
{
  static Double_t root2(sqrt(2.)) ;

  Double_t smin1= (x.min(rangeName) - _mean*_meanSF)/tau;
  Double_t smax1= (x.max(rangeName) - _mean*_meanSF)/tau;
  Double_t c1= sig/(root2*tau);
  Double_t umin1= smin1/(2*c1);
  Double_t umax1= smax1/(2*c1);
  Double_t smin2= (x.min(rangeName) - _mean*_meanSF)/rtau;
  Double_t smax2= (x.max(rangeName) - _mean*_meanSF)/rtau;
  Double_t c2= sig/(root2*rtau);
  Double_t umin2= smin2/(2*c2) ;
  Double_t umax2= smax2/(2*c2) ;

  Double_t eins(1);
  Double_t k(1/tau);
  Double_t term1 = evalCerfInt(sign, tau, -sign*umin1, -sign*umax1, c1);
  Double_t term2 = evalCerfInt(-fsign, rtau, fsign*umin2, fsign*umax2, c2)*fsign*sign;

  // WVE Handle 0/0 numeric divergence
  if (fabs(tau-rtau)<1e-10 && fabs(term1+term2)<1e-10) {
    cout << "epsilon method" << endl ;
    static Double_t epsilon = 1e-4 ;
    return calcSinConvNorm(sign,tau+epsilon,sig,rtau-epsilon,fsign,rangeName) ;
  }
  return (term1+term2)/(eins + k*fsign*sign*rtau) ;
}

// added FMV, 07/24/03
////////////////////////////////////////////////////////////////////////////////

std::complex<Double_t> RooGExpModel::evalCerfInt(Double_t sign, Double_t wt, Double_t tau, Double_t umin, Double_t umax, Double_t c) const
{
  std::complex<Double_t> diff;
  if (_asympInt) {
    diff = std::complex<Double_t>(2,0) ;
  } else {
    diff = std::complex<Double_t>(sign,0.)*(evalCerf(wt,umin,c) - evalCerf(wt,umax,c) + RooMath::erf(umin) - RooMath::erf(umax));
  }
  return std::complex<Double_t>(tau/(1.+wt*wt),0)*std::complex<Double_t>(1,wt)*diff;
}
// added FMV, 08/17/03. Modified FMV, 08/30/03

////////////////////////////////////////////////////////////////////////////////

Double_t RooGExpModel::evalCerfInt(Double_t sign, Double_t tau, Double_t umin, Double_t umax, Double_t c) const
{
  Double_t diff;
  if (_asympInt) {
    diff = 2. ;
  } else {
    if ((umin<-8 && umax>8)||(umax<-8 && umin>8)) {
      // If integral is over >8 sigma, approximate with full integral
      diff = 2. ;
    } else {
      diff = sign*(evalCerfRe(umin,c) - evalCerfRe(umax,c) + RooMath::erf(umin) - RooMath::erf(umax));
    }
  }
  return tau*diff;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooGExpModel::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooGExpModel::generateEvent(Int_t code)
{
  R__ASSERT(code==1) ;
  Double_t xgen ;
  while (true) {
    Double_t xgau = RooRandom::randomGenerator()->Gaus(0,(sigma*ssf));
    Double_t xexp = RooRandom::uniform();
    if (!_flip)
      xgen = xgau + (rlife*rsf)*log(xexp);  // modified, FMV 08/13/03
    else
      xgen = xgau - (rlife*rsf)*log(xexp);

    if (xgen < (x.max() - _mean*_meanSF) && xgen > (x.min() - _mean*_meanSF)) {
      x = xgen + _mean*_meanSF;
      return ;
    }
  }
}
