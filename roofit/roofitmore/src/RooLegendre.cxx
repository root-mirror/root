/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   Nikhef & VU, Gerhard.Raven@nikhef.nl
 *                                                                           *
 * Copyright (c) 2010, Nikhef & VU. All rights reserved.
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooLegendre
    \ingroup Roofit

    Compute the associated Legendre polynomials using ROOT::Math::assoc_legendre().

    Since the Legendre polynomials have a value range of [-1, 1], these cannot be implemented as a PDF.
    They can be used in sums, though, for example using a RooRealSumFunc of RooLegendre plus an offset.
**/

#include "RooLegendre.h"
#include "RooBatchCompute.h"
#include "RooAbsReal.h"

#include "Math/SpecFunc.h"
#include "TMath.h"

#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

ClassImp(RooLegendre);

////////////////////////////////////////////////////////////////////////////////

namespace {
    inline double a(int p, int l, int m) {
        double r = TMath::Factorial(l+m)/TMath::Factorial(m+p)/TMath::Factorial(p)/TMath::Factorial(l-m-2*p);
        r /= pow(2.,m+2*p);
        return p%2==0 ? r : -r ;
    }

    void throwIfNoMathMore() {
#ifndef R__HAS_MATHMORE
      throw std::runtime_error("RooLegendre needs functions from MathMore. It is not available in this root build.");
#endif
    }

    void checkCoeffs(int m1, int l1, int m2, int l2) {
      if (m1 < 0 || m2 < 0) {
        throw std::invalid_argument("RooLegendre: m coefficients need to be >= 0.");
      }
      if (l1 < m1 || l2 < m2) {
        throw std::invalid_argument("RooLegendre: m coefficients need to be smaller than corresponding l.");
      }
    }
}

////////////////////////////////////////////////////////////////////////////////

RooLegendre::RooLegendre() :
  _l1(1),_m1(1),_l2(0),_m2(0)
{
  throwIfNoMathMore();
}

////////////////////////////////////////////////////////////////////////////////
///TODO: for now, we assume that ctheta has a range [-1,1]
/// should map the ctheta range onto this interval, and adjust integrals...

RooLegendre::RooLegendre(const char* name, const char* title, RooAbsReal& ctheta, int l, int m)
 : RooAbsReal(name, title)
 , _ctheta("ctheta", "ctheta", this, ctheta)
 , _l1(l),_m1(m),_l2(0),_m2(0)
{
  checkCoeffs(_m1, _l1, _m2, _l2);

  throwIfNoMathMore();
}

////////////////////////////////////////////////////////////////////////////////

RooLegendre::RooLegendre(const char* name, const char* title, RooAbsReal& ctheta, int l1, int m1, int l2, int m2)
 : RooAbsReal(name, title)
 , _ctheta("ctheta", "ctheta", this, ctheta)
 , _l1(l1),_m1(m1),_l2(l2),_m2(m2)
{
  checkCoeffs(_m1, _l1, _m2, _l2);

  throwIfNoMathMore();
}

////////////////////////////////////////////////////////////////////////////////

RooLegendre::RooLegendre(const RooLegendre& other, const char* name)
    : RooAbsReal(other, name)
    , _ctheta("ctheta", this, other._ctheta)
    , _l1(other._l1), _m1(other._m1)
    , _l2(other._l2), _m2(other._m2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Note: P_0^0 = 1, so P_l^m = P_l^m P_0^0

Double_t RooLegendre::evaluate() const
{
#ifdef R__HAS_MATHMORE
  double r = 1;
  double ctheta = std::max(-1., std::min((double)_ctheta, +1.));
  if (_l1!=0||_m1!=0) r *= ROOT::Math::assoc_legendre(_l1,_m1,ctheta);
  if (_l2!=0||_m2!=0) r *= ROOT::Math::assoc_legendre(_l2,_m2,ctheta);
  if ((_m1+_m2)%2==1) r = -r;
  return r;
#else
  throwIfNoMathMore();
  return 0.;
#endif
}


////////////////////////////////////////////////////////////////////////////////

namespace {
//Author: Emmanouil Michalainas, CERN 26 August 2019

void compute(	size_t batchSize, const int l1, const int m1, const int l2, const int m2,
              double * __restrict output,
              double const * __restrict TH)
{
#ifdef R__HAS_MATHMORE
  double legendre1=1.0, legendreMinus1=1.0;
  if (l1+m1 > 0) {
    legendre1      = ROOT::Math::internal::legendre(l1,m1,1.0);
    legendreMinus1 = ROOT::Math::internal::legendre(l1,m1,-1.0); 
  }
  if (l2+m2 > 0) {
    legendre1      *= ROOT::Math::internal::legendre(l2,m2,1.0); 
    legendreMinus1 *= ROOT::Math::internal::legendre(l2,m2,-1.0); 
  }
  
  for (size_t i=0; i<batchSize; i++) {
    if (TH[i] <= -1.0) {
      output[i] = legendreMinus1;
    } else if (TH[i] >= 1.0) {
      output[i] = legendre1;
    }
    else {
      output[i] = 1.0;
      if (l1+m1 > 0) {
        output[i] *= ROOT::Math::internal::legendre(l1,m1,TH[i]);
      }
      if (l2+m2 > 0) {
        output[i] *= ROOT::Math::internal::legendre(l2,m2,TH[i]);
      }
    }
  }

#else
  (void) batchSize, (void) l1, (void)m1, (void)l2, (void)m2, (void)output, (void)TH;
  throwIfNoMathMore();
#endif
}
};

RooSpan<double> RooLegendre::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {
  RooSpan<const double> cthetaData = _ctheta->getValues(evalData, normSet);
  size_t batchSize = cthetaData.size();
  auto output = evalData.makeBatch(this, batchSize);
  compute(batchSize, _l1, _m1, _l2, _m2, output.data(), cthetaData.data());
  return output;
}


////////////////////////////////////////////////////////////////////////////////

namespace {
  Bool_t fullRange(const RooRealProxy& x ,const char* range)
  {
    return range == 0 || strlen(range) == 0
        ? std::fabs(x.min() + 1.) < 1.e-8 && std::fabs(x.max() - 1.) < 1.e-8
        : std::fabs(x.min(range) + 1.) < 1.e-8 && std::fabs(x.max(range) - 1.) < 1.e-8;
  }
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooLegendre::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // don't support indefinite integrals...
  if (fullRange(_ctheta,rangeName) && matchArgs(allVars, analVars, _ctheta)) return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// this was verified to match mathematica for
/// l1 in [0,2], m1 in [0,l1], l2 in [l1,4], m2 in [0,l2]

Double_t RooLegendre::analyticalIntegral(Int_t code, const char* ) const
{
  R__ASSERT(code==1) ;
  if ( _m1==_m2 )                 return ( _l1 == _l2) ?  TMath::Factorial(_l1+_m2)/TMath::Factorial(_l1-_m1)*double(2)/(2*_l1+1) : 0.;
  if ( (_l1+_l2-_m1-_m2)%2 != 0 ) return 0; // these combinations are odd under x -> -x

  // from B.R. Wong, "On the overlap integral of associated Legendre Polynomials" 1998 J. Phys. A: Math. Gen. 31 1101
  // TODO: update to the result of
  //       H. A. Mavromatis
  //       "A single-sum expression for the overlap integral of two associated Legendre polynomials"
  //       1999 J. Phys. A: Math. Gen. 32 2601
  //       http://iopscience.iop.org/0305-4470/32/13/011/pdf/0305-4470_32_13_011.pdf
  //       For that we need Wigner 3-j, which Lorenzo has added for Root 5.28... (note: check Condon-Shortly convention in this paper!)
  double r=0;
  for (int p1=0; 2*p1 <= _l1-_m1 ;++p1) {
    double a1 = a(p1,_l1,_m1);
    for (int p2=0; 2*p2 <= _l2-_m2 ; ++p2) {
       double a2 = a(p2,_l2,_m2);
       r+= a1*a2*TMath::Gamma( double(_l1+_l2-_m1-_m2-2*p1-2*p2+1)/2 )*TMath::Gamma( double(_m1+_m2+2*p1+2*p2+2)/2 );
    }
  }
  r /= TMath::Gamma( double(_l1+_l2+3)/2 );

  if ((_m1+_m2)%2==1) r = -r;
  return r;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooLegendre::getMaxVal( const RooArgSet& /*vars*/) const {
    if (_m1==0&&_m2==0) return 1;
    // does anyone know the analytical expression for the  max values in case m!=0??
    if (_l1<3&&_l2<3) return 1;
    return 0;
}

namespace {
    inline double maxSingle(int i, int j) {
        R__ASSERT(j<=i);
        //   x0 : 1 (ordinary Legendre)
        if (j==0) return 1;
        R__ASSERT(i<3);
        //   11: 1
        if (i<2) return 1;
        //   21: 3   22: 3
        static const double m2[3] = { 3,3 };
        return m2[j-1];
    }
}
Double_t RooLegendre::maxVal( Int_t /*code*/) const {
    return maxSingle(_l1,_m1)*maxSingle(_l2,_m2);
}
