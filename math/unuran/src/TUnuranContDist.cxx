// @(#)root/unuran:$Id$
// Authors: L. Moneta, J. Leydold Wed Feb 28 2007

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class TUnuranContDist

#include "TUnuranContDist.h"
#include "Math/RichardsonDerivator.h"
#include "Math/WrappedTF1.h"

#include "Math/Integrator.h"

#include "TF1.h"
#include <cassert>
#include <cmath>

ClassImp(TUnuranContDist);

TUnuranContDist::TUnuranContDist(const ROOT::Math::IGenFunction *pdf, const ROOT::Math::IGenFunction *dpdf,
                                 const ROOT::Math::IGenFunction *cdf, bool isLogPdf, bool copyFunc)
   : fPdf(pdf), fDPdf(dpdf), fCdf(cdf), fXmin(1.), fXmax(-1.), fMode(0), fArea(0), fIsLogPdf(isLogPdf), fHasDomain(0),
     fHasMode(0), fHasArea(0), fOwnFunc(copyFunc)
{
   // Constructor from generic function interfaces
   // manage the functions and clone them if flag copyFunc is true
   if (fOwnFunc) {
      if (fPdf) 
         fPdf = fPdf->Clone();
      if (fDPdf)
         fDPdf = fDPdf->Clone();
      if (fCdf)
         fCdf = fCdf->Clone();
   }
}

TUnuranContDist::TUnuranContDist (const ROOT::Math::IGenFunction & pdf, const ROOT::Math::IGenFunction * deriv, bool isLogPdf, bool copyFunc  ) :
   TUnuranContDist(&pdf,deriv, nullptr, isLogPdf, copyFunc)
{}

TUnuranContDist::TUnuranContDist (TF1 * pdf, TF1 * deriv, TF1 * cdf, bool isLogPdf  ) :
   fPdf(  (pdf) ? new ROOT::Math::WrappedTF1 ( *pdf) : nullptr ),
   fDPdf( (deriv) ?  new ROOT::Math::WrappedTF1 ( *deriv) : nullptr ),
   fCdf( (cdf) ?  new ROOT::Math::WrappedTF1 ( *cdf) : nullptr),
   fXmin(1.),
   fXmax(-1.),
   fMode(0),
   fArea(0),
   fIsLogPdf(isLogPdf),
   fHasDomain(0),
   fHasMode(0),
   fHasArea(0),
   fOwnFunc(true)
{
   // Constructor from a TF1 objects
   // function pointers are managed by class
}

TUnuranContDist::TUnuranContDist (TF1 * pdf, TF1 * deriv, bool isLogPdf  ) :
   TUnuranContDist(pdf,deriv, nullptr, isLogPdf)
   {}

TUnuranContDist::TUnuranContDist(const TUnuranContDist & rhs) :
   TUnuranBaseDist(),
   fPdf(nullptr),
   fDPdf(nullptr),
   fCdf(nullptr)
{
   // Implementation of copy constructor
   operator=(rhs);
}

TUnuranContDist & TUnuranContDist::operator = (const TUnuranContDist &rhs)
{
   // Implementation of assignment operator.
   if (this == &rhs) return *this;  // time saving self-test
   fXmin  = rhs.fXmin;
   fXmax  = rhs.fXmax;
   fMode  = rhs.fMode;
   fArea  = rhs.fArea;
   fIsLogPdf  = rhs.fIsLogPdf;
   fHasDomain = rhs.fHasDomain;
   fHasMode   = rhs.fHasMode;
   fHasArea   = rhs.fHasArea;
   fOwnFunc   = rhs.fOwnFunc;
   if (!fOwnFunc) {
      fPdf   = rhs.fPdf;
      fDPdf  = rhs.fDPdf;
      fCdf   = rhs.fCdf;
   }
   else {
      if (fPdf) delete fPdf;
      if (fDPdf) delete fDPdf;
      if (fCdf) delete fCdf;
      fPdf  = (rhs.fPdf)  ? rhs.fPdf->Clone()  : nullptr;
      fDPdf = (rhs.fDPdf) ? rhs.fDPdf->Clone() : nullptr;
      fCdf  = (rhs.fCdf)  ? rhs.fCdf->Clone()  : nullptr;
   }

   return *this;
}

TUnuranContDist::~TUnuranContDist() {
   // destructor implementation
   if (fOwnFunc) {
      if (fPdf) delete fPdf;
      if (fDPdf) delete fDPdf;
      if (fCdf) delete fCdf;
   }
}

void TUnuranContDist::SetCdf(const ROOT::Math::IGenFunction & cdf) {
   //  set cdf distribution using a generic function interface
   fCdf = (fOwnFunc) ? cdf.Clone() : &cdf;
}


void TUnuranContDist::SetCdf(TF1 *  cdf) {
   // set cumulative distribution function from a TF1
   if (!fOwnFunc) {
      // need to manage all functions now
      if (fPdf) fPdf = fPdf->Clone();
      if (fDPdf) fDPdf->Clone();
   }
   else
      if (fCdf) delete fCdf;

   fCdf = (cdf) ? new ROOT::Math::WrappedTF1 ( *cdf) : nullptr;
   fOwnFunc = true;
}

double TUnuranContDist::Pdf ( double x) const {
   // evaluate the pdf of the distribution. Return NaN if pdf is not available
   return (fPdf) ? (*fPdf)(x) : TMath::QuietNaN();
}

double TUnuranContDist::DPdf( double x) const {
   // evaluate the derivative of the pdf
   // if derivative function is not given is evaluated numerically
   // in case a pdf is available, otherwise a NaN is returned
   if (fDPdf) {
      return (*fDPdf)(x);
   }
   if (!fPdf) return TMath::QuietNaN();
   // do numerical derivation using numerical derivation
   ROOT::Math::RichardsonDerivator rd;
   static double gEps = 0.001;
   double h = ( std::abs(x) > 0 ) ?  gEps * std::abs(x) : gEps;
   assert (fPdf != 0);
   return rd.Derivative1( *fPdf, x, h);
}

double TUnuranContDist::Cdf(double x) const {
   // evaluate the integral (cdf)  on the domain
   if (fCdf) {
      return (*fCdf)(x);
   }
   // do numerical integration
   if (!fPdf) return TMath::QuietNaN();
   ROOT::Math::Integrator ig;
   if (fXmin > fXmax) return ig.Integral( *fPdf );
   else
      return ig.Integral( *fPdf, fXmin, fXmax );

}

