// @(#)root/mathcore:$Id$
// Authors: Rene Brun, Anna Kreshuk, Eddy Offermann, Fons Rademakers   29/07/95

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMath
#define ROOT_TMath

#include "TMathBase.h"

#include "TError.h"
#include <algorithm>
#include <limits>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
///
/// TMath
///
/// Encapsulate most frequently used Math functions.
/// NB. The basic functions Min, Max, Abs and Sign are defined
/// in TMathBase.

namespace TMath {

////////////////////////////////////////////////////////////////////////////////
// Fundamental constants

////////////////////////////////////////////////////////////////////////////////
/// \f[ \pi\f]
constexpr Double_t Pi()
{
   return 3.14159265358979323846;
}

////////////////////////////////////////////////////////////////////////////////
/// \f[ 2\pi\f]
constexpr Double_t TwoPi()
{
   return 2.0 * Pi();
}

////////////////////////////////////////////////////////////////////////////////
/// \f[ \frac{\pi}{2} \f]
constexpr Double_t PiOver2()
{
   return Pi() / 2.0;
}

////////////////////////////////////////////////////////////////////////////////
/// \f[ \frac{\pi}{4} \f]
constexpr Double_t PiOver4()
{
   return Pi() / 4.0;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \frac{1.}{\pi}\f$
constexpr Double_t InvPi()
{
   return 1.0 / Pi();
}

////////////////////////////////////////////////////////////////////////////////
/// Conversion from radian to degree:
/// \f[ \frac{180}{\pi} \f]
constexpr Double_t RadToDeg()
{
   return 180.0 / Pi();
}

////////////////////////////////////////////////////////////////////////////////
/// Conversion from degree to radian:
/// \f[ \frac{\pi}{180} \f]
constexpr Double_t DegToRad()
{
   return Pi() / 180.0;
}

////////////////////////////////////////////////////////////////////////////////
/// \f[ \sqrt{2} \f]
constexpr Double_t Sqrt2()
{
   return 1.4142135623730950488016887242097;
}

////////////////////////////////////////////////////////////////////////////////
/// Base of natural log:
///  \f[ e \f]
constexpr Double_t E()
{
   return 2.71828182845904523536;
}

////////////////////////////////////////////////////////////////////////////////
/// Natural log of 10 (to convert log to ln)
constexpr Double_t Ln10()
{
   return 2.30258509299404568402;
}

////////////////////////////////////////////////////////////////////////////////
/// Base-10 log of e  (to convert ln to log)
constexpr Double_t LogE()
{
   return 0.43429448190325182765;
}

////////////////////////////////////////////////////////////////////////////////
/// Velocity of light in \f$ m s^{-1} \f$
constexpr Double_t C()
{
   return 2.99792458e8;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ cm s^{-1} \f$
constexpr Double_t Ccgs()
{
   return 100.0 * C();
}

////////////////////////////////////////////////////////////////////////////////
/// Speed of light uncertainty.
constexpr Double_t CUncertainty()
{
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Gravitational constant in: \f$ m^{3} kg^{-1} s^{-2} \f$
constexpr Double_t G()
{
   // use 2018 value from NIST  (https://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=G)
   return 6.67430e-11;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ cm^{3} g^{-1} s^{-2} \f$
constexpr Double_t Gcgs()
{
   return G() * 1000.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Gravitational constant uncertainty.
constexpr Double_t GUncertainty()
{
   // use 2018 value from NIST
   return 0.00015e-11;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \frac{G}{\hbar C} \f$ in \f$ (GeV/c^{2})^{-2} \f$
constexpr Double_t GhbarC()
{
   // use new value from NIST (2018)
   return 6.70883e-39;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \frac{G}{\hbar C} \f$ uncertainty.
constexpr Double_t GhbarCUncertainty()
{
   // use new value from NIST (2018)
   return 0.00015e-39;
}

////////////////////////////////////////////////////////////////////////////////
/// Standard acceleration of gravity in \f$ m s^{-2} \f$
constexpr Double_t Gn()
{
   return 9.80665;
}

////////////////////////////////////////////////////////////////////////////////
/// Standard acceleration of gravity uncertainty.
constexpr Double_t GnUncertainty()
{
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Planck's constant in \f$ J s \f$
/// \f[ h \f]
constexpr Double_t H()
{
   return 6.62607015e-34;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ erg s \f$
constexpr Double_t Hcgs()
{
   return 1.0e7 * H();
}

////////////////////////////////////////////////////////////////////////////////
/// Planck's constant uncertainty.
constexpr Double_t HUncertainty()
{
   // Planck constant is exact according to 2019 redefinition 
   // (https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units)
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \hbar \f$ in \f$ J s \f$
/// \f[ \hbar = \frac{h}{2\pi} \f]
constexpr Double_t Hbar()
{
   return 1.054571817e-34;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ erg s \f$
constexpr Double_t Hbarcgs()
{
   return 1.0e7 * Hbar();
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \hbar \f$ uncertainty.
constexpr Double_t HbarUncertainty()
{
   // hbar is an exact constant
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ hc \f$ in \f$ J m \f$
constexpr Double_t HC()
{
   return H() * C();
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ erg cm \f$
constexpr Double_t HCcgs()
{
   return Hcgs() * Ccgs();
}

////////////////////////////////////////////////////////////////////////////////
/// Boltzmann's constant in \f$ J K^{-1} \f$
/// \f[ k \f]
constexpr Double_t K()
{
   return 1.380649e-23;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ erg K^{-1} \f$
constexpr Double_t Kcgs()
{
   return 1.0e7 * K();
}

////////////////////////////////////////////////////////////////////////////////
/// Boltzmann's constant uncertainty.
constexpr Double_t KUncertainty()
{
   // constant is exact according to 2019 redefinition 
   // (https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units)
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stefan-Boltzmann constant in \f$ W m^{-2} K^{-4}\f$
/// \f[ \sigma \f]
constexpr Double_t Sigma()
{
   return 5.670373e-8;
}

////////////////////////////////////////////////////////////////////////////////
/// Stefan-Boltzmann constant uncertainty.
constexpr Double_t SigmaUncertainty()
{
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Avogadro constant (Avogadro's Number) in \f$ mol^{-1} \f$
constexpr Double_t Na()
{
   return 6.02214076e+23;
}

////////////////////////////////////////////////////////////////////////////////
/// Avogadro constant (Avogadro's Number) uncertainty.
constexpr Double_t NaUncertainty()
{
   // constant is exact according to 2019 redefinition 
   // (https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units)
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// [Universal gas constant](http://scienceworld.wolfram.com/physics/UniversalGasConstant.html)
/// (\f$ Na K \f$) in \f$ J K^{-1} mol^{-1} \f$
//
constexpr Double_t R()
{
   return K() * Na();
}

////////////////////////////////////////////////////////////////////////////////
/// Universal gas constant uncertainty.
constexpr Double_t RUncertainty()
{
   return R() * ((KUncertainty() / K()) + (NaUncertainty() / Na()));
}

////////////////////////////////////////////////////////////////////////////////
/// [Molecular weight of dry air 1976 US Standard Atmosphere](http://atmos.nmsu.edu/jsdap/encyclopediawork.html)
/// in \f$ kg kmol^{-1} \f$ or \f$ gm mol^{-1} \f$
constexpr Double_t MWair()
{
   return 28.9644;
}

////////////////////////////////////////////////////////////////////////////////
/// [Dry Air Gas Constant (R / MWair)](http://atmos.nmsu.edu/education_and_outreach/encyclopedia/gas_constant.htm)
/// in \f$ J kg^{-1} K^{-1} \f$
constexpr Double_t Rgair()
{
   return (1000.0 * R()) / MWair();
}

////////////////////////////////////////////////////////////////////////////////
/// Euler-Mascheroni Constant.
constexpr Double_t EulerGamma()
{
   return 0.577215664901532860606512090082402431042;
}

////////////////////////////////////////////////////////////////////////////////
/// Elementary charge in \f$ C \f$ .
constexpr Double_t Qe()
{
   return 1.602176634e-19;
}

////////////////////////////////////////////////////////////////////////////////
/// Elementary charge uncertainty.
constexpr Double_t QeUncertainty()
{
   // constant is exact according to 2019 redefinition 
   // (https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units)
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Mathematical Functions

////////////////////////////////////////////////////////////////////////////////
// Trigonometrical Functions

inline Double_t Sin(Double_t);
inline Double_t Cos(Double_t);
inline Double_t Tan(Double_t);
inline Double_t SinH(Double_t);
inline Double_t CosH(Double_t);
inline Double_t TanH(Double_t);
inline Double_t ASin(Double_t);
inline Double_t ACos(Double_t);
inline Double_t ATan(Double_t);
inline Double_t ATan2(Double_t y, Double_t x);
Double_t ASinH(Double_t);
Double_t ACosH(Double_t);
Double_t ATanH(Double_t);
Double_t Hypot(Double_t x, Double_t y);

////////////////////////////////////////////////////////////////////////////////
// Elementary Functions

inline Double_t Ceil(Double_t x);
inline Int_t CeilNint(Double_t x);
inline Double_t Floor(Double_t x);
inline Int_t FloorNint(Double_t x);
template <typename T>
inline Int_t Nint(T x);

inline Double_t Sq(Double_t x);
inline Double_t Sqrt(Double_t x);
inline Double_t Exp(Double_t x);
inline Double_t Ldexp(Double_t x, Int_t exp);
Double_t Factorial(Int_t i);
inline LongDouble_t Power(LongDouble_t x, LongDouble_t y);
inline LongDouble_t Power(LongDouble_t x, Long64_t y);
inline LongDouble_t Power(Long64_t x, Long64_t y);
inline Double_t Power(Double_t x, Double_t y);
inline Double_t Power(Double_t x, Int_t y);
inline Double_t Log(Double_t x);
Double_t Log2(Double_t x);
inline Double_t Log10(Double_t x);
inline Int_t Finite(Double_t x);
inline Int_t Finite(Float_t x);
inline Bool_t IsNaN(Double_t x);
inline Bool_t IsNaN(Float_t x);

inline Double_t QuietNaN();
inline Double_t SignalingNaN();
inline Double_t Infinity();

template <typename T>
struct Limits {
   inline static T Min();
   inline static T Max();
   inline static T Epsilon();
   };

   // Some integer math
   Long_t   Hypot(Long_t x, Long_t y);     // sqrt(px*px + py*py)

   // Comparing floating points
   inline Bool_t AreEqualAbs(Double_t af, Double_t bf, Double_t epsilon) {
      //return kTRUE if absolute difference between af and bf is less than epsilon
      return TMath::Abs(af-bf) < epsilon ||
             TMath::Abs(af - bf) < Limits<Double_t>::Min(); // handle 0 < 0 case

   }
   inline Bool_t AreEqualRel(Double_t af, Double_t bf, Double_t relPrec) {
      //return kTRUE if relative difference between af and bf is less than relPrec
      return TMath::Abs(af - bf) <= 0.5 * relPrec * (TMath::Abs(af) + TMath::Abs(bf)) ||
             TMath::Abs(af - bf) < Limits<Double_t>::Min(); // handle denormals
   }

   /////////////////////////////////////////////////////////////////////////////
   // Array Algorithms

   // Min, Max of an array
   template <typename T> T MinElement(Long64_t n, const T *a);
   template <typename T> T MaxElement(Long64_t n, const T *a);

   // Locate Min, Max element number in an array
   template <typename T> Long64_t  LocMin(Long64_t n, const T *a);
   template <typename Iterator> Iterator LocMin(Iterator first, Iterator last);
   template <typename T> Long64_t  LocMax(Long64_t n, const T *a);
   template <typename Iterator> Iterator LocMax(Iterator first, Iterator last);

   // Hashing
   ULong_t Hash(const void *txt, Int_t ntxt);
   ULong_t Hash(const char *str);

   void BubbleHigh(Int_t Narr, Double_t *arr1, Int_t *arr2);
   void BubbleLow (Int_t Narr, Double_t *arr1, Int_t *arr2);

   Bool_t   Permute(Int_t n, Int_t *a); // Find permutations

   /////////////////////////////////////////////////////////////////////////////
   // Geometrical Functions

   //Sample quantiles
   void      Quantiles(Int_t n, Int_t nprob, Double_t *x, Double_t *quantiles, Double_t *prob,
                       Bool_t isSorted=kTRUE, Int_t *index = 0, Int_t type=7);

   // IsInside
   template <typename T> Bool_t IsInside(T xp, T yp, Int_t np, T *x, T *y);

   // Calculate the Cross Product of two vectors
   template <typename T> T *Cross(const T v1[3],const T v2[3], T out[3]);

   Float_t   Normalize(Float_t v[3]);  // Normalize a vector
   Double_t  Normalize(Double_t v[3]); // Normalize a vector

   //Calculate the Normalized Cross Product of two vectors
   template <typename T> inline T NormCross(const T v1[3],const T v2[3],T out[3]);

   // Calculate a normal vector of a plane
   template <typename T> T *Normal2Plane(const T v1[3],const T v2[3],const T v3[3], T normal[3]);

   /////////////////////////////////////////////////////////////////////////////
   // Polynomial Functions

   Bool_t    RootsCubic(const Double_t coef[4],Double_t &a, Double_t &b, Double_t &c);

   /////////////////////////////////////////////////////////////////////////////
   // Statistic Functions

   Double_t Binomial(Int_t n,Int_t k);  // Calculate the binomial coefficient n over k
   Double_t BinomialI(Double_t p, Int_t n, Int_t k);
   Double_t BreitWigner(Double_t x, Double_t mean=0, Double_t gamma=1);
   Double_t CauchyDist(Double_t x, Double_t t=0, Double_t s=1);
   Double_t ChisquareQuantile(Double_t p, Double_t ndf);
   Double_t FDist(Double_t F, Double_t N, Double_t M);
   Double_t FDistI(Double_t F, Double_t N, Double_t M);
   Double_t Gaus(Double_t x, Double_t mean=0, Double_t sigma=1, Bool_t norm=kFALSE);
   Double_t KolmogorovProb(Double_t z);
   Double_t KolmogorovTest(Int_t na, const Double_t *a, Int_t nb, const Double_t *b, Option_t *option);
   Double_t Landau(Double_t x, Double_t mpv=0, Double_t sigma=1, Bool_t norm=kFALSE);
   Double_t LandauI(Double_t x);
   Double_t LaplaceDist(Double_t x, Double_t alpha=0, Double_t beta=1);
   Double_t LaplaceDistI(Double_t x, Double_t alpha=0, Double_t beta=1);
   Double_t LogNormal(Double_t x, Double_t sigma, Double_t theta=0, Double_t m=1);
   Double_t NormQuantile(Double_t p);
   Double_t Poisson(Double_t x, Double_t par);
   Double_t PoissonI(Double_t x, Double_t par);
   Double_t Prob(Double_t chi2,Int_t ndf);
   Double_t Student(Double_t T, Double_t ndf);
   Double_t StudentI(Double_t T, Double_t ndf);
   Double_t StudentQuantile(Double_t p, Double_t ndf, Bool_t lower_tail=kTRUE);
   Double_t Vavilov(Double_t x, Double_t kappa, Double_t beta2);
   Double_t VavilovI(Double_t x, Double_t kappa, Double_t beta2);
   Double_t Voigt(Double_t x, Double_t sigma, Double_t lg, Int_t r = 4);

   /////////////////////////////////////////////////////////////////////////////
   // Statistics over arrays

   //Mean, Geometric Mean, Median, RMS(sigma)

   template <typename T> Double_t Mean(Long64_t n, const T *a, const Double_t *w=0);
   template <typename Iterator> Double_t Mean(Iterator first, Iterator last);
   template <typename Iterator, typename WeightIterator> Double_t Mean(Iterator first, Iterator last, WeightIterator wfirst);

   template <typename T> Double_t GeomMean(Long64_t n, const T *a);
   template <typename Iterator> Double_t GeomMean(Iterator first, Iterator last);

   template <typename T> Double_t RMS(Long64_t n, const T *a, const Double_t *w=0);
   template <typename Iterator> Double_t RMS(Iterator first, Iterator last);
   template <typename Iterator, typename WeightIterator> Double_t RMS(Iterator first, Iterator last, WeightIterator wfirst);

   template <typename T> Double_t StdDev(Long64_t n, const T *a, const Double_t * w = 0) { return RMS<T>(n,a,w); }
   template <typename Iterator> Double_t StdDev(Iterator first, Iterator last) { return RMS<Iterator>(first,last); }
   template <typename Iterator, typename WeightIterator> Double_t StdDev(Iterator first, Iterator last, WeightIterator wfirst) { return RMS<Iterator,WeightIterator>(first,last,wfirst); }

   template <typename T> Double_t Median(Long64_t n, const T *a,  const Double_t *w=0, Long64_t *work=0);

   //k-th order statistic
   template <class Element, typename Size> Element KOrdStat(Size n, const Element *a, Size k, Size *work = 0);

   /////////////////////////////////////////////////////////////////////////////
   // Special Functions

   Double_t Beta(Double_t p, Double_t q);
   Double_t BetaCf(Double_t x, Double_t a, Double_t b);
   Double_t BetaDist(Double_t x, Double_t p, Double_t q);
   Double_t BetaDistI(Double_t x, Double_t p, Double_t q);
   Double_t BetaIncomplete(Double_t x, Double_t a, Double_t b);

   // Bessel functions
   Double_t BesselI(Int_t n,Double_t x);  /// integer order modified Bessel function I_n(x)
   Double_t BesselK(Int_t n,Double_t x);  /// integer order modified Bessel function K_n(x)
   Double_t BesselI0(Double_t x);         /// modified Bessel function I_0(x)
   Double_t BesselK0(Double_t x);         /// modified Bessel function K_0(x)
   Double_t BesselI1(Double_t x);         /// modified Bessel function I_1(x)
   Double_t BesselK1(Double_t x);         /// modified Bessel function K_1(x)
   Double_t BesselJ0(Double_t x);         /// Bessel function J0(x) for any real x
   Double_t BesselJ1(Double_t x);         /// Bessel function J1(x) for any real x
   Double_t BesselY0(Double_t x);         /// Bessel function Y0(x) for positive x
   Double_t BesselY1(Double_t x);         /// Bessel function Y1(x) for positive x
   Double_t StruveH0(Double_t x);         /// Struve functions of order 0
   Double_t StruveH1(Double_t x);         /// Struve functions of order 1
   Double_t StruveL0(Double_t x);         /// Modified Struve functions of order 0
   Double_t StruveL1(Double_t x);         /// Modified Struve functions of order 1

   Double_t DiLog(Double_t x);
   Double_t Erf(Double_t x);
   Double_t ErfInverse(Double_t x);
   Double_t Erfc(Double_t x);
   Double_t ErfcInverse(Double_t x);
   Double_t Freq(Double_t x);
   Double_t Gamma(Double_t z);
   Double_t Gamma(Double_t a,Double_t x);
   Double_t GammaDist(Double_t x, Double_t gamma, Double_t mu=0, Double_t beta=1);
   Double_t LnGamma(Double_t z);
}

////////////////////////////////////////////////////////////////////////////////
// Trig and other functions

#include <float.h>

#if defined(R__WIN32) && !defined(__CINT__)
#   ifndef finite
#      define finite _finite
#   endif
#endif
#if defined(R__AIX) || defined(R__SOLARIS_CC50) || \
    defined(R__HPUX11) || defined(R__GLIBC) || \
    (defined(R__MACOSX) )
// math functions are defined inline so we have to include them here
#   include <math.h>
#   ifdef R__SOLARIS_CC50
       extern "C" { int finite(double); }
#   endif
// #   if defined(R__GLIBC) && defined(__STRICT_ANSI__)
// #      ifndef finite
// #         define finite __finite
// #      endif
// #      ifndef isnan
// #         define isnan  __isnan
// #      endif
// #   endif
#else
// don't want to include complete <math.h>
extern "C" {
   extern double sin(double);
   extern double cos(double);
   extern double tan(double);
   extern double sinh(double);
   extern double cosh(double);
   extern double tanh(double);
   extern double asin(double);
   extern double acos(double);
   extern double atan(double);
   extern double atan2(double, double);
   extern double sqrt(double);
   extern double exp(double);
   extern double pow(double, double);
   extern double log(double);
   extern double log10(double);
#ifndef R__WIN32
#   if !defined(finite)
       extern int finite(double);
#   endif
#   if !defined(isnan)
       extern int isnan(double);
#   endif
   extern double ldexp(double, int);
   extern double ceil(double);
   extern double floor(double);
#else
   _CRTIMP double ldexp(double, int);
   _CRTIMP double ceil(double);
   _CRTIMP double floor(double);
#endif
}
#endif

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Sin(Double_t x)
   { return sin(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Cos(Double_t x)
   { return cos(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Tan(Double_t x)
   { return tan(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::SinH(Double_t x)
   { return sinh(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::CosH(Double_t x)
   { return cosh(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::TanH(Double_t x)
   { return tanh(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::ASin(Double_t x)
   { 
     return asin(x);
   }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::ACos(Double_t x)
   { 
     return acos(x);
   }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::ATan(Double_t x)
   { return atan(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::ATan2(Double_t y, Double_t x)
   { if (x != 0) return  atan2(y, x);
     if (y == 0) return  0;
     if (y >  0) return  Pi()/2;
     else        return -Pi()/2;
   }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Sq(Double_t x)
   { return x*x; }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Sqrt(Double_t x)
   { return sqrt(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Ceil(Double_t x)
   { return ceil(x); }

////////////////////////////////////////////////////////////////////////////////
inline Int_t TMath::CeilNint(Double_t x)
   { return TMath::Nint(ceil(x)); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Floor(Double_t x)
   { return floor(x); }

////////////////////////////////////////////////////////////////////////////////
inline Int_t TMath::FloorNint(Double_t x)
   { return TMath::Nint(floor(x)); }

////////////////////////////////////////////////////////////////////////////////
/// Round to nearest integer. Rounds half integers to the nearest even integer.
template<typename T>
inline Int_t TMath::Nint(T x)
{
   int i;
   if (x >= 0) {
      i = int(x + 0.5);
      if ( i & 1 && x + 0.5 == T(i) ) i--;
   } else {
      i = int(x - 0.5);
      if ( i & 1 && x - 0.5 == T(i) ) i++;
   }
   return i;
}

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Exp(Double_t x)
   { return exp(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Ldexp(Double_t x, Int_t exp)
   { return ldexp(x, exp); }

////////////////////////////////////////////////////////////////////////////////
inline LongDouble_t TMath::Power(LongDouble_t x, LongDouble_t y)
   { return std::pow(x,y); }

////////////////////////////////////////////////////////////////////////////////
inline LongDouble_t TMath::Power(LongDouble_t x, Long64_t y)
   { return std::pow(x,(LongDouble_t)y); }

////////////////////////////////////////////////////////////////////////////////
inline LongDouble_t TMath::Power(Long64_t x, Long64_t y)
   { return std::pow(x,y); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Power(Double_t x, Double_t y)
   { return pow(x, y); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Power(Double_t x, Int_t y) {
#ifdef R__ANSISTREAM
   return std::pow(x, y);
#else
   return pow(x, (Double_t) y);
#endif
}

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Log(Double_t x)
   { return log(x); }

////////////////////////////////////////////////////////////////////////////////
inline Double_t TMath::Log10(Double_t x)
   { return log10(x); }

////////////////////////////////////////////////////////////////////////////////
/// Check if it is finite with a mask in order to be consistent in presence of
/// fast math.
/// Inspired from the CMSSW FWCore/Utilities package
inline Int_t TMath::Finite(Double_t x)
#if defined(R__FAST_MATH)

{
   const unsigned long long mask = 0x7FF0000000000000LL;
   union { unsigned long long l; double d;} v;
   v.d =x;
   return (v.l&mask)!=mask;
}
#else
#  if defined(R__HPUX11)
   { return isfinite(x); }
#  elif defined(R__MACOSX)
#  ifdef isfinite
   // from math.h
   { return isfinite(x); }
#  else
   // from cmath
   { return std::isfinite(x); }
#  endif
#  else
   { return finite(x); }
#  endif
#endif

////////////////////////////////////////////////////////////////////////////////
/// Check if it is finite with a mask in order to be consistent in presence of
/// fast math.
/// Inspired from the CMSSW FWCore/Utilities package
inline Int_t TMath::Finite(Float_t x)
#if defined(R__FAST_MATH)

{
   const unsigned int mask =  0x7f800000;
   union { unsigned int l; float d;} v;
   v.d =x;
   return (v.l&mask)!=mask;
}
#else
{ return std::isfinite(x); }
#endif

// This namespace provides all the routines necessary for checking if a number
// is a NaN also in presence of optimisations affecting the behaviour of the
// floating point calculations.
// Inspired from the CMSSW FWCore/Utilities package

#if defined (R__FAST_MATH)
namespace ROOT {
namespace Internal {
namespace Math {
// abridged from GNU libc 2.6.1 - in detail from
//   math/math_private.h
//   sysdeps/ieee754/ldbl-96/math_ldbl.h

// part of ths file:
   /*
    * ====================================================
    * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
    *
    * Developed at SunPro, a Sun Microsystems, Inc. business.
    * Permission to use, copy, modify, and distribute this
    * software is freely granted, provided that this notice
    * is preserved.
    * ====================================================
    */

   // A union which permits us to convert between a double and two 32 bit ints.
   typedef union {
      Double_t value;
      struct {
         UInt_t lsw;
         UInt_t msw;
      } parts;
   } ieee_double_shape_type;

#define EXTRACT_WORDS(ix0,ix1,d)                                    \
   do {                                                             \
      ieee_double_shape_type ew_u;                                  \
      ew_u.value = (d);                                             \
      (ix0) = ew_u.parts.msw;                                       \
      (ix1) = ew_u.parts.lsw;                                       \
   } while (0)

   inline Bool_t IsNaN(Double_t x)
   {
      UInt_t hx, lx;

      EXTRACT_WORDS(hx, lx, x);

      lx |= hx & 0xfffff;
      hx &= 0x7ff00000;
      return (hx == 0x7ff00000) && (lx != 0);
   }

   typedef union {
      Float_t value;
      UInt_t word;
   } ieee_float_shape_type;

#define GET_FLOAT_WORD(i,d)                                         \
    do {                                                            \
      ieee_float_shape_type gf_u;                                   \
      gf_u.value = (d);                                             \
      (i) = gf_u.word;                                              \
    } while (0)

   inline Bool_t IsNaN(Float_t x)
   {
      UInt_t wx;
      GET_FLOAT_WORD (wx, x);
      wx &= 0x7fffffff;
      return (Bool_t)(wx > 0x7f800000);
   }
} } } // end NS ROOT::Internal::Math
#endif // End R__FAST_MATH

#if defined(R__FAST_MATH)
   inline Bool_t TMath::IsNaN(Double_t x) { return ROOT::Internal::Math::IsNaN(x); }
   inline Bool_t TMath::IsNaN(Float_t x) { return ROOT::Internal::Math::IsNaN(x); }
#else
   inline Bool_t TMath::IsNaN(Double_t x) { return std::isnan(x); }
   inline Bool_t TMath::IsNaN(Float_t x) { return std::isnan(x); }
#endif

////////////////////////////////////////////////////////////////////////////////
// Wrapper to numeric_limits

////////////////////////////////////////////////////////////////////////////////
/// Returns a quiet NaN as [defined by IEEE 754](http://en.wikipedia.org/wiki/NaN#Quiet_NaN)
inline Double_t TMath::QuietNaN() {

   return std::numeric_limits<Double_t>::quiet_NaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a signaling NaN as defined by IEEE 754](http://en.wikipedia.org/wiki/NaN#Signaling_NaN)
inline Double_t TMath::SignalingNaN() {
   return std::numeric_limits<Double_t>::signaling_NaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an infinity as defined by the IEEE standard
inline Double_t TMath::Infinity() {
   return std::numeric_limits<Double_t>::infinity();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns maximum representation for type T
template<typename T>
inline T TMath::Limits<T>::Min() {
   return (std::numeric_limits<T>::min)();    //N.B. use this signature to avoid class with macro min() on Windows
}

////////////////////////////////////////////////////////////////////////////////
/// Returns minimum double representation
template<typename T>
inline T TMath::Limits<T>::Max() {
   return (std::numeric_limits<T>::max)();  //N.B. use this signature to avoid class with macro max() on Windows
}

////////////////////////////////////////////////////////////////////////////////
/// Returns minimum double representation
template<typename T>
inline T TMath::Limits<T>::Epsilon() {
   return std::numeric_limits<T>::epsilon();
}

////////////////////////////////////////////////////////////////////////////////
// Advanced.

////////////////////////////////////////////////////////////////////////////////
/// Calculate the Normalized Cross Product of two vectors
template <typename T> inline T TMath::NormCross(const T v1[3],const T v2[3],T out[3])
{
   return Normalize(Cross(v1,v2,out));
}

////////////////////////////////////////////////////////////////////////////////
/// Return minimum of array a of length n.
template <typename T>
T TMath::MinElement(Long64_t n, const T *a) {
   return *std::min_element(a,a+n);
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum of array a of length n.
template <typename T>
T TMath::MaxElement(Long64_t n, const T *a) {
   return *std::max_element(a,a+n);
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of array with the minimum element.
/// If more than one element is minimum returns first found.
///
/// Implement here since this one is found to be faster (mainly on 64 bit machines)
/// than stl generic implementation.
/// When performing the comparison,  the STL implementation needs to de-reference both the array iterator
/// and the iterator pointing to the resulting minimum location
template <typename T>
Long64_t TMath::LocMin(Long64_t n, const T *a) {
   if  (n <= 0 || !a) return -1;
   T xmin = a[0];
   Long64_t loc = 0;
   for  (Long64_t i = 1; i < n; i++) {
      if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
      }
   }
   return loc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of array with the minimum element.
/// If more than one element is minimum returns first found.
template <typename Iterator>
Iterator TMath::LocMin(Iterator first, Iterator last) {

   return std::min_element(first, last);
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of array with the maximum element.
/// If more than one element is maximum returns first found.
///
/// Implement here since it is faster (see comment in LocMin function)
template <typename T>
Long64_t TMath::LocMax(Long64_t n, const T *a) {
   if  (n <= 0 || !a) return -1;
   T xmax = a[0];
   Long64_t loc = 0;
   for  (Long64_t i = 1; i < n; i++) {
      if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
      }
   }
   return loc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of array with the maximum element.
/// If more than one element is maximum returns first found.
template <typename Iterator>
Iterator TMath::LocMax(Iterator first, Iterator last)
{

   return std::max_element(first, last);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted mean of an array defined by the iterators.
template <typename Iterator>
Double_t TMath::Mean(Iterator first, Iterator last)
{
   Double_t sum = 0;
   Double_t sumw = 0;
   while ( first != last )
   {
      sum += *first;
      sumw += 1;
      first++;
   }

   return sum/sumw;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted mean of an array defined by the first and
/// last iterators. The w iterator should point to the first element
/// of a vector of weights of the same size as the main array.
template <typename Iterator, typename WeightIterator>
Double_t TMath::Mean(Iterator first, Iterator last, WeightIterator w)
{

   Double_t sum = 0;
   Double_t sumw = 0;
   int i = 0;
   while ( first != last ) {
      if ( *w < 0) {
         ::Error("TMath::Mean","w[%d] = %.4e < 0 ?!",i,*w);
         return 0;
      }
      sum  += (*w) * (*first);
      sumw += (*w) ;
      ++w;
      ++first;
      ++i;
   }
   if (sumw <= 0) {
      ::Error("TMath::Mean","sum of weights == 0 ?!");
      return 0;
   }

   return sum/sumw;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted mean of an array a with length n.
template <typename T>
Double_t TMath::Mean(Long64_t n, const T *a, const Double_t *w)
{
   if (w) {
      return TMath::Mean(a, a+n, w);
   } else {
      return TMath::Mean(a, a+n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the geometric mean of an array defined by the iterators.
/// \f[ GeomMean = (\prod_{i=0}^{n-1} |a[i]|)^{1/n} \f]
template <typename Iterator>
Double_t TMath::GeomMean(Iterator first, Iterator last)
{
   Double_t logsum = 0.;
   Long64_t n = 0;
   while ( first != last ) {
      if (*first == 0) return 0.;
      Double_t absa = (Double_t) TMath::Abs(*first);
      logsum += TMath::Log(absa);
      ++first;
      ++n;
   }

   return TMath::Exp(logsum/n);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the geometric mean of an array a of size n.
/// \f[ GeomMean = (\prod_{i=0}^{n-1} |a[i]|)^{1/n} \f]
template <typename T>
Double_t TMath::GeomMean(Long64_t n, const T *a)
{
   return TMath::GeomMean(a, a+n);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Standard Deviation of an array defined by the iterators.
/// Note that this function returns the sigma(standard deviation) and
/// not the root mean square of the array.
///
/// Use the two pass algorithm, which is slower (! a factor of 2) but much more
/// precise.  Since we have a vector the 2 pass algorithm is still faster than the
/// Welford algorithm. (See also ROOT-5545)
template <typename Iterator>
Double_t TMath::RMS(Iterator first, Iterator last)
{

   Double_t n = 0;

   Double_t tot = 0;
   Double_t mean = TMath::Mean(first,last);
   while ( first != last ) {
      Double_t x = Double_t(*first);
      tot += (x - mean)*(x - mean);
      ++first;
      ++n;
   }
   Double_t rms = (n > 1) ? TMath::Sqrt(tot/(n-1)) : 0.0;
   return rms;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted Standard Deviation of an array defined by the iterators.
/// Note that this function returns the sigma(standard deviation) and
/// not the root mean square of the array.
///
/// As in the unweighted case use the two pass algorithm
template <typename Iterator, typename WeightIterator>
Double_t TMath::RMS(Iterator first, Iterator last, WeightIterator w)
{
   Double_t tot = 0;
   Double_t sumw = 0;
   Double_t sumw2 = 0;
   Double_t mean = TMath::Mean(first,last,w);
   while ( first != last ) {
      Double_t x = Double_t(*first);
      sumw += *w;
      sumw2 += (*w) * (*w);
      tot += (*w) * (x - mean)*(x - mean);
      ++first;
      ++w;
   }
   // use the correction neff/(neff -1) for the unbiased formula
   Double_t rms =  TMath::Sqrt(tot * sumw/ (sumw*sumw - sumw2) );
   return rms;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Standard Deviation of an array a with length n.
/// Note that this function returns the sigma(standard deviation) and
/// not the root mean square of the array.
template <typename T>
Double_t TMath::RMS(Long64_t n, const T *a, const Double_t * w)
{
   return (w) ? TMath::RMS(a, a+n, w) : TMath::RMS(a, a+n);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the Cross Product of two vectors:
///         out = [v1 x v2]
template <typename T> T *TMath::Cross(const T v1[3],const T v2[3], T out[3])
{
   out[0] = v1[1] * v2[2] - v1[2] * v2[1];
   out[1] = v1[2] * v2[0] - v1[0] * v2[2];
   out[2] = v1[0] * v2[1] - v1[1] * v2[0];

   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate a normal vector of a plane.
///
/// \param[in]  p1, p2,p3     3 3D points belonged the plane to define it.
/// \param[out] normal        Pointer to 3D normal vector (normalized)
template <typename T> T * TMath::Normal2Plane(const T p1[3],const T p2[3],const T p3[3], T normal[3])
{
   T v1[3], v2[3];

   v1[0] = p2[0] - p1[0];
   v1[1] = p2[1] - p1[1];
   v1[2] = p2[2] - p1[2];

   v2[0] = p3[0] - p1[0];
   v2[1] = p3[1] - p1[1];
   v2[2] = p3[2] - p1[2];

   NormCross(v1,v2,normal);
   return normal;
}

////////////////////////////////////////////////////////////////////////////////
/// Function which returns kTRUE if point xp,yp lies inside the
/// polygon defined by the np points in arrays x and y, kFALSE otherwise.
/// Note that the polygon may be open or closed.
template <typename T> Bool_t TMath::IsInside(T xp, T yp, Int_t np, T *x, T *y)
{
   Int_t i, j = np-1 ;
   Bool_t oddNodes = kFALSE;

   for (i=0; i<np; i++) {
      if ((y[i]<yp && y[j]>=yp) || (y[j]<yp && y[i]>=yp)) {
         if (x[i]+(yp-y[i])/(y[j]-y[i])*(x[j]-x[i])<xp) {
            oddNodes = !oddNodes;
         }
      }
      j=i;
   }

   return oddNodes;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the median of the array a where each entry i has weight w[i] .
/// Both arrays have a length of at least n . The median is a number obtained
/// from the sorted array a through
///
/// median = (a[jl]+a[jh])/2.  where (using also the sorted index on the array w)
///
/// sum_i=0,jl w[i] <= sumTot/2
/// sum_i=0,jh w[i] >= sumTot/2
/// sumTot = sum_i=0,n w[i]
///
/// If w=0, the algorithm defaults to the median definition where it is
/// a number that divides the sorted sequence into 2 halves.
/// When n is odd or n > 1000, the median is kth element k = (n + 1) / 2.
/// when n is even and n < 1000the median is a mean of the elements k = n/2 and k = n/2 + 1.
///
/// If the weights are supplied (w not 0) all weights must be >= 0
///
/// If work is supplied, it is used to store the sorting index and assumed to be
/// >= n . If work=0, local storage is used, either on the stack if n < kWorkMax
/// or on the heap for n >= kWorkMax .
template <typename T> Double_t TMath::Median(Long64_t n, const T *a,  const Double_t *w, Long64_t *work)
{

   const Int_t kWorkMax = 100;

   if (n <= 0 || !a) return 0;
   Bool_t isAllocated = kFALSE;
   Double_t median;
   Long64_t *ind;
   Long64_t workLocal[kWorkMax];

   if (work) {
      ind = work;
   } else {
      ind = workLocal;
      if (n > kWorkMax) {
         isAllocated = kTRUE;
         ind = new Long64_t[n];
      }
   }

   if (w) {
      Double_t sumTot2 = 0;
      for (Int_t j = 0; j < n; j++) {
         if (w[j] < 0) {
            ::Error("TMath::Median","w[%d] = %.4e < 0 ?!",j,w[j]);
            if (isAllocated)  delete [] ind;
            return 0;
         }
         sumTot2 += w[j];
      }

      sumTot2 /= 2.;

      Sort(n, a, ind, kFALSE);

      Double_t sum = 0.;
      Int_t jl;
      for (jl = 0; jl < n; jl++) {
         sum += w[ind[jl]];
         if (sum >= sumTot2) break;
      }

      Int_t jh;
      sum = 2.*sumTot2;
      for (jh = n-1; jh >= 0; jh--) {
         sum -= w[ind[jh]];
         if (sum <= sumTot2) break;
      }

      median = 0.5*(a[ind[jl]]+a[ind[jh]]);

   } else {

      if (n%2 == 1)
         median = KOrdStat(n, a,n/2, ind);
      else {
         median = 0.5*(KOrdStat(n, a, n/2 -1, ind)+KOrdStat(n, a, n/2, ind));
      }
   }

   if (isAllocated)
      delete [] ind;
   return median;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns k_th order statistic of the array a of size n
/// (k_th smallest element out of n elements).
///
/// C-convention is used for array indexing, so if you want
/// the second smallest element, call KOrdStat(n, a, 1).
///
/// If work is supplied, it is used to store the sorting index and
/// assumed to be >= n. If work=0, local storage is used, either on
/// the stack if n < kWorkMax or on the heap for n >= kWorkMax.
/// Note that the work index array will not contain the sorted indices but
/// all indices of the smaller element in arbitrary order in work[0,...,k-1] and
/// all indices of the larger element in arbitrary order in work[k+1,..,n-1]
/// work[k] will contain instead the index of the returned element.
///
/// Taken from "Numerical Recipes in C++" without the index array
/// implemented by Anna Khreshuk.
///
/// See also the declarations at the top of this file
template <class Element, typename Size>
Element TMath::KOrdStat(Size n, const Element *a, Size k, Size *work)
{

   const Int_t kWorkMax = 100;

   typedef Size Index;

   Bool_t isAllocated = kFALSE;
   Size i, ir, j, l, mid;
   Index arr;
   Index *ind;
   Index workLocal[kWorkMax];
   Index temp;

   if (work) {
      ind = work;
   } else {
      ind = workLocal;
      if (n > kWorkMax) {
         isAllocated = kTRUE;
         ind = new Index[n];
      }
   }

   for (Size ii=0; ii<n; ii++) {
      ind[ii]=ii;
   }
   Size rk = k;
   l=0;
   ir = n-1;
   for(;;) {
      if (ir<=l+1) { //active partition contains 1 or 2 elements
         if (ir == l+1 && a[ind[ir]]<a[ind[l]])
            {temp = ind[l]; ind[l]=ind[ir]; ind[ir]=temp;}
         Element tmp = a[ind[rk]];
         if (isAllocated)
            delete [] ind;
         return tmp;
      } else {
         mid = (l+ir) >> 1; //choose median of left, center and right
         {temp = ind[mid]; ind[mid]=ind[l+1]; ind[l+1]=temp;}//elements as partitioning element arr.
         if (a[ind[l]]>a[ind[ir]])  //also rearrange so that a[l]<=a[l+1]
            {temp = ind[l]; ind[l]=ind[ir]; ind[ir]=temp;}

         if (a[ind[l+1]]>a[ind[ir]])
            {temp=ind[l+1]; ind[l+1]=ind[ir]; ind[ir]=temp;}

         if (a[ind[l]]>a[ind[l+1]])
            {temp = ind[l]; ind[l]=ind[l+1]; ind[l+1]=temp;}

         i=l+1;        //initialize pointers for partitioning
         j=ir;
         arr = ind[l+1];
         for (;;){
            do i++; while (a[ind[i]]<a[arr]);
            do j--; while (a[ind[j]]>a[arr]);
            if (j<i) break;  //pointers crossed, partitioning complete
               {temp=ind[i]; ind[i]=ind[j]; ind[j]=temp;}
         }
         ind[l+1]=ind[j];
         ind[j]=arr;
         if (j>=rk) ir = j-1; //keep active the partition that
         if (j<=rk) l=i;      //contains the k_th element
      }
   }
}

#endif
