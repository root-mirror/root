// @(#)root/mathcore:$Id$
// Authors: L. Moneta    8/2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015 , ROOT MathLib Team                             *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// file for random class
//
//
// Created by: Lorenzo Moneta  : Tue 4 Aug 2015
//
//
#include "Math/RandomFunctions.h"

#include "Math/DistFuncMathCore.h"

#include "TMath.h"

namespace ROOT {
namespace Math {
   

Int_t RandomFunctionsImpl<TRandomEngine>::Binomial(Int_t ntot, Double_t prob)
{
   if (prob < 0 || prob > 1) return 0;
   Int_t n = 0;
   for (Int_t i=0;i<ntot;i++) {
      if (Rndm() > prob) continue;
      n++;
   }
   return n;
}

   ////////////////////////////////////////////////////////////////////////////////
/// Return a number distributed following a BreitWigner function with mean and gamma.

Double_t RandomFunctionsImpl<TRandomEngine>::BreitWigner(Double_t mean, Double_t gamma)
{
   Double_t rval, displ;
   rval = 2*Rndm() - 1;
   displ = 0.5*gamma*TMath::Tan(rval*TMath::PiOver2());

   return (mean+displ);
}

////////////////////////////////////////////////////////////////////////////////
/// Generates random vectors, uniformly distributed over a circle of given radius.
///   Input : r = circle radius
///   Output: x,y a random 2-d vector of length r

void RandomFunctionsImpl<TRandomEngine>::Circle(Double_t &x, Double_t &y, Double_t r)
{
   Double_t phi = Uniform(0,TMath::TwoPi());
   x = r*TMath::Cos(phi);
   y = r*TMath::Sin(phi);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an exponential deviate.
///
///          exp( -t/tau )

Double_t RandomFunctionsImpl<TRandomEngine>::Exp(Double_t tau)
{
   Double_t x = Rndm();              // uniform on ] 0, 1 ]
   Double_t t = -tau * TMath::Log( x ); // convert to exponential distribution
   return t;
}

   

   
double RandomFunctionsImpl<TRandomEngine>::GausBM(double mean, double sigma) {
   double y =  Rndm();
   double z =  Rndm();
   double x = z * 6.28318530717958623;
   double radius = std::sqrt(-2*std::log(y));
   double g = radius * std::sin(x);
   return mean + g * sigma; 
}


   // double GausImpl(TRandomEngine * r, double mean, double sigma) {
   //    double y =  r->Rndm();
   //    double z =  r->Rndm();
   //    double x = z * 6.28318530717958623;
   //    double radius = std::sqrt(-2*std::log(y));
   //    double g = radius * std::sin(x);
   //    return mean + g * sigma; 
   // }

double RandomFunctionsImpl<TRandomEngine>::GausACR(Double_t mean, Double_t sigma)
{
   const Double_t kC1 = 1.448242853;
   const Double_t kC2 = 3.307147487;
   const Double_t kC3 = 1.46754004;
   const Double_t kD1 = 1.036467755;
   const Double_t kD2 = 5.295844968;
   const Double_t kD3 = 3.631288474;
   const Double_t kHm = 0.483941449;
   const Double_t kZm = 0.107981933;
   const Double_t kHp = 4.132731354;
   const Double_t kZp = 18.52161694;
   const Double_t kPhln = 0.4515827053;
   const Double_t kHm1 = 0.516058551;
   const Double_t kHp1 = 3.132731354;
   const Double_t kHzm = 0.375959516;
   const Double_t kHzmp = 0.591923442;
   /*zhm 0.967882898*/

   const Double_t kAs = 0.8853395638;
   const Double_t kBs = 0.2452635696;
   const Double_t kCs = 0.2770276848;
   const Double_t kB  = 0.5029324303;
   const Double_t kX0 = 0.4571828819;
   const Double_t kYm = 0.187308492 ;
   const Double_t kS  = 0.7270572718 ;
   const Double_t kT  = 0.03895759111;

   Double_t result;
   Double_t rn,x,y,z;

   do {
      y = Rndm();

      if (y>kHm1) {
         result = kHp*y-kHp1; break; }

      else if (y<kZm) {
         rn = kZp*y-1;
         result = (rn>0) ? (1+rn) : (-1+rn);
         break;
      }

      else if (y<kHm) {
         rn = Rndm();
         rn = rn-1+rn;
         z = (rn>0) ? 2-rn : -2-rn;
         if ((kC1-y)*(kC3+TMath::Abs(z))<kC2) {
            result = z; break; }
         else {
            x = rn*rn;
            if ((y+kD1)*(kD3+x)<kD2) {
               result = rn; break; }
            else if (kHzmp-y<exp(-(z*z+kPhln)/2)) {
               result = z; break; }
            else if (y+kHzm<exp(-(x+kPhln)/2)) {
               result = rn; break; }
         }
      }

      while (1) {
         x = Rndm();
         y = kYm * Rndm();
         z = kX0 - kS*x - y;
         if (z>0)
            rn = 2+y/x;
         else {
            x = 1-x;
            y = kYm-y;
            rn = -(2+y/x);
         }
         if ((y-kAs+x)*(kCs+x)+kBs<0) {
            result = rn; break; }
         else if (y<x+kT)
            if (rn*rn<4*(kB-log(x))) {
               result = rn; break; }
      }
   } while(0);

   return mean + sigma * result;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a random number following a Landau distribution
/// with location parameter mu and scale parameter sigma:
///      Landau( (x-mu)/sigma )
/// Note that mu is not the mpv(most probable value) of the Landa distribution
/// and sigma is not the standard deviation of the distribution which is not defined.
/// For mu  =0 and sigma=1, the mpv = -0.22278
///
/// The Landau random number generation is implemented using the
/// function landau_quantile(x,sigma), which provides
/// the inverse of the landau cumulative distribution.
/// landau_quantile has been converted from CERNLIB ranlan(G110).

Double_t RandomFunctionsImpl<TRandomEngine>::Landau(Double_t mu, Double_t sigma)
{
   if (sigma <= 0) return 0;
   Double_t x = Rndm();
   Double_t res = mu + ROOT::Math::landau_quantile(x, sigma);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Generates a random integer N according to a Poisson law.
/// Prob(N) = exp(-mean)*mean^N/Factorial(N)
///
/// Use a different procedure according to the mean value.
/// The algorithm is the same used by CLHEP.
/// For lower value (mean < 25) use the rejection method based on
/// the exponential.
/// For higher values use a rejection method comparing with a Lorentzian
/// distribution, as suggested by several authors.
/// This routine since is returning 32 bits integer will not work for values
/// larger than 2*10**9.
/// One should then use the Trandom::PoissonD for such large values.

Int_t RandomFunctionsImpl<TRandomEngine>::Poisson(Double_t mean)
{
   Int_t n;
   if (mean <= 0) return 0;
   if (mean < 25) {
      Double_t expmean = TMath::Exp(-mean);
      Double_t pir = 1;
      n = -1;
      while(1) {
         n++;
         pir *= Rndm();
         if (pir <= expmean) break;
      }
      return n;
   }
   // for large value we use inversion method
   else if (mean < 1E9) {
      Double_t em, t, y;
      Double_t sq, alxm, g;
      Double_t pi = TMath::Pi();

      sq = TMath::Sqrt(2.0*mean);
      alxm = TMath::Log(mean);
      g = mean*alxm - TMath::LnGamma(mean + 1.0);

      do {
         do {
            y = TMath::Tan(pi*Rndm());
            em = sq*y + mean;
         } while( em < 0.0 );

         em = TMath::Floor(em);
         t = 0.9*(1.0 + y*y)* TMath::Exp(em*alxm - TMath::LnGamma(em + 1.0) - g);
      } while( Rndm() > t );

      return static_cast<Int_t> (em);

   }
   else {
      // use Gaussian approximation vor very large values
      n = Int_t(Gaus(0,1)*TMath::Sqrt(mean) + mean +0.5);
      return n;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generates a random number according to a Poisson law.
/// Prob(N) = exp(-mean)*mean^N/Factorial(N)
///
/// This function is a variant of RandomFunctionsImpl<TRandomEngine>::Poisson returning a double
/// instead of an integer.

Double_t RandomFunctionsImpl<TRandomEngine>::PoissonD(Double_t mean)
{
   Int_t n;
   if (mean <= 0) return 0;
   if (mean < 25) {
      Double_t expmean = TMath::Exp(-mean);
      Double_t pir = 1;
      n = -1;
      while(1) {
         n++;
         pir *= Rndm();
         if (pir <= expmean) break;
      }
      return static_cast<Double_t>(n);
   }
   // for large value we use inversion method
   else if (mean < 1E9) {
      Double_t em, t, y;
      Double_t sq, alxm, g;
      Double_t pi = TMath::Pi();

      sq = TMath::Sqrt(2.0*mean);
      alxm = TMath::Log(mean);
      g = mean*alxm - TMath::LnGamma(mean + 1.0);

      do {
         do {
            y = TMath::Tan(pi*Rndm());
            em = sq*y + mean;
         } while( em < 0.0 );

         em = TMath::Floor(em);
         t = 0.9*(1.0 + y*y)* TMath::Exp(em*alxm - TMath::LnGamma(em + 1.0) - g);
      } while( Rndm() > t );

      return em;

   } else {
      // use Gaussian approximation vor very large values
      return Gaus(0,1)*TMath::Sqrt(mean) + mean +0.5;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return 2 numbers distributed following a gaussian with mean=0 and sigma=1.

void RandomFunctionsImpl<TRandomEngine>::Rannor(Double_t &a, Double_t &b)
{
   Double_t r, x, y, z;

   y = Rndm();
   z = Rndm();
   x = z * 6.28318530717958623;
   r = TMath::Sqrt(-2*TMath::Log(y));
   a = r * TMath::Sin(x);
   b = r * TMath::Cos(x);
}
   
////////////////////////////////////////////////////////////////////////////////
/// Generates random vectors, uniformly distributed over the surface
/// of a sphere of given radius.
///   Input : r = sphere radius
///   Output: x,y,z a random 3-d vector of length r
/// Method: (based on algorithm suggested by Knuth and attributed to Robert E Knop)
///         which uses less random numbers than the CERNLIB RN23DIM algorithm

void RandomFunctionsImpl<TRandomEngine>::Sphere(Double_t &x, Double_t &y, Double_t &z, Double_t r)
{
   Double_t a=0,b=0,r2=1;
   while (r2 > 0.25) {
      a  = Rndm() - 0.5;
      b  = Rndm() - 0.5;
      r2 =  a*a + b*b;
   }
   z = r* ( -1. + 8.0 * r2 );

   Double_t scale = 8.0 * r * TMath::Sqrt(0.25 - r2);
   x = a*scale;
   y = b*scale;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a uniform deviate on the interval  (0, x1).

Double_t RandomFunctionsImpl<TRandomEngine>::Uniform(Double_t x1)
{
   Double_t ans = Rndm();
   return x1*ans;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a uniform deviate on the interval (x1, x2).

double RandomFunctionsImpl<TRandomEngine>::Uniform(double a, double b) {
   return (b-a) * Rndm() + a; 
}
   
  
   } // namespace Math
} // namespace ROOT
