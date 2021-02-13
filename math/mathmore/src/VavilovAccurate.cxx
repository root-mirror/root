// @(#)root/mathmore:$Id$
// Authors: B. List 29.4.2010


 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2004 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  * This library is free software; you can redistribute it and/or      *
  * modify it under the terms of the GNU General Public License        *
  * as published by the Free Software Foundation; either version 2     *
  * of the License, or (at your option) any later version.             *
  *                                                                    *
  * This library is distributed in the hope that it will be useful,    *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
  * General Public License for more details.                           *
  *                                                                    *
  * You should have received a copy of the GNU General Public License  *
  * along with this library (see file COPYING); if not, write          *
  * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
  * 330, Boston, MA 02111-1307 USA, or contact the author.             *
  *                                                                    *
  **********************************************************************/

// Implementation file for class VavilovAccurate
//
// Created by: blist  at Thu Apr 29 11:19:00 2010
//
// Last update: Thu Apr 29 11:19:00 2010
//


#include "Math/VavilovAccurate.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/QuantFuncMathCore.h"

#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>


namespace ROOT {
namespace Math {

VavilovAccurate *VavilovAccurate::fgInstance = nullptr;


VavilovAccurate::VavilovAccurate(double kappa, double beta2, double epsilonPM, double epsilon)
{
   Set (kappa, beta2, epsilonPM, epsilon);
}


VavilovAccurate::~VavilovAccurate()
{
   // desctructor (clean up resources)
}

void VavilovAccurate::SetKappaBeta2 (double kappa, double beta2) {
   Set (kappa, beta2);
}

void VavilovAccurate::Set(double kappa, double beta2, double epsilonPM, double epsilon) {
   // Method described in
   // B. Schorr, Programs for the Landau and the Vavilov distributions and the corresponding random numbers,
   // <A HREF="http://dx.doi.org/10.1016/0010-4655(74)90091-5">Computer Phys. Comm. 7 (1974) 215-224</A>.
   fQuantileInit = false;

   fKappa = kappa;
   fBeta2 = beta2;
   fEpsilonPM = epsilonPM;    // epsilon_+ = epsilon_-: determines support (T0, T1)
   fEpsilon = epsilon;

   static const double eu = 0.577215664901532860606;              // Euler's constant
   static const double pi2 = 6.28318530717958647693,              // 2pi
                       rpi = 0.318309886183790671538,             // 1/pi
                       pih = 1.57079632679489661923;              // pi/2
   double h1 = -std::log(fEpsilon)-1.59631259113885503887;        // -ln(fEpsilon) + ln(2/pi**2)
   double deltaEpsilon = 0.001;
   static const double logdeltaEpsilon = -std::log(deltaEpsilon); // 3 ln 10 = -ln(.001);
   double logEpsilonPM = std::log(fEpsilonPM);
   static const double eps = 1e-5;                                // accuracy of root finding for x0

   double xp[9] = {0,
                   9.29,  2.47, 0.89, 0.36, 0.15, 0.07, 0.03, 0.02};
   double xq[7] = {0,
                   0.012, 0.03, 0.08, 0.26, 0.87, 3.83};

   if (kappa < 0.001) {
      std::cerr << "VavilovAccurate::Set: kappa = " << kappa << " - out of range" << std::endl;
      if (kappa < 0.001) kappa = 0.001;
   }
   if (beta2 < 0 || beta2 > 1) {
      std::cerr << "VavilovAccurate::Set: beta2 = " << beta2 << " - out of range" << std::endl;
      if (beta2 < 0) beta2 = -beta2;
      if (beta2 > 1) beta2 = 1;
   }

   // approximation of x_-
   fH[5] = 1-beta2*(1-eu)-logEpsilonPM/kappa;       // eq. 3.9
   fH[6] = beta2;
   fH[7] = 1-beta2;
   double h4 = logEpsilonPM/kappa-(1+beta2*eu);
   double logKappa = std::log(kappa);
   double kappaInv = 1/kappa;
   // Calculate T0 from Eq. (3.6), using x_- = fH[5]
//    double e1h5 = (fH[5] > 40 ) ? 0 : -ROOT::Math::expint (-fH[5]);
//    fT0 = (h4-fH[5]*logKappa-(fH[5]+beta2)*(std::log(fH[5])+e1h5)+std::exp(-fH[5]))/fH[5];
   fT0 = (h4-fH[5]*logKappa-(fH[5]+beta2)*E1plLog(fH[5])+std::exp(-fH[5]))/fH[5];
   int lp = 1;
   while (lp < 9 && kappa < xp[lp]) ++lp;
   int lq = 1;
   while (lq < 7 && kappa >= xq[lq]) ++lq;
   // Solve eq. 3.7 to get fH[0] = x_+
   double delta = 0;
   int ifail = 0;
   do {
      ifail = Rzero(-lp-0.5-delta,lq-7.5+delta,fH[0],eps,1000,&ROOT::Math::VavilovAccurate::G116f2);
      delta += 0.5;
   } while (ifail == 2);

   double q = 1/fH[0];
   // Calculate T1 from Eq. (3.6)
//    double e1h0 = (fH[0] > 40 ) ? 0 : -ROOT::Math::expint (-fH[0]);
//    fT1 = h4*q-logKappa-(1+beta2*q)*(std::log(std::fabs(fH[0]))+e1h0)+std::exp(-fH[0])*q;
   fT1 = h4*q-logKappa-(1+beta2*q)*E1plLog(fH[0])+std::exp(-fH[0])*q;

   fT = fT1-fT0;                         // Eq. (2.5)
   fOmega = pi2/fT;                      // Eq. (2.5)
   fH[1] = kappa*(2+beta2*eu)+h1;
   if(kappa >= 0.07) fH[1] += logdeltaEpsilon;       // reduce fEpsilon by a factor .001 for large kappa
   fH[2] = beta2*kappa;
   fH[3] = kappaInv*fOmega;
   fH[4] = pih*fOmega;

   // Solve log(eq. (4.10)) to get fX0 = N
   ifail = Rzero(5.,MAXTERMS,fX0,eps,1000,&ROOT::Math::VavilovAccurate::G116f1);
//    if (ifail) {
//       std::cerr << "Rzero failed for x0: ifail=" << ifail << ", kappa=" << kappa << ", beta2=" << beta2 << std::endl;
//       std::cerr << "G116f1(" << 5. << ")=" << G116f1(5.) << ", G116f1(" << MAXTERMS << ")=" << G116f1(MAXTERMS)  << std::endl;
//       std::cerr << "fH[0]=" << fH[0] << ", fH[1]=" << fH[1] << ", fH[2]=" << fH[2] << ", fH[3]=" << fH[3] << ", fH[4]=" << fH[4] << std::endl;
//       std::cerr << "fH[5]=" << fH[5] << ", fH[6]=" << fH[6] << ", fH[7]=" << fH[7] << std::endl;
//       std::cerr << "fT0=" << fT0 << ", fT1=" << fT1 << std::endl;
//       std::cerr << "G116f2(" << fH[0] << ")=" << G116f2(fH[0]) << std::endl;
//    }
   if (ifail == 2) {
      fX0 = (G116f1(5) > G116f1(MAXTERMS)) ? MAXTERMS : 5;
   }
   if (fX0 < 5) fX0 = 5;
   else if (fX0 > MAXTERMS) fX0 = MAXTERMS;
   int n = int(fX0+1);
   // logKappa=log(kappa)
   double d = rpi*std::exp(kappa*(1+beta2*(eu-logKappa)));
   fA_pdf[n] = rpi*fOmega;
   fA_cdf[n] = 0;
   q = -1;
   double q2 = 2;
   for (int k = 1; k < n; ++k) {
      int l = n-k;
      double x = fOmega*k;
      double x1 = kappaInv*x;
      double c1 = std::log(x)-ROOT::Math::cosint(x1);
      double c2 = ROOT::Math::sinint(x1);
      double c3 = std::sin(x1);
      double c4 = std::cos(x1);
      double xf1 = kappa*(beta2*c1-c4)-x*c2;
      double xf2 = x*(c1 + fT0) + kappa*(c3+beta2*c2);
      double d1 = q*d*fOmega*std::exp(xf1);
      double s = std::sin(xf2);
      double c = std::cos(xf2);
      fA_pdf[l] = d1*c;
      fB_pdf[l] = -d1*s;
      d1 = q*d*std::exp(xf1)/k;
      fA_cdf[l] = d1*s;
      fB_cdf[l] = d1*c;
      fA_cdf[n] += q2*fA_cdf[l];
      q = -q;
      q2 = -q2;
   }
}

void VavilovAccurate::InitQuantile() const {
   fQuantileInit = true;

   fNQuant = 16;
   // for kappa<0.02: use landau_quantile as first approximation
   if (fKappa < 0.02) return;
   else if (fKappa < 0.05) fNQuant = 32;

   // crude approximation for the median:

   double estmedian = -4.22784335098467134e-01-std::log(fKappa)-fBeta2;
   if (estmedian>1.3) estmedian = 1.3;

   // distribute test values evenly below and above the median
   for (int i = 1; i < fNQuant/2; ++i) {
      double x = fT0 + i*(estmedian-fT0)/(fNQuant/2);
      fQuant[i] = Cdf(x);
      fLambda[i] = x;
   }
   for (int i = fNQuant/2; i < fNQuant-1; ++i) {
      double x = estmedian + (i-fNQuant/2)*(fT1-estmedian)/(fNQuant/2-1);
      fQuant[i] = Cdf(x);
      fLambda[i] = x;
   }

   fQuant[0] = 0;
   fLambda[0] = fT0;
   fQuant[fNQuant-1] = 1;
   fLambda[fNQuant-1] = fT1;

}

double VavilovAccurate::Pdf (double x) const {

   static const double pi = 3.14159265358979323846;       // pi

   int n = int(fX0);
   double f;
   if (x < fT0) {
      f = 0;
   } else if (x <= fT1) {
      double y = x-fT0;
      double u = fOmega*y-pi;
      double cof = 2*cos(u);
      double a1 = 0;
      double a0 = fA_pdf[1];
      double a2 = 0;
      for (int k = 2; k <= n+1; ++k) {
         a2 = a1;
         a1 = a0;
         a0 = fA_pdf[k]+cof*a1-a2;
      }
      double b1 = 0;
      double b0 = fB_pdf[1];
      for (int k = 2; k <= n; ++k) {
         double b2 = b1;
         b1 = b0;
         b0 = fB_pdf[k]+cof*b1-b2;
      }
      f = 0.5*(a0-a2)+b0*sin(u);
   } else {
      f = 0;
   }
   return f;
}

double VavilovAccurate::Pdf (double x, double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Pdf (x);
}

double VavilovAccurate::Cdf (double x) const {

   static const double pi = 3.14159265358979323846;       // pi

   int n = int(fX0);
   double f;
   if (x < fT0) {
      f = 0;
   } else if (x <= fT1) {
      double y = x-fT0;
      double u = fOmega*y-pi;
      double cof = 2*cos(u);
      double a1 = 0;
      double a0 = fA_cdf[1];
      double a2 = 0;
      for (int k = 2; k <= n+1; ++k) {
         a2 = a1;
         a1 = a0;
         a0 = fA_cdf[k]+cof*a1-a2;
      }
      double b1 = 0;
      double b0 = fB_cdf[1];
      for (int k = 2; k <= n; ++k) {
         double b2 = b1;
         b1 = b0;
         b0 = fB_cdf[k]+cof*b1-b2;
      }
      f = 0.5*(a0-a2)+b0*sin(u);
      f += y/fT;
   } else {
      f = 1;
   }
   return f;
}

double VavilovAccurate::Cdf (double x, double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Cdf (x);
}

double VavilovAccurate::Cdf_c (double x) const {

   static const double pi = 3.14159265358979323846;       // pi

   int n = int(fX0);
   double f;
   if (x < fT0) {
      f = 1;
   } else if (x <= fT1) {
      double y = fT1-x;
      double u = fOmega*y-pi;
      double cof = 2*cos(u);
      double a1 = 0;
      double a0 = fA_cdf[1];
      double a2 = 0;
      for (int k = 2; k <= n+1; ++k) {
         a2 = a1;
         a1 = a0;
         a0 = fA_cdf[k]+cof*a1-a2;
      }
      double b1 = 0;
      double b0 = fB_cdf[1];
      for (int k = 2; k <= n; ++k) {
         double b2 = b1;
         b1 = b0;
         b0 = fB_cdf[k]+cof*b1-b2;
      }
      f = -0.5*(a0-a2)+b0*sin(u);
      f += y/fT;
   } else {
      f = 0;
   }
   return f;
}

double VavilovAccurate::Cdf_c (double x, double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Cdf_c (x);
}

double VavilovAccurate::Quantile (double z) const {
   if (z < 0 || z > 1) return std::numeric_limits<double>::signaling_NaN();

   if (!fQuantileInit) InitQuantile();

   double x;
   if (fKappa < 0.02) {
      x = ROOT::Math::landau_quantile (z*(1-2*fEpsilonPM) + fEpsilonPM);
      if (x < fT0+5*fEpsilon) x = fT0+5*fEpsilon;
      else if (x > fT1-10*fEpsilon) x = fT1-10*fEpsilon;
   }
   else {
      // yes, I know what a binary search is, but linear search is faster for small n!
      int i = 1;
      while (z > fQuant[i]) ++i;
      assert (i < fNQuant);

      assert (i >= 1);
      assert (i < fNQuant);

      // initial solution
      double f = (z-fQuant[i-1])/(fQuant[i]-fQuant[i-1]);
      assert (f >= 0);
      assert (f <= 1);
      assert (fQuant[i] > fQuant[i-1]);

      x = f*fLambda[i] + (1-f)*fLambda[i-1];
   }
   if (fabs(x-fT0) < fEpsilon || fabs(x-fT1) < fEpsilon) return x;

   assert (x > fT0 && x < fT1);
   double dx;
   int n = 0;
   do {
      ++n;
      double y = Cdf(x)-z;
      double y1 = Pdf (x);
      dx =  - y/y1;
      x = x + dx;
      // protect against shooting beyond the support
      if (x < fT0) x = 0.5*(fT0+x-dx);
      else if (x > fT1) x = 0.5*(fT1+x-dx);
      assert (x > fT0 && x < fT1);
   } while (fabs(dx) > fEpsilon && n < 100);
   return x;
}

double VavilovAccurate::Quantile (double z, double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Quantile (z);
}

double VavilovAccurate::Quantile_c (double z) const {
   if (z < 0 || z > 1) return std::numeric_limits<double>::signaling_NaN();

   if (!fQuantileInit) InitQuantile();

   double z1 = 1-z;

   double x;
   if (fKappa < 0.02) {
      x = ROOT::Math::landau_quantile (z1*(1-2*fEpsilonPM) + fEpsilonPM);
      if (x < fT0+5*fEpsilon) x = fT0+5*fEpsilon;
      else if (x > fT1-10*fEpsilon) x = fT1-10*fEpsilon;
   }
   else {
      // yes, I know what a binary search is, but linear search is faster for small n!
      int i = 1;
      while (z1 > fQuant[i]) ++i;
      assert (i < fNQuant);

//       int i0=0, i1=fNQuant, i;
//       for (int it = 0; it < LOG2fNQuant; ++it) {
//         i = (i0+i1)/2;
//         if (z > fQuant[i]) i0 = i;
//         else i1 = i;
//       }
//       assert (i1-i0 == 1);

      assert (i >= 1);
      assert (i < fNQuant);

      // initial solution
      double f = (z1-fQuant[i-1])/(fQuant[i]-fQuant[i-1]);
      assert (f >= 0);
      assert (f <= 1);
      assert (fQuant[i] > fQuant[i-1]);

      x = f*fLambda[i] + (1-f)*fLambda[i-1];
   }
   if (fabs(x-fT0) < fEpsilon || fabs(x-fT1) < fEpsilon) return x;

   assert (x > fT0 && x < fT1);
   double dx;
   int n = 0;
   do {
      ++n;
      double y = Cdf_c(x)-z;
      double y1 = -Pdf (x);
      dx =  - y/y1;
      x = x + dx;
      // protect against shooting beyond the support
      if (x < fT0) x = 0.5*(fT0+x-dx);
      else if (x > fT1) x = 0.5*(fT1+x-dx);
      assert (x > fT0 && x < fT1);
   } while (fabs(dx) > fEpsilon && n < 100);
   return x;
}

double VavilovAccurate::Quantile_c (double z, double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Quantile_c (z);
}

VavilovAccurate *VavilovAccurate::GetInstance() {
   if (!fgInstance) fgInstance = new VavilovAccurate (1, 1);
   return fgInstance;
}

VavilovAccurate *VavilovAccurate::GetInstance(double kappa, double beta2) {
   if (!fgInstance) fgInstance = new VavilovAccurate (kappa, beta2);
   else if (kappa != fgInstance->fKappa || beta2 != fgInstance->fBeta2) fgInstance->Set (kappa, beta2);
   return fgInstance;
}

double vavilov_accurate_pdf (double x, double kappa, double beta2) {
   VavilovAccurate *vavilov = VavilovAccurate::GetInstance (kappa, beta2);
   return vavilov->Pdf (x);
}

double vavilov_accurate_cdf_c (double x, double kappa, double beta2) {
   VavilovAccurate *vavilov = VavilovAccurate::GetInstance (kappa, beta2);
   return vavilov->Cdf_c (x);
}

double vavilov_accurate_cdf (double x, double kappa, double beta2) {
   VavilovAccurate *vavilov = VavilovAccurate::GetInstance (kappa, beta2);
   return vavilov->Cdf (x);
}

double vavilov_accurate_quantile (double z, double kappa, double beta2) {
   VavilovAccurate *vavilov = VavilovAccurate::GetInstance (kappa, beta2);
   return vavilov->Quantile (z);
}

double vavilov_accurate_quantile_c (double z, double kappa, double beta2) {
   VavilovAccurate *vavilov = VavilovAccurate::GetInstance (kappa, beta2);
   return vavilov->Quantile_c (z);
}

double VavilovAccurate::G116f1 (double x) const {
   // fH[1] = kappa*(2+beta2*eu) -ln(fEpsilon) + ln(2/pi**2)
   // fH[2] = beta2*kappa
   // fH[3] = omwga/kappa
   // fH[4] = pi/2 *fOmega
   // log of Eq. (4.10)
   return fH[1]+fH[2]*std::log(fH[3]*x)-fH[4]*x;
}

double VavilovAccurate::G116f2 (double x) const {
   // fH[5] = 1-beta2*(1-eu)-logEpsilonPM/kappa;       // eq. 3.9
   // fH[6] = beta2;
   // fH[7] = 1-beta2;
   // Eq. 3.7 of Schorr
//   return fH[5]-x+fH[6]*(std::log(std::fabs(x))-ROOT::Math::expint (-x))-fH[7]*std::exp(-x);
   return fH[5]-x+fH[6]*E1plLog(x)-fH[7]*std::exp(-x);
}

int VavilovAccurate::Rzero (double a, double b, double& x0,
                     double eps, int mxf, double (VavilovAccurate::*f)(double)const) const {

   double xa, xb, fa, fb, r;

   if (a <= b) {
      xa = a;
      xb = b;
   } else {
      xa = b;
      xb = a;
   }
   fa = (this->*f)(xa);
   fb = (this->*f)(xb);

   if(fa*fb > 0) {
      r = -2*(xb-xa);
      x0 = 0;
//      std::cerr << "VavilovAccurate::Rzero: fail=2, f(" << a << ")=" << (this->*f) (a)
//                << ", f(" << b << ")=" << (this->*f) (b) << std::endl;
      return 2;
   }
   int mc = 0;

   bool recalcF12 = true;
   bool recalcFab = true;
   bool fail      = false;


   double x1=0, x2=0, f1=0, f2=0, fx=0, ee=0;
   do {
      if (recalcF12) {
         x0 = 0.5*(xa+xb);
         r = x0-xa;
         ee = eps*(std::abs(x0)+1);
         if(r <= ee) break;
         f1 = fa;
         x1 = xa;
         f2 = fb;
         x2 = xb;
      }
      if (recalcFab) {
         fx = (this->*f)(x0);
         ++mc;
         if(mc > mxf) {
            fail = true;
            break;
         }
         if(fx*fa > 0) {
            xa = x0;
            fa = fx;
         } else {
            xb = x0;
            fb = fx;
         }
      }
      recalcF12 = true;
      recalcFab = true;

      double u1 = f1-f2;
      double u2 = x1-x2;
      double u3 = f2-fx;
      double u4 = x2-x0;
      if(u2 == 0 || u4 == 0) continue;
      double f3 = fx;
      double x3 = x0;
      u1 = u1/u2;
      u2 = u3/u4;
      double ca = u1-u2;
      double cb = (x1+x2)*u2-(x2+x0)*u1;
      double cc = (x1-x0)*f1-x1*(ca*x1+cb);
      if(ca == 0) {
         if(cb == 0) continue;
         x0 = -cc/cb;
      } else {
         u3 = cb/(2*ca);
         u4 = u3*u3-cc/ca;
         if(u4 < 0) continue;
         x0 = -u3 + (x0+u3 >= 0 ? +1 : -1)*std::sqrt(u4);
      }
      if(x0 < xa || x0 > xb) continue;

      recalcF12 = false;

      r = std::abs(x0-x3) < std::abs(x0-x2) ? std::abs(x0-x3) : std::abs(x0-x2);
      ee = eps*(std::abs(x0)+1);
      if (r > ee) {
         f1 = f2;
         x1 = x2;
         f2 = f3;
         x2 = x3;
         continue;
      }

      recalcFab = false;

      fx = (this->*f) (x0);
      if (fx == 0) break;
      double xx, ff;
      if (fx*fa < 0) {
         xx = x0-ee;
         if (xx <= xa) break;
         ff = (this->*f)(xx);
         fb = ff;
         xb = xx;
      } else {
         xx = x0+ee;
         if(xx >= xb) break;
         ff = (this->*f)(xx);
         fa = ff;
         xa = xx;
      }
      if (fx*ff <= 0) break;
      mc += 2;
      if (mc > mxf) {
         fail = true;
         break;
      }
      f1 = f3;
      x1 = x3;
      f2 = fx;
      x2 = x0;
      x0 = xx;
      fx = ff;
   }
   while (true);

   if (fail) {
      r = -0.5*std::abs(xb-xa);
      x0 = 0;
      std::cerr << "VavilovAccurate::Rzero: fail=" << fail << ", f(" << a << ")=" << (this->*f) (a)
                << ", f(" << b << ")=" << (this->*f) (b) << std::endl;
      return 1;
   }

   r = ee;
   return 0;
}

// Calculates log(|x|)+E_1(x)
double VavilovAccurate::E1plLog (double x) {
   static const double eu = 0.577215664901532860606;      // Euler's constant
   double absx = std::fabs(x);
   if (absx < 1E-4) {
      return (x-0.25*x)*x-eu;
   }
   else if (x > 35) {
      return log (x);
   }
   else if (x < -50) {
      return -ROOT::Math::expint (-x);
   }
   return log(absx) -ROOT::Math::expint (-x);
}

double VavilovAccurate::GetLambdaMin() const {
   return fT0;
}

double VavilovAccurate::GetLambdaMax() const {
   return fT1;
}

double VavilovAccurate::GetKappa()     const {
   return fKappa;
}

double VavilovAccurate::GetBeta2()     const {
   return fBeta2;
}

double VavilovAccurate::Mode() const {
   double x = -4.22784335098467134e-01-std::log(fKappa)-fBeta2;
   if (x>-0.223172) x = -0.223172;
   double eps = 0.01;
   double dx;

   do {
      double p0 = Pdf (x - eps);
      double p1 = Pdf (x);
      double p2 = Pdf (x + eps);
      double y1 = 0.5*(p2-p0)/eps;
      double y2 = (p2-2*p1+p0)/(eps*eps);
      dx = - y1/y2;
      x += dx;
      if (fabs(dx) < eps) eps = 0.1*fabs(dx);
   } while (fabs(dx) > 1E-5);
   return x;
}

double VavilovAccurate::Mode(double kappa, double beta2) {
   if (kappa != fKappa || beta2 != fBeta2) Set (kappa, beta2);
   return Mode();
}

double VavilovAccurate::GetEpsilonPM() const {
   return fEpsilonPM;
}

double VavilovAccurate::GetEpsilon()   const {
   return fEpsilon;
}

double VavilovAccurate::GetNTerms()    const {
   return fX0;
}



} // namespace Math
} // namespace ROOT
