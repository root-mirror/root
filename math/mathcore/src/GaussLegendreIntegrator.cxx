// @(#)root/mathcore:$Id$
// Authors: David Gonzalez Maline    01/2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Math/GaussLegendreIntegrator.h"
#include "Math/Error.h"
#include "Math/IFunction.h"
#include "Math/IFunctionfwd.h"
#include "Math/IntegratorOptions.h"
#include <cmath>
#include <string.h>
#include <algorithm>

namespace ROOT {
namespace Math {

   GaussLegendreIntegrator::GaussLegendreIntegrator(int num, double eps) :
      GaussIntegrator(eps, eps)
{
   // Basic contructor
   fNum = num;
   fX = 0;
   fW = 0;

   CalcGaussLegendreSamplingPoints();
}

GaussLegendreIntegrator::~GaussLegendreIntegrator()
{
   // Default Destructor


   delete [] fX;
   delete [] fW;
}

void GaussLegendreIntegrator::SetNumberPoints(int num)
{
   // Set the number of points used in the calculation of the integral

   fNum = num;
   CalcGaussLegendreSamplingPoints();
}

void GaussLegendreIntegrator::GetWeightVectors(double *x, double *w) const
{
   // Returns the arrays x and w.

   std::copy(fX,fX+fNum, x);
   std::copy(fW,fW+fNum, w);
}


double GaussLegendreIntegrator::DoIntegral(double a, double b, const IGenFunction* function)
{
   // Gauss-Legendre integral, see CalcGaussLegendreSamplingPoints.

   if (fNum<=0 || fX == 0 || fW == 0)
      return 0;

   fUsedOnce = true;

   const double a0 = (b + a)/2;
   const double b0 = (b - a)/2;

   double xx[1];

   double result = 0.0;
   for (int i=0; i<fNum; i++)
   {
      xx[0] = a0 + b0*fX[i];
      result += fW[i] * (*function)(xx);
   }

   fLastResult = result*b0;
   return fLastResult;
}


void GaussLegendreIntegrator::SetRelTolerance (double eps)
{
   // Set the desired relative Error.
   fEpsRel = eps;
   CalcGaussLegendreSamplingPoints();
}

void GaussLegendreIntegrator::SetAbsTolerance (double)
{   MATH_WARN_MSG("ROOT::Math::GaussLegendreIntegrator", "There is no Absolute Tolerance!");  }



void GaussLegendreIntegrator::CalcGaussLegendreSamplingPoints()
{
   // Given the number of sampling points this routine fills the
   // arrays x and w.

   if (fNum<=0 || fEpsRel<=0)
      return;

   if ( fX )
      delete [] fX;

   if ( fW )
      delete [] fW;

   fX = new double[fNum];
   fW = new double[fNum];

   // The roots of symmetric is the interval, so we only have to find half of them
   const unsigned int m = (fNum+1)/2;

   double z, pp, p1,p2, p3;

   // Loop over the desired roots
   for (unsigned int i=0; i<m; i++) {
      z = std::cos(3.14159265358979323846*(i+0.75)/(fNum+0.5));

      // Starting with the above approximation to the i-th root, we enter
      // the main loop of refinement by Newton's method
      do {
         p1=1.0;
         p2=0.0;

         // Loop up the recurrence relation to get the Legendre
         // polynomial evaluated at z
         for (int j=0; j<fNum; j++)
         {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1.0);
         }
         // p1 is now the desired Legendre polynomial. We next compute pp, its
         // derivative, by a standard relation involving also p2, the polynomial
         // of one lower order
         pp = fNum*(z*p1-p2)/(z*z-1.0);
         // Newton's method
         z -= p1/pp;

      } while (std::fabs(p1/pp) > fEpsRel);

      // Put root and its symmetric counterpart
      fX[i]       = -z;
      fX[fNum-i-1] =  z;

      // Compute the weight and put its symmetric counterpart
      fW[i]       = 2.0/((1.0-z*z)*pp*pp);
      fW[fNum-i-1] = fW[i];
   }
}

ROOT::Math::IntegratorOneDimOptions  GaussLegendreIntegrator::Options() const {
   ROOT::Math::IntegratorOneDimOptions opt;
   opt.SetAbsTolerance(0);
   opt.SetRelTolerance(fEpsRel);
   opt.SetWKSize(0);
   opt.SetNPoints(fNum);
   opt.SetIntegrator("GaussLegendre");
   return opt;
}

void GaussLegendreIntegrator::SetOptions(const ROOT::Math::IntegratorOneDimOptions & opt)
{
   //   set integration options
//    std::cout << "fEpsilon = " << fEpsilon << std::endl;
//    std::cout << opt.RelTolerance() << " abs " << opt.AbsTolerance() << std::endl;
   //double tol = opt.RelTolerance(); fEpsilon = tol;
   fEpsRel = opt.RelTolerance();
//    std::cout << "fEpsilon = " << fEpsilon << std::endl;
   fNum = opt.NPoints();
   if (fNum <= 7)  MATH_WARN_MSGVAL("GaussLegendreIntegrator::SetOptions","setting a low number of points ",fNum);
   CalcGaussLegendreSamplingPoints();
}

} // end namespace Math
} // end namespace ROOT
