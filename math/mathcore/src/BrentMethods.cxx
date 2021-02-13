// @(#)root/mathcore:$Id$
// Authors: David Gonzalez Maline    01/2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Math/BrentMethods.h"
#include "Math/IFunction.h"
#include "Math/IFunctionfwd.h"

#include <cmath>
#include <algorithm>

#include "Math/Util.h"
#include "Math/Error.h"

#include <iostream>

namespace ROOT {
namespace Math {

namespace BrentMethods {

   double MinimStep(const IGenFunction* function, int type, double &xmin, double &xmax, double fy, int npx,bool logStep)
{
   //   Grid search implementation, used to bracket the minimum and later
   //   use Brent's method with the bracketed interval
   //   The step of the search is set to (xmax-xmin)/npx
   //   type: 0,1 returns MinimumX
   //         2,3 returns MaximumX
   //         4-returns best X corresponding to fy

   if (logStep) {
      xmin = std::log(xmin);
      xmax = std::log(xmax);
   }


   if (npx < 2) return 0.5*(xmax-xmin); // no bracketing - return just mid-point
   double dx = (xmax-xmin)/(npx-1);
   double xxmin = (logStep) ? std::exp(xmin) : xmin;
   double yymin;
   double xmiddle = 0;
   // found an interval containing root (for type=4).
   //  For type !=4 an interval is always found
   bool foundInterval = (type != 4); 
   if (type < 2)
      yymin = (*function)(xxmin);
   else if (type < 4)
      yymin = -(*function)(xxmin);
   else { 
      yymin = (*function)(xxmin)-fy;

      // case root is at the interval boundaries
      if (yymin==0) {
         xmin = xxmin; 
         xmax = xxmin; 
         return xxmin; 
      }
   }
   for (int i=1; i<=npx-1; i++) {
      double x = xmin + i*dx;
      if (logStep) x = std::exp(x);
      double y = 0;
      if (type < 2)
         y = (*function)(x);
      else if (type < 4)
         y = -(*function)(x);
      else
         y = (*function)(x)-fy;
      
      if ( type < 4 && y < yymin) {
         xxmin = x;
         yymin = y;
      }
      // when looking for root break at first instance
      if (type == 4 ) {
         // if root is at interval boundaries
         if (y == 0) {
            xmin = x;
            xmax = x;
            xmiddle = x;
            foundInterval = true;
            break; 
         }
         // found good interval if sign product is negative or equal zero to
         if (std::copysign(1.,y)*std::copysign(1.,yymin) < 0 ) {
            xmin = xxmin; // previous value
            xmax = x;  // current value
            xmiddle = 0.5*(xmax+xmin);
            foundInterval = true; 
            break; 
         }
         // continue bracketing
         xxmin = x;
         yymin = y;
      }
   }
   
   if (type < 4 ) {
      if (logStep) {
         xmin = std::exp(xmin);
         xmax = std::exp(xmax);
      }
      xmin = std::max(xmin,xxmin-dx);
      xmax = std::min(xmax,xxmin+dx);
      xmiddle = std::min(xxmin,xmax);
   } 
   
   //std::cout << "bracketing result " << xmin << "  " << xmax  << "middle value " <<  xmiddle << std::endl;
   if (!foundInterval) {
      MATH_INFO_MSG("BrentMethods::MinimStep", "Grid search failed to find a root in the  interval ");
      std::string msg = "xmin = ";
      msg += ROOT::Math::Util::ToString(xmin);
      msg += std::string(" xmax = ");
      msg += ROOT::Math::Util::ToString(xmax);
      msg += std::string(" npts = ");
      msg += ROOT::Math::Util::ToString(npx);
      MATH_INFO_MSG("BrentMethods::MinimStep", msg.c_str());
      xmin = 1;
      xmax = 0; 
   }

   return xmiddle;
}
   double MinimBrent(const IGenFunction* function, int type, double &xmin, double &xmax, double xmiddle, double fy,  bool &ok, int &niter, double epsabs, double epsrel, int itermax)
{
   //Finds a minimum of a function, if the function is unimodal  between xmin and xmax
   //This method uses a combination of golden section search and parabolic interpolation
   //Details about convergence and properties of this algorithm can be
   //found in the book by R.P.Brent "Algorithms for Minimization Without Derivatives"
   //or in the "Numerical Recipes", chapter 10.2
   // convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
   //
   //type: 0-returns MinimumX
   //      1-returns Minimum
   //      2-returns MaximumX
   //      3-returns Maximum
   //      4-returns X corresponding to fy
   //if ok=true the method has converged

   const double c = 3.81966011250105097e-01; //  (3.-std::sqrt(5.))/2.  (comes from golden section)
   double u, v, w, x, fv, fu, fw, fx, e, p, q, r, t2, d=0, m, tol;
   v = w = x = xmiddle;
   e=0;

   double a=xmin;
   double b=xmax;
   if (type < 2)
      fv = fw = fx = (*function)(x);
   else if (type < 4)
      fv = fw = fx = -(*function)(x);
   else
      fv = fw = fx = std::fabs((*function)(x)-fy);

   for (int i=0; i<itermax; i++){
      m=0.5*(a + b);
      tol = epsrel*(std::fabs(x))+epsabs;
      t2 = 2*tol;

      if (std::fabs(x-m) <= (t2-0.5*(b-a))) {
         //converged, return x
         ok=true;
         niter = i-1;
         if (type==1)
            return fx;
         else if (type==3)
            return -fx;
         else
            return x;
      }

      if (std::fabs(e)>tol){
         //fit parabola
         r = (x-w)*(fx-fv);
         q = (x-v)*(fx-fw);
         p = (x-v)*q - (x-w)*r;
         q = 2*(q-r);
         if (q>0) p=-p;
         else q=-q;
         r=e;   // Deltax before last
         e=d;   // last delta x
         // current Deltax = p/q
         // take a parabolic  step only if:
         // Deltax < 0.5* (DeltaX before last) && Deltax > a && Deltax < b
         // (a BUG in testing this condition is fixed 11/3/2010 (with revision  32544)
         if (std::fabs(p) >= std::fabs(0.5*q*r) || p <= q*(a-x) || p >= q*(b-x)) {
            //  condition fails - do not take parabolic step
            e=(x>=m ? a-x : b-x);
            d = c*e;
         } else {
            // take a parabolic interpolation step
            d = p/q;
            u = x+d;
            if (u-a < t2 || b-u < t2)
               //d=TMath::Sign(tol, m-x);
               d=(m-x >= 0) ? std::fabs(tol) : -std::fabs(tol);
         }
      } else {
         e=(x>=m ? a-x : b-x);
         d = c*e;
      }
      u = (std::fabs(d)>=tol ? x+d : x+ ((d >= 0) ? std::fabs(tol) : -std::fabs(tol)) );
      if (type < 2)
         fu = (*function)(u);
      else if (type < 4)
         fu = -(*function)(u);
      else
         fu = std::fabs((*function)(u)-fy);
      //update a, b, v, w and x
      if (fu<=fx){
         if (u<x) b=x;
         else a=x;
         v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
      } else {
         if (u<x) a=u;
         else b=u;
         if (fu<=fw || w==x){
            v=w; fv=fw; w=u; fw=fu;
         }
         else if (fu<=fv || v==x || v==w){
            v=u; fv=fu;
         }
      }
   }
   //didn't converge
   ok = false;
   xmin = a;
   xmax = b;
   niter = itermax;
   return x;

}

} // end namespace BrentMethods
} // end namespace Math
} // ned namespace ROOT
