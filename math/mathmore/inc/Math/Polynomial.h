// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

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

// Header file for class Polynomial
//
// Created by: Lorenzo Moneta  at Wed Nov 10 17:46:19 2004
//
// Last update: Wed Nov 10 17:46:19 2004
//
#ifndef ROOT_Math_Polynomial
#define ROOT_Math_Polynomial

#include <complex>
#include <vector>

#include "Math/ParamFunction.h"

// #ifdef _WIN32
// #pragma warning(disable : 4250)
// #endif

namespace ROOT {
namespace Math {

//_____________________________________________________________________________________
  /**
     Parametric Function class describing polynomials of order n.

     <em>P(x) = p[0] + p[1]*x + p[2]*x**2 + ....... + p[n]*x**n</em>

     The class implements also the derivatives, \a dP(x)/dx and the \a dP(x)/dp[i].

     The class provides also the method to find the roots of the polynomial.
     It uses analytical methods up to quartic polynomials.

     Implements both the Parameteric function interface and the gradient interface
     since it provides the analytical gradient with respect to x


     @ingroup ParamFunc
  */

class Polynomial : public ParamFunction<IParamGradFunction>,
                   public IGradientOneDim
{


public:

 typedef  ParamFunction<IParamGradFunction> ParFunc;
   /**
      Construct a Polynomial function of order n.
      The number of Parameters is n+1.
   */

   Polynomial(unsigned int n = 0);

   /**
      Construct a Polynomial of degree  1 : a*x + b
   */
   Polynomial(double a, double b);

   /**
      Construct a Polynomial of degree  2 : a*x**2 + b*x + c
   */
   Polynomial(double a, double b, double c);

   /**
      Construct a Polynomial of degree  3 : a*x**3 + b*x**2 + c*x + d
   */
   Polynomial(double a, double b, double c, double d);

   /**
      Construct a Polynomial of degree  4 : a*x**4 + b*x**3 + c*x**2 + dx  + e
   */
   Polynomial(double a, double b, double c, double d, double e);


   virtual ~Polynomial() {}

   // use default copy-ctor and assignment operators



//   using ParamFunction::operator();


   /**
      Find the polynomial roots.
      For n <= 4, the roots are found analytically while for larger order an iterative numerical method is used
      The numerical method used is from GSL (see <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_6.html#SEC53" )
   */
   const std::vector<std::complex <double> > & FindRoots();


   /**
      Find the only the real polynomial roots.
      For n <= 4, the roots are found analytically while for larger order an iterative numerical method is used
      The numerical method used is from GSL (see <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_6.html#SEC53" )
   */
   std::vector<double > FindRealRoots();


   /**
      Find the polynomial roots using always an iterative numerical methods
      The numerical method used is from GSL (see <A HREF="http://www.gnu.org/software/gsl/manual/gsl-ref_6.html#SEC53" )
   */
   const std::vector<std::complex <double> > & FindNumRoots();

   /**
      Order of Polynomial
   */
   unsigned int Order() const { return fOrder; }


   IGenFunction * Clone() const;

   /**
       Optimized method to evaluate at the same time the function value and derivative at a point x.
       Implement the interface specified bby ROOT::Math::IGradientOneDim.
       In the case of polynomial there is no advantage to compute both at the same time
   */
   void FdF (double x, double & f, double & df) const {
      f = (*this)(x);
      df = Derivative(x);
   }


private:

   double DoEvalPar ( double x, const double * p ) const ;

   double DoDerivative (double x) const ;

   double DoParameterDerivative(double x, const double * p, unsigned int ipar) const;


   // cache order = number of params - 1)
   unsigned int fOrder;

   // cache Parameters for Gradient
   mutable std::vector<double> fDerived_params;

   // roots

   std::vector< std::complex < double > > fRoots;

};

} // namespace Math
} // namespace ROOT


#endif /* ROOT_Math_Polynomial */
