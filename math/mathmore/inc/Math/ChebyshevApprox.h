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

// Header file for class ChebyshevApprox
//
// Created by: moneta  at Thu Dec  2 14:51:15 2004
//
// Last update: Thu Dec  2 14:51:15 2004
//
#ifndef ROOT_Math_ChebyshevApprox
#define ROOT_Math_ChebyshevApprox



/**
   @defgroup FuncApprox Function Approximation (ChebyshevApprox)
   Numerical algorithm from the \ref MathMore library and implemented using the
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/">GSL</A> library
 
   @ingroup NumAlgo
 */


#include "Math/IFunctionfwd.h"

#include "Math/GSLFunctionAdapter.h"

#include <utility>
#include <cstddef>


namespace ROOT {
namespace Math {

class GSLChebSeries;
class GSLFunctionWrapper;

//____________________________________________________________________________
/**
   Class describing a Chebyshev series which can be used to approximate a
   function in a defined range [a,b] using Chebyshev polynomials.
   It uses the algorithm from
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Chebyshev-Approximations.html">GSL</A>

   This class does not support copying
   @ingroup FuncApprox
 */


class ChebyshevApprox {

public:


   /**
      Construct a Chebyshev series approximation to a Function f in range [a,b];
      constructor based on functions of type IGenFunction
   */

  ChebyshevApprox(const ROOT::Math::IGenFunction & f, double a, double b, size_t n);

   /**
      Construct a Chebyshev series approximation to a Function f in range [a,b];
      constructor based on free functions with gsl_function type signature
   */
   ChebyshevApprox(GSLFuncPointer f, void *p, double a, double b, size_t n);

   // destructor
   virtual ~ChebyshevApprox();


private:

   /**
      construct a Chebyshev series or order n
      The series must be initialized from a function
   */
   ChebyshevApprox(size_t n);

// usually copying is non trivial, so we make this unaccessible
   ChebyshevApprox(const ChebyshevApprox &);
   ChebyshevApprox & operator = (const ChebyshevApprox &);

public:

   /**
       Evaluate the series at a given point x
   */
   double operator() ( double x) const;

   /**
      Evaluate the series at a given point x estimating both the series result and its absolute error.
      The error estimate is made from the first neglected term in the series.
      A pair containing result and error is returned
   */
   std::pair<double, double>  EvalErr( double x) const;

   /**
      Evaluate the series at a given point, to (at most) the given order n
   */
   double operator() ( double x, size_t n) const;

   /**
      evaluate the series at a given point x to the given order n,
      estimating both the series result and its absolute error.
      The error estimate is made from the first neglected term in the series.
      A pair containing result and error is returned
   */
   std::pair<double, double>  EvalErr( double x, size_t n) const;

   /**
      Compute the derivative of the series and return a pointer to a new Chebyshev series with the
      derivatives coefficients. The returned pointer must be managed by the user.
   */
   //TO DO: implement copying to return by value
   ChebyshevApprox * Deriv();

   /**
      Compute the integral of the series and return a pointer to a new Chebyshev series with the
      integral coefficients. The lower limit of the integration is the left range value a.
      The returned pointer must be managed by the user
   */
   //TO DO: implement copying to return by value
   ChebyshevApprox * Integral();

protected:

   /**
       Initialize series passing function and range
   */
   void Initialize( GSLFuncPointer f, void * params, double a, double b);

private:

   size_t fOrder;

   GSLChebSeries * fSeries;
   GSLFunctionWrapper * fFunction;     // pointer to function

};

} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_ChebyshevApprox */
