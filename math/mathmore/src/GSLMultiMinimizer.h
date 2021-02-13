// @(#)root/mathmore:$Id$
// Author: L. Moneta Tue Dec 19 14:09:15 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
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

// Header file for class GSLMultiMinimizer

#ifndef ROOT_Math_GSLMultiMinimizer
#define ROOT_Math_GSLMultiMinimizer

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_version.h"
#include "GSLMultiMinFunctionWrapper.h"

#include "Math/Error.h"

#include "Math/IFunction.h"

#include <cassert>
#include <string>

namespace ROOT {

   namespace Math {


/**
   GSLMultiMinimizer class , for minimizing multi-dimensional function
   using derivatives

   @ingroup MultiMin

*/
class GSLMultiMinimizer {

public:

   /**
      Default constructor
   */
   GSLMultiMinimizer (ROOT::Math::EGSLMinimizerType type)  :
      fMinimizer(0),
      fType(0),
      fVec(0)
   {
      switch(type)
      {
      case ROOT::Math::kConjugateFR :
         fType = gsl_multimin_fdfminimizer_conjugate_fr;
         break;
      case ROOT::Math::kConjugatePR :
         fType = gsl_multimin_fdfminimizer_conjugate_pr;
         break;
      case ROOT::Math::kVectorBFGS :
         fType = gsl_multimin_fdfminimizer_vector_bfgs;
         break;
      case ROOT::Math::kVectorBFGS2 :
#if (GSL_MAJOR_VERSION > 1) || ((GSL_MAJOR_VERSION == 1) && (GSL_MINOR_VERSION >= 9))
         // bfgs2 is available only for v>= 1.9
         fType = gsl_multimin_fdfminimizer_vector_bfgs2;
#else
         MATH_INFO_MSG("GSLMultiMinimizer","minimizer BFSG2 does not exist with this GSL version , use BFGS");
         fType = gsl_multimin_fdfminimizer_vector_bfgs;
#endif
         break;
      case ROOT::Math::kSteepestDescent:
         fType = gsl_multimin_fdfminimizer_steepest_descent;
         break;
      default:
         fType = gsl_multimin_fdfminimizer_conjugate_fr;
         break;
      }

   }

   /**
      Destructor
   */
   ~GSLMultiMinimizer ()  {
      if (fMinimizer != 0 ) gsl_multimin_fdfminimizer_free(fMinimizer);
      // can free vector (is copied inside)
      if (fVec != 0) gsl_vector_free(fVec);
   }

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   GSLMultiMinimizer(const GSLMultiMinimizer &) {}

   /**
      Assignment operator
   */
   GSLMultiMinimizer & operator = (const GSLMultiMinimizer & rhs)  {
      if (this == &rhs) return *this;  // time saving self-test
      return *this;
   }

public:

   /**
      set the function to be minimize the initial minimizer parameters,
      step size and tolerance in the line search
    */
   int Set(const ROOT::Math::IMultiGradFunction & func, const double * x, double stepSize, double tol) {
      // create function wrapper
      fFunc.SetFunction(func);
      // create minimizer object (free previous one if already existing)
      unsigned int ndim = func.NDim();
      CreateMinimizer( ndim );
      // set initial values
      if (fVec != 0) gsl_vector_free(fVec);
      fVec = gsl_vector_alloc( ndim );
      std::copy(x,x+ndim, fVec->data);
      assert(fMinimizer != 0);
      return gsl_multimin_fdfminimizer_set(fMinimizer, fFunc.GetFunc(), fVec, stepSize, tol);
   }

   /// create the minimizer from the type and size
   void CreateMinimizer(unsigned int n) {
      if (fMinimizer) gsl_multimin_fdfminimizer_free(fMinimizer);
      fMinimizer = gsl_multimin_fdfminimizer_alloc(fType, n);
   }

   std::string Name() const {
      if (fMinimizer == 0) return "undefined";
      return std::string(gsl_multimin_fdfminimizer_name(fMinimizer) );
   }

   int Iterate() {
      if (fMinimizer == 0) return -1;
      return gsl_multimin_fdfminimizer_iterate(fMinimizer);
   }

   /// x values at the minimum
   double * X() const {
      if (fMinimizer == 0) return 0;
      gsl_vector * x =  gsl_multimin_fdfminimizer_x(fMinimizer);
      return x->data;
   }

   /// function value at the minimum
   double Minimum() const {
      if (fMinimizer == 0) return 0;
      return gsl_multimin_fdfminimizer_minimum(fMinimizer);
   }

   /// gradient value at the minimum
   double * Gradient() const {
      if (fMinimizer == 0) return 0;
      gsl_vector * g =  gsl_multimin_fdfminimizer_gradient(fMinimizer);
      return g->data;
   }

   /// restart minimization from current point
   int Restart() {
      if (fMinimizer == 0) return -1;
      return gsl_multimin_fdfminimizer_restart(fMinimizer);
   }

   /// test gradient (ask from minimizer gradient vector)
   int TestGradient(double absTol) const {
      if (fMinimizer == 0) return -1;
      gsl_vector * g =  gsl_multimin_fdfminimizer_gradient(fMinimizer);
      return gsl_multimin_test_gradient( g, absTol);
   }

   /// test gradient (require a vector gradient)
   int TestGradient(const double * g, double absTol) const {
      if (fVec == 0 ) return -1;
      unsigned int n = fVec->size;
      if (n == 0 ) return -1;
      std::copy(g,g+n, fVec->data);
      return gsl_multimin_test_gradient( fVec, absTol);
   }


private:

   gsl_multimin_fdfminimizer * fMinimizer;
   GSLMultiMinDerivFunctionWrapper fFunc;
   const gsl_multimin_fdfminimizer_type * fType;
   // cached vector to avoid re-allocating every time a new one
   mutable gsl_vector * fVec;

};

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_GSLMultiMinimizer */
