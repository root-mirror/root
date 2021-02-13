// @(#)root/mathmore:$Id$
// Author: L. Moneta Wed Oct 18 11:48:00 2006

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


// Header file for class GSLMinimizer

#ifndef ROOT_Math_GSLMinimizer
#define ROOT_Math_GSLMinimizer

#include "Math/Minimizer.h"

#include "Math/IFunctionfwd.h"

#include "Math/IParamFunctionfwd.h"

#include "Math/BasicMinimizer.h"


namespace ROOT {

namespace Math {


   /**
      enumeration specifying the types of GSL minimizers
      @ingroup MultiMin
   */
   enum EGSLMinimizerType {
      kConjugateFR,
      kConjugatePR,
      kVectorBFGS,
      kVectorBFGS2,
      kSteepestDescent
   };


   class GSLMultiMinimizer;

   class MinimTransformFunction;


//_____________________________________________________________________________________
/**
   GSLMinimizer class.
   Implementation of the ROOT::Math::Minimizer interface using the GSL multi-dimensional
   minimization algorithms.

   See <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html">GSL doc</A>
   from more info on the GSL minimization algorithms.

   The class implements the ROOT::Math::Minimizer interface and can be instantiated using the
   ROOT plugin manager (plugin name is "GSLMultiMin"). The varius minimization algorithms
   (conjugatefr, conjugatepr, bfgs, etc..) can be passed as enumerations and also as a string.
   The default algorithm is conjugatefr (Fletcher-Reeves conjugate gradient algorithm).

   @ingroup MultiMin
*/
class GSLMinimizer : public ROOT::Math::BasicMinimizer {

public:

   /**
      Default constructor
   */
   GSLMinimizer (ROOT::Math::EGSLMinimizerType type = ROOT::Math::kConjugateFR  );

   /**
      Constructor with a string giving name of algorithm
    */
   GSLMinimizer (const char *  type  );

   /**
      Destructor
   */
   virtual ~GSLMinimizer ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   GSLMinimizer(const GSLMinimizer &) : BasicMinimizer() {}

   /**
      Assignment operator
   */
   GSLMinimizer & operator = (const GSLMinimizer & rhs) {
      if (this == &rhs) return *this;  // time saving self-test
      return *this;
   }

public:

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func);

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func) { BasicMinimizer::SetFunction(func);}

   /// method to perform the minimization
   virtual  bool Minimize();


   /// return expected distance reached from the minimum
   virtual double Edm() const { return 0; } // not impl. }


   /// return pointer to gradient values at the minimum
   virtual const double *  MinGradient() const;

   /// number of function calls to reach the minimum
   virtual unsigned int NCalls() const;


   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return false; }

   /// return errors at the minimum
   virtual const double * Errors() const {
      return 0;
   }

   /** return covariance matrices elements
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */
   virtual double CovMatrix(unsigned int , unsigned int ) const { return 0; }




protected:

private:


   ROOT::Math::GSLMultiMinimizer * fGSLMultiMin;

   double fLSTolerance;  // Line Search Tolerance

};

   } // end namespace Fit

} // end namespace ROOT



#endif /* ROOT_Math_GSLMinimizer */
