// @(#)root/mathmore:$Id$
// Author: L. Moneta Wed Dec 20 17:16:32 2006

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

// Header file for class GSLNLSMinimizer

#ifndef ROOT_Math_GSLNLSMinimizer
#define ROOT_Math_GSLNLSMinimizer



#include "Math/BasicMinimizer.h"

#include "Math/IFunctionfwd.h"

#include "Math/IParamFunctionfwd.h"

#include "Math/FitMethodFunction.h"

#include "Math/MinimTransformVariable.h"

#include <vector>

namespace ROOT {

   namespace Math {

      class GSLMultiFit;


//________________________________________________________________________________
/**
    LSResidualFunc class description.
    Internal class used for accessing the residuals of the Least Square function
    and their derivates which are estimated numerically using GSL numerical derivation.
    The class contains a pointer to the fit method function and an index specifying
    the i-th residual and wraps it in a multi-dim gradient function interface
    ROOT::Math::IGradientFunctionMultiDim.
    The class is used by ROOT::Math::GSLNLSMinimizer (GSL non linear least square fitter)

    @ingroup MultiMin
*/
class LSResidualFunc : public IMultiGradFunction {
public:

   //default ctor (required by CINT)
   LSResidualFunc() : fIndex(0), fChi2(0)
   {}


   LSResidualFunc(const ROOT::Math::FitMethodFunction & func, unsigned int i) :
      fIndex(i),
      fChi2(&func),
      fX2(std::vector<double>(func.NDim() ) )
   {}


   // copy ctor
   LSResidualFunc(const LSResidualFunc & rhs) :
      IMultiGenFunction(),
      IMultiGradFunction()
   {
      operator=(rhs);
   }

   // assignment
   LSResidualFunc & operator= (const LSResidualFunc & rhs)
   {
      fIndex = rhs.fIndex;
      fChi2 = rhs.fChi2;
      fX2 = rhs.fX2;
      return *this;
   }

   IMultiGenFunction * Clone() const {
      return new LSResidualFunc(*fChi2,fIndex);
   }

   unsigned int NDim() const { return fChi2->NDim(); }

   void Gradient( const double * x, double * g) const {
      double f0 = 0;
      FdF(x,f0,g);
   }

   void FdF (const double * x, double & f, double * g) const {
      unsigned int n = NDim();
      std::copy(x,x+n,fX2.begin());
      const double kEps = 1.0E-4;
      f = DoEval(x);
      for (unsigned int i = 0; i < n; ++i) {
         fX2[i] += kEps;
         g[i] =  ( DoEval(&fX2.front()) - f )/kEps;
         fX2[i] = x[i];
      }
   }


private:

   double DoEval (const double * x) const {
      return fChi2->DataElement(x, fIndex);
   }

   double DoDerivative(const double * x, unsigned int icoord) const {
      //return  ROOT::Math::Derivator::Eval(*this, x, icoord, 1E-8);
      std::copy(x,x+NDim(),fX2.begin());
      const double kEps = 1.0E-4;
      fX2[icoord] += kEps;
      return ( DoEval(&fX2.front()) - DoEval(x) )/kEps;
   }

   unsigned int fIndex;
   const ROOT::Math::FitMethodFunction * fChi2;
   mutable std::vector<double> fX2;  // cached vector
};


//_____________________________________________________________________________________________________
/**
   GSLNLSMinimizer class for Non Linear Least Square fitting
   It Uses the Levemberg-Marquardt algorithm from
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html">
   GSL Non Linear Least Square fitting</A>.

   @ingroup MultiMin
*/
class GSLNLSMinimizer : public  ROOT::Math::BasicMinimizer {

public:

   /**
      Default constructor
   */
   GSLNLSMinimizer (int type = 0);

   /**
      Destructor (no operations)
   */
   ~GSLNLSMinimizer ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   GSLNLSMinimizer(const GSLNLSMinimizer &) : ROOT::Math::BasicMinimizer() {}

   /**
      Assignment operator
   */
   GSLNLSMinimizer & operator = (const GSLNLSMinimizer & rhs)  {
      if (this == &rhs) return *this;  // time saving self-test
      return *this;
   }

public:

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func);

   /// set gradient the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func);


   /// method to perform the minimization
   virtual  bool Minimize();


   /// return expected distance reached from the minimum
   virtual double Edm() const { return fEdm; } // not impl. }


   /// return pointer to gradient values at the minimum
   virtual const double *  MinGradient() const;

   /// number of function calls to reach the minimum
   virtual unsigned int NCalls() const { return (fChi2Func) ? fChi2Func->NCalls() : 0; }

   /// number of free variables (real dimension of the problem)
   /// this is <= Function().NDim() which is the total
//   virtual unsigned int NFree() const { return fNFree; }

   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return true; }

   /// return errors at the minimum
   virtual const double * Errors() const { return (fErrors.size() > 0) ? &fErrors.front() : 0; }
//  {
//       static std::vector<double> err;
//       err.resize(fDim);
//       return &err.front();
//    }

   /** return covariance matrices elements
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */
   virtual double CovMatrix(unsigned int , unsigned int ) const;

   /// return covariance matrix status
   virtual int CovMatrixStatus() const;

protected:


private:

   unsigned int fNFree;      // dimension of the internal function to be minimized
   unsigned int fSize;        // number of fit points (residuals)

   ROOT::Math::GSLMultiFit * fGSLMultiFit;        // pointer to GSL multi fit solver
   const ROOT::Math::FitMethodFunction * fChi2Func;      // pointer to Least square function

   double fEdm;                                   // edm value
   double fLSTolerance;                           // Line Search Tolerance
   std::vector<double> fErrors;
   std::vector<double> fCovMatrix;              //  cov matrix (stored as cov[ i * dim + j]
   std::vector<LSResidualFunc> fResiduals;   //! transient Vector of the residual functions



};

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_GSLNLSMinimizer */
