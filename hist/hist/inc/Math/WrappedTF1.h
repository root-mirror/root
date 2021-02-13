// @(#)root/mathmore:$Id$
// Author: L. Moneta Wed Sep  6 09:52:26 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class WrappedTF1

#ifndef ROOT_Math_WrappedTF1
#define ROOT_Math_WrappedTF1


#include "Math/IParamFunction.h"

#include "TF1.h"
#include <string>
#include <vector>

namespace ROOT {

   namespace Math {


      /**
         Class to Wrap a ROOT Function class (like TF1)  in a IParamFunction interface
         of one dimensions to be used in the ROOT::Math numerical algorithms
         The wrapper does not own bby default the TF1 pointer, so it assumes it exists during the wrapper lifetime

         The class from ROOT version 6.03  does not contain anymore a copy of the parameters. The parameters are
         stored in the TF1 class.


         @ingroup CppFunctions
      */
      class WrappedTF1 : public ROOT::Math::IParamGradFunction, public ROOT::Math::IGradientOneDim {

      public:

         typedef  ROOT::Math::IGradientOneDim     IGrad;
         typedef  ROOT::Math::IParamGradFunction  BaseGradFunc;
         typedef  ROOT::Math::IParamGradFunction::BaseFunc BaseFunc;


         /**
            constructor from a TF1 function pointer.
         */
         WrappedTF1(TF1 &f);

         /**
            Destructor (no operations). TF1 Function pointer is not owned
         */
         virtual ~WrappedTF1() {}

         /**
            Copy constructor
         */
         WrappedTF1(const WrappedTF1 &rhs);

         /**
            Assignment operator
         */
         WrappedTF1 &operator = (const WrappedTF1 &rhs);

         /** @name interface inherited from IFunction */

         /**
             Clone the wrapper but not the original function
         */
         ROOT::Math::IGenFunction *Clone() const
         {
            return  new WrappedTF1(*this);
         }


         /** @name interface inherited from IParamFunction */

         /// get the parameter values (return values cachen inside, those inside TF1 might be different)
         const double *Parameters() const
         {
            //return  (fParams.size() > 0) ? &fParams.front() : 0;
            return fFunc->GetParameters();
         }

         /// set parameter values
         /// need to call also SetParameters in TF1 in ace some other operations (re-normalizations) are needed
         void SetParameters(const double *p)
         {
            //std::copy(p,p+fParams.size(),fParams.begin());
            fFunc->SetParameters(p);
         }

         /// return number of parameters
         unsigned int NPar() const
         {
            //return fParams.size();
            return fFunc->GetNpar();
         }

         /// return parameter name (this is stored in TF1)
         std::string ParameterName(unsigned int i) const
         {
            return std::string(fFunc->GetParName(i));
         }


         using BaseGradFunc::operator();

         /// evaluate the derivative of the function with respect to the parameters
         void  ParameterGradient(double x, const double *par, double *grad) const;

         /// calculate function and derivative at same time (required by IGradient interface)
         void FdF(double x, double &f, double &deriv) const
         {
            f = DoEval(x);
            deriv = DoDerivative(x);
         }

         /// precision value used for calculating the derivative step-size
         /// h = eps * |x|. The default is 0.001, give a smaller in case function changes rapidly
         static void SetDerivPrecision(double eps);

         /// get precision value used for calculating the derivative step-size
         static double GetDerivPrecision();

      private:


         /// evaluate function passing coordinates x and vector of parameters
         double DoEvalPar(double x, const double *p) const
         {
            fX[0] = x;
            if (fFunc->GetMethodCall()) fFunc->InitArgs(fX, p);  // needed for interpreted functions
            return fFunc->EvalPar(fX, p);
         }

         /// evaluate function using the cached parameter values (of TF1)
         /// re-implement for better efficiency
         double DoEval(double x) const
         {
            // no need to call InitArg for interpreted functions (done in ctor)
            // use EvalPar since it is much more efficient than Eval
            fX[0] = x;
            //const double * p = (fParams.size() > 0) ? &fParams.front() : 0;
            return fFunc->EvalPar(fX, 0);
         }

         /// return the function derivatives w.r.t. x
         double DoDerivative(double  x) const;

         /// evaluate the derivative of the function with respect to the parameters
         double  DoParameterDerivative(double x, const double *p, unsigned int ipar) const;

         bool fLinear;                 // flag for linear functions
         bool fPolynomial;             // flag for polynomial functions
         TF1 *fFunc;                   // pointer to ROOT function
         mutable double fX[1];         //! cached vector for x value (needed for TF1::EvalPar signature)
         //std::vector<double> fParams;  //  cached vector with parameter values

      };

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_WrappedTF1 */
