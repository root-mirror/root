// @(#)root/tmva $Id$
// Author: Simon Pfreundschuh

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////
// Generic tests for the derivatives and gradiens of acitvation,    //
// loss and regularization functions. Each function generates a     //
// random 10 x 10 matrix and uses a central finite difference and   //
// to numerically compute the derivative of the function            //
// w.r.t. this element. The result is compared to the result        //
// obtained by the corresponding analytic derivative implemented by //
// the evaluateDerivative(...), evaluateGradients(...),             //
// addRegularizationGradients(...) functions.                       //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/Net.h"
#include "Utility.h"

using namespace TMVA::DNN;

//______________________________________________________________________________
//
//  Activation Functions
//______________________________________________________________________________

/*! Generic function that numerically computes the derivative of a matrix
 *  function f and the analytical solution given by df the function signatures
 *  are assumed to be
 *  - void f(Matrix_t &X)
 *  - void df(Matrix_t &Y, const Matrix_t &X) -> derivative of f at X(i,j) is
 *  The function f is supposed to apply the corresponding mathematical function
 *  to each element in the provided matrix X. The function df is expected to
 *  set each element in Y to the derivative of the corresponding mathematical
 *  function evaluated at the corresponding element in X.
 */
template<typename Architecture, typename F, typename dF>
    auto testDerivatives(F f, dF df,
                         typename Architecture::Scalar_t dx)
    -> typename Architecture::Scalar_t
{
   using Scalar_t   = typename Architecture::Scalar_t;
   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;

   Scalar_t maximum_error = 0.0;

   for (size_t i = 0; i < 100; i++)
   {
      Matrix_t mX(10,10), mY(10,10);
      randomMatrix(mY);
      
      Tensor_t X(mX,2); 
      Tensor_t Y(mY,2); 

      df(X, Y);
      Scalar_t dy = X(0,0);

      // copy from 
      Architecture::Copy(X, Y);
      X(0,0) += dx;
      f(X);
      Scalar_t y1 = X(0,0);
      Architecture::Copy(X, Y);
      X(0,0) -= dx;
      f(X);
      Scalar_t y0 = X(0,0);
      Scalar_t dy_num = (y1 - y0) / (2.0 * dx);
      Scalar_t error = relativeError(dy_num, dy);
      maximum_error = std::max(maximum_error, error);
   }

   return maximum_error;
}

/*! Test derivatives of all activation functions and return the maximum relative
 *  error. Prints the result for each function to the stdout. */
//______________________________________________________________________________
template<typename Architecture>
auto testActivationFunctionDerivatives(bool useFastTanh = false)
    -> typename Architecture::Scalar_t
{
   using Scalar_t   = typename Architecture::Scalar_t;
   //using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;

   // Test only differentiable activation functions.
   std::vector<EActivationFunction> EActivationFunctions
   = {EActivationFunction::kIdentity,
      EActivationFunction::kSigmoid,
      EActivationFunction::kTanh,
      EActivationFunction::kSoftSign,
      EActivationFunction::kGauss};

   Scalar_t error, maximum_error;
   maximum_error = 0.0;

   for (auto & af : EActivationFunctions)
   {
      auto f = [&af](Tensor_t &X) { evaluate<Architecture>(X, af); };
      auto df = [& af](Tensor_t &X, const Tensor_t &Y)
      {
         evaluateDerivative<Architecture>(X, af, Y);
      };

      auto h = std::sqrt(std::numeric_limits<Scalar_t>::epsilon());
      // in case of tanh and using VDT h (derivative step size  must be much larger)
      if (useFastTanh && af == EActivationFunction::kTanh ) h = 1.E-3; 

      error = testDerivatives<Architecture>(f, df, h);

      std::cout << "Testing " << static_cast<int>(af) << ": ";
      std::cout << "Maximum Relative Error = " << error << std::endl;

      maximum_error = std::max(maximum_error, error);
   }

   return maximum_error;
}

//______________________________________________________________________________
//
//  Loss functions.
//______________________________________________________________________________

/*! Similar to testDerivatives only that here the mathematical function is
 *  expected to be a matrix functional, i.e. to be mapping a matrix to a
 *  scalar value. The scalar value is supposed to be computed by the provided
 *  function object f, while the function object is just like above. */
template<typename Architecture, typename F, typename dF>
    auto testGradients(F f, dF df,
                       typename Architecture::Scalar_t dx)
    -> typename Architecture::Scalar_t
{
    using Scalar_t   = typename Architecture::Scalar_t;
    using Matrix_t = typename Architecture::Matrix_t;

    Scalar_t maximum_error = 0.0;

    for (size_t i = 0; i < 100; i++)
    {
       Matrix_t X(10, 10), Y(10, 10), Z(10, 10), W(10, 10);
       randomMatrix(X);
       randomMatrix(Y);
       randomMatrix(W);

       df(Z, Y, X, W);
       Scalar_t dy = Z(0, 0);

       X(0, 0) += dx;
       Scalar_t y1 = f(Y, X, W);
       X(0, 0) -= 2.0 * dx;
       Scalar_t y0 = f(Y, X, W);
       Scalar_t dy_num = (y1 - y0) / (2.0 * dx);

       Scalar_t error = relativeError(dy_num, dy);
       maximum_error = std::max(maximum_error, error);
    }

    return maximum_error;
}

/*! Test gradients of all loss function for the given architecture type and
 *  return the maximum relative error. Prints results for each function to
 *  standard out. */
//______________________________________________________________________________
template<typename Architecture>
auto testLossFunctionGradients()
    -> typename Architecture::Scalar_t
{
    using Scalar_t   = typename Architecture::Scalar_t;
    using Matrix_t = typename Architecture::Matrix_t;

    std::vector<ELossFunction> LossFunctions
        = {ELossFunction::kMeanSquaredError,
           ELossFunction::kCrossEntropy,
           ELossFunction::kSoftmaxCrossEntropy};

    Scalar_t error, maximum_error;
    maximum_error = 0.0;

    for (auto & lf : LossFunctions)
    {
       auto f = [lf](const Matrix_t &Y, const Matrix_t &Z, const Matrix_t &W) {
          return evaluate<Architecture>(lf, Y, Z, W);
       };
       auto df = [&lf](Matrix_t &X, const Matrix_t &Y, const Matrix_t &Z, const Matrix_t &W) {
          evaluateGradients<Architecture>(X, lf, Y, Z, W);
       };

       auto h = 100.0 * std::sqrt(std::numeric_limits<Scalar_t>::epsilon());
       error = testGradients<Architecture>(f, df, h);

       std::cout << "Testing " << static_cast<char>(lf) << ": ";
       std::cout << "Maximum Relative Error = " << error << std::endl;

       maximum_error = std::max(maximum_error, error);
    }

    return maximum_error;
}

//______________________________________________________________________________
//
//  Regularization.
//______________________________________________________________________________

/*! Test the computation of gradients for all differentiable regularization types,
 *  which is so far only L2 and no regularization and print the results to standard
 *  out */
template<typename Architecture>
auto testRegularizationGradients()
    -> typename Architecture::Scalar_t
{
    using Scalar_t   = typename Architecture::Scalar_t;
    using Matrix_t = typename Architecture::Matrix_t;

    std::vector<ERegularization> Regularizations
        = {ERegularization::kNone,
           ERegularization::kL2};

    Scalar_t error, maximum_error;
    maximum_error = 0.0;

    for (auto & r : Regularizations)
    {
       auto f = [r](const Matrix_t &, const Matrix_t &Y, const Matrix_t & /*W*/) {
          return regularization<Architecture>(Y, r);
       };
       auto df = [&r](Matrix_t &X, const Matrix_t &, const Matrix_t &Y, const Matrix_t & /*W*/) {
          applyMatrix(X, [](double) { return 0.0; });
          addRegularizationGradients<Architecture>(X, Y, (Scalar_t)1.0, r);
       };

       error = testGradients<Architecture>(f, df, 1.0);

       std::cout << "Testing " << static_cast<char>(r) << ": ";
       std::cout << "Maximum Relative Error = " << error << std::endl;

       maximum_error = std::max(maximum_error, error);
    }

    return maximum_error;
}
