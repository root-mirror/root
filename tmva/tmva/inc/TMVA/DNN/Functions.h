// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 20/06/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////
// Contains function enums for activation and output functions, as //
// well as generic evaluation functions, that delegate the call to //
// the corresponding evaluation kernel.                            //
/////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_FUNCTIONS
#define TMVA_DNN_FUNCTIONS

namespace TMVA
{
namespace DNN
{
//______________________________________________________________________________
//
//  Enum Definitions
//______________________________________________________________________________

/*! Enum that represents layer activation functions. */
enum class EActivationFunction
{
   kIdentity = 0,
   kRelu     = 1,
   kSigmoid  = 2,
   kTanh     = 3,
   kSymmRelu = 4,
   kSoftSign = 5,
   kGauss    = 6,
   kFastTanh = 7

};

/*! Enum that represents output functions */
enum class EOutputFunction
{
   kIdentity = 'I',
   kSigmoid  = 'S',
   kSoftmax  = 'M'
};

/*! Enum that represents objective functions for the net, i.e. functions
*  that take the output from the last layer in the net together with the
*  truths and return the objective function values that is to be minimized
*  in the training process. */
enum class ELossFunction
{
    kCrossEntropy        = 'C',
    kMeanSquaredError    = 'R',
    kSoftmaxCrossEntropy = 'S'
};

/*! Enum representing the regularization type applied for a given layer */
enum class ERegularization
{
    kNone = '0',
    kL1   = '1',
    kL2   = '2'
    };

/* Enum represnting the initialization method used for this layer. */
enum class EInitialization {
    kGauss    = 'G',
    kUniform  = 'U',
    kIdentity = 'I',
    kZero = 'Z',
    kGlorotNormal = 'X',
    kGlorotUniform = 'F',
};

/// Enum representing the optimizer used for training.
enum class EOptimizer {
   kSGD = 0,
   kAdam = 1,
   kAdagrad = 2,
   kRMSProp = 3,
   kAdadelta = 4,
};

//______________________________________________________________________________
//
//  Activation Functions
//______________________________________________________________________________

/*! Apply the given activation function to each value in the given
*  tensor A. */
template<typename Architecture_t>
inline void evaluate(typename Architecture_t::Tensor_t &A,
                    EActivationFunction f)
{
    switch(f)
    {
    case EActivationFunction::kIdentity : break;
    case EActivationFunction::kRelu :     Architecture_t::Relu(A);
        break;
    case EActivationFunction::kSigmoid  :  Architecture_t::Sigmoid(A);
        break;
    case EActivationFunction::kTanh     :  Architecture_t::Tanh(A);
        break;
    case EActivationFunction::kSymmRelu :  Architecture_t::SymmetricRelu(A);
        break;
    case EActivationFunction::kSoftSign :  Architecture_t::SoftSign(A);
        break;
    case EActivationFunction::kGauss    :  Architecture_t::Gauss(A);
        break;
    case EActivationFunction::kFastTanh :  Architecture_t::FastTanh(A);
        break;
    }
}

/*! Compute the first partial derivative of the activation function for
*  the values given in tensor A and write the results into B. */
//______________________________________________________________________________
template<typename Architecture_t>
inline void evaluateDerivative(typename Architecture_t::Tensor_t & B,
                                EActivationFunction f,
                                const typename Architecture_t::Tensor_t & A)
{
    switch(f)
    {
    case EActivationFunction::kIdentity : Architecture_t::IdentityDerivative(B, A);
        break;
    case EActivationFunction::kRelu     : Architecture_t::ReluDerivative(B, A);
        break;
    case EActivationFunction::kSigmoid  : Architecture_t::SigmoidDerivative(B, A);
        break;
    case EActivationFunction::kTanh     : Architecture_t::TanhDerivative(B, A);
        break;
    case EActivationFunction::kSymmRelu : Architecture_t::SymmetricReluDerivative(B, A);
        break;
    case EActivationFunction::kSoftSign : Architecture_t::SoftSignDerivative(B, A);
        break;
    case EActivationFunction::kGauss    : Architecture_t::GaussDerivative(B, A);
        break;
    case EActivationFunction::kFastTanh : Architecture_t::FastTanhDerivative(B, A);
        break;
    }
}

// matrix version of the function (for backward comp.)
template<typename Architecture_t>
inline void evaluateMatrix( typename Architecture_t::Matrix_t &A,
                        EActivationFunction f)
{
    typename Architecture_t::Tensor_t t(A);
    evaluate<Architecture_t>(t,f);
}

template<typename Architecture_t>
inline void evaluateDerivativeMatrix( typename Architecture_t::Matrix_t &B,
                        EActivationFunction f,
                        const typename Architecture_t::Matrix_t & A)
{
    typename Architecture_t::Tensor_t t(B);
    evaluateDerivative<Architecture_t>(t,f, typename Architecture_t::Tensor_t(A));
}
//______________________________________________________________________________
//
//  Output Functions
//______________________________________________________________________________

/*! Apply the given output function to each value in the given
*  tensor A. */
template<typename Architecture_t>
inline void evaluate(typename Architecture_t::Matrix_t &A,
                    EOutputFunction f,
                    const typename Architecture_t::Matrix_t &X)
{
    switch(f)
    {
    case EOutputFunction::kIdentity : Architecture_t::Copy(A, X);
                                      break;
    case EOutputFunction::kSigmoid  : Architecture_t::Sigmoid(A, X);
                                      break;
    case EOutputFunction::kSoftmax  : Architecture_t::Softmax(A, X);
                                      break;
    }
}

//______________________________________________________________________________
//
//  Loss Functions
//______________________________________________________________________________

/*! Compute the value of the objective function f for given activations
*  of the ouput layer and the truth Y. */
template <typename Architecture_t>
inline auto evaluate(ELossFunction f, const typename Architecture_t::Matrix_t &Y,
                     const typename Architecture_t::Matrix_t &output, const typename Architecture_t::Matrix_t &weights)
   -> decltype(Architecture_t::CrossEntropy(Y, output, weights))
{
    switch(f)
    {
    case ELossFunction::kCrossEntropy: return Architecture_t::CrossEntropy(Y, output, weights);
    case ELossFunction::kMeanSquaredError: return Architecture_t::MeanSquaredError(Y, output, weights);
    case ELossFunction::kSoftmaxCrossEntropy: return Architecture_t::SoftmaxCrossEntropy(Y, output, weights);
    }
    return 0.0;
}

/*! Compute the gradient of the given output function f for given activations
*  output of the output layer and truth Y and write the results into dY. */
//______________________________________________________________________________
template <typename Architecture_t>
inline void evaluateGradients(typename Architecture_t::Matrix_t &dY, ELossFunction f,
                              const typename Architecture_t::Matrix_t &Y,
                              const typename Architecture_t::Matrix_t &output,
                              const typename Architecture_t::Matrix_t &weights)
{
    switch(f)
    {
    case ELossFunction::kCrossEntropy: Architecture_t::CrossEntropyGradients(dY, Y, output, weights); break;
    case ELossFunction::kMeanSquaredError: Architecture_t::MeanSquaredErrorGradients(dY, Y, output, weights); break;
    case ELossFunction::kSoftmaxCrossEntropy :
       Architecture_t::SoftmaxCrossEntropyGradients(dY, Y, output, weights);
       break;
    }
}


//______________________________________________________________________________
//
// Regularization
//______________________________________________________________________________

/*! Evaluate the regularization functional for a given weight matrix. */
template<typename Architecture_t>
inline auto regularization(const typename Architecture_t::Matrix_t &A,
                    ERegularization R)
-> decltype(Architecture_t::L1Regularization(A))
{
    switch(R)
    {
    case ERegularization::kNone :
        return 0.0;
    case ERegularization::kL1 :
        return Architecture_t::L1Regularization(A);
    case ERegularization::kL2 :
        return Architecture_t::L2Regularization(A);
    }
    return 0.0;
}

/*! Add the regularization gradient corresponding to weight matrix W, to
*  the matrix A. */
//______________________________________________________________________________
template<typename Architecture_t>
inline void addRegularizationGradients(typename Architecture_t::Matrix_t &A,
                                       const typename Architecture_t::Matrix_t &W,
                                       typename Architecture_t::Scalar_t weightDecay,
                                       ERegularization R)
{
    switch(R)
    {
    case ERegularization::kNone :
        break;
    case ERegularization::kL1 :
        Architecture_t::AddL1RegularizationGradients(A, W, weightDecay);
        break;
    case ERegularization::kL2 :
        Architecture_t::AddL2RegularizationGradients(A, W, weightDecay);
        break;
    }
}

//______________________________________________________________________________
//
// Initialization
//______________________________________________________________________________

template<typename Architecture_t>
inline void initialize(typename Architecture_t::Matrix_t & A,
                       EInitialization m)
{
   switch(m) {
   case EInitialization::kGauss    : Architecture_t::InitializeGauss(A);
       break;
   case EInitialization::kUniform  : Architecture_t::InitializeUniform(A);
       break;
   case EInitialization::kIdentity : Architecture_t::InitializeIdentity(A);
       break;
   case EInitialization::kZero     : Architecture_t::InitializeZero(A);
       break;
   case EInitialization::kGlorotNormal    : Architecture_t::InitializeGlorotNormal(A);
       break;
   case EInitialization::kGlorotUniform  : Architecture_t::InitializeGlorotUniform(A);
       break;
   }
}

} // namespace DNN
} // namespace TMVA

#endif
