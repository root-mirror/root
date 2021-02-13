// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 19/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 ///////////////////////////////////////////////////////////////////
 // Implementation of the activation functions for multi-threaded //
 // CPU architectures using Roots TThreadExecutor and BLAS.            //
 ///////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cpu.h"
#include <math.h>

#ifdef R__HAS_VDT
#include "vdt/tanh.h"
#endif


namespace TMVA
{
namespace DNN
{

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::ActivationFunctionForward(Tensor_t & X, EActivationFunction activFunct,
                                             const ActivationDescriptor_t /* activationDescr */,
                                             const double /* coef */, const AFloat /*alpha */, const AFloat /*beta*/)
{
   // scaling and translation is not yet implemented
   TMVA::DNN::evaluate<TCpu<AFloat>>( X, activFunct);
}
//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::ActivationFunctionBackward(Tensor_t & dX, const Tensor_t & /* Y */,
                                                const Tensor_t & dY, const Tensor_t & X,
                                                EActivationFunction activFunct,
                                                const ActivationDescriptor_t /* activationDescr */,
                                                const AFloat /* alpha */, const AFloat /* beta */)
{
   // scaling and translation not yet implemented
   // output tensor (Y) could also be used to speed up derivative calculation
   // compute dx = f'(x)
   TMVA::DNN::evaluateDerivative<TCpu<AFloat>>(dX, activFunct, X);
    // Compute element-wise product.  dx = f'(x) * dY
   Hadamard(dX, dY);
}
//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::IdentityDerivative(TCpuTensor<AFloat> & B,
                                      const TCpuTensor<AFloat> &/*A*/)
{
   auto f = [](AFloat) {return 1.0;};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::Relu(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return (x < 0.0) ? 0.0 : x;};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::ReluDerivative(TCpuTensor<AFloat> & B,
                                               const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {return (x < 0.0) ? 0.0 : 1.0;};
   B.MapFrom(f, A);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::Sigmoid(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return 1.0 / (1.0 + exp(-x));};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SigmoidDerivative(TCpuTensor<AFloat> & B,
                                     const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {
      AFloat sig = 1.0 / (1.0 + exp(-x));
      return sig * (1.0 - sig);
   };
   B.MapFrom(f, A);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::Tanh(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return tanh(x);};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::TanhDerivative(TCpuTensor<AFloat> & B,
                                  const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {
      AFloat t = tanh(x);
      return 1 - t * t;
   };
   B.MapFrom(f, A);
}

#ifdef R__HAS_VDT
//______________________________________________________________________________
template <>
void TCpu<float>::FastTanh(TCpuTensor<float> &B)
{
   auto f = [](float x) { return vdt::fast_tanhf(x); };
   B.Map(f);
}
template <>
void TCpu<double>::FastTanh(TCpuTensor<double> &B)
{
   auto f = [](double x) { return vdt::fast_tanh(x); };
   B.Map(f);
}

//______________________________________________________________________________
template <>
void TCpu<float>::FastTanhDerivative(TCpuTensor<float> &B, const TCpuTensor<float> &A)
{
   auto f = [](float x) {
      double t = vdt::fast_tanhf(x);
      return 1 - t * t;
   };
   B.MapFrom(f, A);
}
template <>
void TCpu<double>::FastTanhDerivative(TCpuTensor<double> &B, const TCpuTensor<double> &A)
{
   auto f = [](double x) {
      double t = vdt::fast_tanh(x);
      return 1 - t * t;
   };
   B.MapFrom(f, A);
}

#else   // when VDT is not available
//______________________________________________________________________________
template <typename AFloat>
void TCpu<AFloat>::FastTanh(TCpuTensor<AFloat> &B)
{
   TCpu<AFloat>::Tanh(B);
}

//______________________________________________________________________________
template <typename AFloat>
void TCpu<AFloat>::FastTanhDerivative(TCpuTensor<AFloat> &B, const TCpuTensor<AFloat> &A)
{
   TCpu<AFloat>::TanhDerivative(B, A);
}
#endif

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SymmetricRelu(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return fabs(x);};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SymmetricReluDerivative(TCpuTensor<AFloat> & B,
                                           const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {
      return (x < 0.0) ? -1.0 : 1.0;
   };
   B.MapFrom(f, A);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SoftSign(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return x / (1 + fabs(x));};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SoftSignDerivative(TCpuTensor<AFloat> & B,
                                      const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {
      x = 1.0 + fabs(x);
      x = 1.0 / (x * x);
      return x;
   };
   B.MapFrom(f, A);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::Gauss(TCpuTensor<AFloat> & B)
{
   auto f = [](AFloat x) {return exp(- x * x);};
   B.Map(f);
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::GaussDerivative(TCpuTensor<AFloat> & B,
                                   const TCpuTensor<AFloat> &A)
{
   auto f = [](AFloat x) {return - 2.0 * x * exp(- x * x);};
   B.MapFrom(f, A);
}

} // namespace DNN
} // namespace TMVA
