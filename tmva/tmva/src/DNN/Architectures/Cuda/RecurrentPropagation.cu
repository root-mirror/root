// @(#)root/tmva/tmva/dnn:$Id$
// Author: Saurav Shekhar 23/06/17

/*************************************************************************
 * Copyright (C) 2017, Saurav Shekhar                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 //////////////////////////////////////////////////////////////////
 // Implementation of the functions required for the forward and //
 // backward propagation of activations through a neural network //
 // for CUDA architectures.                                      //
 //////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cuda.h"
#include "TMVA/DNN/Architectures/Cuda/Device.h"
#include "Kernels.cuh"

namespace TMVA
{
namespace DNN
{

//____________________________________________________________________________
template<typename AFloat>
TCudaMatrix<AFloat> &  TCuda<AFloat>::RecurrentLayerBackward(TCudaMatrix<AFloat> & state_gradients_backward, // BxH
                                           TCudaMatrix<AFloat> & input_weight_gradients,
                                           TCudaMatrix<AFloat> & state_weight_gradients,
                                           TCudaMatrix<AFloat> & bias_gradients,
                                           TCudaMatrix<AFloat> & df, //DxH
                                           const TCudaMatrix<AFloat> & state, // BxH
                                           const TCudaMatrix<AFloat> & weights_input, // HxD
                                           const TCudaMatrix<AFloat> & weights_state, // HxH
                                           const TCudaMatrix<AFloat> & input,  // BxD
                                           TCudaMatrix<AFloat> & input_gradient)
{

   Fatal("TCuda::RecurrentLayerBackward", "Recurrent layers are not supported in the native Cuda architecture!!!");

#if 0
   /// LM: This needs to be fixed !

   // Compute element-wise product.
   TCuda<AFloat>::Hadamard(df, state_gradients_backward); // B x H

   // Input gradients.
   if (input_gradient.GetNoElements() > 0) {
      TCuda<AFloat>::Multiply(input_gradient, df, weights_input);
   }

   // State gradients.
   if (state_gradients_backward.GetNoElements() > 0) {
      TCuda<AFloat>::Multiply(state_gradients_backward, df, weights_state);
   }

   // Weights gradients
   if (input_weight_gradients.GetNoElements() > 0) {
      TCudaMatrix<AFloat> tmp(input_weight_gradients);
      TCuda<AFloat>::TransposeMultiply(input_weight_gradients, df, input); // H x B . B x D
      TCuda<AFloat>::ScaleAdd(input_weight_gradients, tmp, 1);
   }
   if (state_weight_gradients.GetNoElements() > 0) {
      TCudaMatrix<AFloat> tmp(state_weight_gradients);
      TCuda<AFloat>::TransposeMultiply(state_weight_gradients, df, state); // H x B . B x H
      TCuda<AFloat>::ScaleAdd(state_weight_gradients, tmp, 1);
   }

   // Bias gradients.
   if (bias_gradients.GetNoElements() > 0) {
      TCuda<AFloat>::SumColumns(bias_gradients, df);
   }
#endif

   return input_gradient;


}

} // namespace DNN
} // namespace TMVA
