// @(#)root/tmva/tmva/dnn:$Id$ 
// Author: Saurav Shekhar 23/06/17

/*************************************************************************
 * Copyright (C) 2017, Saurav Shekhar                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////
// Implementation of the functions required for the forward and    //
// backward propagation of activations through a recurrent neural  //
// network in the reference implementation.                        //
/////////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Reference.h"
#include <vector>

namespace TMVA {
namespace DNN  {

//______________________________________________________________________________
template<typename Scalar_t>
auto TReference<Scalar_t>::RecurrentLayerBackward(TMatrixT<Scalar_t> & state_gradients_backward, // BxH
                                                  TMatrixT<Scalar_t> & input_weight_gradients,
                                                  TMatrixT<Scalar_t> & state_weight_gradients,
                                                  TMatrixT<Scalar_t> & bias_gradients,
                                                  TMatrixT<Scalar_t> & df, //DxH
                                                  const TMatrixT<Scalar_t> & state, // BxH
                                                  const TMatrixT<Scalar_t> & weights_input, // HxD 
                                                  const TMatrixT<Scalar_t> & weights_state, // HxH
                                                  const TMatrixT<Scalar_t> & input,  // BxD
                                                  TMatrixT<Scalar_t> & input_gradient)
-> Matrix_t &
{
   // Compute element-wise product.
   for (size_t i = 0; i < (size_t) df.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t) df.GetNcols(); j++) {
         df(i,j) *= state_gradients_backward(i,j);      // B x H
      }
   }
   
   // Input gradients.
   if (input_gradient.GetNoElements() > 0) {
      input_gradient.Mult(df, weights_input);     // B x H . H x D = B x D
   }
   // State gradients
   if (state_gradients_backward.GetNoElements() > 0) {
      state_gradients_backward.MultT(df, weights_state);  // B x H . H x H = B x H
   }
   
   // Weights gradients.
   if (input_weight_gradients.GetNoElements() > 0) {
      TMatrixT<Scalar_t> tmp(input_weight_gradients);
      input_weight_gradients.TMult(df, input);             // H x B . B x D
      input_weight_gradients.Plus(tmp, input_weight_gradients);
   }
   if (state_weight_gradients.GetNoElements() > 0) {
      TMatrixT<Scalar_t> tmp(state_weight_gradients);
      state_weight_gradients.TMult(df, state);             // H x B . B x H
      state_weight_gradients.Plus(tmp, state_weight_gradients);
   }
   
   // Bias gradients. B x H -> H x 1
   if (bias_gradients.GetNoElements() > 0) {
      for (size_t j = 0; j < (size_t) df.GetNcols(); j++) {
         Scalar_t sum = 0.0;
         for (size_t i = 0; i < (size_t) df.GetNrows(); i++) {
            sum += df(i,j);
         }
         bias_gradients(j,0) = sum;
      }
   }
   return input_gradient;
}


} // namespace DNN
} // namespace TMVA
