// @(#)root/tmva/tmva/dnn:$Id$ 
// Authors: Surya S Dwivedi 15/07/2019, Saurav Shekhar 23/06/17
/*************************************************************************
 * Copyright (C) 2019, Surya S Dwivedi, Saurav Shekhar                   *
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

namespace TMVA 
{
namespace DNN  
{
  
//______________________________________________________________________________
template<typename Scalar_t>
auto TReference<Scalar_t>::RecurrentLayerBackward(TMatrixT<Scalar_t> & state_gradients_backward, // BxH
                                                  TMatrixT<Scalar_t> & input_weight_gradients,
                                                  TMatrixT<Scalar_t> & state_weight_gradients,
                                                  TMatrixT<Scalar_t> & bias_gradients,
                                                  TMatrixT<Scalar_t> & df, //BxH
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
      state_gradients_backward.Mult(df, weights_state);  // B x H . H x H = B x H
   }
   
   // Weights gradients.
   if (input_weight_gradients.GetNoElements() > 0) {
      TMatrixT<Scalar_t> tmp(input_weight_gradients);
      input_weight_gradients.TMult(df, input);             // H x B . B x D
      input_weight_gradients += tmp;
   }
   if (state_weight_gradients.GetNoElements() > 0) {
      TMatrixT<Scalar_t> tmp(state_weight_gradients);
      state_weight_gradients.TMult(df, state);             // H x B . B x H
      state_weight_gradients += tmp;
   }
   
   // Bias gradients. B x H -> H x 1
   if (bias_gradients.GetNoElements() > 0) {
      // this loops on state size
      for (size_t j = 0; j < (size_t) df.GetNcols(); j++) {
         Scalar_t sum = 0.0;
         // this loops on batch size summing all gradient contributions in a batch
         for (size_t i = 0; i < (size_t) df.GetNrows(); i++) {
            sum += df(i,j);
         }
         bias_gradients(j,0) += sum;
      }
   }
   
   return input_gradient;
}

 
//______________________________________________________________________________
template <typename Scalar_t>
auto TReference<Scalar_t>::LSTMLayerBackward(TMatrixT<Scalar_t> & state_gradients_backward,
                                             TMatrixT<Scalar_t> & cell_gradients_backward,
                                             TMatrixT<Scalar_t> & input_weight_gradients,
                                             TMatrixT<Scalar_t> & forget_weight_gradients,
                                             TMatrixT<Scalar_t> & candidate_weight_gradients,
                                             TMatrixT<Scalar_t> & output_weight_gradients,
                                             TMatrixT<Scalar_t> & input_state_weight_gradients,
                                             TMatrixT<Scalar_t> & forget_state_weight_gradients,
                                             TMatrixT<Scalar_t> & candidate_state_weight_gradients,
                                             TMatrixT<Scalar_t> & output_state_weight_gradients,
                                             TMatrixT<Scalar_t> & input_bias_gradients,
                                             TMatrixT<Scalar_t> & forget_bias_gradients,
                                             TMatrixT<Scalar_t> & candidate_bias_gradients,
                                             TMatrixT<Scalar_t> & output_bias_gradients,
                                             TMatrixT<Scalar_t> & di,
                                             TMatrixT<Scalar_t> & df,
                                             TMatrixT<Scalar_t> & dc,
                                             TMatrixT<Scalar_t> & dout,
                                             const TMatrixT<Scalar_t> & precStateActivations,
                                             const TMatrixT<Scalar_t> & precCellActivations,
                                             const TMatrixT<Scalar_t> & fInput,
                                             const TMatrixT<Scalar_t> & fForget,
                                             const TMatrixT<Scalar_t> & fCandidate,
                                             const TMatrixT<Scalar_t> & fOutput,
                                             const TMatrixT<Scalar_t> & weights_input,
                                             const TMatrixT<Scalar_t> & weights_forget,
                                             const TMatrixT<Scalar_t> & weights_candidate,
                                             const TMatrixT<Scalar_t> & weights_output,
                                             const TMatrixT<Scalar_t> & weights_input_state,
                                             const TMatrixT<Scalar_t> & weights_forget_state,
                                             const TMatrixT<Scalar_t> & weights_candidate_state,
                                             const TMatrixT<Scalar_t> & weights_output_state,
                                             const TMatrixT<Scalar_t> & input,
                                             TMatrixT<Scalar_t> & input_gradient,
                                             TMatrixT<Scalar_t> & cell_gradient,
                                             TMatrixT<Scalar_t> & cell_tanh)
-> Matrix_t & 
{
   // cell gradient
   Hadamard(cell_gradient, fOutput); 
   Hadamard(cell_gradient, state_gradients_backward);
   cell_gradient += cell_gradients_backward;
   cell_gradients_backward = cell_gradient;
   Hadamard(cell_gradients_backward, fForget);

   // candidate gradient
   TMatrixT<Scalar_t> candidate_gradient(cell_gradient);
   Hadamard(candidate_gradient, fInput);
   Hadamard(candidate_gradient, dc);
 
   // input gate gradient
   TMatrixT<Scalar_t> input_gate_gradient(cell_gradient);
   Hadamard(input_gate_gradient, fCandidate);
   Hadamard(input_gate_gradient, di);

   // forget gradient
   TMatrixT<Scalar_t> forget_gradient(cell_gradient);
   Hadamard(forget_gradient, precCellActivations);
   Hadamard(forget_gradient, df);

   // output gradient
   TMatrixT<Scalar_t> output_gradient(cell_tanh);
   Hadamard(output_gradient, state_gradients_backward);
   Hadamard(output_gradient, dout);

   // input gradient
   TMatrixT<Scalar_t> tmpInp(input_gradient);
   tmpInp.Mult(input_gate_gradient, weights_input);
   input_gradient = tmpInp;
   tmpInp.Mult(forget_gradient, weights_forget);
   input_gradient += tmpInp;
   tmpInp.Mult(candidate_gradient, weights_candidate);
   input_gradient += tmpInp;
   tmpInp.Mult(output_gradient, weights_output);
   input_gradient += tmpInp;

   // state gradient backwards
   TMatrixT<Scalar_t> tmpState(state_gradients_backward);
   tmpState.Mult(input_gate_gradient, weights_input_state);
   state_gradients_backward = tmpState;
   tmpState.Mult(forget_gradient, weights_forget_state);
   state_gradients_backward += tmpState;
   tmpState.Mult(candidate_gradient, weights_candidate_state);
   state_gradients_backward += tmpState;
   tmpState.Mult(output_gradient, weights_output_state);
   state_gradients_backward += tmpState;

   //input weight gradients
   TMatrixT<Scalar_t> tmp(input_weight_gradients);
   input_weight_gradients.TMult(input_gate_gradient, input);             
   input_weight_gradients += tmp;
   tmp = forget_weight_gradients;
   forget_weight_gradients.TMult(forget_gradient, input);
   forget_weight_gradients += tmp;
   tmp = candidate_weight_gradients;
   candidate_weight_gradients.TMult(candidate_gradient, input);
   candidate_weight_gradients += tmp;
   tmp = output_weight_gradients;
   output_weight_gradients.TMult(output_gradient, input);
   output_weight_gradients += tmp;

   // state weight gradients
   TMatrixT<Scalar_t> tmp1(input_state_weight_gradients);
   input_state_weight_gradients.TMult(input_gate_gradient, precStateActivations);        
   input_state_weight_gradients += tmp1;   
   tmp1 = forget_state_weight_gradients;
   forget_state_weight_gradients.TMult(forget_gradient, precStateActivations);
   forget_state_weight_gradients += tmp1;
   tmp1 = candidate_state_weight_gradients;
   candidate_state_weight_gradients.TMult(candidate_gradient, precStateActivations);
   candidate_state_weight_gradients += tmp1;
   tmp1 = output_state_weight_gradients;
   output_state_weight_gradients.TMult(output_gradient, precStateActivations);
   output_state_weight_gradients += tmp1;

   // bias gradients
   for (size_t j = 0; j < (size_t) df.GetNcols(); j++) {
      Scalar_t sum_inp = 0.0, sum_forget = 0.0, sum_candidate = 0.0, sum_out = 0.0;
      // this loops on batch size summing all gradient contributions in a batch
      for (size_t i = 0; i < (size_t) df.GetNrows(); i++) {
         sum_inp += input_gate_gradient(i,j);
         sum_forget += forget_gradient(i,j);
         sum_candidate += candidate_gradient(i,j);
         sum_out += output_gradient(i,j);
      }
      input_bias_gradients(j,0) += sum_inp;
      forget_bias_gradients(j,0) += sum_forget;
      candidate_bias_gradients(j,0) += sum_candidate;
      output_bias_gradients(j,0) += sum_out;
   }
   
   return input_gradient;
}



//______________________________________________________________________________
template <typename Scalar_t>
auto TReference<Scalar_t>::GRULayerBackward(TMatrixT<Scalar_t> & state_gradients_backward,
                                            TMatrixT<Scalar_t> & reset_weight_gradients,
                                            TMatrixT<Scalar_t> & update_weight_gradients,
                                            TMatrixT<Scalar_t> & candidate_weight_gradients,
                                            TMatrixT<Scalar_t> & reset_state_weight_gradients,
                                            TMatrixT<Scalar_t> & update_state_weight_gradients,
                                            TMatrixT<Scalar_t> & candidate_state_weight_gradients,
                                            TMatrixT<Scalar_t> & reset_bias_gradients,
                                            TMatrixT<Scalar_t> & update_bias_gradients,
                                            TMatrixT<Scalar_t> & candidate_bias_gradients,
                                            TMatrixT<Scalar_t> & dr,
                                            TMatrixT<Scalar_t> & du,
                                            TMatrixT<Scalar_t> & dc,
                                            const TMatrixT<Scalar_t> & precStateActivations,
                                            const TMatrixT<Scalar_t> & fReset,
                                            const TMatrixT<Scalar_t> & fUpdate,
                                            const TMatrixT<Scalar_t> & fCandidate,
                                            const TMatrixT<Scalar_t> & weights_reset,
                                            const TMatrixT<Scalar_t> & weights_update,
                                            const TMatrixT<Scalar_t> & weights_candidate,
                                            const TMatrixT<Scalar_t> & weights_reset_state,
                                            const TMatrixT<Scalar_t> & weights_update_state,
                                            const TMatrixT<Scalar_t> & weights_candidate_state,
                                            const TMatrixT<Scalar_t> & input,
                                            TMatrixT<Scalar_t> & input_gradient)
-> Matrix_t &
{
   // reset gradient
   TMatrixT<Scalar_t> reset_gradient(fUpdate);
   for (size_t j = 0; j < (size_t) reset_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) reset_gradient.GetNrows(); i++) {
         reset_gradient(i,j) = 1 - reset_gradient(i,j);
      }
   }
   Hadamard(reset_gradient, dc);
   Hadamard(reset_gradient, state_gradients_backward);
   TMatrixT<Scalar_t> tmpMul(precStateActivations);
   tmpMul.Mult(reset_gradient, weights_candidate_state);
   Hadamard(tmpMul, precStateActivations);
   Hadamard(tmpMul, dr);
   reset_gradient = tmpMul;
   
   // update gradient
   TMatrixT<Scalar_t> update_gradient(precStateActivations); // H X 1
   for (size_t j = 0; j < (size_t) update_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) update_gradient.GetNrows(); i++) {
         update_gradient(i,j) = update_gradient(i,j) - fCandidate(i,j);
      }
   }
   Hadamard(update_gradient, du);
   Hadamard(update_gradient, state_gradients_backward);

   // candidate gradient
   TMatrixT<Scalar_t> candidate_gradient(fUpdate);
   for (size_t j = 0; j < (size_t) candidate_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) candidate_gradient.GetNrows(); i++) {
         candidate_gradient(i,j) = 1 - candidate_gradient(i,j);
      }
   }
   Hadamard(candidate_gradient, dc);
   Hadamard(candidate_gradient, state_gradients_backward);

   // calculating state_gradient_backwards term by term
   // term 1
   TMatrixT<Scalar_t> temp(state_gradients_backward);
   TMatrixT<Scalar_t> term(fUpdate); // H X 1
   Hadamard(term, temp);
   state_gradients_backward = term;

   //term 2
   term = precStateActivations;
   Hadamard(term, du);
   Hadamard(term, temp);
   TMatrixT<Scalar_t> var(precStateActivations);
   var.Mult(term, weights_update_state);
   term = var;
   state_gradients_backward += term;

   // term 3
   term = fCandidate;
   for (size_t j = 0; j < (size_t) term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) term.GetNrows(); i++) {
         term(i,j) = - term(i,j);
      }
   }
   Hadamard(term, du);
   Hadamard(term, temp);
   var.Mult(term, weights_update_state);
   term = var;
   state_gradients_backward += term;   

   // term 4
   term = fUpdate;
   for (size_t j = 0; j < (size_t) term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) term.GetNrows(); i++) {
         term(i,j) = 1 - term(i,j);
      }
   }
   Hadamard(term, dc);
   Hadamard(term, temp);
   var.Mult(term, weights_candidate_state);
   Hadamard(var, fReset);
   term = var;
   state_gradients_backward += term;

   // term 5
   term = fUpdate;
   for (size_t j = 0; j < (size_t) term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) term.GetNrows(); i++) {
         term(i,j) = 1 - term(i,j);
      }
   }
   Hadamard(term, dc);
   Hadamard(term, temp);
   var.Mult(term, weights_candidate_state);
   Hadamard(var, precStateActivations);
   Hadamard(var, dr);
   term.Mult(var, weights_reset_state);
   state_gradients_backward += term;

   // input gradients
   TMatrixT<Scalar_t> tmpInp(input_gradient);
   tmpInp.Mult(reset_gradient, weights_reset);
   input_gradient = tmpInp;
   tmpInp.Mult(update_gradient, weights_update);
   input_gradient += tmpInp;
   tmpInp.Mult(candidate_gradient, weights_candidate);
   input_gradient += tmpInp;
   
   //input weight gradients
   TMatrixT<Scalar_t> tmp(reset_weight_gradients);
   reset_weight_gradients.TMult(reset_gradient, input);             
   reset_weight_gradients += tmp;
   tmp = update_weight_gradients;
   update_weight_gradients.TMult(update_gradient, input);
   update_weight_gradients += tmp;
   tmp = candidate_weight_gradients;
   candidate_weight_gradients.TMult(candidate_gradient, input);
   candidate_weight_gradients += tmp;

   // state weight gradients
   TMatrixT<Scalar_t> tmp1(reset_state_weight_gradients);
   reset_state_weight_gradients.TMult(reset_gradient, precStateActivations);        
   reset_state_weight_gradients += tmp1;  
   tmp1 = update_state_weight_gradients;
   update_state_weight_gradients.TMult(update_gradient, precStateActivations);
   update_state_weight_gradients += tmp1;
   tmp1 = candidate_state_weight_gradients;
   TMatrixT<Scalar_t> tmp2(fReset);
   Hadamard(tmp2, precStateActivations);
   candidate_state_weight_gradients.TMult(candidate_gradient, tmp2);
   candidate_state_weight_gradients += tmp1;

   // bias gradients
   for (size_t j = 0; j < (size_t) du.GetNcols(); j++) {
      Scalar_t sum_reset = 0.0, sum_update = 0.0, sum_candidate = 0.0;
      // this loops on batch size summing all gradient contributions in a batch
      for (size_t i = 0; i < (size_t) du.GetNrows(); i++) {
         sum_reset += reset_gradient(i,j);
         sum_update += update_gradient(i,j);
         sum_candidate += candidate_gradient(i,j);
      }
      reset_bias_gradients(j,0) += sum_reset;
      update_bias_gradients(j,0) += sum_update;
      candidate_bias_gradients(j,0) += sum_candidate;
   }

   return input_gradient;
}

} // namespace DNN
} // namespace TMVA
