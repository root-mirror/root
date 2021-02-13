// @(#)root/tmva/tmva/dnn:$Id$
// Authors: Surya S Dwivedi 01/08/2019, Saurav Shekhar 23/06/17
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
// network in the TCpu architecture                                //
/////////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cpu.h"
#include "TMVA/DNN/Architectures/Cpu/Blas.h"

namespace TMVA
{
namespace DNN
{

template<typename AFloat>
auto TCpu<AFloat>::RecurrentLayerBackward(TCpuMatrix<AFloat> & state_gradients_backward, // BxH
                                          TCpuMatrix<AFloat> & input_weight_gradients,
                                          TCpuMatrix<AFloat> & state_weight_gradients,
                                          TCpuMatrix<AFloat> & bias_gradients,
                                          TCpuMatrix<AFloat> & df, //BxH
                                          const TCpuMatrix<AFloat> & state, // BxH
                                          const TCpuMatrix<AFloat> & weights_input, // HxD
                                          const TCpuMatrix<AFloat> & weights_state, // HxH
                                          const TCpuMatrix<AFloat> & input,  // BxD
                                          TCpuMatrix<AFloat> & input_gradient)
-> TCpuMatrix<AFloat> &
{
   // std::cout << "Recurrent Propo" << std::endl;
   // TMVA_DNN_PrintTCpuMatrix(df,"DF");
   // TMVA_DNN_PrintTCpuMatrix(state_gradients_backward,"State grad");
   // TMVA_DNN_PrintTCpuMatrix(input_weight_gradients,"input w grad");
   // TMVA_DNN_PrintTCpuMatrix(state,"state");
   // TMVA_DNN_PrintTCpuMatrix(input,"input");

   // Compute element-wise product.
   //Hadamard(df, state_gradients_backward);  // B x H

   // Input gradients.
   if (input_gradient.GetNoElements() > 0) {
      Multiply(input_gradient, df, weights_input);
   }

   // State gradients.
   if (state_gradients_backward.GetNoElements() > 0) {
      Multiply(state_gradients_backward, df, weights_state);
   }

   // compute the gradients
   // Perform the operation in place by readding the result on the same gradient matrix
   // e.g. W += D * X

   // Weights gradients
   if (input_weight_gradients.GetNoElements() > 0) {
      TransposeMultiply(input_weight_gradients, df, input, 1. , 1.); // H x B . B x D
   }

   if (state_weight_gradients.GetNoElements() > 0) {
      TransposeMultiply(state_weight_gradients, df, state, 1. , 1. ); // H x B . B x H
   }

   // Bias gradients.
   if (bias_gradients.GetNoElements() > 0) {
      SumColumns(bias_gradients, df, 1., 1.);  // could be probably do all here
   }

   return input_gradient;
}

//______________________________________________________________________________
template <typename Scalar_t>
auto inline TCpu<Scalar_t>::LSTMLayerBackward(TCpuMatrix<Scalar_t> & state_gradients_backward,
												          TCpuMatrix<Scalar_t> & cell_gradients_backward,
												          TCpuMatrix<Scalar_t> & input_weight_gradients,
												          TCpuMatrix<Scalar_t> & forget_weight_gradients,
												          TCpuMatrix<Scalar_t> & candidate_weight_gradients,
												          TCpuMatrix<Scalar_t> & output_weight_gradients,
												          TCpuMatrix<Scalar_t> & input_state_weight_gradients,
												          TCpuMatrix<Scalar_t> & forget_state_weight_gradients,
												          TCpuMatrix<Scalar_t> & candidate_state_weight_gradients,
												          TCpuMatrix<Scalar_t> & output_state_weight_gradients,
												          TCpuMatrix<Scalar_t> & input_bias_gradients,
												          TCpuMatrix<Scalar_t> & forget_bias_gradients,
												          TCpuMatrix<Scalar_t> & candidate_bias_gradients,
												          TCpuMatrix<Scalar_t> & output_bias_gradients,
												          TCpuMatrix<Scalar_t> & di,
												          TCpuMatrix<Scalar_t> & df,
												          TCpuMatrix<Scalar_t> & dc,
												          TCpuMatrix<Scalar_t> & dout,
												          const TCpuMatrix<Scalar_t> & precStateActivations,
												          const TCpuMatrix<Scalar_t> & precCellActivations,
												          const TCpuMatrix<Scalar_t> & fInput,
												          const TCpuMatrix<Scalar_t> & fForget,
												          const TCpuMatrix<Scalar_t> & fCandidate,
												          const TCpuMatrix<Scalar_t> & fOutput,
												          const TCpuMatrix<Scalar_t> & weights_input,
												          const TCpuMatrix<Scalar_t> & weights_forget,
												          const TCpuMatrix<Scalar_t> & weights_candidate,
												          const TCpuMatrix<Scalar_t> & weights_output,
												          const TCpuMatrix<Scalar_t> & weights_input_state,
												          const TCpuMatrix<Scalar_t> & weights_forget_state,
												          const TCpuMatrix<Scalar_t> & weights_candidate_state,
												          const TCpuMatrix<Scalar_t> & weights_output_state,
												          const TCpuMatrix<Scalar_t> & input,
												          TCpuMatrix<Scalar_t> & input_gradient,
                                              TCpuMatrix<Scalar_t> & cell_gradient,
                                              TCpuMatrix<Scalar_t> & cell_tanh)
-> TCpuMatrix<Scalar_t> &
{
   //some temporary varibales used later
   TCpuMatrix<Scalar_t> tmpInp(input_gradient.GetNrows(), input_gradient.GetNcols());
   TCpuMatrix<Scalar_t> tmpState(state_gradients_backward.GetNrows(), state_gradients_backward.GetNcols());

   TCpuMatrix<Scalar_t> input_gate_gradient(fInput.GetNrows(), fInput.GetNcols());
   TCpuMatrix<Scalar_t> forget_gradient(fForget.GetNrows(), fForget.GetNcols());
   TCpuMatrix<Scalar_t> candidate_gradient(fCandidate.GetNrows(), fCandidate.GetNcols());
   TCpuMatrix<Scalar_t> output_gradient(fOutput.GetNrows(), fOutput.GetNcols());

   // cell gradient
   Hadamard(cell_gradient, fOutput);
   Hadamard(cell_gradient, state_gradients_backward);
   ScaleAdd(cell_gradient, cell_gradients_backward);
   Copy(cell_gradients_backward, cell_gradient);
   Hadamard(cell_gradients_backward, fForget);

   // candidate gradient
   Copy(candidate_gradient, cell_gradient);
   Hadamard(candidate_gradient, fInput);
   Hadamard(candidate_gradient, dc);

   // input gate gradient
   Copy(input_gate_gradient, cell_gradient);
   Hadamard(input_gate_gradient, fCandidate);
   Hadamard(input_gate_gradient, di);

   // forget gradient
   Copy(forget_gradient, cell_gradient);
   Hadamard(forget_gradient, precCellActivations);
   Hadamard(forget_gradient, df);

   // output grdient
   Copy(output_gradient, cell_tanh);
   Hadamard(output_gradient, state_gradients_backward);
   Hadamard(output_gradient, dout);

   // input gradient
   Multiply(tmpInp, input_gate_gradient, weights_input);
   Copy(input_gradient, tmpInp);
   Multiply(tmpInp, forget_gradient, weights_forget);
   ScaleAdd(input_gradient, tmpInp);
   Multiply(tmpInp, candidate_gradient, weights_candidate);
   ScaleAdd(input_gradient, tmpInp);
   Multiply(tmpInp, output_gradient, weights_output);
   ScaleAdd(input_gradient, tmpInp);

   // state gradient backwards
   Multiply(tmpState, input_gate_gradient, weights_input_state);
   Copy(state_gradients_backward, tmpState);
   Multiply(tmpState, forget_gradient, weights_forget_state);
   ScaleAdd(state_gradients_backward, tmpState);
   Multiply(tmpState, candidate_gradient, weights_candidate_state);
   ScaleAdd(state_gradients_backward, tmpState);
   Multiply(tmpState, output_gradient, weights_output_state);
   ScaleAdd(state_gradients_backward, tmpState);

   // input weight gradient
   TransposeMultiply(input_weight_gradients, input_gate_gradient, input, 1. , 1.); // H x B . B x D
   TransposeMultiply(forget_weight_gradients, forget_gradient, input, 1. , 1.);
   TransposeMultiply(candidate_weight_gradients, candidate_gradient, input, 1. , 1.);
   TransposeMultiply(output_weight_gradients, output_gradient, input, 1. , 1.);

   // state weight gradients
   TransposeMultiply(input_state_weight_gradients, input_gate_gradient, precStateActivations, 1. , 1. ); // H x B . B x H
   TransposeMultiply(forget_state_weight_gradients, forget_gradient, precStateActivations, 1. , 1. );
   TransposeMultiply(candidate_state_weight_gradients, candidate_gradient, precStateActivations, 1. , 1. );
   TransposeMultiply(output_state_weight_gradients, output_gradient, precStateActivations, 1. , 1. );

   // bias gradient
   SumColumns(input_bias_gradients, input_gate_gradient, 1., 1.);
   SumColumns(forget_bias_gradients, forget_gradient, 1., 1.);
   SumColumns(candidate_bias_gradients, candidate_gradient, 1., 1.);
   SumColumns(output_bias_gradients, output_gradient, 1., 1.);

   return input_gradient;
}


//______________________________________________________________________________
template <typename Scalar_t>
auto inline TCpu<Scalar_t>::GRULayerBackward(TCpuMatrix<Scalar_t> & state_gradients_backward,
                                             TCpuMatrix<Scalar_t> & reset_weight_gradients,
                                             TCpuMatrix<Scalar_t> & update_weight_gradients,
                                             TCpuMatrix<Scalar_t> & candidate_weight_gradients,
                                             TCpuMatrix<Scalar_t> & reset_state_weight_gradients,
                                             TCpuMatrix<Scalar_t> & update_state_weight_gradients,
                                             TCpuMatrix<Scalar_t> & candidate_state_weight_gradients,
                                             TCpuMatrix<Scalar_t> & reset_bias_gradients,
                                             TCpuMatrix<Scalar_t> & update_bias_gradients,
                                             TCpuMatrix<Scalar_t> & candidate_bias_gradients,
                                             TCpuMatrix<Scalar_t> & dr,
                                             TCpuMatrix<Scalar_t> & du,
                                             TCpuMatrix<Scalar_t> & dc,
                                             const TCpuMatrix<Scalar_t> & precStateActivations,
                                             const TCpuMatrix<Scalar_t> & fReset,
                                             const TCpuMatrix<Scalar_t> & fUpdate,
                                             const TCpuMatrix<Scalar_t> & fCandidate,
                                             const TCpuMatrix<Scalar_t> & weights_reset,
                                             const TCpuMatrix<Scalar_t> & weights_update,
                                             const TCpuMatrix<Scalar_t> & weights_candidate,
                                             const TCpuMatrix<Scalar_t> & weights_reset_state,
                                             const TCpuMatrix<Scalar_t> & weights_update_state,
                                             const TCpuMatrix<Scalar_t> & weights_candidate_state,
                                             const TCpuMatrix<Scalar_t> & input,
                                             TCpuMatrix<Scalar_t> & input_gradient,
                                             bool resetGateAfter)
-> TCpuMatrix<Scalar_t> &
{
   // reset gradient
   int r = fUpdate.GetNrows(), c = fUpdate.GetNcols();
   TCpuMatrix<Scalar_t> reset_gradient(r, c);
   Copy(reset_gradient, fUpdate);
   for (size_t j = 0; j < (size_t)reset_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)reset_gradient.GetNrows(); i++) {
         reset_gradient(i, j) = 1 - reset_gradient(i, j);
      }
   }
   Hadamard(reset_gradient, dc);
   Hadamard(reset_gradient, state_gradients_backward);
   TCpuMatrix<Scalar_t> tmpMul(r, c);

   if (!resetGateAfter) {
      // case resetGateAfter is false    U * ( r * h)
      // dr = h * (UT * dy)
      Multiply(tmpMul, reset_gradient, weights_candidate_state);
      Hadamard(tmpMul, precStateActivations);
   } else {
      // case true :   r * ( U * h) -->  dr = dy * (U * h)
      MultiplyTranspose(tmpMul, precStateActivations, weights_candidate_state);
      Hadamard(tmpMul, reset_gradient);
   }
   Hadamard(tmpMul, dr);
   Copy(reset_gradient, tmpMul);

   // update gradient
   TCpuMatrix<Scalar_t> update_gradient(r, c); // H X 1
   Copy(update_gradient, precStateActivations);
   for (size_t j = 0; j < (size_t)update_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)update_gradient.GetNrows(); i++) {
         update_gradient(i, j) = update_gradient(i, j) - fCandidate(i, j);
      }
   }
   Hadamard(update_gradient, du);
   Hadamard(update_gradient, state_gradients_backward);

   // candidate gradient
   TCpuMatrix<Scalar_t> candidate_gradient(r, c);
   Copy(candidate_gradient, fUpdate);
   for (size_t j = 0; j < (size_t)candidate_gradient.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)candidate_gradient.GetNrows(); i++) {
         candidate_gradient(i, j) = 1 - candidate_gradient(i, j);
      }
   }
   Hadamard(candidate_gradient, dc);
   Hadamard(candidate_gradient, state_gradients_backward);

   // calculating state gradient backwards term by term
   // term 1
   TCpuMatrix<Scalar_t> temp(r, c);
   Copy(temp, state_gradients_backward);
   TCpuMatrix<Scalar_t> term(r, c); // H X 1
   Copy(term, fUpdate);
   Hadamard(term, temp);
   Copy(state_gradients_backward, term);

   // term 2
   Copy(term, precStateActivations);
   Hadamard(term, du);
   Hadamard(term, temp);
   TCpuMatrix<Scalar_t> var(r, c);
   Multiply(var, term, weights_update_state);
   Copy(term, var);
   ScaleAdd(state_gradients_backward, term);

   // term 3
   Copy(term, fCandidate);
   for (size_t j = 0; j < (size_t)term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)term.GetNrows(); i++) {
         term(i, j) = -term(i, j);
      }
   }
   Hadamard(term, du);
   Hadamard(term, temp);
   Multiply(var, term, weights_update_state);
   Copy(term, var);
   ScaleAdd(state_gradients_backward, term);

   // term 4
   Copy(term, fUpdate);
   for (size_t j = 0; j < (size_t)term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)term.GetNrows(); i++) {
         term(i, j) = 1 - term(i, j);
      }
   }
   Hadamard(term, dc);
   Hadamard(term, temp);

   if (!resetGateAfter) {
      // case resetGateAfter is false   : U * ( r * h)
      // dh = r * (UT * dy)
      Multiply(var, term, weights_candidate_state);
      Hadamard(var, fReset);
   } else {
      // case resetGateAfter = true
      // dh = UT * ( r * dy )
      Hadamard(term, fReset);
      Multiply(var, term, weights_candidate_state);
   }
   //
   Copy(term, var);
   ScaleAdd(state_gradients_backward, term);

   // term 5
   Copy(term, fUpdate);
   for (size_t j = 0; j < (size_t)term.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t)term.GetNrows(); i++) {
         term(i, j) = 1 - term(i, j);
      }
   }
   // here we re-compute dr (probably we could be more eficient)
   Hadamard(term, dc);
   Hadamard(term, temp);
   if (!resetGateAfter) {
      // case reset gate after = false
      // recompute dr/dh (as above for dr): // dr = h * (UT * dy)
      Multiply(var, term, weights_candidate_state);
      Hadamard(var, precStateActivations);
   } else {
      // case = true  dr = dy * (U * h)
      MultiplyTranspose(var, precStateActivations, weights_candidate_state);
      Hadamard(var, term);
   }
   Hadamard(var, dr);
   Multiply(term, var, weights_reset_state);
   ScaleAdd(state_gradients_backward, term);

   // input gradients
   TCpuMatrix<Scalar_t> tmpInp(input_gradient.GetNrows(), input_gradient.GetNcols());
   Multiply(tmpInp, reset_gradient, weights_reset);
   Copy(input_gradient, tmpInp);
   Multiply(tmpInp, update_gradient, weights_update);
   ScaleAdd(input_gradient, tmpInp);
   Multiply(tmpInp, candidate_gradient, weights_candidate);
   ScaleAdd(input_gradient, tmpInp);

   // input weight gradients
   TransposeMultiply(reset_weight_gradients, reset_gradient, input, 1., 1.); // H x B . B x D
   TransposeMultiply(update_weight_gradients, update_gradient, input, 1., 1.);
   TransposeMultiply(candidate_weight_gradients, candidate_gradient, input, 1., 1.);

   // state weight gradients
   TransposeMultiply(reset_state_weight_gradients, reset_gradient, precStateActivations, 1., 1.); // H x B . B x H
   TransposeMultiply(update_state_weight_gradients, update_gradient, precStateActivations, 1., 1.);
   TCpuMatrix<Scalar_t> tempvar(r, c);

   // candidate weight gradients
   // impl case reseyGateAfter = false
   if (!resetGateAfter) {
      // dU = ( h * r) * dy
      Copy(tempvar, precStateActivations);
      Hadamard(tempvar, fReset);
      TransposeMultiply(candidate_state_weight_gradients, candidate_gradient, tempvar, 1., 1.);
   } else {
      // case resetAfter=true
      // dU  = h * ( r * dy)
      Copy(tempvar, candidate_gradient);
      Hadamard(tempvar, fReset);
      TransposeMultiply(candidate_state_weight_gradients, tempvar, precStateActivations, 1., 1.);
   }

   // bias gradients
   SumColumns(reset_bias_gradients, reset_gradient, 1., 1.); // could be probably do all here
   SumColumns(update_bias_gradients, update_gradient, 1., 1.);
   SumColumns(candidate_bias_gradients, candidate_gradient, 1., 1.);

   return input_gradient;
}

} // namespace DNN
} // namespace TMVA
