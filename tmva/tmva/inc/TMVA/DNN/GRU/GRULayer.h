// @(#)root/tmva/tmva/dnn/gru:$Id$
// Author: Surya S Dwivedi 03/07/19

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class : BasicGRULayer                                                         *
 *                                                                                *
 * Description:                                                                   *
 *       NeuralNetwork                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *       Surya S Dwivedi  <surya2191997@gmail.com> - IIT Kharagpur, India         *
 *                                                                                *
 * Copyright (c) 2005-2019:                                                       *
 * All rights reserved.                                                           *
 *       CERN, Switzerland                                                        *
 *                                                                                *
 * For the licensing terms see $ROOTSYS/LICENSE.                                  *
 * For the list of contributors see $ROOTSYS/README/CREDITS.                      *
 **********************************************************************************/

//#pragma once

//////////////////////////////////////////////////////////////////////
// This class implements the GRU layer. GRU is a variant of vanilla 
// RNN which is capable of learning long range dependencies. 
//////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_GRU_LAYER
#define TMVA_DNN_GRU_LAYER

#include <cmath>
#include <iostream>
#include <vector>

#include "TMatrix.h"
#include "TMVA/DNN/Functions.h"

namespace TMVA
{
namespace DNN
{
namespace GRU
{

//______________________________________________________________________________
//
// Basic GRU Layer
//______________________________________________________________________________

/** \class BasicGRULayer
      Generic implementation
*/
template<typename Architecture_t>
      class TBasicGRULayer : public VGeneralLayer<Architecture_t>
{

public:

   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using Tensor_t = std::vector<Matrix_t>;

  

private:

   size_t fStateSize;                           ///< Hidden state size for GRU
   size_t fTimeSteps;                           ///< Timesteps for GRU

   bool fRememberState;                         ///< Remember state in next pass

   DNN::EActivationFunction fF1;                ///< Activation function: sigmoid
   DNN::EActivationFunction fF2;                ///< Activaton function: tanh

   Matrix_t fResetValue;                        ///< Computed reset gate values
   Matrix_t fUpdateValue;                       ///< Computed forget gate values
   Matrix_t fCandidateValue;                    ///< Computed candidate values
   Matrix_t fState;                             ///< Hidden state of GRU


   Matrix_t &fWeightsResetGate;                 ///< Reset Gate weights for input, fWeights[0]
   Matrix_t &fWeightsResetGateState;            ///< Input Gate weights for prev state, fWeights[1]
   Matrix_t &fResetGateBias;                    ///< Input Gate bias

   Matrix_t &fWeightsUpdateGate;                ///< Update Gate weights for input, fWeights[2]
   Matrix_t &fWeightsUpdateGateState;           ///< Update Gate weights for prev state, fWeights[3]
   Matrix_t &fUpdateGateBias;                   ///< Update Gate bias

   Matrix_t &fWeightsCandidate;                 ///< Candidate Gate weights for input, fWeights[4]
   Matrix_t &fWeightsCandidateState;            ///< Candidate Gate weights for prev state, fWeights[5]
   Matrix_t &fCandidateBias;                    ///< Candidate Gate bias


   std::vector<Matrix_t> reset_gate_value;      ///< Reset gate value for every time step
   std::vector<Matrix_t> update_gate_value;     ///< Update gate value for every time step
   std::vector<Matrix_t> candidate_gate_value;  ///< Candidate gate value for every time step
  
   std::vector<Matrix_t> fDerivativesReset;     ///< First fDerivatives of the activations reset gate
   std::vector<Matrix_t> fDerivativesUpdate;    ///< First fDerivatives of the activations update gate
   std::vector<Matrix_t> fDerivativesCandidate; ///< First fDerivatives of the activations candidate gate
   
   Matrix_t &fWeightsResetGradients;            ///< Gradients w.r.t the reset gate - input weights
   Matrix_t &fWeightsResetStateGradients;       ///< Gradients w.r.t the reset gate - hidden state weights
   Matrix_t &fResetBiasGradients;               ///< Gradients w.r.t the reset gate - bias weights
   Matrix_t &fWeightsUpdateGradients;           ///< Gradients w.r.t the update gate - input weights
   Matrix_t &fWeightsUpdateStateGradients;      ///< Gradients w.r.t the update gate - hidden state weights
   Matrix_t &fUpdateBiasGradients;              ///< Gradients w.r.t the update gate - bias weights
   Matrix_t &fWeightsCandidateGradients;        ///< Gradients w.r.t the candidate gate - input weights
   Matrix_t &fWeightsCandidateStateGradients;   ///< Gradients w.r.t the candidate gate - hidden state weights
   Matrix_t &fCandidateBiasGradients;           ///< Gradients w.r.t the candidate gate - bias weights

public:

   /*! Constructor */
   TBasicGRULayer(size_t batchSize, size_t stateSize, size_t inputSize,
                   size_t timeSteps, bool rememberState = false,
                   DNN::EActivationFunction f1 = DNN::EActivationFunction::kSigmoid,
                   DNN::EActivationFunction f2 = DNN::EActivationFunction::kTanh,
                   bool training = true, DNN::EInitialization fA = DNN::EInitialization::kZero);

   /*! Copy Constructor */
   TBasicGRULayer(const TBasicGRULayer &);

   /*! Initialize the hidden state and cell state method. */
   void InitState(DNN::EInitialization m = DNN::EInitialization::kZero);
    
   /*! Computes the next hidden state 
    *  and next cell state with given input matrix. */
   void Forward(Tensor_t &input, bool isTraining = true);

   /*! Forward for a single cell (time unit) */
   void CellForward(Matrix_t &updateGateValues, Matrix_t &candidateValues);

   /*! Backpropagates the error. Must only be called directly at the corresponding
    *  call to Forward(...). */
   void Backward(Tensor_t &gradients_backward,
                 const Tensor_t &activations_backward,
                 std::vector<Matrix_t> &inp1,
                 std::vector<Matrix_t> &inp2);
    
   /* Updates weights and biases, given the learning rate */
   void Update(const Scalar_t learningRate);
    
   /*! Backward for a single time unit
    *  a the corresponding call to Forward(...). */
   Matrix_t & CellBackward(Matrix_t & state_gradients_backward,
                           const Matrix_t & precStateActivations,
                           const Matrix_t & reset_gate, const Matrix_t & update_gate,   
                           const Matrix_t & candidate_gate, 
                           const Matrix_t & input, Matrix_t & input_gradient,
                           Matrix_t &dr, Matrix_t &du, Matrix_t &dc, size_t t);
    
   /*! Decides the values we'll update (NN with Sigmoid) */
   void ResetGate(const Matrix_t &input, Matrix_t &di);
    
   /*! Forgets the past values (NN with Sigmoid) */
   void UpdateGate(const Matrix_t &input, Matrix_t &df);
    
   /*! Decides the new candidate values (NN with Tanh) */
   void CandidateValue(const Matrix_t &input, Matrix_t &dc);
       
   /*! Prints the info about the layer */
   void Print() const;
    
   /*! Writes the information and the weights about the layer in an XML node. */
   void AddWeightsXMLTo(void *parent) override;
    
   /*! Read the information and the weights about the layer from XML node. */
   void ReadWeightsFromXML(void *parent) override;
    
   /*! Getters */
   size_t GetInputSize()               const { return this->GetInputWidth(); }
   size_t GetTimeSteps()               const { return fTimeSteps; }
   size_t GetStateSize()               const { return fStateSize; }

   inline bool DoesRememberState()       const { return fRememberState; }
   
   inline DNN::EActivationFunction     GetActivationFunctionF1()        const { return fF1; }
   inline DNN::EActivationFunction     GetActivationFunctionF2()        const { return fF2; }

   const Matrix_t                    & GetResetGateValue()                const { return fResetValue; }
   Matrix_t                          & GetResetGateValue()                      { return fResetValue; }
   const Matrix_t                    & GetCandidateValue()                const { return fCandidateValue; }
   Matrix_t                          & GetCandidateValue()                      { return fCandidateValue; }
   const Matrix_t                    & GetUpdateGateValue()               const { return fUpdateValue; }
   Matrix_t                          & GetUpdateGateValue()                     { return fUpdateValue; }

   const Matrix_t                    & GetState()                   const { return fState; }
   Matrix_t                          & GetState()                         { return fState; }

   const Matrix_t                    & GetWeightsResetGate()              const { return fWeightsResetGate; }
   Matrix_t                          & GetWeightsResetGate()                    { return fWeightsResetGate; }
   const Matrix_t                    & GetWeightsCandidate()              const { return fWeightsCandidate; }
   Matrix_t                          & GetWeightsCandidate()                    { return fWeightsCandidate; }
   const Matrix_t                    & GetWeightsUpdateGate()             const { return fWeightsUpdateGate; }
   Matrix_t                          & GetWeightsUpdateGate()                   { return fWeightsUpdateGate; }
   
   const Matrix_t                    & GetWeightsResetGateState()         const { return fWeightsResetGateState; }
   Matrix_t                          & GetWeightsResetGateState()               { return fWeightsResetGateState; }
   const Matrix_t                    & GetWeightsUpdateGateState()        const { return fWeightsUpdateGateState; }
   Matrix_t                          & GetWeightsUpdateGateState()              { return fWeightsUpdateGateState; }
   const Matrix_t                    & GetWeightsCandidateState()         const { return fWeightsCandidateState; }
   Matrix_t                          & GetWeightsCandidateState()               { return fWeightsCandidateState; }

   const std::vector<Matrix_t>       & GetDerivativesReset()              const { return fDerivativesReset; }
   std::vector<Matrix_t>             & GetDerivativesReset()                    { return fDerivativesReset; }
   const Matrix_t                    & GetResetDerivativesAt(size_t i)    const { return fDerivativesReset[i]; }
   Matrix_t                          & GetResetDerivativesAt(size_t i)           { return fDerivativesReset[i]; }
   const std::vector<Matrix_t>       & GetDerivativesUpdate()              const { return fDerivativesUpdate; }
   std::vector<Matrix_t>             & GetDerivativesUpdate()                    { return fDerivativesUpdate; }
   const Matrix_t                    & GetUpdateDerivativesAt(size_t i)    const { return fDerivativesUpdate[i]; }
   Matrix_t                          & GetUpdateDerivativesAt(size_t i)          { return fDerivativesUpdate[i]; }
   const std::vector<Matrix_t>       & GetDerivativesCandidate()           const { return fDerivativesCandidate; }
   std::vector<Matrix_t>             & GetDerivativesCandidate()                 { return fDerivativesCandidate; }
   const Matrix_t                    & GetCandidateDerivativesAt(size_t i) const { return fDerivativesCandidate[i]; }
   Matrix_t                          & GetCandidateDerivativesAt(size_t i)       { return fDerivativesCandidate[i]; }
   
   const std::vector<Matrix_t>       & GetResetGateTensor()              const { return reset_gate_value; }
   std::vector<Matrix_t>             & GetResetGateTensor()                    { return reset_gate_value; }
   const Matrix_t                    & GetResetGateTensorAt(size_t i)    const { return reset_gate_value[i]; }
   Matrix_t                          & GetResetGateTensorAt(size_t i)           { return reset_gate_value[i]; }
   const std::vector<Matrix_t>       & GetUpdateGateTensor()              const { return update_gate_value; }
   std::vector<Matrix_t>             & GetUpdateGateTensor()                    { return update_gate_value; }
   const Matrix_t                    & GetUpdateGateTensorAt(size_t i)    const { return update_gate_value[i]; }
   Matrix_t                          & GetUpdateGateTensorAt(size_t i)          { return update_gate_value[i]; }
   const std::vector<Matrix_t>       & GetCandidateGateTensor()           const { return candidate_gate_value; }
   std::vector<Matrix_t>             & GetCandidateGateTensor()                 { return candidate_gate_value; }
   const Matrix_t                    & GetCandidateGateTensorAt(size_t i) const { return candidate_gate_value[i]; }
   Matrix_t                          & GetCandidateGateTensorAt(size_t i)       { return candidate_gate_value[i]; }
  
   

   const Matrix_t                   & GetResetGateBias()         const { return fResetGateBias; }
   Matrix_t                         & GetResetGateBias()               { return fResetGateBias; }
   const Matrix_t                   & GetUpdateGateBias()        const { return fUpdateGateBias; }
   Matrix_t                         & GetUpdateGateBias()              { return fUpdateGateBias; }
   const Matrix_t                   & GetCandidateBias()         const { return fCandidateBias; }
   Matrix_t                         & GetCandidateBias()               { return fCandidateBias; }
   
   const Matrix_t                   & GetWeightsResetGradients()        const { return fWeightsResetGradients; }
   Matrix_t                         & GetWeightsResetGradients()              { return fWeightsResetGradients; }
   const Matrix_t                   & GetWeightsResetStateGradients()   const { return fWeightsResetStateGradients; }
   Matrix_t                         & GetWeightsResetStateGradients()         { return fWeightsResetStateGradients; }
   const Matrix_t                   & GetResetBiasGradients()           const { return fResetBiasGradients; }
   Matrix_t                         & GetResetBiasGradients()                 { return fResetBiasGradients; }
   const Matrix_t                   & GetWeightsUpdateGradients()      const { return fWeightsUpdateGradients; }
   Matrix_t                         & GetWeightsUpdateGradients()            { return fWeightsUpdateGradients; }
   const Matrix_t                   & GetWeigthsUpdateStateGradients() const { return fWeightsUpdateStateGradients; }
   Matrix_t                         & GetWeightsUpdateStateGradients()       { return fWeightsUpdateStateGradients; }
   const Matrix_t                   & GetUpdateBiasGradients()         const { return fUpdateBiasGradients; }
   Matrix_t                         & GetUpdateBiasGradients()               { return fUpdateBiasGradients; }
   const Matrix_t                   & GetWeightsCandidateGradients()      const { return fWeightsCandidateGradients; }
   Matrix_t                         & GetWeightsCandidateGradients()            { return fWeightsCandidateGradients; }
   const Matrix_t                   & GetWeightsCandidateStateGradients() const { return fWeightsCandidateStateGradients; }
   Matrix_t                         & GetWeightsCandidateStateGradients()       { return fWeightsCandidateStateGradients; }
   const Matrix_t                   & GetCandidateBiasGradients()         const { return fCandidateBiasGradients; }
   Matrix_t                         & GetCandidateBiasGradients()               { return fCandidateBiasGradients; }

};

//______________________________________________________________________________
//
// Basic GRU-Layer Implementation
//______________________________________________________________________________

template <typename Architecture_t>
TBasicGRULayer<Architecture_t>::TBasicGRULayer(size_t batchSize, size_t stateSize, size_t inputSize, size_t timeSteps,
                                                 bool rememberState, DNN::EActivationFunction f1, DNN::EActivationFunction f2,
                                                 bool /* training */, DNN::EInitialization fA)
   : VGeneralLayer<Architecture_t>(batchSize, 1, timeSteps, inputSize, 1, timeSteps, stateSize, 6,
                                   {stateSize, stateSize, stateSize, stateSize, stateSize, stateSize}, 
                                   {inputSize, stateSize, inputSize, stateSize, inputSize, stateSize}, 
                                   3, {stateSize, stateSize, stateSize}, {1, 1, 1}, batchSize, timeSteps, stateSize, fA),
   fStateSize(stateSize),
   fTimeSteps(timeSteps),
   fRememberState(rememberState),
   fF1(f1),
   fF2(f2),
   fResetValue(batchSize, stateSize),
   fCandidateValue(batchSize, stateSize),
   fUpdateValue(batchSize, stateSize),
   fState(batchSize, stateSize),
   fWeightsResetGate(this->GetWeightsAt(0)),
   fWeightsResetGateState(this->GetWeightsAt(1)),
   fResetGateBias(this->GetBiasesAt(0)),
   fWeightsUpdateGate(this->GetWeightsAt(2)),
   fWeightsUpdateGateState(this->GetWeightsAt(3)),
   fUpdateGateBias(this->GetBiasesAt(1)),
   fWeightsCandidate(this->GetWeightsAt(4)),
   fWeightsCandidateState(this->GetWeightsAt(5)),
   fCandidateBias(this->GetBiasesAt(2)),
   fWeightsResetGradients(this->GetWeightGradientsAt(0)),
   fWeightsResetStateGradients(this->GetWeightGradientsAt(1)),
   fResetBiasGradients(this->GetBiasGradientsAt(0)),
   fWeightsUpdateGradients(this->GetWeightGradientsAt(2)),
   fWeightsUpdateStateGradients(this->GetWeightGradientsAt(3)),
   fUpdateBiasGradients(this->GetBiasGradientsAt(1)),
   fWeightsCandidateGradients(this->GetWeightGradientsAt(4)),
   fWeightsCandidateStateGradients(this->GetWeightGradientsAt(5)),
   fCandidateBiasGradients(this->GetBiasGradientsAt(2))
{
   for (size_t i = 0; i < timeSteps; ++i) {
      fDerivativesReset.emplace_back(batchSize, stateSize);
      fDerivativesUpdate.emplace_back(batchSize, stateSize);
      fDerivativesCandidate.emplace_back(batchSize, stateSize);
      reset_gate_value.emplace_back(batchSize, stateSize);
      update_gate_value.emplace_back(batchSize, stateSize);
      candidate_gate_value.emplace_back(batchSize, stateSize);
   }
}

 //______________________________________________________________________________
template <typename Architecture_t>
TBasicGRULayer<Architecture_t>::TBasicGRULayer(const TBasicGRULayer &layer)
   : VGeneralLayer<Architecture_t>(layer), 
      fStateSize(layer.fStateSize),
      fTimeSteps(layer.fTimeSteps),
      fRememberState(layer.fRememberState),
      fF1(layer.GetActivationFunctionF1()),
      fF2(layer.GetActivationFunctionF2()),
      fResetValue(layer.GetBatchSize(), layer.GetStateSize()),
      fCandidateValue(layer.GetBatchSize(), layer.GetStateSize()),
      fUpdateValue(layer.GetBatchSize(), layer.GetStateSize()),
      fState(layer.GetBatchSize(), layer.GetStateSize()),
      fWeightsResetGate(this->GetWeightsAt(0)),
      fWeightsResetGateState(this->GetWeightsAt(1)),
      fResetGateBias(this->GetBiasesAt(0)),
      fWeightsUpdateGate(this->GetWeightsAt(2)),
      fWeightsUpdateGateState(this->GetWeightsAt(3)),
      fUpdateGateBias(this->GetBiasesAt(1)),
      fWeightsCandidate(this->GetWeightsAt(4)),
      fWeightsCandidateState(this->GetWeightsAt(5)),
      fCandidateBias(this->GetBiasesAt(2)),
      fWeightsResetGradients(this->GetWeightGradientsAt(0)),
      fWeightsResetStateGradients(this->GetWeightGradientsAt(1)),
      fResetBiasGradients(this->GetBiasGradientsAt(0)),
      fWeightsUpdateGradients(this->GetWeightGradientsAt(2)),
      fWeightsUpdateStateGradients(this->GetWeightGradientsAt(3)),
      fUpdateBiasGradients(this->GetBiasGradientsAt(1)),
      fWeightsCandidateGradients(this->GetWeightGradientsAt(4)),
      fWeightsCandidateStateGradients(this->GetWeightGradientsAt(5)),
      fCandidateBiasGradients(this->GetBiasGradientsAt(2))
{
   for (size_t i = 0; i < fTimeSteps; ++i) {
      fDerivativesReset.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesReset[i], layer.GetResetDerivativesAt(i));
        
      fDerivativesUpdate.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesUpdate[i], layer.GetUpdateDerivativesAt(i));
        
      fDerivativesCandidate.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesCandidate[i], layer.GetCandidateDerivativesAt(i));
        
      reset_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(reset_gate_value[i], layer.GetResetGateTensorAt(i));

      update_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(update_gate_value[i], layer.GetUpdateGateTensorAt(i));

      candidate_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(candidate_gate_value[i], layer.GetCandidateGateTensorAt(i));
   }
    
   // Gradient matrices not copied
   Architecture_t::Copy(fState, layer.GetState());

   // Copy each gate values.
   Architecture_t::Copy(fResetValue, layer.GetResetGateValue());
   Architecture_t::Copy(fCandidateValue, layer.GetCandidateValue());
   Architecture_t::Copy(fUpdateValue, layer.GetUpdateGateValue());
}
 
//______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::ResetGate(const Matrix_t &input, Matrix_t &dr)
-> void
{
   /*! Computes reset gate values according to equation:
    *  input = act(W_input . input + W_state . state + bias)
    *  activation function: sigmoid. */
   const DNN::EActivationFunction fRst = this->GetActivationFunctionF1();
   Matrix_t tmpState(fResetValue.GetNrows(), fResetValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsResetGateState);
   Architecture_t::MultiplyTranspose(fResetValue, input, fWeightsResetGate);
   Architecture_t::ScaleAdd(fResetValue, tmpState);
   Architecture_t::AddRowWise(fResetValue, fResetGateBias);
   DNN::evaluateDerivative<Architecture_t>(dr, fRst, fResetValue);
   DNN::evaluate<Architecture_t>(fResetValue, fRst);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::UpdateGate(const Matrix_t &input, Matrix_t &du)
-> void
{
   /*! Computes update gate values according to equation:
    *  forget = act(W_input . input + W_state . state + bias)
    *  activation function: sigmoid. */
   const DNN::EActivationFunction fUpd = this->GetActivationFunctionF1();
   Matrix_t tmpState(fUpdateValue.GetNrows(), fUpdateValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsUpdateGateState);
   Architecture_t::MultiplyTranspose(fUpdateValue, input, fWeightsUpdateGate);
   Architecture_t::ScaleAdd(fUpdateValue, tmpState);
   Architecture_t::AddRowWise(fUpdateValue, fUpdateGateBias);
   DNN::evaluateDerivative<Architecture_t>(du, fUpd, fUpdateValue);
   DNN::evaluate<Architecture_t>(fUpdateValue, fUpd);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::CandidateValue(const Matrix_t &input, Matrix_t &dc)
-> void
{
   /*! candidate_value = act(W_input . input + W_state . (reset*state) + bias)
    *  activation function = tanh. */
   const DNN::EActivationFunction fCan = this->GetActivationFunctionF2();
   Matrix_t tmpState(fResetValue);
   Architecture_t::Hadamard(tmpState, fState);
   Matrix_t tmp(fCandidateValue.GetNrows(), fCandidateValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmp, tmpState, fWeightsCandidateState);
   Architecture_t::MultiplyTranspose(fCandidateValue, input, fWeightsCandidate);
   Architecture_t::ScaleAdd(fCandidateValue, tmp);
   Architecture_t::AddRowWise(fCandidateValue, fCandidateBias);
   DNN::evaluateDerivative<Architecture_t>(dc, fCan, fCandidateValue);
   DNN::evaluate<Architecture_t>(fCandidateValue, fCan);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::Forward(Tensor_t &input, bool /* isTraining = true */)
-> void
{
   // D : input size
   // H : state size
   // T : time size
   // B : batch size

   Tensor_t arrInput;
   for (size_t t = 0; t < fTimeSteps; ++t) {
      arrInput.emplace_back(this->GetBatchSize(), this->GetInputWidth()); // T x B x D
   }
   Architecture_t::Rearrange(arrInput, input); // B x T x D

   Tensor_t arrOutput;
   for (size_t t = 0; t < fTimeSteps;++t) {
      arrOutput.emplace_back(this->GetBatchSize(), fStateSize); // T x B x H 
   }
  
   if (!this->fRememberState) {
      InitState(DNN::EInitialization::kZero);
   }

   /*! Pass each gate values to CellForward() to calculate
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < fTimeSteps; ++t) {
      /* Feed forward network: value of each gate being computed at each timestep t. */
      ResetGate(arrInput[t], fDerivativesReset[t]);
      UpdateGate(arrInput[t], fDerivativesUpdate[t]);
      CandidateValue(arrInput[t], fDerivativesCandidate[t]);
   
      Architecture_t::Copy(this->GetResetGateTensorAt(t), fResetValue);
      Architecture_t::Copy(this->GetUpdateGateTensorAt(t), fUpdateValue);
      Architecture_t::Copy(this->GetCandidateGateTensorAt(t), fCandidateValue);
       
      CellForward(fUpdateValue, fCandidateValue);
      Architecture_t::Copy(arrOutput[t], fState);
   }

   Architecture_t::Rearrange(this->GetOutput(), arrOutput);  // B x T x D
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::CellForward(Matrix_t &updateGateValues, Matrix_t &candidateValues)
-> void
{ 
   Matrix_t tmp(updateGateValues); // H X 1
   for (size_t j = 0; j < (size_t) tmp.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) tmp.GetNrows(); i++) {
         tmp(i,j) = 1 - tmp(i,j);
      }
   }

   // Update state
   Architecture_t::Hadamard(fState, tmp);
   Architecture_t::Hadamard(updateGateValues, candidateValues);
   Architecture_t::ScaleAdd(fState, updateGateValues);
}

 //____________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::Backward(Tensor_t &gradients_backward,           // B x T x D
                                                      const Tensor_t &activations_backward,   // B x T x D
                                                      std::vector<Matrix_t> & /*inp1*/,
                                                      std::vector<Matrix_t> & /*inp2*/)
-> void
{
   // gradients_backward is activationGradients of layer before it, which is input layer.
   // Currently, gradients_backward is for input(x) and not for state.
   // For the state it can be:
   Matrix_t state_gradients_backward(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(state_gradients_backward, DNN::EInitialization::kZero); // B x H

   // if dummy is false gradients_backward will be written back on the matrix
   bool dummy = false;
   if (gradients_backward.size() == 0 || gradients_backward[0].GetNrows() == 0 || gradients_backward[0].GetNcols() == 0) {
      dummy = true;
   }

   Tensor_t arr_gradients_backward;
   for (size_t t = 0; t < fTimeSteps; ++t) {
      arr_gradients_backward.emplace_back(this->GetBatchSize(), this->GetInputSize()); // T x B x D
   }
   
   //Architecture_t::Rearrange(arr_gradients_backward, gradients_backward); // B x T x D
   // activations_backward is input.
   Tensor_t arr_activations_backward;
   for (size_t t = 0; t < fTimeSteps; ++t) {
      arr_activations_backward.emplace_back(this->GetBatchSize(), this->GetInputSize()); // T x B x D
   }
   Architecture_t::Rearrange(arr_activations_backward, activations_backward); // B x T x D

   /*! For backpropagation, we need to calculate loss. For loss, output must be known.
    *  We obtain outputs during forward propagation and place the results in arr_output tensor. */
   Tensor_t arr_output;
   for (size_t t = 0; t < fTimeSteps; ++t) {
      arr_output.emplace_back(this->GetBatchSize(), fStateSize); // B x H
   }
   Architecture_t::Rearrange(arr_output, this->GetOutput());

   Matrix_t initState(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(initState, DNN::EInitialization::kZero); // B x H

   // This will take partial derivative of state[t] w.r.t state[t-1]
   Tensor_t arr_actgradients;
   for (size_t t = 0; t < fTimeSteps; ++t) {
      arr_actgradients.emplace_back(this->GetBatchSize(), fStateSize);
   }
   Architecture_t::Rearrange(arr_actgradients, this->GetActivationGradients());

   /*! There are total 8 different weight matrices and 4 bias vectors.
    *  Re-initialize them with zero because it should have some value. (can't be garbage values) */

   // Reset Gate.
   fWeightsResetGradients.Zero();
   fWeightsResetStateGradients.Zero();
   fResetBiasGradients.Zero();

   // Update Gate.
   fWeightsUpdateGradients.Zero();
   fWeightsUpdateStateGradients.Zero();
   fUpdateBiasGradients.Zero();

   // Candidate Gate.
   fWeightsCandidateGradients.Zero();
   fWeightsCandidateStateGradients.Zero();
   fCandidateBiasGradients.Zero();


   for (size_t t = fTimeSteps; t > 0; t--) {
      // Store the sum of gradients obtained at each timestep during backward pass.
      Architecture_t::ScaleAdd(state_gradients_backward, arr_actgradients[t-1]);
      if (t > 1) {
         const Matrix_t &prevStateActivations = arr_output[t-2];
         // During forward propagation, each gate value calculates their gradients.
         CellBackward(state_gradients_backward, prevStateActivations, 
                      this->GetResetGateTensorAt(t-1), this->GetUpdateGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1),   
                      arr_activations_backward[t-1], arr_gradients_backward[t-1],
                      fDerivativesReset[t-1], fDerivativesUpdate[t-1],
                      fDerivativesCandidate[t-1], t-1);
      } else {
         const Matrix_t &prevStateActivations = initState;
         CellBackward(state_gradients_backward, prevStateActivations, 
                      this->GetResetGateTensorAt(t-1), this->GetUpdateGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1),   
                      arr_activations_backward[t-1], arr_gradients_backward[t-1],
                      fDerivativesReset[t-1], fDerivativesUpdate[t-1],
                      fDerivativesCandidate[t-1], t-1);
        }
   }

   if (!dummy) {
      Architecture_t::Rearrange(gradients_backward, arr_gradients_backward );
   }

}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::CellBackward(Matrix_t & state_gradients_backward,
                                                          const Matrix_t & precStateActivations,  
                                                          const Matrix_t & reset_gate, const Matrix_t & update_gate,   
                                                          const Matrix_t & candidate_gate, 
                                                          const Matrix_t & input, Matrix_t & input_gradient,
                                                          Matrix_t &dr, Matrix_t &du, Matrix_t &dc,
                                                          size_t t)
-> Matrix_t &
{   
   /*! Call here GRULayerBackward() to pass parameters i.e. gradient
    *  values obtained from each gate during forward propagation. */
   return Architecture_t::GRULayerBackward(state_gradients_backward,
                                            fWeightsResetGradients, fWeightsUpdateGradients, fWeightsCandidateGradients,
                                            fWeightsResetStateGradients, fWeightsUpdateStateGradients,
                                            fWeightsCandidateStateGradients, fResetBiasGradients, fUpdateBiasGradients,
                                            fCandidateBiasGradients, dr, du, dc, 
                                            precStateActivations,
                                            reset_gate, update_gate, candidate_gate, 
                                            fWeightsResetGate, fWeightsUpdateGate, fWeightsCandidate,
                                            fWeightsResetGateState, fWeightsUpdateGateState, fWeightsCandidateState,
                                            input, input_gradient);
}
	
 //______________________________________________________________________________
template <typename Architecture_t>
auto TBasicGRULayer<Architecture_t>::InitState(DNN::EInitialization /* m */) 
-> void
{
   DNN::initialize<Architecture_t>(this->GetState(),  DNN::EInitialization::kZero);
}

 //______________________________________________________________________________
template<typename Architecture_t>
auto TBasicGRULayer<Architecture_t>::Print() const
-> void
{
   std::cout << " GRU Layer: \t ";
   std::cout << " (NInput = " << this->GetInputSize();  // input size 
   std::cout << ", NState = " << this->GetStateSize();  // hidden state size
   std::cout << ", NTime  = " << this->GetTimeSteps() << " )";  // time size
   std::cout << "\tOutput = ( " << this->GetOutput().size() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " )\n";
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::AddWeightsXMLTo(void *parent) 
-> void
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "GRULayer");
    
   // Write all other info like outputSize, cellSize, inputSize, timeSteps, rememberState
   gTools().xmlengine().NewAttr(layerxml, 0, "OutputSize", gTools().StringFromInt(this->GetStateSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "InputSize", gTools().StringFromInt(this->GetInputSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "TimeSteps", gTools().StringFromInt(this->GetTimeSteps()));
   gTools().xmlengine().NewAttr(layerxml, 0, "RememberState", gTools().StringFromInt(this->DoesRememberState()));
   
   // write weights and bias matrices
   this->WriteMatrixToXML(layerxml, "ResetWeights", this->GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "ResetStateWeights", this->GetWeightsAt(1));
   this->WriteMatrixToXML(layerxml, "ResetBiases", this->GetBiasesAt(0));
   this->WriteMatrixToXML(layerxml, "UpdateWeights", this->GetWeightsAt(2));
   this->WriteMatrixToXML(layerxml, "UpdateStateWeights", this->GetWeightsAt(3));
   this->WriteMatrixToXML(layerxml, "UpdateBiases", this->GetBiasesAt(1));
   this->WriteMatrixToXML(layerxml, "CandidateWeights", this->GetWeightsAt(4));
   this->WriteMatrixToXML(layerxml, "CandidateStateWeights", this->GetWeightsAt(5));
   this->WriteMatrixToXML(layerxml, "CandidateBiases", this->GetBiasesAt(2));
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::ReadWeightsFromXML(void *parent) 
-> void
{
	// Read weights and biases
   this->ReadMatrixXML(parent, "ResetWeights", this->GetWeightsAt(0));
   this->ReadMatrixXML(parent, "ResetStateWeights", this->GetWeightsAt(1));
   this->ReadMatrixXML(parent, "ResetBiases", this->GetBiasesAt(0));
   this->ReadMatrixXML(parent, "UpdateWeights", this->GetWeightsAt(2));
   this->ReadMatrixXML(parent, "UpdateStateWeights", this->GetWeightsAt(3));
   this->ReadMatrixXML(parent, "UpdateBiases", this->GetBiasesAt(1));
   this->ReadMatrixXML(parent, "CandidateWeights", this->GetWeightsAt(4));
   this->ReadMatrixXML(parent, "CandidateStateWeights", this->GetWeightsAt(5));
   this->ReadMatrixXML(parent, "CandidateBiases", this->GetBiasesAt(2));
}

} // namespace RNN
} // namespace DNN
} // namespace TMVA

#endif // GRU_LAYER_H