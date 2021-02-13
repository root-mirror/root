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
namespace RNN
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
   using Tensor_t = typename Architecture_t::Tensor_t;

   using LayerDescriptor_t = typename Architecture_t::RecurrentDescriptor_t;
   using WeightsDescriptor_t = typename Architecture_t::FilterDescriptor_t;
   using TensorDescriptor_t = typename Architecture_t::TensorDescriptor_t;
   using HelperDescriptor_t = typename Architecture_t::DropoutDescriptor_t;

   using RNNWorkspace_t = typename Architecture_t::RNNWorkspace_t;
   using RNNDescriptors_t = typename Architecture_t::RNNDescriptors_t;

private:

   size_t fStateSize;                           ///< Hidden state size for GRU
   size_t fTimeSteps;                           ///< Timesteps for GRU

   bool fRememberState;                         ///< Remember state in next pass
   bool fReturnSequence = false;                ///< Return in output full sequence or just last element
   bool fResetGateAfter = false;                ///< GRU variant to Apply the reset gate multiplication afterwards (used by cuDNN)

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

   Matrix_t fCell;                              ///< EMpty matrix for GRU

   // Tensor representing all weights (used by cuDNN)
   Tensor_t fWeightsTensor;         ///< Tensor for all weights
   Tensor_t fWeightGradientsTensor; ///< Tensor for all weight gradients

   // tensors used internally for the forward and backward pass
   Tensor_t fX;  ///<  cached input tensor as T x B x I
   Tensor_t fY;  ///<  cached output tensor as T x B x S
   Tensor_t fDx; ///< cached   gradient on the input (output of backward)   as T x B x I
   Tensor_t fDy; ///< cached  activation gradient (input of backward)   as T x B x S

   TDescriptors *fDescriptors = nullptr; ///< Keeps all the RNN descriptors
   TWorkspace *fWorkspace = nullptr;     // workspace needed for GPU computation (CudNN)

public:

   /*! Constructor */
   TBasicGRULayer(size_t batchSize, size_t stateSize, size_t inputSize,
                   size_t timeSteps, bool rememberState = false, bool returnSequence = false,
                   bool resetGateAfter = false,
                   DNN::EActivationFunction f1 = DNN::EActivationFunction::kSigmoid,
                   DNN::EActivationFunction f2 = DNN::EActivationFunction::kTanh,
                   bool training = true, DNN::EInitialization fA = DNN::EInitialization::kZero);

   /*! Copy Constructor */
   TBasicGRULayer(const TBasicGRULayer &);

   /*! Initialize the weights according to the given initialization
    **  method. */
   virtual void Initialize();

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
                 const Tensor_t &activations_backward);

   /* Updates weights and biases, given the learning rate */
   void Update(const Scalar_t learningRate);

   /*! Backward for a single time unit
    *  a the corresponding call to Forward(...). */
   Matrix_t & CellBackward(Matrix_t & state_gradients_backward,
                           const Matrix_t & precStateActivations,
                           const Matrix_t & reset_gate, const Matrix_t & update_gate,
                           const Matrix_t & candidate_gate,
                           const Matrix_t & input, Matrix_t & input_gradient,
                           Matrix_t &dr, Matrix_t &du, Matrix_t &dc);

   /*! Decides the values we'll update (NN with Sigmoid) */
   void ResetGate(const Matrix_t &input, Matrix_t &di);

   /*! Forgets the past values (NN with Sigmoid) */
   void UpdateGate(const Matrix_t &input, Matrix_t &df);

   /*! Decides the new candidate values (NN with Tanh) */
   void CandidateValue(const Matrix_t &input, Matrix_t &dc);

   /*! Prints the info about the layer */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   void ReadWeightsFromXML(void *parent);

   /*! Getters */
   size_t GetInputSize()               const { return this->GetInputWidth(); }
   size_t GetTimeSteps()               const { return fTimeSteps; }
   size_t GetStateSize()               const { return fStateSize; }

   inline bool DoesRememberState()       const { return fRememberState; }
   inline bool DoesReturnSequence() const { return fReturnSequence; }

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
   const Matrix_t                    &GetCell()                     const { return fCell; }
   Matrix_t                          & GetCell()                          { return  fCell; }

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

   Tensor_t &GetWeightsTensor() { return fWeightsTensor; }
   const Tensor_t &GetWeightsTensor() const { return fWeightsTensor; }
   Tensor_t &GetWeightGradientsTensor() { return fWeightGradientsTensor; }
   const Tensor_t &GetWeightGradientsTensor() const { return fWeightGradientsTensor; }

   Tensor_t &GetX() { return fX; }
   Tensor_t &GetY() { return fY; }
   Tensor_t &GetDX() { return fDx; }
   Tensor_t &GetDY() { return fDy; }
};


//______________________________________________________________________________
//
// Basic GRU-Layer Implementation
//______________________________________________________________________________

template <typename Architecture_t>
TBasicGRULayer<Architecture_t>::TBasicGRULayer(size_t batchSize, size_t stateSize, size_t inputSize, size_t timeSteps,
                                               bool rememberState, bool returnSequence, bool resetGateAfter, DNN::EActivationFunction f1,
                                               DNN::EActivationFunction f2, bool /* training */,
                                               DNN::EInitialization fA)
   : VGeneralLayer<Architecture_t>(batchSize, 1, timeSteps, inputSize, 1, (returnSequence) ? timeSteps : 1, stateSize,
                                   6, {stateSize, stateSize, stateSize, stateSize, stateSize, stateSize},
                                   {inputSize, inputSize, inputSize, stateSize, stateSize, stateSize}, 3,
                                   {stateSize, stateSize, stateSize}, {1, 1, 1}, batchSize,
                                   (returnSequence) ? timeSteps : 1, stateSize, fA),
     fStateSize(stateSize), fTimeSteps(timeSteps), fRememberState(rememberState), fReturnSequence(returnSequence), fResetGateAfter(resetGateAfter),
     fF1(f1), fF2(f2), fResetValue(batchSize, stateSize), fUpdateValue(batchSize, stateSize),
     fCandidateValue(batchSize, stateSize), fState(batchSize, stateSize), fWeightsResetGate(this->GetWeightsAt(0)),
     fWeightsResetGateState(this->GetWeightsAt(3)), fResetGateBias(this->GetBiasesAt(0)),
     fWeightsUpdateGate(this->GetWeightsAt(1)), fWeightsUpdateGateState(this->GetWeightsAt(4)),
     fUpdateGateBias(this->GetBiasesAt(1)), fWeightsCandidate(this->GetWeightsAt(2)),
     fWeightsCandidateState(this->GetWeightsAt(5)), fCandidateBias(this->GetBiasesAt(2)),
     fWeightsResetGradients(this->GetWeightGradientsAt(0)), fWeightsResetStateGradients(this->GetWeightGradientsAt(3)),
     fResetBiasGradients(this->GetBiasGradientsAt(0)), fWeightsUpdateGradients(this->GetWeightGradientsAt(1)),
     fWeightsUpdateStateGradients(this->GetWeightGradientsAt(4)), fUpdateBiasGradients(this->GetBiasGradientsAt(1)),
     fWeightsCandidateGradients(this->GetWeightGradientsAt(2)),
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
   Architecture_t::InitializeGRUTensors(this);
}

 //______________________________________________________________________________
template <typename Architecture_t>
TBasicGRULayer<Architecture_t>::TBasicGRULayer(const TBasicGRULayer &layer)
   : VGeneralLayer<Architecture_t>(layer),
      fStateSize(layer.fStateSize),
      fTimeSteps(layer.fTimeSteps),
      fRememberState(layer.fRememberState),
      fReturnSequence(layer.fReturnSequence),
      fResetGateAfter(layer.fResetGateAfter),
      fF1(layer.GetActivationFunctionF1()),
      fF2(layer.GetActivationFunctionF2()),
      fResetValue(layer.GetBatchSize(), layer.GetStateSize()),
      fUpdateValue(layer.GetBatchSize(), layer.GetStateSize()),
      fCandidateValue(layer.GetBatchSize(), layer.GetStateSize()),
      fState(layer.GetBatchSize(), layer.GetStateSize()),
      fWeightsResetGate(this->GetWeightsAt(0)),
      fWeightsResetGateState(this->GetWeightsAt(3)),
      fResetGateBias(this->GetBiasesAt(0)),
      fWeightsUpdateGate(this->GetWeightsAt(1)),
      fWeightsUpdateGateState(this->GetWeightsAt(4)),
      fUpdateGateBias(this->GetBiasesAt(1)),
      fWeightsCandidate(this->GetWeightsAt(2)),
      fWeightsCandidateState(this->GetWeightsAt(5)),
      fCandidateBias(this->GetBiasesAt(2)),
      fWeightsResetGradients(this->GetWeightGradientsAt(0)),
      fWeightsResetStateGradients(this->GetWeightGradientsAt(3)),
      fResetBiasGradients(this->GetBiasGradientsAt(0)),
      fWeightsUpdateGradients(this->GetWeightGradientsAt(1)),
      fWeightsUpdateStateGradients(this->GetWeightGradientsAt(4)),
      fUpdateBiasGradients(this->GetBiasGradientsAt(1)),
      fWeightsCandidateGradients(this->GetWeightGradientsAt(2)),
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

   Architecture_t::InitializeGRUTensors(this);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TBasicGRULayer<Architecture_t>::Initialize()
{
   VGeneralLayer<Architecture_t>::Initialize();

   Architecture_t::InitializeGRUDescriptors(fDescriptors, this);
   Architecture_t::InitializeGRUWorkspace(fWorkspace, fDescriptors, this);

   //cuDNN only supports resetGate after
   if (Architecture_t::IsCudnn())
      fResetGateAfter = true;
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
   DNN::evaluateDerivativeMatrix<Architecture_t>(dr, fRst, fResetValue);
   DNN::evaluateMatrix<Architecture_t>(fResetValue, fRst);
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
   DNN::evaluateDerivativeMatrix<Architecture_t>(du, fUpd, fUpdateValue);
   DNN::evaluateMatrix<Architecture_t>(fUpdateValue, fUpd);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::CandidateValue(const Matrix_t &input, Matrix_t &dc)
-> void
{
   /*!
        vanilla GRU:
        candidate_value = act(W_input . input + W_state . (reset*state) + bias)

        but CuDNN uses reset_after variant that is faster (with bias mode = input)
        (apply reset gate multiplication after matrix multiplication)
        candidate_value = act(W_input . input + reset * (W_state . state) + bias

        activation function = tanh.

    */

   const DNN::EActivationFunction fCan = this->GetActivationFunctionF2();
   Matrix_t tmp(fCandidateValue.GetNrows(), fCandidateValue.GetNcols());
   if (!fResetGateAfter) {
      Matrix_t tmpState(fResetValue); // I think here tmpState uses fResetValue buffer
      Architecture_t::Hadamard(tmpState, fState);
      Architecture_t::MultiplyTranspose(tmp, tmpState, fWeightsCandidateState);
   } else {
      // variant GRU used in cuDNN slightly faster
      Architecture_t::MultiplyTranspose(tmp, fState, fWeightsCandidateState);
      Architecture_t::Hadamard(tmp, fResetValue);
   }
   Architecture_t::MultiplyTranspose(fCandidateValue, input, fWeightsCandidate);
   Architecture_t::ScaleAdd(fCandidateValue, tmp);
   Architecture_t::AddRowWise(fCandidateValue, fCandidateBias);
   DNN::evaluateDerivativeMatrix<Architecture_t>(dc, fCan, fCandidateValue);
   DNN::evaluateMatrix<Architecture_t>(fCandidateValue, fCan);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::Forward(Tensor_t &input, bool isTraining )
-> void
{
   // for Cudnn
   if (Architecture_t::IsCudnn()) {

      // input size is stride[1] of input tensor that is B x T x inputSize
      assert(input.GetStrides()[1] == this->GetInputSize());

      Tensor_t &x = this->fX;
      Tensor_t &y = this->fY;
      Architecture_t::Rearrange(x, input);

      const auto &weights = this->GetWeightsAt(0);

      auto &hx = this->fState;
      auto &cx = this->fCell;
      // use same for hy and cy
      auto &hy = this->fState;
      auto &cy = this->fCell;

      auto rnnDesc = static_cast<RNNDescriptors_t &>(*fDescriptors);
      auto rnnWork = static_cast<RNNWorkspace_t &>(*fWorkspace);

      Architecture_t::RNNForward(x, hx, cx, weights, y, hy, cy, rnnDesc, rnnWork, isTraining);

      if (fReturnSequence) {
         Architecture_t::Rearrange(this->GetOutput(), y); // swap B and T from y to Output
      } else {
         // tmp is a reference to y (full cudnn output)
         Tensor_t tmp = (y.At(y.GetShape()[0] - 1)).Reshape({y.GetShape()[1], 1, y.GetShape()[2]});
         Architecture_t::Copy(this->GetOutput(), tmp);
      }

      return;
   }

   // D : input size
   // H : state size
   // T : time size
   // B : batch size

   Tensor_t arrInput ( fTimeSteps, this->GetBatchSize(), this->GetInputWidth());
   // for (size_t t = 0; t < fTimeSteps; ++t) {
   //    arrInput.emplace_back(this->GetBatchSize(), this->GetInputWidth()); // T x B x D
   // }
   Architecture_t::Rearrange(arrInput, input); // B x T x D

   Tensor_t arrOutput ( fTimeSteps, this->GetBatchSize(), fStateSize );
   // for (size_t t = 0; t < fTimeSteps;++t) {
   //    arrOutput.emplace_back(this->GetBatchSize(), fStateSize); // T x B x H
   // }

   if (!this->fRememberState) {
      InitState(DNN::EInitialization::kZero);
   }

   /*! Pass each gate values to CellForward() to calculate
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < fTimeSteps; ++t) {
      /* Feed forward network: value of each gate being computed at each timestep t. */
      ResetGate(arrInput[t], fDerivativesReset[t]);
      Architecture_t::Copy(this->GetResetGateTensorAt(t), fResetValue);
      UpdateGate(arrInput[t], fDerivativesUpdate[t]);
      Architecture_t::Copy(this->GetUpdateGateTensorAt(t), fUpdateValue);

      CandidateValue(arrInput[t], fDerivativesCandidate[t]);
      Architecture_t::Copy(this->GetCandidateGateTensorAt(t), fCandidateValue);


      CellForward(fUpdateValue, fCandidateValue);

      // Architecture_t::PrintTensor(Tensor_t(fState), "state output");

      Matrix_t arrOutputMt = arrOutput[t];
      Architecture_t::Copy(arrOutputMt, fState);
   }

   if (fReturnSequence)
      Architecture_t::Rearrange(this->GetOutput(), arrOutput); // B x T x D
   else {
      // get T[end[]]
      Tensor_t tmp = arrOutput.At(fTimeSteps - 1); // take last time step
      // shape of tmp is  for CPU (columnwise) B x D ,   need to reshape to  make a B x D x 1
      //  and transpose it to 1 x D x B  (this is how output is expected in columnmajor format)
      tmp = tmp.Reshape({tmp.GetShape()[0], tmp.GetShape()[1], 1});
      assert(tmp.GetSize() == this->GetOutput().GetSize());
      assert(tmp.GetShape()[0] == this->GetOutput().GetShape()[2]); // B is last dim in output and first in tmp
      Architecture_t::Rearrange(this->GetOutput(), tmp);
      // keep array output
      fY = arrOutput;
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::CellForward(Matrix_t &updateGateValues, Matrix_t &candidateValues)
-> void
{
   Architecture_t::Hadamard(fState, updateGateValues);

   // this will reuse content of updateGateValues
   Matrix_t tmp(updateGateValues); // H X 1
   for (size_t j = 0; j < (size_t) tmp.GetNcols(); j++) {
      for (size_t i = 0; i < (size_t) tmp.GetNrows(); i++) {
         tmp(i,j) = 1 - tmp(i,j);
      }
   }

   // Update state
   Architecture_t::Hadamard(candidateValues, tmp);
   Architecture_t::ScaleAdd(fState, candidateValues);
}

//____________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::Backward(Tensor_t &gradients_backward,           // B x T x D
                                                      const Tensor_t &activations_backward)   // B x T x D
-> void
{
   // BACKWARD for CUDNN
   if (Architecture_t::IsCudnn()) {

      Tensor_t &x = this->fX;
      Tensor_t &y = this->fY;
      Tensor_t &dx = this->fDx;
      Tensor_t &dy = this->fDy;

      // input size is stride[1] of input tensor that is B x T x inputSize
      assert(activations_backward.GetStrides()[1] == this->GetInputSize());


      Architecture_t::Rearrange(x, activations_backward);

      if (!fReturnSequence) {

         // Architecture_t::InitializeZero(dy);
         Architecture_t::InitializeZero(dy);

         // Tensor_t tmp1 = y.At(y.GetShape()[0] - 1).Reshape({y.GetShape()[1], 1, y.GetShape()[2]});
         Tensor_t tmp2 = dy.At(dy.GetShape()[0] - 1).Reshape({dy.GetShape()[1], 1, dy.GetShape()[2]});

         // Architecture_t::Copy(tmp1, this->GetOutput());
         Architecture_t::Copy(tmp2, this->GetActivationGradients());
      } else {
         Architecture_t::Rearrange(y, this->GetOutput());
         Architecture_t::Rearrange(dy, this->GetActivationGradients());
      }

      // Architecture_t::PrintTensor(this->GetOutput(), "output before bwd");

      // for cudnn Matrix_t and Tensor_t are same type
      const auto &weights = this->GetWeightsTensor();
      auto &weightGradients = this->GetWeightGradientsTensor();

      // note that cudnnRNNBackwardWeights accumulate the weight gradients.
      // We need then to initialize the tensor to zero every time
      Architecture_t::InitializeZero(weightGradients);

      // hx is fState
      auto &hx = this->GetState();
      auto &cx = this->GetCell();
      // use same for hy and cy
      auto &dhy = hx;
      auto &dcy = cx;
      auto &dhx = hx;
      auto &dcx = cx;

      auto rnnDesc = static_cast<RNNDescriptors_t &>(*fDescriptors);
      auto rnnWork = static_cast<RNNWorkspace_t &>(*fWorkspace);

      Architecture_t::RNNBackward(x, hx, cx, y, dy, dhy, dcy, weights, dx, dhx, dcx, weightGradients, rnnDesc, rnnWork);

      // Architecture_t::PrintTensor(this->GetOutput(), "output after bwd");

      if (gradients_backward.GetSize() != 0)
         Architecture_t::Rearrange(gradients_backward, dx);

      return;
   }

   // gradients_backward is activationGradients of layer before it, which is input layer.
   // Currently, gradients_backward is for input(x) and not for state.
   // For the state it can be:
   Matrix_t state_gradients_backward(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(state_gradients_backward, DNN::EInitialization::kZero); // B x H

   // if dummy is false gradients_backward will be written back on the matrix
   bool dummy = false;
   if (gradients_backward.GetSize() == 0 || gradients_backward[0].GetNrows() == 0 || gradients_backward[0].GetNcols() == 0) {
      dummy = true;
   }

   Tensor_t arr_gradients_backward ( fTimeSteps, this->GetBatchSize(), this->GetInputSize());


   //Architecture_t::Rearrange(arr_gradients_backward, gradients_backward); // B x T x D
   // activations_backward is input.
   Tensor_t arr_activations_backward ( fTimeSteps, this->GetBatchSize(), this->GetInputSize());

   Architecture_t::Rearrange(arr_activations_backward, activations_backward); // B x T x D

   /*! For backpropagation, we need to calculate loss. For loss, output must be known.
    *  We obtain outputs during forward propagation and place the results in arr_output tensor. */
   Tensor_t arr_output ( fTimeSteps, this->GetBatchSize(), fStateSize);

   Matrix_t initState(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(initState, DNN::EInitialization::kZero); // B x H

   // This will take partial derivative of state[t] w.r.t state[t-1]
   Tensor_t arr_actgradients ( fTimeSteps, this->GetBatchSize(), fStateSize);

   if (fReturnSequence) {
      Architecture_t::Rearrange(arr_output, this->GetOutput());
      Architecture_t::Rearrange(arr_actgradients, this->GetActivationGradients());
   } else {
      //
      arr_output = fY;
      Architecture_t::InitializeZero(arr_actgradients);
      // need to reshape to pad a time dimension = 1 (note here is columnmajor tensors)
      Tensor_t tmp_grad = arr_actgradients.At(fTimeSteps - 1).Reshape({this->GetBatchSize(), fStateSize, 1});
      assert(tmp_grad.GetSize() == this->GetActivationGradients().GetSize());
      assert(tmp_grad.GetShape()[0] ==
             this->GetActivationGradients().GetShape()[2]); // B in tmp is [0] and [2] in input act. gradients

      Architecture_t::Rearrange(tmp_grad, this->GetActivationGradients());
   }

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
         Matrix_t dx = arr_gradients_backward[t-1];
         // During forward propagation, each gate value calculates their gradients.
         CellBackward(state_gradients_backward, prevStateActivations,
                      this->GetResetGateTensorAt(t-1), this->GetUpdateGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1),
                      arr_activations_backward[t-1], dx ,
                      fDerivativesReset[t-1], fDerivativesUpdate[t-1],
                      fDerivativesCandidate[t-1]);
      } else {
         const Matrix_t &prevStateActivations = initState;
         Matrix_t dx = arr_gradients_backward[t-1];
         CellBackward(state_gradients_backward, prevStateActivations,
                      this->GetResetGateTensorAt(t-1), this->GetUpdateGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1),
                      arr_activations_backward[t-1], dx ,
                      fDerivativesReset[t-1], fDerivativesUpdate[t-1],
                      fDerivativesCandidate[t-1]);
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
                                                          Matrix_t &dr, Matrix_t &du, Matrix_t &dc)
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
                                           input, input_gradient, fResetGateAfter);
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
   std::cout << "\tOutput = ( " << this->GetOutput().GetFirstSize() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " )\n";
}

//______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicGRULayer<Architecture_t>::AddWeightsXMLTo(void *parent)
-> void
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "GRULayer");

   // Write all other info like outputSize, cellSize, inputSize, timeSteps, rememberState
   gTools().xmlengine().NewAttr(layerxml, 0, "StateSize", gTools().StringFromInt(this->GetStateSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "InputSize", gTools().StringFromInt(this->GetInputSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "TimeSteps", gTools().StringFromInt(this->GetTimeSteps()));
   gTools().xmlengine().NewAttr(layerxml, 0, "RememberState", gTools().StringFromInt(this->DoesRememberState()));
   gTools().xmlengine().NewAttr(layerxml, 0, "ReturnSequence", gTools().StringFromInt(this->DoesReturnSequence()));
   gTools().xmlengine().NewAttr(layerxml, 0, "ResetGateAfter", gTools().StringFromInt(this->fResetGateAfter));

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

} // namespace GRU
} // namespace DNN
} // namespace TMVA

#endif // GRU_LAYER_H
