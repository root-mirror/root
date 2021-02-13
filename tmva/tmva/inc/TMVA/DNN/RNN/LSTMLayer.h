// @(#)root/tmva/tmva/dnn/lstm:$Id$
// Author: Surya S Dwivedi 27/05/19

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class : BasicLSTMLayer                                                         *
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
// This class implements the LSTM layer. LSTM is a variant of vanilla
// RNN which is capable of learning long range dependencies.
//////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_LSTM_LAYER
#define TMVA_DNN_LSTM_LAYER

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
// Basic LSTM Layer
//______________________________________________________________________________

/** \class BasicLSTMLayer
      Generic implementation
*/
template<typename Architecture_t>
      class TBasicLSTMLayer : public VGeneralLayer<Architecture_t>
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

   size_t fStateSize;                           ///< Hidden state size for LSTM
   size_t fCellSize;                            ///< Cell state size of LSTM
   size_t fTimeSteps;                           ///< Timesteps for LSTM

   bool fRememberState;                         ///< Remember state in next pass
   bool fReturnSequence = false;                ///< Return in output full sequence or just last element

   DNN::EActivationFunction fF1;                ///< Activation function: sigmoid
   DNN::EActivationFunction fF2;                ///< Activaton function: tanh

   Matrix_t fInputValue;                        ///< Computed input gate values
   Matrix_t fCandidateValue;                    ///< Computed candidate values
   Matrix_t fForgetValue;                       ///< Computed forget gate values
   Matrix_t fOutputValue;                       ///< Computed output gate values
   Matrix_t fState;                             ///< Hidden state of LSTM
   Matrix_t fCell;                              ///< Cell state of LSTM

   Matrix_t &fWeightsInputGate;                 ///< Input Gate weights for input, fWeights[0]
   Matrix_t &fWeightsInputGateState;            ///< Input Gate weights for prev state, fWeights[1]
   Matrix_t &fInputGateBias;                    ///< Input Gate bias

   Matrix_t &fWeightsForgetGate;                ///< Forget Gate weights for input, fWeights[2]
   Matrix_t &fWeightsForgetGateState;           ///< Forget Gate weights for prev state, fWeights[3]
   Matrix_t &fForgetGateBias;                   ///< Forget Gate bias

   Matrix_t &fWeightsCandidate;                 ///< Candidate Gate weights for input, fWeights[4]
   Matrix_t &fWeightsCandidateState;            ///< Candidate Gate weights for prev state, fWeights[5]
   Matrix_t &fCandidateBias;                    ///< Candidate Gate bias

   Matrix_t &fWeightsOutputGate;                ///< Output Gate weights for input, fWeights[6]
   Matrix_t &fWeightsOutputGateState;           ///< Output Gate weights for prev state, fWeights[7]
   Matrix_t &fOutputGateBias;                   ///< Output Gate bias

   std::vector<Matrix_t> input_gate_value;      ///< input gate value for every time step
   std::vector<Matrix_t> forget_gate_value;     ///< forget gate value for every time step
   std::vector<Matrix_t> candidate_gate_value;  ///< candidate gate value for every time step
   std::vector<Matrix_t> output_gate_value;     ///< output gate value for every time step
   std::vector<Matrix_t> cell_value;            ///< cell value for every time step
   std::vector<Matrix_t> fDerivativesInput;     ///< First fDerivatives of the activations input gate
   std::vector<Matrix_t> fDerivativesForget;    ///< First fDerivatives of the activations forget gate
   std::vector<Matrix_t> fDerivativesCandidate; ///< First fDerivatives of the activations candidate gate
   std::vector<Matrix_t> fDerivativesOutput;    ///< First fDerivatives of the activations output gate

   Matrix_t &fWeightsInputGradients;            ///< Gradients w.r.t the input gate - input weights
   Matrix_t &fWeightsInputStateGradients;       ///< Gradients w.r.t the input gate - hidden state weights
   Matrix_t &fInputBiasGradients;               ///< Gradients w.r.t the input gate - bias weights
   Matrix_t &fWeightsForgetGradients;           ///< Gradients w.r.t the forget gate - input weights
   Matrix_t &fWeightsForgetStateGradients;      ///< Gradients w.r.t the forget gate - hidden state weights
   Matrix_t &fForgetBiasGradients;              ///< Gradients w.r.t the forget gate - bias weights
   Matrix_t &fWeightsCandidateGradients;        ///< Gradients w.r.t the candidate gate - input weights
   Matrix_t &fWeightsCandidateStateGradients;   ///< Gradients w.r.t the candidate gate - hidden state weights
   Matrix_t &fCandidateBiasGradients;           ///< Gradients w.r.t the candidate gate - bias weights
   Matrix_t &fWeightsOutputGradients;           ///< Gradients w.r.t the output gate - input weights
   Matrix_t &fWeightsOutputStateGradients;      ///< Gradients w.r.t the output gate - hidden state weights
   Matrix_t &fOutputBiasGradients;              ///< Gradients w.r.t the output gate - bias weights

   // Tensor representing all weights (used by cuDNN)
   Tensor_t fWeightsTensor;                     ///< Tensor for all weights
   Tensor_t fWeightGradientsTensor;             ///< Tensor for all weight gradients

   // tensors used internally for the forward and backward pass
   Tensor_t fX;  ///<  cached input tensor as T x B x I
   Tensor_t fY;  ///<  cached output tensor as T x B x S
   Tensor_t fDx; ///< cached   gradient on the input (output of backward)   as T x B x I
   Tensor_t fDy; ///< cached  activation gradient (input of backward)   as T x B x S

   TDescriptors *fDescriptors = nullptr; ///< Keeps all the RNN descriptors
   TWorkspace *fWorkspace = nullptr;     // workspace needed for GPU computation (CudNN)

public:

   /*! Constructor */
   TBasicLSTMLayer(size_t batchSize, size_t stateSize, size_t inputSize, size_t timeSteps, bool rememberState = false,
                   bool returnSequence = false,
                   DNN::EActivationFunction f1 = DNN::EActivationFunction::kSigmoid,
                   DNN::EActivationFunction f2 = DNN::EActivationFunction::kTanh, bool training = true,
                   DNN::EInitialization fA = DNN::EInitialization::kZero);

   /*! Copy Constructor */
   TBasicLSTMLayer(const TBasicLSTMLayer &);

   /*! Initialize the weights according to the given initialization
    **  method. */
   virtual void Initialize();

   /*! Initialize the hidden state and cell state method. */
   void InitState(DNN::EInitialization m = DNN::EInitialization::kZero);

   /*! Computes the next hidden state
    *  and next cell state with given input matrix. */
   void Forward(Tensor_t &input, bool isTraining = true);

   /*! Forward for a single cell (time unit) */
   void CellForward(Matrix_t &inputGateValues, const Matrix_t &forgetGateValues,
                  const Matrix_t &candidateValues, const Matrix_t &outputGateValues);

   /*! Backpropagates the error. Must only be called directly at the corresponding
    *  call to Forward(...). */
   void Backward(Tensor_t &gradients_backward,
                 const Tensor_t &activations_backward);

   /* Updates weights and biases, given the learning rate */
   void Update(const Scalar_t learningRate);

   /*! Backward for a single time unit
    *  a the corresponding call to Forward(...). */
   Matrix_t & CellBackward(Matrix_t & state_gradients_backward,
                           Matrix_t & cell_gradients_backward,
                           const Matrix_t & precStateActivations, const Matrix_t & precCellActivations,
                           const Matrix_t & input_gate, const Matrix_t & forget_gate,
                           const Matrix_t & candidate_gate, const Matrix_t & output_gate,
                           const Matrix_t & input, Matrix_t & input_gradient,
                           Matrix_t &di, Matrix_t &df, Matrix_t &dc, Matrix_t &dout, size_t t);

   /*! Decides the values we'll update (NN with Sigmoid) */
   void InputGate(const Matrix_t &input, Matrix_t &di);

   /*! Forgets the past values (NN with Sigmoid) */
   void ForgetGate(const Matrix_t &input, Matrix_t &df);

   /*! Decides the new candidate values (NN with Tanh) */
   void CandidateValue(const Matrix_t &input, Matrix_t &dc);

   /*! Computes output values (NN with Sigmoid) */
   void OutputGate(const Matrix_t &input, Matrix_t &dout);

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
   size_t GetCellSize()                const { return fCellSize; }

   inline bool DoesRememberState()       const { return fRememberState; }
   inline bool DoesReturnSequence()      const { return fReturnSequence; }

   inline DNN::EActivationFunction     GetActivationFunctionF1()        const { return fF1; }
   inline DNN::EActivationFunction     GetActivationFunctionF2()        const { return fF2; }

   const Matrix_t                    & GetInputGateValue()                const { return fInputValue; }
   Matrix_t                          & GetInputGateValue()                      { return fInputValue; }
   const Matrix_t                    & GetCandidateValue()                const { return fCandidateValue; }
   Matrix_t                          & GetCandidateValue()                      { return fCandidateValue; }
   const Matrix_t                    & GetForgetGateValue()               const { return fForgetValue; }
   Matrix_t                          & GetForgetGateValue()                     { return fForgetValue; }
   const Matrix_t                    & GetOutputGateValue()               const { return fOutputValue; }
   Matrix_t                          & GetOutputGateValue()                     { return fOutputValue; }

   const Matrix_t                    & GetState()                   const { return fState; }
   Matrix_t                          & GetState()                         { return fState; }
   const Matrix_t                    & GetCell()                    const { return fCell; }
   Matrix_t                          & GetCell()                          { return fCell; }

   const Matrix_t                    & GetWeightsInputGate()              const { return fWeightsInputGate; }
   Matrix_t                          & GetWeightsInputGate()                    { return fWeightsInputGate; }
   const Matrix_t                    & GetWeightsCandidate()              const { return fWeightsCandidate; }
   Matrix_t                          & GetWeightsCandidate()                    { return fWeightsCandidate; }
   const Matrix_t                    & GetWeightsForgetGate()             const { return fWeightsForgetGate; }
   Matrix_t                          & GetWeightsForgetGate()                   { return fWeightsForgetGate; }
   const Matrix_t                    & GetWeightsOutputGate()             const { return fWeightsOutputGate; }
   Matrix_t                          & GetWeightsOutputGate()                   { return fWeightsOutputGate; }
   const Matrix_t                    & GetWeightsInputGateState()         const { return fWeightsInputGateState; }
   Matrix_t                          & GetWeightsInputGateState()               { return fWeightsInputGateState; }
   const Matrix_t                    & GetWeightsForgetGateState()        const { return fWeightsForgetGateState; }
   Matrix_t                          & GetWeightsForgetGateState()              { return fWeightsForgetGateState; }
   const Matrix_t                    & GetWeightsCandidateState()         const { return fWeightsCandidateState; }
   Matrix_t                          & GetWeightsCandidateState()               { return fWeightsCandidateState; }
   const Matrix_t                    & GetWeightsOutputGateState()        const { return fWeightsOutputGateState; }
   Matrix_t                          & GetWeightsOutputGateState()              { return fWeightsOutputGateState; }

   const std::vector<Matrix_t>       & GetDerivativesInput()              const { return fDerivativesInput; }
   std::vector<Matrix_t>             & GetDerivativesInput()                    { return fDerivativesInput; }
   const Matrix_t                    & GetInputDerivativesAt(size_t i)    const { return fDerivativesInput[i]; }
   Matrix_t                          & GetInputDerivativesAt(size_t i)           { return fDerivativesInput[i]; }
   const std::vector<Matrix_t>       & GetDerivativesForget()              const { return fDerivativesForget; }
   std::vector<Matrix_t>             & GetDerivativesForget()                    { return fDerivativesForget; }
   const Matrix_t                    & GetForgetDerivativesAt(size_t i)    const { return fDerivativesForget[i]; }
   Matrix_t                          & GetForgetDerivativesAt(size_t i)          { return fDerivativesForget[i]; }
   const std::vector<Matrix_t>       & GetDerivativesCandidate()           const { return fDerivativesCandidate; }
   std::vector<Matrix_t>             & GetDerivativesCandidate()                 { return fDerivativesCandidate; }
   const Matrix_t                    & GetCandidateDerivativesAt(size_t i) const { return fDerivativesCandidate[i]; }
   Matrix_t                          & GetCandidateDerivativesAt(size_t i)       { return fDerivativesCandidate[i]; }
   const std::vector<Matrix_t>       & GetDerivativesOutput()              const { return fDerivativesOutput; }
   std::vector<Matrix_t>             & GetDerivativesOutput()                    { return fDerivativesOutput; }
   const Matrix_t                    & GetOutputDerivativesAt(size_t i)    const { return fDerivativesOutput[i]; }
   Matrix_t                          & GetOutputDerivativesAt(size_t i)          { return fDerivativesOutput[i]; }

   const std::vector<Matrix_t>       & GetInputGateTensor()              const { return input_gate_value; }
   std::vector<Matrix_t>             & GetInputGateTensor()                    { return input_gate_value; }
   const Matrix_t                    & GetInputGateTensorAt(size_t i)    const { return input_gate_value[i]; }
   Matrix_t                          & GetInputGateTensorAt(size_t i)           { return input_gate_value[i]; }
   const std::vector<Matrix_t>       & GetForgetGateTensor()              const { return forget_gate_value; }
   std::vector<Matrix_t>             & GetForgetGateTensor()                    { return forget_gate_value; }
   const Matrix_t                    & GetForgetGateTensorAt(size_t i)    const { return forget_gate_value[i]; }
   Matrix_t                          & GetForgetGateTensorAt(size_t i)          { return forget_gate_value[i]; }
   const std::vector<Matrix_t>       & GetCandidateGateTensor()           const { return candidate_gate_value; }
   std::vector<Matrix_t>             & GetCandidateGateTensor()                 { return candidate_gate_value; }
   const Matrix_t                    & GetCandidateGateTensorAt(size_t i) const { return candidate_gate_value[i]; }
   Matrix_t                          & GetCandidateGateTensorAt(size_t i)       { return candidate_gate_value[i]; }
   const std::vector<Matrix_t>       & GetOutputGateTensor()              const { return output_gate_value; }
   std::vector<Matrix_t>             & GetOutputGateTensor()                    { return output_gate_value; }
   const Matrix_t                    & GetOutputGateTensorAt(size_t i)    const { return output_gate_value[i]; }
   Matrix_t                          & GetOutputGateTensorAt(size_t i)          { return output_gate_value[i]; }
   const std::vector<Matrix_t>       & GetCellTensor()                    const { return cell_value; }
   std::vector<Matrix_t>             & GetCellTensor()                          { return cell_value; }
   const Matrix_t                    & GetCellTensorAt(size_t i)          const { return cell_value[i]; }
   Matrix_t                          & GetCellTensorAt(size_t i)                { return cell_value[i]; }

   const Matrix_t                   & GetInputGateBias()         const { return fInputGateBias; }
   Matrix_t                         & GetInputGateBias()               { return fInputGateBias; }
   const Matrix_t                   & GetForgetGateBias()        const { return fForgetGateBias; }
   Matrix_t                         & GetForgetGateBias()              { return fForgetGateBias; }
   const Matrix_t                   & GetCandidateBias()         const { return fCandidateBias; }
   Matrix_t                         & GetCandidateBias()               { return fCandidateBias; }
   const Matrix_t                   & GetOutputGateBias()        const { return fOutputGateBias; }
   Matrix_t                         & GetOutputGateBias()              { return fOutputGateBias; }
   const Matrix_t                   & GetWeightsInputGradients()        const { return fWeightsInputGradients; }
   Matrix_t                         & GetWeightsInputGradients()              { return fWeightsInputGradients; }
   const Matrix_t                   & GetWeightsInputStateGradients()   const { return fWeightsInputStateGradients; }
   Matrix_t                         & GetWeightsInputStateGradients()         { return fWeightsInputStateGradients; }
   const Matrix_t                   & GetInputBiasGradients()           const { return fInputBiasGradients; }
   Matrix_t                         & GetInputBiasGradients()                 { return fInputBiasGradients; }
   const Matrix_t                   & GetWeightsForgetGradients()      const { return fWeightsForgetGradients; }
   Matrix_t                         & GetWeightsForgetGradients()            { return fWeightsForgetGradients; }
   const Matrix_t                   & GetWeigthsForgetStateGradients() const { return fWeightsForgetStateGradients; }
   Matrix_t                         & GetWeightsForgetStateGradients()       { return fWeightsForgetStateGradients; }
   const Matrix_t                   & GetForgetBiasGradients()         const { return fForgetBiasGradients; }
   Matrix_t                         & GetForgetBiasGradients()               { return fForgetBiasGradients; }
   const Matrix_t                   & GetWeightsCandidateGradients()      const { return fWeightsCandidateGradients; }
   Matrix_t                         & GetWeightsCandidateGradients()            { return fWeightsCandidateGradients; }
   const Matrix_t                   & GetWeightsCandidateStateGradients() const { return fWeightsCandidateStateGradients; }
   Matrix_t                         & GetWeightsCandidateStateGradients()       { return fWeightsCandidateStateGradients; }
   const Matrix_t                   & GetCandidateBiasGradients()         const { return fCandidateBiasGradients; }
   Matrix_t                         & GetCandidateBiasGradients()               { return fCandidateBiasGradients; }
   const Matrix_t                   & GetWeightsOutputGradients()        const { return fWeightsOutputGradients; }
   Matrix_t                         & GetWeightsOutputGradients()              { return fWeightsOutputGradients; }
   const Matrix_t                   & GetWeightsOutputStateGradients()   const { return fWeightsOutputStateGradients; }
   Matrix_t                         & GetWeightsOutputStateGradients()         { return fWeightsOutputStateGradients; }
   const Matrix_t                   & GetOutputBiasGradients()           const { return fOutputBiasGradients; }
   Matrix_t                         & GetOutputBiasGradients()                 { return fOutputBiasGradients; }

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
// Basic LSTM-Layer Implementation
//______________________________________________________________________________

template <typename Architecture_t>
TBasicLSTMLayer<Architecture_t>::TBasicLSTMLayer(size_t batchSize, size_t stateSize, size_t inputSize, size_t timeSteps,
                                                 bool rememberState, bool returnSequence, DNN::EActivationFunction f1,
                                                 DNN::EActivationFunction f2, bool /* training */,
                                                 DNN::EInitialization fA)
   : VGeneralLayer<Architecture_t>(
        batchSize, 1, timeSteps, inputSize, 1, (returnSequence) ? timeSteps : 1, stateSize, 8,
        {stateSize, stateSize, stateSize, stateSize, stateSize, stateSize, stateSize, stateSize},
        {inputSize, inputSize, inputSize, inputSize, stateSize, stateSize, stateSize, stateSize}, 4,
        {stateSize, stateSize, stateSize, stateSize}, {1, 1, 1, 1}, batchSize, (returnSequence) ? timeSteps : 1,
        stateSize, fA),
     fStateSize(stateSize), fCellSize(stateSize), fTimeSteps(timeSteps), fRememberState(rememberState),
     fReturnSequence(returnSequence), fF1(f1), fF2(f2), fInputValue(batchSize, stateSize),
     fCandidateValue(batchSize, stateSize), fForgetValue(batchSize, stateSize), fOutputValue(batchSize, stateSize),
     fState(batchSize, stateSize), fCell(batchSize, stateSize), fWeightsInputGate(this->GetWeightsAt(0)),
     fWeightsInputGateState(this->GetWeightsAt(4)), fInputGateBias(this->GetBiasesAt(0)),
     fWeightsForgetGate(this->GetWeightsAt(1)), fWeightsForgetGateState(this->GetWeightsAt(5)),
     fForgetGateBias(this->GetBiasesAt(1)), fWeightsCandidate(this->GetWeightsAt(2)),
     fWeightsCandidateState(this->GetWeightsAt(6)), fCandidateBias(this->GetBiasesAt(2)),
     fWeightsOutputGate(this->GetWeightsAt(3)), fWeightsOutputGateState(this->GetWeightsAt(7)),
     fOutputGateBias(this->GetBiasesAt(3)), fWeightsInputGradients(this->GetWeightGradientsAt(0)),
     fWeightsInputStateGradients(this->GetWeightGradientsAt(4)), fInputBiasGradients(this->GetBiasGradientsAt(0)),
     fWeightsForgetGradients(this->GetWeightGradientsAt(1)),
     fWeightsForgetStateGradients(this->GetWeightGradientsAt(5)), fForgetBiasGradients(this->GetBiasGradientsAt(1)),
     fWeightsCandidateGradients(this->GetWeightGradientsAt(2)),
     fWeightsCandidateStateGradients(this->GetWeightGradientsAt(6)),
     fCandidateBiasGradients(this->GetBiasGradientsAt(2)), fWeightsOutputGradients(this->GetWeightGradientsAt(3)),
     fWeightsOutputStateGradients(this->GetWeightGradientsAt(7)), fOutputBiasGradients(this->GetBiasGradientsAt(3))
{
   for (size_t i = 0; i < timeSteps; ++i) {
      fDerivativesInput.emplace_back(batchSize, stateSize);
      fDerivativesForget.emplace_back(batchSize, stateSize);
      fDerivativesCandidate.emplace_back(batchSize, stateSize);
      fDerivativesOutput.emplace_back(batchSize, stateSize);
      input_gate_value.emplace_back(batchSize, stateSize);
      forget_gate_value.emplace_back(batchSize, stateSize);
      candidate_gate_value.emplace_back(batchSize, stateSize);
      output_gate_value.emplace_back(batchSize, stateSize);
      cell_value.emplace_back(batchSize, stateSize);
   }
   Architecture_t::InitializeLSTMTensors(this);
}

 //______________________________________________________________________________
template <typename Architecture_t>
TBasicLSTMLayer<Architecture_t>::TBasicLSTMLayer(const TBasicLSTMLayer &layer)
   : VGeneralLayer<Architecture_t>(layer),
      fStateSize(layer.fStateSize),
      fCellSize(layer.fCellSize),
      fTimeSteps(layer.fTimeSteps),
      fRememberState(layer.fRememberState),
      fReturnSequence(layer.fReturnSequence),
      fF1(layer.GetActivationFunctionF1()),
      fF2(layer.GetActivationFunctionF2()),
      fInputValue(layer.GetBatchSize(), layer.GetStateSize()),
      fCandidateValue(layer.GetBatchSize(), layer.GetStateSize()),
      fForgetValue(layer.GetBatchSize(), layer.GetStateSize()),
      fOutputValue(layer.GetBatchSize(), layer.GetStateSize()),
      fState(layer.GetBatchSize(), layer.GetStateSize()),
      fCell(layer.GetBatchSize(), layer.GetCellSize()),
      fWeightsInputGate(this->GetWeightsAt(0)),
      fWeightsInputGateState(this->GetWeightsAt(4)),
      fInputGateBias(this->GetBiasesAt(0)),
      fWeightsForgetGate(this->GetWeightsAt(1)),
      fWeightsForgetGateState(this->GetWeightsAt(5)),
      fForgetGateBias(this->GetBiasesAt(1)),
      fWeightsCandidate(this->GetWeightsAt(2)),
      fWeightsCandidateState(this->GetWeightsAt(6)),
      fCandidateBias(this->GetBiasesAt(2)),
      fWeightsOutputGate(this->GetWeightsAt(3)),
      fWeightsOutputGateState(this->GetWeightsAt(7)),
      fOutputGateBias(this->GetBiasesAt(3)),
      fWeightsInputGradients(this->GetWeightGradientsAt(0)),
      fWeightsInputStateGradients(this->GetWeightGradientsAt(4)),
      fInputBiasGradients(this->GetBiasGradientsAt(0)),
      fWeightsForgetGradients(this->GetWeightGradientsAt(1)),
      fWeightsForgetStateGradients(this->GetWeightGradientsAt(5)),
      fForgetBiasGradients(this->GetBiasGradientsAt(1)),
      fWeightsCandidateGradients(this->GetWeightGradientsAt(2)),
      fWeightsCandidateStateGradients(this->GetWeightGradientsAt(6)),
      fCandidateBiasGradients(this->GetBiasGradientsAt(2)),
      fWeightsOutputGradients(this->GetWeightGradientsAt(3)),
      fWeightsOutputStateGradients(this->GetWeightGradientsAt(7)),
      fOutputBiasGradients(this->GetBiasGradientsAt(3))
{
   for (size_t i = 0; i < fTimeSteps; ++i) {
      fDerivativesInput.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesInput[i], layer.GetInputDerivativesAt(i));

      fDerivativesForget.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesForget[i], layer.GetForgetDerivativesAt(i));

      fDerivativesCandidate.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesCandidate[i], layer.GetCandidateDerivativesAt(i));

      fDerivativesOutput.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(fDerivativesOutput[i], layer.GetOutputDerivativesAt(i));

      input_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(input_gate_value[i], layer.GetInputGateTensorAt(i));

      forget_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(forget_gate_value[i], layer.GetForgetGateTensorAt(i));

      candidate_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(candidate_gate_value[i], layer.GetCandidateGateTensorAt(i));

      output_gate_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(output_gate_value[i], layer.GetOutputGateTensorAt(i));

      cell_value.emplace_back(layer.GetBatchSize(), layer.GetStateSize());
      Architecture_t::Copy(cell_value[i], layer.GetCellTensorAt(i));
   }

   // Gradient matrices not copied
   Architecture_t::Copy(fState, layer.GetState());
   Architecture_t::Copy(fCell, layer.GetCell());

   // Copy each gate values.
   Architecture_t::Copy(fInputValue, layer.GetInputGateValue());
   Architecture_t::Copy(fCandidateValue, layer.GetCandidateValue());
   Architecture_t::Copy(fForgetValue, layer.GetForgetGateValue());
   Architecture_t::Copy(fOutputValue, layer.GetOutputGateValue());

   Architecture_t::InitializeLSTMTensors(this);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TBasicLSTMLayer<Architecture_t>::Initialize()
{
   VGeneralLayer<Architecture_t>::Initialize();

   Architecture_t::InitializeLSTMDescriptors(fDescriptors, this);
   Architecture_t::InitializeLSTMWorkspace(fWorkspace, fDescriptors, this);
}

//______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::InputGate(const Matrix_t &input, Matrix_t &di)
-> void
{
   /*! Computes input gate values according to equation:
    *  input = act(W_input . input + W_state . state + bias)
    *  activation function: sigmoid. */
   const DNN::EActivationFunction fInp = this->GetActivationFunctionF1();
   Matrix_t tmpState(fInputValue.GetNrows(), fInputValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsInputGateState);
   Architecture_t::MultiplyTranspose(fInputValue, input, fWeightsInputGate);
   Architecture_t::ScaleAdd(fInputValue, tmpState);
   Architecture_t::AddRowWise(fInputValue, fInputGateBias);
   DNN::evaluateDerivativeMatrix<Architecture_t>(di, fInp, fInputValue);
   DNN::evaluateMatrix<Architecture_t>(fInputValue, fInp);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::ForgetGate(const Matrix_t &input, Matrix_t &df)
-> void
{
   /*! Computes forget gate values according to equation:
    *  forget = act(W_input . input + W_state . state + bias)
    *  activation function: sigmoid. */
   const DNN::EActivationFunction fFor = this->GetActivationFunctionF1();
   Matrix_t tmpState(fForgetValue.GetNrows(), fForgetValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsForgetGateState);
   Architecture_t::MultiplyTranspose(fForgetValue, input, fWeightsForgetGate);
   Architecture_t::ScaleAdd(fForgetValue, tmpState);
   Architecture_t::AddRowWise(fForgetValue, fForgetGateBias);
   DNN::evaluateDerivativeMatrix<Architecture_t>(df, fFor, fForgetValue);
   DNN::evaluateMatrix<Architecture_t>(fForgetValue, fFor);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::CandidateValue(const Matrix_t &input, Matrix_t &dc)
-> void
{
   /*! Candidate value will be used to scale input gate values followed by Hadamard product.
    *  candidate_value = act(W_input . input + W_state . state + bias)
    *  activation function = tanh. */
   const DNN::EActivationFunction fCan = this->GetActivationFunctionF2();
   Matrix_t tmpState(fCandidateValue.GetNrows(), fCandidateValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsCandidateState);
   Architecture_t::MultiplyTranspose(fCandidateValue, input, fWeightsCandidate);
   Architecture_t::ScaleAdd(fCandidateValue, tmpState);
   Architecture_t::AddRowWise(fCandidateValue, fCandidateBias);
   DNN::evaluateDerivativeMatrix<Architecture_t>(dc, fCan, fCandidateValue);
   DNN::evaluateMatrix<Architecture_t>(fCandidateValue, fCan);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::OutputGate(const Matrix_t &input, Matrix_t &dout)
-> void
{
   /*! Output gate values will be used to calculate next hidden state and output values.
    *  output = act(W_input . input + W_state . state + bias)
    *  activation function = sigmoid. */
   const DNN::EActivationFunction fOut = this->GetActivationFunctionF1();
   Matrix_t tmpState(fOutputValue.GetNrows(), fOutputValue.GetNcols());
   Architecture_t::MultiplyTranspose(tmpState, fState, fWeightsOutputGateState);
   Architecture_t::MultiplyTranspose(fOutputValue, input, fWeightsOutputGate);
   Architecture_t::ScaleAdd(fOutputValue, tmpState);
   Architecture_t::AddRowWise(fOutputValue, fOutputGateBias);
   DNN::evaluateDerivativeMatrix<Architecture_t>(dout, fOut, fOutputValue);
   DNN::evaluateMatrix<Architecture_t>(fOutputValue, fOut);
}



 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::Forward(Tensor_t &input, bool  isTraining )
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
      // Tensor_t cx({1}); // not used for normal RNN
      // Tensor_t cy({1}); // not used for normal RNN

      // hx is fState - tensor are of right shape
      auto &hx = this->fState;
      //auto &cx = this->fCell;
      auto &cx = this->fCell; // pass an empty cell state
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

   // Standard CPU implementation

   // D : input size
   // H : state size
   // T : time size
   // B : batch size

   Tensor_t arrInput( fTimeSteps, this->GetBatchSize(), this->GetInputWidth());
   //Tensor_t &arrInput = this->GetX();

   Architecture_t::Rearrange(arrInput, input); // B x T x D

   Tensor_t arrOutput ( fTimeSteps, this->GetBatchSize(), fStateSize);


   if (!this->fRememberState) {
      InitState(DNN::EInitialization::kZero);
   }

   /*! Pass each gate values to CellForward() to calculate
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < fTimeSteps; ++t) {
      /* Feed forward network: value of each gate being computed at each timestep t. */
      Matrix_t arrInputMt = arrInput[t];
      InputGate(arrInputMt, fDerivativesInput[t]);
      ForgetGate(arrInputMt, fDerivativesForget[t]);
      CandidateValue(arrInputMt, fDerivativesCandidate[t]);
      OutputGate(arrInputMt, fDerivativesOutput[t]);

      Architecture_t::Copy(this->GetInputGateTensorAt(t), fInputValue);
      Architecture_t::Copy(this->GetForgetGateTensorAt(t), fForgetValue);
      Architecture_t::Copy(this->GetCandidateGateTensorAt(t), fCandidateValue);
      Architecture_t::Copy(this->GetOutputGateTensorAt(t), fOutputValue);

      CellForward(fInputValue, fForgetValue, fCandidateValue, fOutputValue);
      Matrix_t arrOutputMt = arrOutput[t];
      Architecture_t::Copy(arrOutputMt, fState);
      Architecture_t::Copy(this->GetCellTensorAt(t), fCell);
   }

   // check if full output needs to be returned
   if (fReturnSequence)
      Architecture_t::Rearrange(this->GetOutput(), arrOutput); // B x T x D
   else {
      // get T[end[]]
      Tensor_t tmp = arrOutput.At(fTimeSteps - 1); // take last time step
      // shape of tmp is  for CPU (columnwise) B x D ,   need to reshape to  make a B x D x 1
      //  and transpose it to 1 x D x B  (this is how output is expected in columnmajor format)
      tmp = tmp.Reshape( {tmp.GetShape()[0], tmp.GetShape()[1], 1});
      assert(tmp.GetSize() == this->GetOutput().GetSize());
      assert( tmp.GetShape()[0] == this->GetOutput().GetShape()[2]);  // B is last dim in output and first in tmp
      Architecture_t::Rearrange(this->GetOutput(), tmp);
      // keep array output
      fY = arrOutput;
   }
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::CellForward(Matrix_t &inputGateValues, const Matrix_t &forgetGateValues,
                                                         const Matrix_t &candidateValues, const Matrix_t &outputGateValues)
-> void
{

   // Update cell state.
   Architecture_t::Hadamard(fCell, forgetGateValues);
   Architecture_t::Hadamard(inputGateValues, candidateValues);
   Architecture_t::ScaleAdd(fCell, inputGateValues);

   Matrix_t cache(fCell.GetNrows(), fCell.GetNcols());
   Architecture_t::Copy(cache, fCell);

   // Update hidden state.
   const DNN::EActivationFunction fAT = this->GetActivationFunctionF2();
   DNN::evaluateMatrix<Architecture_t>(cache, fAT);

   /*! The Hadamard product of output_gate_value . tanh(cell_state)
    *  will be copied to next hidden state (passed to next LSTM cell)
    *  and we will update our outputGateValues also. */
   Architecture_t::Copy(fState, cache);
   Architecture_t::Hadamard(fState, outputGateValues);
}

 //____________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::Backward(Tensor_t &gradients_backward,           // B x T x D
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
         // dy is a tensor of shape (rowmajor for Cudnn): T x B x S
         // and this->ActivatuonGradients is  B x (T=1) x S
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
      //auto &cx = this->GetCell();
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
   // CPU implementation

   // gradients_backward is activationGradients of layer before it, which is input layer.
   // Currently, gradients_backward is for input(x) and not for state.
   // For the state it can be:
   Matrix_t state_gradients_backward(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(state_gradients_backward, DNN::EInitialization::kZero); // B x H


   Matrix_t cell_gradients_backward(this->GetBatchSize(), fStateSize); // B x H
   DNN::initialize<Architecture_t>(cell_gradients_backward, DNN::EInitialization::kZero); // B x H

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
   Tensor_t arr_output (  fTimeSteps, this->GetBatchSize(), fStateSize);

   Matrix_t initState(this->GetBatchSize(), fCellSize); // B x H
   DNN::initialize<Architecture_t>(initState, DNN::EInitialization::kZero); // B x H

   // This will take partial derivative of state[t] w.r.t state[t-1]

   Tensor_t arr_actgradients(fTimeSteps, this->GetBatchSize(), fStateSize);

   if (fReturnSequence) {
      Architecture_t::Rearrange(arr_output, this->GetOutput());
      Architecture_t::Rearrange(arr_actgradients, this->GetActivationGradients());
   } else {
      // here for CPU need to transpose the input activatuon gradients into the right format
      arr_output = fY;
      Architecture_t::InitializeZero(arr_actgradients);
      // need to reshape to pad a time dimension = 1 (note here is columnmajor tensors)
      Tensor_t tmp_grad = arr_actgradients.At(fTimeSteps - 1).Reshape( {this->GetBatchSize(), fStateSize, 1});
      assert(tmp_grad.GetSize() == this->GetActivationGradients().GetSize());
      assert(tmp_grad.GetShape()[0] == this->GetActivationGradients().GetShape()[2]);  // B in tmp is [0] and [2] in input act. gradients

      Architecture_t::Rearrange(tmp_grad, this->GetActivationGradients());
   }

   /*! There are total 8 different weight matrices and 4 bias vectors.
    *  Re-initialize them with zero because it should have some value. (can't be garbage values) */

   // Input Gate.
   fWeightsInputGradients.Zero();
   fWeightsInputStateGradients.Zero();
   fInputBiasGradients.Zero();

   // Forget Gate.
   fWeightsForgetGradients.Zero();
   fWeightsForgetStateGradients.Zero();
   fForgetBiasGradients.Zero();

   // Candidate Gate.
   fWeightsCandidateGradients.Zero();
   fWeightsCandidateStateGradients.Zero();
   fCandidateBiasGradients.Zero();

   // Output Gate.
   fWeightsOutputGradients.Zero();
   fWeightsOutputStateGradients.Zero();
   fOutputBiasGradients.Zero();


   for (size_t t = fTimeSteps; t > 0; t--) {
      // Store the sum of gradients obtained at each timestep during backward pass.
      Architecture_t::ScaleAdd(state_gradients_backward, arr_actgradients[t-1]);
      if (t > 1) {
         const Matrix_t &prevStateActivations = arr_output[t-2];
         const Matrix_t &prevCellActivations = this->GetCellTensorAt(t-2);
         // During forward propagation, each gate value calculates their gradients.
         Matrix_t dx = arr_gradients_backward[t-1];
         CellBackward(state_gradients_backward, cell_gradients_backward,
         	          prevStateActivations, prevCellActivations,
                      this->GetInputGateTensorAt(t-1), this->GetForgetGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1), this->GetOutputGateTensorAt(t-1),
                      arr_activations_backward[t-1], dx,
                      fDerivativesInput[t-1], fDerivativesForget[t-1],
                      fDerivativesCandidate[t-1], fDerivativesOutput[t-1], t-1);
      } else {
         const Matrix_t &prevStateActivations = initState;
         const Matrix_t &prevCellActivations = initState;
         Matrix_t dx = arr_gradients_backward[t-1];
         CellBackward(state_gradients_backward, cell_gradients_backward,
         	          prevStateActivations, prevCellActivations,
                      this->GetInputGateTensorAt(t-1), this->GetForgetGateTensorAt(t-1),
                      this->GetCandidateGateTensorAt(t-1), this->GetOutputGateTensorAt(t-1),
                      arr_activations_backward[t-1], dx,
                      fDerivativesInput[t-1], fDerivativesForget[t-1],
                      fDerivativesCandidate[t-1], fDerivativesOutput[t-1], t-1);
        }
   }

   if (!dummy) {
      Architecture_t::Rearrange(gradients_backward, arr_gradients_backward );
   }

}


 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::CellBackward(Matrix_t & state_gradients_backward,
                                                          Matrix_t & cell_gradients_backward,
                                                          const Matrix_t & precStateActivations, const Matrix_t & precCellActivations,
                                                          const Matrix_t & input_gate, const Matrix_t & forget_gate,
                                                          const Matrix_t & candidate_gate, const Matrix_t & output_gate,
                                                          const Matrix_t & input, Matrix_t & input_gradient,
                                                          Matrix_t &di, Matrix_t &df, Matrix_t &dc, Matrix_t &dout,
                                                          size_t t)
-> Matrix_t &
{
   /*! Call here LSTMLayerBackward() to pass parameters i.e. gradient
    *  values obtained from each gate during forward propagation. */


   // cell gradient for current time step
   const DNN::EActivationFunction fAT = this->GetActivationFunctionF2();
   Matrix_t cell_gradient(this->GetCellTensorAt(t).GetNrows(), this->GetCellTensorAt(t).GetNcols());
   DNN::evaluateDerivativeMatrix<Architecture_t>(cell_gradient, fAT, this->GetCellTensorAt(t));

   // cell tanh value for current time step
   Matrix_t cell_tanh(this->GetCellTensorAt(t).GetNrows(), this->GetCellTensorAt(t).GetNcols());
   Architecture_t::Copy(cell_tanh, this->GetCellTensorAt(t));
   DNN::evaluateMatrix<Architecture_t>(cell_tanh, fAT);

   return Architecture_t::LSTMLayerBackward(state_gradients_backward, cell_gradients_backward,
                                            fWeightsInputGradients, fWeightsForgetGradients, fWeightsCandidateGradients,
                                            fWeightsOutputGradients, fWeightsInputStateGradients, fWeightsForgetStateGradients,
                                            fWeightsCandidateStateGradients, fWeightsOutputStateGradients, fInputBiasGradients, fForgetBiasGradients,
                                            fCandidateBiasGradients, fOutputBiasGradients, di, df, dc, dout,
                                            precStateActivations, precCellActivations,
                                            input_gate, forget_gate, candidate_gate, output_gate,
                                            fWeightsInputGate, fWeightsForgetGate, fWeightsCandidate, fWeightsOutputGate,
                                            fWeightsInputGateState, fWeightsForgetGateState, fWeightsCandidateState,
                                            fWeightsOutputGateState, input, input_gradient,
                                            cell_gradient, cell_tanh);
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto TBasicLSTMLayer<Architecture_t>::InitState(DNN::EInitialization /* m */)
-> void
{
   DNN::initialize<Architecture_t>(this->GetState(),  DNN::EInitialization::kZero);
   DNN::initialize<Architecture_t>(this->GetCell(),  DNN::EInitialization::kZero);
}

 //______________________________________________________________________________
template<typename Architecture_t>
auto TBasicLSTMLayer<Architecture_t>::Print() const
-> void
{
   std::cout << " LSTM Layer: \t ";
   std::cout << " (NInput = " << this->GetInputSize();  // input size
   std::cout << ", NState = " << this->GetStateSize();  // hidden state size
   std::cout << ", NTime  = " << this->GetTimeSteps() << " )";  // time size
   std::cout << "\tOutput = ( " << this->GetOutput().GetFirstSize() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " )\n";
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
-> void
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "LSTMLayer");

   // Write all other info like outputSize, cellSize, inputSize, timeSteps, rememberState
   gTools().xmlengine().NewAttr(layerxml, 0, "StateSize", gTools().StringFromInt(this->GetStateSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "CellSize", gTools().StringFromInt(this->GetCellSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "InputSize", gTools().StringFromInt(this->GetInputSize()));
   gTools().xmlengine().NewAttr(layerxml, 0, "TimeSteps", gTools().StringFromInt(this->GetTimeSteps()));
   gTools().xmlengine().NewAttr(layerxml, 0, "RememberState", gTools().StringFromInt(this->DoesRememberState()));
   gTools().xmlengine().NewAttr(layerxml, 0, "ReturnSequence", gTools().StringFromInt(this->DoesReturnSequence()));

   // write weights and bias matrices
   this->WriteMatrixToXML(layerxml, "InputWeights", this->GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "InputStateWeights", this->GetWeightsAt(1));
   this->WriteMatrixToXML(layerxml, "InputBiases", this->GetBiasesAt(0));
   this->WriteMatrixToXML(layerxml, "ForgetWeights", this->GetWeightsAt(2));
   this->WriteMatrixToXML(layerxml, "ForgetStateWeights", this->GetWeightsAt(3));
   this->WriteMatrixToXML(layerxml, "ForgetBiases", this->GetBiasesAt(1));
   this->WriteMatrixToXML(layerxml, "CandidateWeights", this->GetWeightsAt(4));
   this->WriteMatrixToXML(layerxml, "CandidateStateWeights", this->GetWeightsAt(5));
   this->WriteMatrixToXML(layerxml, "CandidateBiases", this->GetBiasesAt(2));
   this->WriteMatrixToXML(layerxml, "OuputWeights", this->GetWeightsAt(6));
   this->WriteMatrixToXML(layerxml, "OutputStateWeights", this->GetWeightsAt(7));
   this->WriteMatrixToXML(layerxml, "OutputBiases", this->GetBiasesAt(3));
}

 //______________________________________________________________________________
template <typename Architecture_t>
auto inline TBasicLSTMLayer<Architecture_t>::ReadWeightsFromXML(void *parent)
-> void
{
	// Read weights and biases
   this->ReadMatrixXML(parent, "InputWeights", this->GetWeightsAt(0));
   this->ReadMatrixXML(parent, "InputStateWeights", this->GetWeightsAt(1));
   this->ReadMatrixXML(parent, "InputBiases", this->GetBiasesAt(0));
   this->ReadMatrixXML(parent, "ForgetWeights", this->GetWeightsAt(2));
   this->ReadMatrixXML(parent, "ForgetStateWeights", this->GetWeightsAt(3));
   this->ReadMatrixXML(parent, "ForgetBiases", this->GetBiasesAt(1));
   this->ReadMatrixXML(parent, "CandidateWeights", this->GetWeightsAt(4));
   this->ReadMatrixXML(parent, "CandidateStateWeights", this->GetWeightsAt(5));
   this->ReadMatrixXML(parent, "CandidateBiases", this->GetBiasesAt(2));
   this->ReadMatrixXML(parent, "OuputWeights", this->GetWeightsAt(6));
   this->ReadMatrixXML(parent, "OutputStateWeights", this->GetWeightsAt(7));
   this->ReadMatrixXML(parent, "OutputBiases", this->GetBiasesAt(3));
}

} // namespace LSTM
} // namespace DNN
} // namespace TMVA

#endif // LSTM_LAYER_H