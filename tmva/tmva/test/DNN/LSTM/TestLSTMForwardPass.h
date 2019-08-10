// @(#)root/tmva $Id$
// Author: Surya S Dwivedi 07/06/2019

/*************************************************************************
 * Copyright (C) 2019, Surya S Dwivedi                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Generic tests of the LSTM-Layer Forward Pass                   //
////////////////////////////////////////////////////////////////////

#ifndef TMVA_TEST_DNN_TEST_LSTM_TEST_LSTM_FWDPASS_H
#define TMVA_TEST_DNN_TEST_LSTM_TEST_LSTM_FWDPASS_H

#include <iostream>
#include <vector>

#include "../Utility.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::LSTM;

//______________________________________________________________________________
/* Prints out Tensor, printTensor1(A, matrix) */
template <typename Architecture>
auto printTensor1(const std::vector<typename Architecture::Matrix_t> &A, const std::string & name = "matrix")
-> void
{
   std::cout << name << "\n";
   for (size_t l = 0; l < A.size(); ++l) {
      for (size_t i = 0; i < (size_t) A[l].GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) A[l].GetNcols(); ++j) {
            std::cout << A[l](i, j) << " ";
         }
         std::cout << "\n";
      }
      std::cout << "********\n";
  } 
}

//______________________________________________________________________________
/* Prints out Matrix, printMatrix1(A, matrix) */
template <typename Architecture>
auto printMatrix1(const typename Architecture::Matrix_t &A, const std::string name = "matrix")
-> void
{
   std::cout << name << "\n";
   for (size_t i = 0; i < (size_t) A.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) A.GetNcols(); ++j) {
         std::cout << A(i, j) << " ";
      }
      std::cout << "\n";
   }
   std::cout << "********\n";
}

double sigmoid(double x) { 
   return 1 /( 1 + exp(-x)); 
}

/*! Generic sample test for forward propagation in LSTM network. */
//______________________________________________________________________________
template <typename Architecture>
auto testForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize)
-> Double_t
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = std::vector<Matrix_t>;
   using LSTMLayer_t = TBasicLSTMLayer<Architecture>;
   using Net_t = TDeepNet<Architecture>;


   //______________________________________________________________________________
   /* Input Gate: Numerical example. 
    * Reference: https://medium.com/@aidangomez/let-s-do-this-f9b699de31d9 
    * TODO: Numerical example for other gates to verify forward pass values and 
    * backward pass values. */
   //______________________________________________________________________________

   // Defining inputs.
   std::vector<TMatrixT<Double_t>> XRef(timeSteps, TMatrixT<Double_t>(batchSize, inputSize));  // T x B x D
   Tensor_t XArch, arr_XArch;


   for (size_t i = 0; i < batchSize; ++i) {
      arr_XArch.emplace_back(timeSteps, inputSize); // B x T x D
   }
   
   for (size_t i = 0; i < timeSteps; ++i) {
      randomMatrix(XRef[i]);
      XArch.emplace_back(XRef[i]);
   }

   Architecture::Rearrange(arr_XArch, XArch); // B x T x D

   Net_t lstm(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   LSTMLayer_t* layer = lstm.AddBasicLSTMLayer(stateSize, inputSize, timeSteps);

   layer->Initialize();

   /*! unpack weights for each gate. */
   TMatrixT<Double_t> weightsInput = layer->GetWeightsInputGate();         // H x D
   TMatrixT<Double_t> weightsCandidate = layer->GetWeightsCandidate();     // H x D
   TMatrixT<Double_t> weightsForget = layer->GetWeightsForgetGate();       // H x D
   TMatrixT<Double_t> weightsOutput = layer->GetWeightsOutputGate();       // H x D
   TMatrixT<Double_t> weightsInputState = layer->GetWeightsInputGateState();       // H x H
   TMatrixT<Double_t> weightsCandidateState = layer->GetWeightsCandidateState();   // H x H
   TMatrixT<Double_t> weightsForgetState = layer->GetWeightsForgetGateState();     // H x H
   TMatrixT<Double_t> weightsOutputState = layer->GetWeightsOutputGateState();     // H x H
   TMatrixT<Double_t> inputBiases = layer->GetInputGateBias();                     // H x 1
   TMatrixT<Double_t> candidateBiases = layer->GetCandidateBias();                 // H x 1
   TMatrixT<Double_t> forgetBiases = layer->GetForgetGateBias();                   // H x 1
   TMatrixT<Double_t> outputBiases = layer->GetOutputGateBias();                   // H x 1
 
   /*! Get previous hidden state and previous cell state. */
   TMatrixT<Double_t> hiddenState = layer->GetState();       // B x H
   TMatrixT<Double_t> cellState = layer->GetCell();          // B x H
  

   /*! Get each gate values. */
   TMatrixT<Double_t> inputGate = layer->GetInputGateValue();            // B x H
   TMatrixT<Double_t> candidateValue = layer->GetCandidateValue();       // B x H
   TMatrixT<Double_t> forgetGate = layer->GetForgetGateValue();          // B x H
   TMatrixT<Double_t> outputGate = layer->GetOutputGateValue();          // B x H

   /*! Temporary Matrices. */
   TMatrixT<Double_t> inputTmp(batchSize, stateSize);
   TMatrixT<Double_t> candidateTmp(batchSize, stateSize);
   TMatrixT<Double_t> forgetTmp(batchSize, stateSize);
   TMatrixT<Double_t> outputTmp(batchSize, stateSize);
   TMatrixT<Double_t> Tmp(batchSize, stateSize);

   lstm.Forward(arr_XArch);

   Tensor_t outputArch = layer->GetOutput();

   Tensor_t arr_outputArch;
   for (size_t t = 0; t < timeSteps; ++t) {
      arr_outputArch.emplace_back(batchSize, stateSize); // T x B x H
   }

   Architecture::Rearrange(arr_outputArch, outputArch); // B x T x H

   Double_t maximumError = 0.0;

   /*! Element-wise matrix multiplication of previous hidden
    *  state and weights of previous state followed by computing
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < timeSteps; ++t) {
      inputTmp.MultT(hiddenState, weightsInputState);
      inputGate.MultT(XRef[t], weightsInput);
      inputGate += inputTmp;

      candidateTmp.MultT(hiddenState, weightsCandidateState);
      candidateValue.MultT(XRef[t], weightsCandidate);
      candidateValue += candidateTmp;

      forgetTmp.MultT(hiddenState, weightsForgetState);
      forgetGate.MultT(XRef[t], weightsForget);
      forgetGate += forgetTmp;

      outputTmp.MultT(hiddenState, weightsOutputState);
      outputGate.MultT(XRef[t], weightsOutput);
      outputGate += outputTmp;

      /*! Adding bias in each gate. */
      for (size_t i = 0; i < (size_t) inputGate.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) inputGate.GetNcols(); ++j) {
            inputGate(i, j) += inputBiases(j, 0);
         }
      }
    
      for (size_t i = 0; i < (size_t) candidateValue.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) candidateValue.GetNcols(); ++j) {
            candidateValue(i, j) += candidateBiases(j, 0);
        }
      }

      for (size_t i = 0; i < (size_t) forgetGate.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) forgetGate.GetNcols(); ++j) {
            forgetGate(i, j) += forgetBiases(j, 0);
         }
      }
    
      for (size_t i = 0; i < (size_t) outputGate.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) outputGate.GetNcols(); ++j) {
           outputGate(i, j) += outputBiases(j, 0);
         }
      }

      /*! Apply activation function to each computed gate values. */
      applyMatrix(inputGate, [](double i) { return sigmoid(i); });
      applyMatrix(candidateValue, [](double c) { return tanh(c); });
      applyMatrix(forgetGate, [](double f) { return sigmoid(f); });
      applyMatrix(outputGate, [](double o) { return sigmoid(o); });
    
      /*! Computing next cell state and next hidden state. */

    
      Architecture::Hadamard(inputGate, candidateValue);
      Architecture::Hadamard(forgetGate, cellState);
      cellState = inputGate + forgetGate;

      Tmp = cellState;
      applyMatrix(Tmp, [](double y) { return tanh(y); }); 
      Architecture::Hadamard(outputGate, Tmp);
      hiddenState = outputGate;

      TMatrixT<Double_t> output = arr_outputArch[t]; 
      Double_t error = maximumRelativeError(output, hiddenState);
      std::cout << "Time " << t << " Error: " << error << "\n";

      maximumError = std::max(error, maximumError);
   }

   return maximumError;
}

#endif // TMVA_TEST_DNN_TEST_RNN_TEST_LSTM_FWDPASS_H
