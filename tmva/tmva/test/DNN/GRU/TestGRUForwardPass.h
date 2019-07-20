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
// Generic tests of the GRU-Layer Forward Pass                   //
////////////////////////////////////////////////////////////////////

#ifndef TMVA_TEST_DNN_TEST_GRU_TEST_GRU_FWDPASS_H
#define TMVA_TEST_DNN_TEST_GRU_TEST_GRU_FWDPASS_H

#include <iostream>
#include <vector>

#include "../Utility.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::GRU;

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

/*! Generic sample test for forward propagation in GRU network. */
//______________________________________________________________________________
template <typename Architecture>
auto testForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize)
-> Double_t
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = std::vector<Matrix_t>;
   using GRULayer_t = TBasicGRULayer<Architecture>;
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

   Net_t gru(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   GRULayer_t* layer = gru.AddBasicGRULayer(stateSize, inputSize, timeSteps);

   layer->Initialize();

   /*! unpack weights for each gate. */
   TMatrixT<Double_t> weightsReset = layer->GetWeightsResetGate();         // H x D
   TMatrixT<Double_t> weightsCandidate = layer->GetWeightsCandidate();     // H x D
   TMatrixT<Double_t> weightsUpdate = layer->GetWeightsUpdateGate();       // H x D
   TMatrixT<Double_t> weightsResetState = layer->GetWeightsResetGateState();       // H x H
   TMatrixT<Double_t> weightsCandidateState = layer->GetWeightsCandidateState();   // H x H
   TMatrixT<Double_t> weightsUpdateState = layer->GetWeightsUpdateGateState();     // H x H
   TMatrixT<Double_t> resetBiases = layer->GetResetGateBias();                     // H x 1
   TMatrixT<Double_t> candidateBiases = layer->GetCandidateBias();                 // H x 1
   TMatrixT<Double_t> updateBiases = layer->GetUpdateGateBias();                   // H x 1
 
   /*! Get previous hidden state and previous cell state. */
   TMatrixT<Double_t> hiddenState = layer->GetState();       // B x H
  
   /*! Get each gate values. */
   TMatrixT<Double_t> resetGate = layer->GetResetGateValue();            // B x H
   TMatrixT<Double_t> candidateValue = layer->GetCandidateValue();       // B x H
   TMatrixT<Double_t> updateGate = layer->GetUpdateGateValue();          // B x H

   /*! Temporary Matrices. */
   TMatrixT<Double_t> resetTmp(batchSize, stateSize);
   TMatrixT<Double_t> candidateTmp(batchSize, stateSize);
   TMatrixT<Double_t> updateTmp(batchSize, stateSize);
   TMatrixT<Double_t> Tmp(batchSize, stateSize);

   gru.Forward(arr_XArch);

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
      resetTmp.MultT(hiddenState, weightsResetState);
      resetGate.MultT(XRef[t], weightsReset);
      resetGate += resetTmp;

      updateTmp.MultT(hiddenState, weightsUpdateState);
      updateGate.MultT(XRef[t], weightsUpdate);
      updateGate += updateTmp;

      /*! Adding bias in each gate. */
      for (size_t i = 0; i < (size_t) resetGate.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) resetGate.GetNcols(); ++j) {
            resetGate(i, j) += resetBiases(j, 0);
         }
      }
    
      for (size_t i = 0; i < (size_t) updateGate.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) updateGate.GetNcols(); ++j) {
            updateGate(i, j) += updateBiases(j, 0);
         }
      }
    

      /*! Apply activation function to each computed gate values. */
      applyMatrix(resetGate, [](double i) { return sigmoid(i); });
      applyMatrix(updateGate, [](double f) { return sigmoid(f); });

      Architecture::Hadamard(resetGate, hiddenState);
      candidateTmp.MultT(resetGate, weightsCandidateState);
      candidateValue.MultT(XRef[t], weightsCandidate);
      candidateValue += candidateTmp;

      for (size_t i = 0; i < (size_t) candidateValue.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t) candidateValue.GetNcols(); ++j) {
            candidateValue(i, j) += candidateBiases(j, 0);
        }
      }

      applyMatrix(candidateValue, [](double c) { return tanh(c); });
     
    
      /*! Computing next cell state and next hidden state. */

      Matrix_t tmp(updateGate); // H X 1
      for (size_t j = 0; j < (size_t) tmp.GetNcols(); j++) {
         for (size_t i = 0; i < (size_t) tmp.GetNrows(); i++) {
            tmp(i,j) = 1 - tmp(i,j);
         }
      }

      // Update state
      Architecture::Hadamard(hiddenState, tmp);
      Architecture::Hadamard(updateGate, candidateValue);
      Architecture::ScaleAdd(hiddenState, updateGate);


      TMatrixT<Double_t> output = arr_outputArch[t]; 
      Double_t error = maximumRelativeError(output, hiddenState);
      std::cout << "Time " << t << " Error: " << error << "\n";

      maximumError = std::max(error, maximumError);
   }

   return maximumError;
}

#endif // TMVA_TEST_DNN_TEST_RNN_TEST_GRU_FWDPASS_H
