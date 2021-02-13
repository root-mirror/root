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
#include <string>

#include "../Utility.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

//______________________________________________________________________________
/* Prints out Tensor, printTensor1(A, matrix) */
template <typename Architecture>
auto printTensor1(const typename Architecture::Tensor_t &A, const std::string & name = "matrix")
-> void
{
   std::cout << name << "\n";
   for (size_t l = 0; l < A.GetFirstSize(); ++l) {
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
auto testForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize, bool debug = false,
                     double tol = 1.E-5) -> Bool_t
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using LSTMLayer_t = TBasicLSTMLayer<Architecture>;
   using Net_t = TDeepNet<Architecture>;

   std::cout << "Testing LSTM Forward for  bs = " << batchSize << " timeSteps = " << timeSteps
             << " inputSize = " << inputSize << " outputSize = " << stateSize << std::endl;

   //______________________________________________________________________________
   /* Input Gate: Numerical example.
    * Reference: https://medium.com/@aidangomez/let-s-do-this-f9b699de31d9
    * TODO: Numerical example for other gates to verify forward pass values and
    * backward pass values. */
   //______________________________________________________________________________

   // Defining inputs.
   Tensor_t XRef = Architecture::CreateTensor( timeSteps,batchSize, inputSize);
   Tensor_t arr_XArch = Architecture::CreateTensor( batchSize, timeSteps,inputSize);
   Tensor_t XArch = Architecture::CreateTensor( timeSteps,batchSize, inputSize);


   for (size_t i = 0; i < timeSteps; ++i) {
      Matrix_t m = XRef[i];
      randomMatrix(m);
   }

   Architecture::Copy(XArch, XRef); // Copy from XRef to XArch
   Architecture::Rearrange( arr_XArch, XRef);  // rearrange to B x T x D

   bool outputFullSequence = (timeSteps > 1) ? true : false;

   Net_t lstm(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   LSTMLayer_t* layer = lstm.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, outputFullSequence); // output the full sequence

   layer->Initialize();
   layer->Print();

   /*! unpack weights for each gate. */
   Matrix_t weightsInput = layer->GetWeightsInputGate();         // H x D
   Matrix_t weightsCandidate = layer->GetWeightsCandidate();     // H x D
   Matrix_t weightsForget = layer->GetWeightsForgetGate();       // H x D
   Matrix_t weightsOutput = layer->GetWeightsOutputGate();       // H x D
   Matrix_t weightsInputState = layer->GetWeightsInputGateState();       // H x H
   Matrix_t weightsCandidateState = layer->GetWeightsCandidateState();   // H x H
   Matrix_t weightsForgetState = layer->GetWeightsForgetGateState();     // H x H
   Matrix_t weightsOutputState = layer->GetWeightsOutputGateState();     // H x H
   Matrix_t inputBiases = layer->GetInputGateBias();                     // H x 1
   Matrix_t candidateBiases = layer->GetCandidateBias();                 // H x 1
   Matrix_t forgetBiases = layer->GetForgetGateBias();                   // H x 1
   Matrix_t outputBiases = layer->GetOutputGateBias();                   // H x 1

   /*! Get previous hidden state and previous cell state. */
   Matrix_t hiddenState(batchSize, stateSize);       // B x H
   Architecture::Copy(hiddenState, layer->GetState());
   Matrix_t cellState(batchSize, stateSize);          // B x H
   Architecture::Copy(cellState, layer->GetCell());

   /*! Get each gate values. */
   Matrix_t inputGate = layer->GetInputGateValue();            // B x H
   Matrix_t candidateValue = layer->GetCandidateValue();       // B x H
   Matrix_t forgetGate = layer->GetForgetGateValue();          // B x H
   Matrix_t outputGate = layer->GetOutputGateValue();          // B x H

   /*! Temporary Matrices. */
   Matrix_t inputTmp(batchSize, stateSize);
   Matrix_t candidateTmp(batchSize, stateSize);
   Matrix_t forgetTmp(batchSize, stateSize);
   Matrix_t outputTmp(batchSize, stateSize);
   Matrix_t Tmp(batchSize, stateSize);
   lstm.Forward(arr_XArch);

   Tensor_t outputArch = layer->GetOutput();

   Tensor_t arr_outputArch = Architecture::CreateTensor( timeSteps, batchSize, stateSize);


   Architecture::Rearrange(arr_outputArch, outputArch); // B x T x H

   Double_t maximumError = 0.0;
   Architecture::InitializeZero(hiddenState);
   Architecture::InitializeZero(cellState);

   /*! Element-wise matrix multiplication of previous hidden
    *  state and weights of previous state followed by computing
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < timeSteps; ++t) {

      Architecture::MultiplyTranspose(inputTmp, hiddenState, weightsInputState);
      Architecture::MultiplyTranspose(inputGate, XRef[t], weightsInput);
      Architecture::ScaleAdd(inputGate, inputTmp);

      Architecture::MultiplyTranspose(candidateTmp, hiddenState, weightsCandidateState);
      Architecture::MultiplyTranspose(candidateValue, XRef[t], weightsCandidate);
      Architecture::ScaleAdd(candidateValue, candidateTmp);

      Architecture::MultiplyTranspose(forgetTmp, hiddenState, weightsForgetState);
      Architecture::MultiplyTranspose(forgetGate, XRef[t], weightsForget);
      Architecture::ScaleAdd(forgetGate, forgetTmp);

      Architecture::MultiplyTranspose(outputTmp, hiddenState, weightsOutputState);
      Architecture::MultiplyTranspose(outputGate, XRef[t], weightsOutput);
      Architecture::ScaleAdd(outputGate, outputTmp);

      /*! Adding bias in each gate. */
      Architecture::AddRowWise(inputGate, inputBiases);
      Architecture::AddRowWise(candidateValue, candidateBiases);
      Architecture::AddRowWise(forgetGate, forgetBiases);
      Architecture::AddRowWise(outputGate, outputBiases);


      /*! Apply activation function to each computed gate values. */
      applyMatrix(inputGate, [](double i) { return sigmoid(i); });
      applyMatrix(candidateValue, [](double c) { return tanh(c); });
      applyMatrix(forgetGate, [](double f) { return sigmoid(f); });
      applyMatrix(outputGate, [](double o) { return sigmoid(o); });

      /*! Computing next cell state and next hidden state. */
      Architecture::Hadamard(inputGate, candidateValue);
      Architecture::Hadamard(forgetGate, cellState);
      Architecture::Copy(cellState, inputGate);
      Architecture::ScaleAdd(cellState, forgetGate);

      Architecture::Copy(Tmp, cellState);
      applyMatrix(Tmp, [](double y) { return tanh(y); });
      Architecture::Hadamard(outputGate, Tmp);
      Architecture::Copy(hiddenState, outputGate);

      if (debug) {
         Architecture::PrintTensor(Tensor_t(hiddenState), "Expected output at current time");
         Architecture::PrintTensor(arr_outputArch[t], "LSTM re-arranged forward output");
         Architecture::PrintTensor(outputArch, "LSTM raw forward output");
      }
      Matrix_t output = arr_outputArch[t];
      Double_t error = maximumRelativeError(output, hiddenState);
      std::cout << "Time " << t << " Error: " << error << "\n";

      maximumError = std::max(error, maximumError);
   }

   if (maximumError > tol)
      std::cout << "ERROR: - LSTM Forward pass test failed !  - Max deviation is ";
   else
      std::cout << " Test LSTM forward passed ! -   Max deviation is ";

   std::cout << maximumError << std::endl;

   return maximumError < tol;

}
//______________________________________________________________________________
template <typename Arch1, typename Arch2>
auto CompareForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize, bool outputFull = true,
                        bool useFixedInput = false, bool debug = false,
                        double tol = 1.E-5) -> Bool_t
{
   using Scalar1 = typename Arch1::Scalar_t;
   using Matrix1 = typename Arch1::Matrix_t;
   using Tensor1 = typename Arch1::Tensor_t;
   using Net1 = TDeepNet<Arch1>;

   using Tensor2 = typename Arch2::Tensor_t;
   using Net2 = TDeepNet<Arch2>;

   std::cout << "Compare LSTM Forward output results for  bs = " << batchSize << " timeSteps = " << timeSteps
             << " inputSize = " << inputSize << " outputSize = " << stateSize;
   std::cout << " using Arch1 = " << typeid(Arch1).name() << " and Arch2 " << typeid(Arch2).name() << std::endl;

   // Defining inputs.
   Tensor1 XArch1 = Arch1::CreateTensor(batchSize, timeSteps, inputSize);
   Tensor2 XArch2 = Arch2::CreateTensor(batchSize, timeSteps, inputSize);

   if (useFixedInput) { // shuld use t = 2 input = 2 bs = 1
      Scalar1 xinput[] = {-0.1, 0.5, -0.5, 0.9, -0.3, 1.0};
      R__ASSERT(batchSize == 1);
      R__ASSERT(Arch1::GetTensorLayout() == TMVA::Experimental::MemoryLayout::ColumnMajor);
      // assume Arch1 is column major
      XArch1 = Tensor1(xinput, {timeSteps, inputSize, 1}, TMVA::Experimental::MemoryLayout::ColumnMajor);
   } else {
      for (size_t i = 0; i < batchSize; ++i) {
         Matrix1 m = XArch1[i];
         randomMatrix(m);
      }
   }

   //copy diff arch does not work for this tensors shape
   // if different layout one is transpose
   //Arch2::CopyDiffArch(XArch2, XArch1);

   //Architecture::Rearrange(arr_XArch, XRef); // rearrange to B x T x D
   for (size_t i = 0; i < batchSize; ++i) {
      auto mx1 = XArch1.At(i).GetMatrix();
      TMatrixT<typename Arch1::Scalar_t> tmat = mx1;
      // ned to transpose if different layout
      if (Arch2::GetTensorLayout() != Arch1::GetTensorLayout())
         tmat.T();
      auto mx2 = XArch2.At(i);
      Arch2::Copy(mx2, Tensor2(tmat));
   }

      Net1 lstm1(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
                 (useFixedInput) ? EInitialization::kIdentity : EInitialization::kGauss);
      auto layer1 =
         lstm1.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, outputFull); // output the full sequence

      layer1->Initialize();
      layer1->Print();
      layer1->Forward(XArch1);

      Tensor1 output1 = layer1->GetOutput();

      Net2 lstm2(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
                 EInitialization::kGauss);
      auto layer2 =
         lstm2.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, outputFull); // output the full sequence

      // need to initialize for setting up cudnn descriptors
      layer2->Initialize();
      layer2->CopyParameters(*layer1);

      if (debug) {
         // print weights
         std::cout << "LSTM weights  " << std::endl;
         for (size_t i = 0; i < layer1->GetWeights().size(); ++i) {
            std::cout << "Weights component - " << i << std::endl;
            Arch1::PrintTensor(Tensor1(layer1->GetWeightsAt(i)), "weights Arch1");
            Arch2::PrintTensor(Tensor2(layer2->GetWeightsAt(i)), "weights Arch2");
         }
      }

      layer2->Forward(XArch2);

      Tensor2 output2 = layer2->GetOutput();

      if (debug) {
         Arch1::PrintTensor(XArch1, "Input of Architecture 1");
         Arch2::PrintTensor(XArch2, "Input of Architecture 2");
         Arch1::PrintTensor(output1, "Output of Architecture 1");
         Arch2::PrintTensor(output2, "Output of Architecture 2");
      }

      // shape output CPU (Column major ) is  T x S x B
      // shape output GPU (row major ) is B x T x S
      // so I need to transpose

      assert(batchSize == output1.GetFirstSize());
      assert(batchSize == output2.GetFirstSize());

      Double_t max_error = 0;
      for (size_t i = 0; i < batchSize; ++i) {
         auto m1 = output1.At(i).GetMatrix();
         auto m2 = output2.At(i).GetMatrix();
         TMatrixT<typename Arch1::Scalar_t> mat1 = m1;
         TMatrixT<typename Arch2::Scalar_t> mat2 = m2;
         // need to transpose if different layout
         if (Arch2::GetTensorLayout() != Arch1::GetTensorLayout())
            mat2.T();
         if (debug) {
            mat1.Print();
            mat2.Print();
         }

         Double_t error = maximumRelativeError(mat1, mat2);
         if (error > max_error)
            max_error = error;
      }

      if (max_error > tol)
         std::cout << "ERROR: - LSTM Forward pass test failed !  - Max deviation is ";
      else
         std::cout << " Test LSTM forward passed ! -   Max deviation is ";
      std::cout << max_error << std::endl;

      return max_error < tol;
   }

#endif // TMVA_TEST_DNN_TEST_RNN_TEST_LSTM_FWDPASS_H
