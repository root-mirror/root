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

/*! Generic sample test for forward propagation in GRU network. */
//______________________________________________________________________________
template <typename Architecture>
auto testForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize, bool debug = false, double tol = 1.E-5)
-> Bool_t
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using GRULayer_t = TBasicGRULayer<Architecture>;
   using Net_t = TDeepNet<Architecture>;


   //______________________________________________________________________________
   /* Input Gate: Numerical example.
    * Reference: https://medium.com/@aidangomez/let-s-do-this-f9b699de31d9
    * TODO: Numerical example for other gates to verify forward pass values and
    * backward pass values. */
   //______________________________________________________________________________

   // Defining inputs.
   Tensor_t XRef;
   Tensor_t XArch, arr_XArch;

   std::cout << "Testing GRU Forward for  bs = " << batchSize << " timeSteps = " << timeSteps
             << " inputSize = " << inputSize << " outputSize = " << stateSize << std::endl;

   XRef = Architecture::CreateTensor( timeSteps,batchSize, inputSize);
   XArch = Architecture::CreateTensor( batchSize, timeSteps,inputSize);
   arr_XArch = Architecture::CreateTensor( timeSteps,batchSize, inputSize);


   for (size_t i = 0; i < timeSteps; ++i) {
      Matrix_t m = XRef[i];
      randomMatrix(m);
   }

   Architecture::Copy(arr_XArch, XRef); // B x T x D
   Architecture::Rearrange( XArch, XRef);  // Copy from XRef to XArch

   Net_t gru(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   // make a GRU returning the full sequence
   bool outputFullSequence = (timeSteps > 1) ? true : false;
   GRULayer_t* layer = gru.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, outputFullSequence);

   layer->Initialize();
   layer->Print();

   /*! unpack weights for each gate. */
   Matrix_t weightsReset = layer->GetWeightsResetGate();         // H x D
   Matrix_t weightsCandidate = layer->GetWeightsCandidate();     // H x D
   Matrix_t weightsUpdate = layer->GetWeightsUpdateGate();       // H x D
   Matrix_t weightsResetState = layer->GetWeightsResetGateState();       // H x H
   Matrix_t weightsCandidateState = layer->GetWeightsCandidateState();   // H x H
   Matrix_t weightsUpdateState = layer->GetWeightsUpdateGateState();     // H x H
   Matrix_t resetBiases = layer->GetResetGateBias();                     // H x 1
   Matrix_t candidateBiases = layer->GetCandidateBias();                 // H x 1
   Matrix_t updateBiases = layer->GetUpdateGateBias();                   // H x 1

   /*! Get previous hidden state and previous cell state. */
   Matrix_t hiddenState(batchSize, stateSize);       // B x H
   Architecture::Copy(hiddenState, layer->GetState());

   /*! Get each gate values. */
   Matrix_t resetGate = layer->GetResetGateValue();            // B x H
   Matrix_t candidateValue = layer->GetCandidateValue();       // B x H
   Matrix_t updateGate = layer->GetUpdateGateValue();          // B x H

   /*! Temporary Matrices. */
   Matrix_t resetTmp(batchSize, stateSize);
   Matrix_t candidateTmp(batchSize, stateSize);
   Matrix_t updateTmp(batchSize, stateSize);
   Matrix_t tmp(batchSize, stateSize);

   gru.Forward(XArch);

   Tensor_t outputArch = layer->GetOutput();

   Tensor_t arr_outputArch = Architecture::CreateTensor( timeSteps,batchSize, stateSize);

   //Architecture::PrintTensor(outputArch,"GRU output");

   Architecture::Rearrange(arr_outputArch, outputArch); // B x T x H

   if (debug)
      Architecture::PrintTensor(arr_outputArch,"GRU FOrward output");

   Double_t maximumError = 0.0;

   Architecture::InitializeZero(hiddenState);

   /*! Element-wise matrix multiplication of previous hidden
    *  state and weights of previous state followed by computing
    *  next hidden state and next cell state. */
   for (size_t t = 0; t < timeSteps; ++t) {

      Architecture::MultiplyTranspose(resetTmp, hiddenState, weightsResetState);
      Architecture::MultiplyTranspose(resetGate, XRef[t], weightsReset);
      Architecture::ScaleAdd(resetGate, resetTmp);

      Architecture::MultiplyTranspose(updateTmp, hiddenState, weightsUpdateState);
      Architecture::MultiplyTranspose(updateGate, XRef[t], weightsUpdate);
      Architecture::ScaleAdd(updateGate, updateTmp);

      Architecture::AddRowWise(resetGate, resetBiases);
      Architecture::AddRowWise(updateGate, updateBiases);

      /*! Apply activation function to each computed gate values. */
      applyMatrix(resetGate, [](double i) { return sigmoid(i); });
      applyMatrix(updateGate, [](double f) { return sigmoid(f); });

      if (debug)
         Architecture::PrintTensor(Tensor_t(resetGate), "reset gate");

      Architecture::Hadamard(resetGate, hiddenState);

      Architecture::MultiplyTranspose(candidateTmp, resetGate, weightsCandidateState);
      Architecture::MultiplyTranspose(candidateValue, XRef[t], weightsCandidate);
      Architecture::ScaleAdd(candidateValue, candidateTmp);

      Architecture::AddRowWise(candidateValue, candidateBiases);
      applyMatrix(candidateValue, [](double c) { return tanh(c); });

      /*! Computing next cell state and next hidden state. */

      if (debug) {
         Architecture::PrintTensor(Tensor_t(updateGate), "update gate");
         Architecture::PrintTensor(Tensor_t(candidateValue), "candidate value");
      }

      Architecture::Copy(tmp, updateGate); // H X 1
      for (size_t j = 0; j < (size_t) tmp.GetNcols(); j++) {
         for (size_t i = 0; i < (size_t) tmp.GetNrows(); i++) {
            tmp(i,j) = 1 - tmp(i,j);
         }
      }

      if (debug) Architecture::PrintTensor(Tensor_t(tmp), "tmp");

      // Update state
      Architecture::Hadamard(updateGate, hiddenState);
      Architecture::Hadamard(candidateValue, tmp);
      Architecture::ScaleAdd(candidateValue, updateGate);
      Architecture::Copy(hiddenState, candidateValue);


      TMatrixT<Double_t> output = arr_outputArch[t];


      if (debug) {
         std::cout << "GRU output at current time " << std::endl;
         output.Print();
         Architecture::PrintTensor(Tensor_t(hiddenState),"expected output at current time");
      }

      Double_t error = maximumRelativeError(output, hiddenState);
      std::cout << "Time " << t << " Delta: " << error << "\n";

      maximumError = std::max(error, maximumError);
   }

   if (maximumError > tol)
      std::cout << "ERROR: - GRU Forward pass test failed !  - Max dev is ";
   else
      std::cout << " Test GRU forward passed ! -   Max dev is ";

   std::cout << maximumError << std::endl;
   return maximumError < tol;
}

//______________________________________________________________________________
template <typename Arch1, typename Arch2>
auto CompareForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize, size_t inputSize, bool outputFull = true,
                        bool useFixedInput = false, bool debug = false, double tol = 1.E-5) -> Bool_t
{
   using Scalar1 = typename Arch1::Scalar_t;
   using Matrix1 = typename Arch1::Matrix_t;
   using Tensor1 = typename Arch1::Tensor_t;
   using Net1 = TDeepNet<Arch1>;

   using Tensor2 = typename Arch2::Tensor_t;
   using Net2 = TDeepNet<Arch2>;

   std::cout << "Compare GRU Forward output results for  bs = " << batchSize << " timeSteps = " << timeSteps
             << " inputSize = " << inputSize << " outputSize = " << stateSize;
   std::cout << " using Arch1 = " << typeid(Arch1).name() << " and Arch2 " << typeid(Arch2).name() << std::endl;

   // Defining inputs.
   Tensor1 XArch1 = Arch1::CreateTensor(batchSize, timeSteps, inputSize);
   Tensor2 XArch2 = Arch2::CreateTensor(batchSize, timeSteps, inputSize);

   if (useFixedInput) { // shuld use t = 2 input = 2 bs = 1
      Scalar1 xinput[] = {-0.8, 0.5, -0.5, 0.9, -0.3, 1.0};
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

   // copy diff arch does not work for this tensors shape
   // if different layout one is transpose
   // Arch2::CopyDiffArch(XArch2, XArch1);

   // Architecture::Rearrange(arr_XArch, XRef); // rearrange to B x T x D
   for (size_t i = 0; i < batchSize; ++i) {
      auto mx1 = XArch1.At(i).GetMatrix();
      TMatrixT<typename Arch1::Scalar_t> tmat = mx1;
      // ned to transpose if different layout
      if (Arch2::GetTensorLayout() != Arch1::GetTensorLayout())
         tmat.T();
      auto mx2 = XArch2.At(i);
      Arch2::Copy(mx2, Tensor2(tmat));
   }

   Net1 gru1(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
              EInitialization::kGauss);
   auto layer1 =
      gru1.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, outputFull, true); // output the full sequence and use resetgateAfter

   layer1->Initialize();

   // in case of fixed input set all GRU weights equal to 1.
   if (useFixedInput) {
      for (size_t i = 0; i < layer1->GetWeights().size(); ++i) {
         auto &w = layer1->GetWeightsAt(i);
         float _w0 = - float(w.GetNcols());
         for (size_t j = 0; j < w.GetNrows(); ++j) {
            for (size_t k = 0; k < w.GetNcols(); ++k) {
               if (_w0 == 0)
                  _w0 = _w0 + 1;
               w(j, k) = _w0;
               _w0 = _w0 + 1;
            }
         }
      }
   }

   layer1->Print();
   layer1->Forward(XArch1);

   Tensor1 output1 = layer1->GetOutput();

   Net2 gru2(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
              EInitialization::kGauss);
   auto layer2 =
      gru2.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, outputFull); // output the full sequence

   // need to initialize for setting up cudnn descriptors
   layer2->Initialize();
   layer2->CopyParameters(*layer1);

   if (debug) {
      // print weights
      std::cout << "GRU weights  " << std::endl;
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
      std::cout << "ERROR: - GRU Forward pass test failed !  - Max deviation is ";
   else
      std::cout << " Test GRU forward passed ! -   Max deviation is ";
   std::cout << max_error << std::endl;

   return max_error < tol;
}

#endif // TMVA_TEST_DNN_TEST_RNN_TEST_GRU_FWDPASS_H
