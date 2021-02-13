// @(#)root/tmva $Id$
// Author: Saurav Shekhar

/*************************************************************************
 * Copyright (C) 2017, Saurav Shekhar                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Generic tests of the RNNLayer Forward pass                     //
////////////////////////////////////////////////////////////////////

#ifndef TMVA_TEST_DNN_TEST_RNN_TEST_FWDPASS_H
#define TMVA_TEST_DNN_TEST_RNN_TEST_FWDPASS_H

#include <iostream>
#include <vector>
#include <string>

#include "../Utility.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

// template <typename Architecture>
// auto printTensor1(const std::vector<typename Architecture::Matrix_t> &A, const std::string name = "matrix")
// -> void
// {
//   std::cout << name << "\n";
//   for (size_t l = 0; l < A.size(); ++l) {
//       for (size_t i = 0; i < A[l].GetNrows(); ++i) {
//         for (size_t j = 0; j < A[l].GetNcols(); ++j) {
//             std::cout << A[l](i, j) << " ";
//         }
//         std::cout << "\n";
//       }
//       std::cout << "********\n";
//   }
// }

// template <typename Architecture>
// auto printMatrix1(const typename Architecture::Matrix_t &A, const std::string name = "matrix")
// -> void
// {
//   std::cout << name << "\n";
//   for (size_t i = 0; i < A.GetNrows(); ++i) {
//     for (size_t j = 0; j < A.GetNcols(); ++j) {
//         std::cout << A(i, j) << " ";
//     }
//     std::cout << "\n";
//   }
//   std::cout << "********\n";
// }


/*! Generate a DeepNet, test forward pass */
//______________________________________________________________________________
template <typename Architecture>
auto testForwardPass(size_t timeSteps, size_t batchSize, size_t stateSize,
                               size_t inputSize)
-> Double_t
{
   //using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;
   using RNNLayer_t = TBasicRNNLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;

   // std::vector<TMatrixT<Double_t>> XRef(timeSteps, TMatrixT<Double_t>(batchSize, inputSize));    // T x B x D
   Tensor_t XArch (timeSteps, batchSize, inputSize);    // T x B x D
   Tensor_t XRef = XArch;
   Tensor_t arr_XArch(batchSize, timeSteps, inputSize); // B x T x D

   //randomBatch(XArch);
   double x_input[] = {-1, 1, -2, 2, -3, 3, -0.5, 0.5, -0.8, 0.9, -2, 1.5, -0.2, 0.1, -0.5, 0.4, -1, 1.};
   TMatrixD Input(3, 6, x_input);
   for (size_t i = 0; i < batchSize; ++i) {
      for (size_t l = 0; l < timeSteps; ++l) {
         for (size_t m = 0; m < inputSize; ++m) {
            XArch(i, l, m) = Input(l, m);
         }
      }
   }

   Architecture::Rearrange(arr_XArch, XArch);     // B x T x D

   Net_t rnn(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
             EInitialization::kGlorotUniform);
   RNNLayer_t* layer = rnn.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, true);

   layer->Initialize();

   auto weightsInput = layer->GetWeightsInput();  // H x D
   TMatrixD weightsState = layer->GetWeightsState();  // H x H
   TMatrixD biases = layer->GetBiasesAt(0);           // H x 1
   TMatrixD state = layer->GetState();                // B x H
   TMatrixD tmp(batchSize, stateSize);

   rnn.Forward(arr_XArch);
   Tensor_t outputArch = layer->GetOutput();    // B x T x H
   Tensor_t arr_outputArch (timeSteps, batchSize, stateSize); // T x B x H

   Architecture::Rearrange(arr_outputArch, outputArch);

   Double_t maximumError = 0.0;
   for (size_t t = 0; t < timeSteps; ++t) {
      tmp.MultT(state, weightsState);
      state.MultT(XRef.At(t).GetMatrix(), weightsInput);
      state += tmp;
      // adding bias
      for (size_t i = 0; i < (size_t) state.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t) state.GetNcols(); j++) {
            state(i,j) += biases(j,0);
         }
      }
      // activation fn
      applyMatrix(state, [](double x){return tanh(x);});
      TMatrixD output = arr_outputArch.At(t).GetMatrix();

      Double_t error = maximumRelativeError(output, state);
      std::cout << "Time " << t << " Error: " << error << "\n";
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
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

   std::cout << "Compare rnn Forward output results for  bs = " << batchSize << " timeSteps = " << timeSteps
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

   Net1 rnn1(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
              (useFixedInput) ? EInitialization::kIdentity : EInitialization::kGauss);
   auto layer1 =
      rnn1.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, outputFull); // output the full sequence

   layer1->Initialize();
   layer1->Print();
   layer1->Forward(XArch1);

   Tensor1 output1 = layer1->GetOutput();

   Net2 rnn2(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
              EInitialization::kGauss);
   auto layer2 =
      rnn2.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, outputFull); // output the full sequence

   // need to initialize for setting up cudnn descriptors
   layer2->Initialize();
   layer2->CopyParameters(*layer1);

   if (debug) {
      // print weights
      std::cout << "rnn weights  " << std::endl;
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
      std::cout << "ERROR: - rnn Forward pass test failed !  - Max deviation is ";
   else
      std::cout << " Test rnn forward passed ! -   Max deviation is ";
   std::cout << max_error << std::endl;

   return max_error < tol;
}

#endif
