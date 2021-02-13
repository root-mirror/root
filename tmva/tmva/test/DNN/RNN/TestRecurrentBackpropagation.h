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
// Generic tests of the RNNLayer Backward pass                    //
////////////////////////////////////////////////////////////////////

#ifndef TMVA_TEST_DNN_TEST_RNN_TEST_BWDPASS_H
#define TMVA_TEST_DNN_TEST_RNN_TEST_BWDPASS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "../Utility.h"
#include "Math/Functor.h"
#include "Math/RichardsonDerivator.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

template <typename Architecture>
auto printTensor(const typename Architecture::Tensor_t &A, const std::string name = "matrix")
-> void
{
   Architecture::PrintTensor(A,name);
   std::cout << "********\n";
}

template <typename Architecture>
auto printMatrix(const typename Architecture::Matrix_t &A, const std::string name = "matrix")
-> void
{
  Architecture::PrintTensor(typename Architecture::Tensor_t(A),name);
  std::cout << "********\n";
}

/*! Compute the loss of the net as a function of the weight at index (i,j) in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_weight(TDeepNet<Architecture> &net, typename Architecture::Tensor_t & X,
                         const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                         size_t k, size_t i, size_t j, typename Architecture::Scalar_t xvalue) ->
   typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;

    Scalar_t prev_value = net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j);
    net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j) = xvalue;
    Scalar_t res = net.Loss(X, Y, W, false, false);
    net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j) = prev_value;

    //std::cout << "compute loss for weight  " << std::setprecision(12) << xvalue << "  " << prev_value << " result " << res << std::setprecision(6) << std::endl;
    return res;
}

/*! Compute the loss of the net as a function of the weight at index i in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_bias(TDeepNet<Architecture> &net, typename Architecture::Tensor_t & X,
                       const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                       size_t k, size_t i, typename Architecture::Scalar_t xvalue) -> typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;

    Scalar_t prev_value = net.GetLayerAt(l)->GetBiasesAt(k).operator()(i,0);
    net.GetLayerAt(l)->GetBiasesAt(k).operator()(i,0) = xvalue;
    Scalar_t res = net.Loss(X, Y, W, false, false);
    net.GetLayerAt(l)->GetBiasesAt(k).operator()(i,0) = prev_value;
    return res;
}

/*! Generate a DeepNet, test backward pass */
//______________________________________________________________________________
template <typename Architecture>
bool testRecurrentBackpropagation(size_t timeSteps, size_t batchSize, size_t stateSize,
                                  size_t inputSize, typename Architecture::Scalar_t dx = 1.E-5,
                                  std::vector<bool> options = {}, bool debug = false)

{
   bool failed = false;
   const int nOpts = 5;  // size of options
   if (options.size() <  nOpts) options.resize(nOpts);
   bool randomInput = !options[0];
   bool addDenseLayer = options[1];
   bool addExtraRNN = options[2];
   bool returnLastSequence = options[3];
   bool generateWeights = options[4];

   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;
   using RNNLayer_t = TBasicRNNLayer<Architecture>;
   //using DenseLayer_t = TDenseLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;
   using Scalar_t = typename Architecture::Scalar_t;

   if (debug) std::cout << std::endl;
   std::cout
      << "******************************************************************************************************\n";
   std::cout << "Testing Weight Backprop using RNN with batchsize = " << batchSize << " input = " << inputSize
             << " state = " << stateSize << " time = " << timeSteps;
   if (randomInput)
      std::cout << " using a random input";
   else
      std::cout << " with a fixed input";
   if (addDenseLayer)
      std::cout << " and a dense layer";
   if (addExtraRNN)
      std::cout << " and an extra RNN";
   if (returnLastSequence)
      std::cout << " and full output";
   std::cout << std::endl;
   std::cout
      << "******************************************************************************************************\n";
   if (debug) std::cout << std::endl;

   //std::vector<Matrix_t<Double_t>> XRef(batchSize, Matrix_t<Double_t>(timeSteps, inputSize));    // B x T x D
   Tensor_t XArch = Architecture::CreateTensor ( batchSize, timeSteps, inputSize); // B x T x D

   // for random input (default)
   if (randomInput) {
   for (size_t i = 0; i < batchSize; ++i) {
         for (size_t l = 0; l < (size_t) timeSteps; ++l) {
            for (size_t m = 0; m < (size_t) inputSize; ++m) {
               //XArch[i](l, m) = i + l + m;
               XArch(i, l, m) = gRandom->Uniform(-1,1);
            }
         }
      }
   }
   else {
      R__ASSERT(inputSize <= 6);
      R__ASSERT(timeSteps <= 3);
      R__ASSERT(batchSize <= 1);
      double x_input[] = {-1,   1,  -2,  2, -3,  3 ,
                          -0.5, 0.5,-0.8,0.9, -2, 1.5,
                          -0.2, 0.1,-0.5,0.4, -1, 1.};

      TMatrixD Input(3,6,x_input);
      for (size_t i = 0; i < batchSize; ++i) {
         for (size_t l = 0; l < timeSteps; ++l) {
            for (size_t m = 0; m < inputSize; ++m) {
               XArch(i,l,m) = Input(l,m);
            }
         }
      }
   }
   if (debug) printTensor<Architecture>(XArch,"Input Tensor");

   size_t outputSize = (returnLastSequence) ? timeSteps*stateSize : stateSize;
   if (addDenseLayer) outputSize = 1;

   Matrix_t Y(batchSize, outputSize), weights(batchSize, 1);
   //randomMatrix(Y);
   for (size_t i = 0; i < (size_t) Y.GetNrows(); ++i) {
     for (size_t j = 0; j < (size_t) Y.GetNcols(); ++j) {
        Y(i, j) = gRandom->Integer(2); //1;// (i + j)/2.0 - 0.75;
      }
   }
   fillMatrix(weights, 1.0);

   // return full output sequence of first RNN in case of an extra one or in case we required
   bool returnFirstSequence = addExtraRNN || returnLastSequence;
   TMVA::DNN::EActivationFunction afunc = TMVA::DNN::EActivationFunction::kTanh;

   Net_t rnn(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
             EInitialization::kGlorotUniform);
   RNNLayer_t* rnnlayer = rnn.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, returnFirstSequence, afunc);  // don't use tanh in test due to limited vdt precision
   RNNLayer_t *rnnlayer2 = nullptr;
   if (addExtraRNN)
      rnnlayer2 = rnn.AddBasicRNNLayer(stateSize, stateSize, timeSteps, false, returnLastSequence, // do not return state at end
                           TMVA::DNN::EActivationFunction::kTanh);
   //add reshape layer to flatten output
   size_t outputRNNSize = (returnLastSequence) ? timeSteps * stateSize : stateSize;
   rnn.AddReshapeLayer(1, 1, outputRNNSize, true);


   // DenseLayer_t * dlayer1 = nullptr;
   // DenseLayer_t * dlayer2 = nullptr;
   if (addDenseLayer) {
      //dlayer1 = rnn.AddDenseLayer(10, TMVA::DNN::EActivationFunction::kSigmoid);
      rnn.AddDenseLayer(1, TMVA::DNN::EActivationFunction::kTanh);
   }

   rnn.Initialize();

   if (debug) rnn.Print();

   // if (! Architecture::IsCudnn()) {
   auto &wi = rnnlayer->GetWeightsInput();
   auto &wh = rnnlayer->GetWeightsState();
   auto &b = rnnlayer->GetBiasesState();
   if (generateWeights) {
      gRandom->SetSeed(123);
      for (size_t i = 0; i < stateSize; ++i) {
         for (size_t j = 0; j < inputSize; ++j) {
            wi(i, j) = gRandom->Uniform(-1, 1);
         }
         wi(i, i) = 1.;
      }

      for (size_t i = 0; i < stateSize; ++i) {
         for (size_t j = 0; j < stateSize; ++j) {
            wh(i, j) = gRandom->Uniform(-1, 1);
         }
         wh(i, i) = 0.5;
      }

      for (size_t i = 0; i < (size_t)b.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t)b.GetNcols(); ++j) {
            b(i, j) = gRandom->Uniform(-0.5, 0.5);
         }
      }
   }
   if (debug) {
      printTensor<Architecture>(wi, "Input weights");
      printTensor<Architecture>(wh, "State weights");
      printTensor<Architecture>(b, "State Bias");
      if (addExtraRNN) {
         assert(rnnlayer2);
         printTensor<Architecture>(rnnlayer2->GetWeightsInput(), "Input weights RNN 2");
         printTensor<Architecture>(rnnlayer2->GetWeightsState(), "State weights RNN 2");
         printTensor<Architecture>(rnnlayer2->GetBiasesState(), "State Bias RNN 2");
      }
   }
//}
// else {
// auto &wi = rnnlayer->GetWeightsAt(0);
// if (debug)
//    printTensor<Architecture>(wi, "RNN weights");
// }

   if (debug) {
      if (addDenseLayer && rnn.GetLayers().back()->GetWeights().size() > 0)
         printTensor<Architecture>(rnn.GetLayers().back()->GetWeightsAt(0), "weight last layer (DENSE)");

      printTensor<Architecture>(XArch, "input Tensor");
   }

   rnn.Forward(XArch, true);

   if (debug)
      printTensor<Architecture>(rnnlayer->GetOutput(), "output RNN layer after forward");

   rnn.Backward(XArch, Y, weights);



   if (debug) {
      for (size_t l = 0;  l < rnn.GetDepth(); l++) {
         auto &out = rnn.GetLayerAt(l)->GetOutput();
         rnn.GetLayerAt(l)->Print();
         printTensor<Architecture>(out, "output layer");
      }
   }


   std::string maxerrorType;

   ROOT::Math::RichardsonDerivator deriv;

   for (size_t l = rnn.GetDepth() - 1; (int)l >= 0; l--) {
      auto layer = rnn.GetLayerAt(l);
      if (debug) {
         std::cout << "\n\n************************************* \n";
         std::cout << "\tTesting weight gradients:      layer: " << l << " / " << rnn.GetDepth();
         std::cout << std::flush;
         std::cout << std::endl;
         layer->Print();
         std::cout << "************************************* \n\n";
      }

      if (debug) {
         auto &dyi = layer->GetActivationGradients();
         printTensor<Architecture>(dyi, "dy for layer");
      }

      if (layer->GetWeights().size() == 0)
         continue;

      // Weights Input, k = 0
      auto &Wi = layer->GetWeightsAt(0);
      auto &dWi = layer->GetWeightGradientsAt(0);

      if (debug)
         printTensor<Architecture>(dWi, "wgradient");

      Scalar_t maximum_error = 0.0;

      for (size_t i = 0; i < (size_t)Wi.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t)Wi.GetNcols(); ++j) {
            auto f = [&rnn, &XArch, &Y, &weights, l, i, j](Scalar_t x) {
               return evaluate_net_weight(rnn, XArch, Y, weights, l, 0, i, j, x);
            };
            ROOT::Math::Functor1D func(f);
            double dy = deriv.Derivative1(func, Wi(i, j), 1.E-5);
            Scalar_t dy_ref = dWi(i, j);

            // Compute the relative error if dy != 0.
            Scalar_t error;
            std::string errorType;
            if (std::fabs(dy_ref) > 1e-15) {
               error = std::fabs((dy - dy_ref) / dy_ref);
               errorType = "relative";
            } else {
               error = std::fabs(dy - dy_ref);
               errorType = "absolute";
            }

            if (debug)
               std::cout << "Weight-input gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref
                         << std::endl;

            if (error >= maximum_error) {
               maximum_error = error;
               maxerrorType = errorType;
            }
         }
      }

      std::cout << "\rTesting weight input gradients:      ";
      std::cout << "maximum error (" << maxerrorType << "): " << print_error(maximum_error) << std::endl;
      if (maximum_error > 1.E-2) {
         std::cerr << "\033[31m Error \033[39m in weight input gradients" << std::endl;
         failed = true;
      }

      // for Cudnn all weights are collapsed in one
      // if (Architecture::IsCudnn())
      //    continue;

      // if other layers (not RNN) continue
      if (layer->GetWeights().size() == 1)
         continue;


      /// testing weight state gradient

      // Weights State, k = 1
      maximum_error = 0;
      auto &Ws = layer->GetWeightsAt(1);
      auto &dWs = layer->GetWeightGradientsAt(1);
      for (size_t i = 0; i < (size_t)Ws.GetNrows(); ++i) {
         for (size_t j = 0; j < (size_t)Ws.GetNcols(); ++j) {
            auto f = [&rnn, &XArch, &Y, &weights, l, i, j](Scalar_t x) {
               return evaluate_net_weight(rnn, XArch, Y, weights, l, 1, i, j, x);
            };
            ROOT::Math::Functor1D func(f);
            double dy = deriv.Derivative1(func, Ws(i, j), dx);
            Scalar_t dy_ref = dWs(i, j);

            // Compute the relative error if dy != 0.
            Scalar_t error;
            std::string errorType;
            if (std::fabs(dy_ref) > 1e-15) {
               error = std::fabs((dy - dy_ref) / dy_ref);
               errorType = "relative";
            } else {
               error = std::fabs(dy - dy_ref);
               errorType = "absolute";
            }

            if (error >= maximum_error) {
               maximum_error = error;
               maxerrorType = errorType;
            }
            if (debug)
               std::cout << "Weight-state gradient (" << i << "," << j << ") : (num, BP) " << dy << " , " << dy_ref
                         << std::endl;
         }
      }

      std::cout << "\rTesting weight state gradients:      ";
      std::cout << "maximum error (" << maxerrorType << "): " << print_error(maximum_error) << std::endl;
      if (maximum_error > 1.E-2) {
         std::cerr << "\033[31m Error \033[39m in weight state gradients" << std::endl;
         failed = true;
      }

      // testing bias gradients
      maximum_error = 0;
      auto &B = layer->GetBiasesAt(0);
      auto &dB = layer->GetBiasGradientsAt(0);
      for (size_t i = 0; i < (size_t)B.GetNrows(); ++i) {
         auto f = [&rnn, &XArch, &Y, &weights, l, i](Scalar_t x) {
            return evaluate_net_bias(rnn, XArch, Y, weights, l, 0, i, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, B(i, 0), 1.E-5);
         Scalar_t dy_ref = dB(i, 0);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-15) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if (error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
      }
      if (debug) std::cout << "Bias gradient (" << i << ") : (Num, BP) " << dy << " , " << dy_ref << std::endl;
   }

   std::cout << "\rTesting bias gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in bias state gradients" << std::endl;
      failed = true;
   }
   }


   //return std::max(maximum_error, smaximum_error);
   return failed;
}

/*! Generate a DeepNet, test backward pass */
//______________________________________________________________________________
template <typename Architecture>
auto testRecurrentBackpropagationBiases(size_t timeSteps, size_t batchSize, size_t stateSize,
                                        size_t inputSize, typename Architecture::Scalar_t dx)
-> Double_t
{
   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;
   using RNNLayer_t = TBasicRNNLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;
   using Scalar_t = typename Architecture::Scalar_t;

   //std::vector<TMatrixT<Double_t>> XRef(batchSize, TMatrixT<Double_t>(timeSteps, inputSize));    // T x B x D
   //Tensor_t XRef( batchSize, timeSteps, inputSize);
   Tensor_t XArch( batchSize, timeSteps, inputSize);
   //for (size_t i = 0; i < batchSize; ++i) XArch.emplace_back(timeSteps, inputSize); // B x T x D
   randomBatch(XArch);


   Matrix_t Y(batchSize, stateSize), weights(batchSize, 1);
   randomMatrix(Y);
   fillMatrix(weights, 1.0);

   std::cout << "Testing Bias Backprop using RNN with batchsize = " << batchSize << " input = " << inputSize << " state = " << stateSize << " time = " << timeSteps << std::endl;

   Net_t rnn(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   RNNLayer_t* layer = rnn.AddBasicRNNLayer(stateSize, inputSize, timeSteps);
   rnn.AddReshapeLayer(1, timeSteps, stateSize, true);

   rnn.Initialize();
   rnn.Forward(XArch);
   rnn.Backward(XArch, Y, weights);

   Scalar_t maximum_error = 0.0;
   std::string merrorType;

   auto &B = layer->GetBiasesAt(0);
   auto &dB = layer->GetBiasGradientsAt(0);
   for (Int_t i = 0; i < B.GetNrows(); ++i) {
      auto f = [&rnn, &XArch, &Y, &weights, i](Scalar_t x) {
          return evaluate_net_bias(rnn, XArch, Y, weights, 0, 0, i, x);
      };
      Scalar_t dy = finiteDifference(f, dx) / (2.0 * dx);
      Scalar_t dy_ref = dB(i, 0);

      // Compute the relative error if dy != 0.
      Scalar_t error;
      std::string errorType;
      if (std::fabs(dy_ref) > 1e-15) {
         error = std::fabs((dy - dy_ref) / dy_ref);
         errorType = "relative";
      } else {
         error = std::fabs(dy - dy_ref);
         errorType = "absolute";
      }

      if ( error >= maximum_error) {
            maximum_error = error;
            merrorType = errorType;
      }
      std::cout << "Bias gradient (" << i << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      //maximum_error = std::max(error, maximum_error);
   }

   std::cout << "\rTesting bias gradients:      ";
   std::cout << "maximum error (" << merrorType << "): "  << print_error(maximum_error) << std::endl;

   return maximum_error;
}

#endif
