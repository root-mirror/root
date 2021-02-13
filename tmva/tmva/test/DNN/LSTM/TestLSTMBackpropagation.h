// @(#)root/tmva $Id$
// Author: Surya S Dwivedi

/*************************************************************************
 * Copyright (C) 2019, Surya S Dwivedi                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Generic tests of the LSTMLayer Backward pass                    //
////////////////////////////////////////////////////////////////////

#ifndef TMVA_TEST_DNN_TEST_LSTM_TEST_BWDPASS_H
#define TMVA_TEST_DNN_TEST_LSTM_TEST_BWDPASS_H

#include <iostream>
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

template <typename Architecture>
auto printMatrix(const typename Architecture::Matrix_t &A, const std::string name = "matrix")
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
    //std::cout << "compute loss for weight  " << xvalue << "  " << prev_value << " result " << res << std::endl;
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
bool testLSTMBackpropagation(size_t timeSteps, size_t batchSize, size_t stateSize,
                                  size_t inputSize, typename Architecture::Scalar_t dx = 1.E-5,
                                  std::vector<bool> options = {}, bool debug = false)

{
   bool failed = false;
   const int nOpts = 4; // size of options
   if (options.size() < nOpts)
      options.resize(nOpts);
   bool randomInput = !options[0];
   bool addDenseLayer = options[1];
   bool addExtraLSTM = options[2];
   bool returnLastSequence = options[3];

   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;
   using LSTMLayer_t = TBasicLSTMLayer<Architecture>;
   using DenseLayer_t = TDenseLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;
   using Scalar_t = typename Architecture::Scalar_t;

   if (debug)
      std::cout << std::endl;
   std::cout
      << "******************************************************************************************************\n";
   std::cout << "Testing Weight Backprop using LSTM with batchsize = " << batchSize << " input = " << inputSize
             << " state = " << stateSize << " time = " << timeSteps;
   if (randomInput)
      std::cout << " using a random input";
   else
      std::cout << " with a fixed input";
   if (addDenseLayer)
      std::cout << " and a dense layer";
   if (addExtraLSTM)
      std::cout << " and an extra LSTM";
   if (returnLastSequence)
      std::cout << " and full output";
   std::cout << std::endl;
   std::cout
      << "******************************************************************************************************\n";
   if (debug)
      std::cout << std::endl;

   Tensor_t XArch = Architecture::CreateTensor(batchSize,timeSteps, inputSize);

   // for random input (default)
   if (randomInput) {
   for (size_t i = 0; i < batchSize; ++i) {
         auto  mat = XArch[i];
         for (size_t l = 0; l < (size_t) XArch[i].GetNrows(); ++l) {
            for (size_t m = 0; m < (size_t) XArch[i].GetNcols(); ++m) {
               mat(l, m) = gRandom->Uniform(-1,1);
               //XArch[i](0, 0) = 0.5;
               //XArch[i](1, 0) = 0.5;
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
         auto mat = XArch[i];
         // time 0
         for (size_t l = 0; l < timeSteps; ++l) {
            for (size_t m = 0; m < inputSize; ++m) {
               mat(l,m) = Input(l,m);
            }
         }
      }
      gRandom->SetSeed(1); // for weights initizialization
   }
   if (debug) printTensor<Architecture>(XArch,"input");

   size_t outputSize = (returnLastSequence) ? timeSteps * stateSize : stateSize;
   if (addDenseLayer) outputSize = 1;

   Matrix_t Y(batchSize, outputSize), weights(batchSize, 1);
   //randomMatrix(Y);
   for (size_t i = 0; i < (size_t) Y.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Y.GetNcols(); ++j) {
         Y(i, j) = gRandom->Integer(2);
      }
   }
    if (debug) printTensor<Architecture>(Y,"ground truth ");
   fillMatrix(weights, 1.0);

   bool returnFirstSequence = addExtraLSTM || returnLastSequence;
   Net_t lstm(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   LSTMLayer_t* layer = lstm.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, returnFirstSequence );

   // input size second RNN is = stateSize
   if (addExtraLSTM) lstm.AddBasicLSTMLayer(stateSize, stateSize, timeSteps, false, returnLastSequence);

   size_t outputLSTMSize = (returnLastSequence) ? timeSteps * stateSize : stateSize;
   lstm.AddReshapeLayer(1, 1, outputLSTMSize, true);

   DenseLayer_t * dlayer1 = nullptr;
   DenseLayer_t * dlayer2 = nullptr;
   if (addDenseLayer) {
      dlayer1 = lstm.AddDenseLayer(10, TMVA::DNN::EActivationFunction::kTanh);
      dlayer2 = lstm.AddDenseLayer(1, TMVA::DNN::EActivationFunction::kIdentity);
   }


   lstm.Initialize();

   if (debug)
      lstm.Print();

   if (debug)
      Architecture::PrintTensor(XArch,"input");

   auto lstmLayer = lstm.GetLayerAt(0);

   if (debug) {
      for (size_t i = 0; i < lstmLayer->GetWeights().size(); ++i)
         Architecture::PrintTensor(Tensor_t(lstmLayer->GetWeightsAt(i)),
                                   std::string(TString::Format("LSTM weights - %d",(int) i).Data()));
   }

   lstm.Forward(XArch, true);

   if (debug)
      Architecture::PrintTensor(lstmLayer->GetOutput(), "output forward");

   lstm.Backward(XArch, Y, weights);

   std::cout << "after backward ! " << std::endl;

   if (debug) {
      Architecture::PrintTensor(lstmLayer->GetActivationGradients(), "dy");
      // Architecture::PrintTensor(lstmLayer->GetActivationGradients,"dx");

      for (size_t i = 0; i < lstmLayer->GetWeights().size(); ++i)
         Architecture::PrintTensor(Tensor_t(lstmLayer->GetWeightGradientsAt(i)),
                                   std::string(TString::Format("LSTM weight Gradients - %d", (int) i).Data()));
   }

   if (debug)  {
      auto & out = layer->GetOutput();
      printTensor<Architecture>(out,"output");
      if (dlayer1) {
         auto & out2 = dlayer1->GetOutput();
         printTensor<Architecture>(out2,"dense layer1 output");
         auto & out3 = dlayer2->GetOutput();
         printTensor<Architecture>(out3,"dense layer2 output");
      }
   }


   Scalar_t maximum_error = 0.0;
   std::string maxerrorType;

   ROOT::Math::RichardsonDerivator deriv;


   // Weights Input Gate, k = 0

   auto & Wi = layer->GetWeightsAt(0);
   auto & dWi = layer->GetWeightGradientsAt(0);
   for (size_t i = 0; i < (size_t) Wi.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wi.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 0, i, j, x);
         };


         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wi(i,j), 1.E-3);
         Scalar_t dy_ref = dWi(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if (debug) std::cout << "Weights input gate gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;

         if (error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
      }
   }

   std::cout << "\rTesting weights input gate gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights input gate gradients" << std::endl;
      failed = true;
   }

   // Weights Input Gate State, k = 1
   maximum_error = 0;


   auto & Wis = layer->GetWeightsAt(1);
   auto & dWis = layer->GetWeightGradientsAt(1);
   for (size_t i = 0; i < (size_t) Wis.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wis.GetNcols(); ++j) {

         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 1, i, j, x);
         };

         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wis(i,j), dx);
         Scalar_t dy_ref = dWis(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights input gate-state gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights input gate-state gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights input gate-state gradients" << std::endl;
      failed = true;
   }

    // Weights Forget Gate State, k = 2
   maximum_error = 0;
   auto & Wf = layer->GetWeightsAt(2);
   auto & dWf = layer->GetWeightGradientsAt(2);
   for (size_t i = 0; i < (size_t) Wf.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wf.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 2, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wf(i,j), dx);
         Scalar_t dy_ref = dWf(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights forget gate gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights forget gate gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights forget gate gradients" << std::endl;
      failed = true;
   }

    // Weights Forget Gate State, k = 3
   maximum_error = 0;
   auto & Wfs = layer->GetWeightsAt(3);
   auto & dWfs = layer->GetWeightGradientsAt(3);
   for (size_t i = 0; i < (size_t) Wfs.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wfs.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 3, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wfs(i,j), dx);
         Scalar_t dy_ref = dWfs(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights forget gate-state gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights forget gate-state gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights forget gate-state gradients" << std::endl;
      failed = true;
   }


   // Weights Candidate Gate k = 4
   maximum_error = 0;
   auto & Wc = layer->GetWeightsAt(4);
   auto & dWc = layer->GetWeightGradientsAt(4);
   for (size_t i = 0; i < (size_t) Wc.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wc.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 4, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wc(i,j), dx);
         Scalar_t dy_ref = dWc(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights candidate gate gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights candidate gate gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights candidate gate gradients" << std::endl;
      failed = true;
   }


   // Weights Candidate Gate State, k = 5
   maximum_error = 0;
   auto & Wcs = layer->GetWeightsAt(5);
   auto & dWcs = layer->GetWeightGradientsAt(5);
   for (size_t i = 0; i < (size_t) Wcs.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wcs.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 5, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wcs(i,j), dx);
         Scalar_t dy_ref = dWcs(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights candidate gate-state gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights candidate gate-state gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights candidate gate-state gradients" << std::endl;
      failed = true;
   }


   // Weights Output Gate, k = 6
   maximum_error = 0;
   auto & Wo = layer->GetWeightsAt(6);
   auto & dWo = layer->GetWeightGradientsAt(6);
   for (size_t i = 0; i < (size_t) Wo.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wo.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 6, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wo(i,j), dx);
         Scalar_t dy_ref = dWo(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights output gate gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights output gate gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights output gate gradients" << std::endl;
      failed = true;
   }

   // Weights Output Gate State, k = 7
   maximum_error = 0;
   auto & Wos = layer->GetWeightsAt(7);
   auto & dWos = layer->GetWeightGradientsAt(7);
   for (size_t i = 0; i < (size_t) Wos.GetNrows(); ++i) {
      for (size_t j = 0; j < (size_t) Wos.GetNcols(); ++j) {
         auto f = [&lstm, &XArch, &Y, &weights, i, j](Scalar_t x) {
             return evaluate_net_weight(lstm, XArch, Y, weights, 0, 7, i, j, x);
         };
         ROOT::Math::Functor1D func(f);
         double dy = deriv.Derivative1(func, Wos(i,j), dx);
         Scalar_t dy_ref = dWos(i, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         std::string errorType;
         if (std::fabs(dy_ref) > 1e-7) {
            error = std::fabs((dy - dy_ref) / dy_ref);
            errorType = "relative";
         } else {
            error = std::fabs(dy - dy_ref);
            errorType = "absolute";
         }

         if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
         }
         if (debug) std::cout << "Weights output gate-state gradient (" << i << "," << j << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
      }
   }

   std::cout << "\rTesting weights output gate-state gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in weights output gate-state gradients" << std::endl;
      failed = true;
   }


   // testing input gate bias gradients
   maximum_error = 0;
   auto &Bi = layer->GetBiasesAt(0);
   auto &dBi = layer->GetBiasGradientsAt(0);
   for (size_t i = 0;  i < (size_t) Bi.GetNrows(); ++i) {
      auto f = [&lstm, &XArch, &Y, &weights, i](Scalar_t x) {
          return evaluate_net_bias(lstm, XArch, Y, weights, 0, 0, i, x);
      };
      ROOT::Math::Functor1D func(f);
      double dy = deriv.Derivative1(func, Bi(i,0), 1.E-5);
      Scalar_t dy_ref = dBi(i, 0);

      // Compute the relative error if dy != 0.
      Scalar_t error;
      std::string errorType;
      if (std::fabs(dy_ref) > 1e-7) {
         error = std::fabs((dy - dy_ref) / dy_ref);
         errorType = "relative";
      } else {
         error = std::fabs(dy - dy_ref);
         errorType = "absolute";
      }

      if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
      }
      if (debug) std::cout << "input gate bias gradient (" << i << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
   }

   std::cout << "\rTesting input gate bias gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in input gate bias gradients" << std::endl;
      failed = true;
   }

   // testing foget gate bias gradients
   maximum_error = 0;
   auto & Bf = layer->GetBiasesAt(1);
   auto & dBf = layer->GetBiasGradientsAt(1);
   for (size_t i = 0;  i < (size_t) Bf.GetNrows(); ++i) {
      auto f = [&lstm, &XArch, &Y, &weights, i](Scalar_t x) {
          return evaluate_net_bias(lstm, XArch, Y, weights, 0, 1, i, x);
      };
      ROOT::Math::Functor1D func(f);
      double dy = deriv.Derivative1(func, Bf(i,0), 1.E-5);
      Scalar_t dy_ref = dBf(i, 0);

      // Compute the relative error if dy != 0.
      Scalar_t error;
      std::string errorType;
      if (std::fabs(dy_ref) > 1e-7) {
         error = std::fabs((dy - dy_ref) / dy_ref);
         errorType = "relative";
      } else {
         error = std::fabs(dy - dy_ref);
         errorType = "absolute";
      }

      if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
      }
      if (debug) std::cout << "Forget gate bias gradient (" << i << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
   }

   std::cout << "\rTesting forget gate bias gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in forget gate bias gradients" << std::endl;
      failed = true;
   }

   // testing candidate gate bias gradients
   maximum_error = 0;
   auto & Bc = layer->GetBiasesAt(2);
   auto & dBc = layer->GetBiasGradientsAt(2);
   for (size_t i = 0;  i < (size_t) Bc.GetNrows(); ++i) {
      auto f = [&lstm, &XArch, &Y, &weights, i](Scalar_t x) {
          return evaluate_net_bias(lstm, XArch, Y, weights, 0, 2, i, x);
      };
      ROOT::Math::Functor1D func(f);
      double dy = deriv.Derivative1(func, Bc(i,0), 1.E-5);
      Scalar_t dy_ref = dBc(i, 0);

      // Compute the relative error if dy != 0.
      Scalar_t error;
      std::string errorType;
      if (std::fabs(dy_ref) > 1e-7) {
         error = std::fabs((dy - dy_ref) / dy_ref);
         errorType = "relative";
      } else {
         error = std::fabs(dy - dy_ref);
         errorType = "absolute";
      }

      if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
      }
      if (debug) std::cout << "Candidate gate bias gradient (" << i << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
   }

   std::cout << "\rTesting candidate gate bias gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in candidate gate bias gradients" << std::endl;
      failed = true;
   }

   // testing output gate bias gradients
   maximum_error = 0;
   auto & Bo = layer->GetBiasesAt(3);
   auto & dBo = layer->GetBiasGradientsAt(3);
   for (size_t i = 0;  i < (size_t) Bo.GetNrows(); ++i) {
      auto f = [&lstm, &XArch, &Y, &weights, i](Scalar_t x) {
          return evaluate_net_bias(lstm, XArch, Y, weights, 0, 3, i, x);
      };
      ROOT::Math::Functor1D func(f);
      double dy = deriv.Derivative1(func, Bo(i,0), 1.E-5);
      Scalar_t dy_ref = dBo(i, 0);

      // Compute the relative error if dy != 0.
      Scalar_t error;
      std::string errorType;
      if (std::fabs(dy_ref) > 1e-7) {
         error = std::fabs((dy - dy_ref) / dy_ref);
         errorType = "relative";
      } else {
         error = std::fabs(dy - dy_ref);
         errorType = "absolute";
      }

      if ( error >= maximum_error) {
            maximum_error = error;
            maxerrorType = errorType;
      }
      if (debug) std::cout << "Output gate bias gradient (" << i << ") : (comp, ref) " << dy << " , " << dy_ref << std::endl;
   }

   std::cout << "\rTesting output gate bias gradients:      ";
   std::cout << "maximum error (" << maxerrorType << "): "  << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-2) {
      std::cerr << "\033[31m Error \033[39m in output gate bias gradients" << std::endl;
      failed = true;
   }


   //return std::max(maximum_error, smaximum_error);
   return failed;
}

#endif
