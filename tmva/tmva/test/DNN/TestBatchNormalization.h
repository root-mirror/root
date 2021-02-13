// @(#)root/tmva $Id$
// Author: Simon Pfreundschuh

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh, Saurav Shekhar                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Generic tests of the backpropagation algorithm.                //
//                                                                //
// All tests randomly generate a net with identity activation     //
// functions, i.e.  which is completely linear and then tests the //
// computed gradients for each layer using numerical              //
// derivation. The restriction to linear nets is to avoid the     //
// required division by the finite difference interval used to    //
// approximate the numerical derivatives, which would otherwise   //
// cause precision loss.                                          //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"
#include "TMVA/DNN/BatchNormLayer.h"
#include "Utility.h"
#include "TMath.h"
#include "Math/RichardsonDerivator.h"
#include "Math/Functor.h"


using namespace TMVA::DNN;

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
    Scalar_t res = net.Loss(X, Y, W,  true, false);  // set inTraining=true when computing Loss
    net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j) = prev_value;
    //std::cout << "compute loss for weight  " << xvalue << "  " << prev_value << " result " << res << std::endl;
    return res;
}

template <typename Architecture>
double check_gradients(TDeepNet<Architecture> &net, typename Architecture::Tensor_t & X,
typename Architecture::Matrix_t &Y, typename Architecture::Matrix_t & weights, double dx_eps) {

   using Scalar_t = typename Architecture::Scalar_t;
   double maximum_error = 0.0;

   // Compute derivatives for all weights using finite differences and
   // compare to result obtained from backpropagation.
   ROOT::Math::RichardsonDerivator deriv;
   for (size_t il = 0; il < net.GetLayers().size(); il++) {
      int l = net.GetLayers().size() - il - 1;
      // if (l < 1) continue;
      auto layer = net.GetLayerAt(l);
      for (size_t k = 0; k < layer->GetWeights().size(); k++) {
         // if (k != 1 ) continue;
         std::cout << "\rTesting weight gradients   for    layer " << l << std::endl;
         std::cout << std::flush;
         auto &dW = layer->GetWeightGradientsAt(k);
         std::cout << "weight gradient for layer " << l << " component " << k << std::endl;
         dW.Print();
         auto &W = layer->GetWeightsAt(k);
         std::cout << "weights for layer " << l << std::endl;
         W.Print();

         size_t nw1 = W.GetNrows();
         size_t nw2 = W.GetNcols();
         for (size_t i = 0; i < nw1; ++i) {
            for (size_t j = 0; j < nw2; j++) {
               auto f = [&net, &X, &Y, &weights, l, k, i, j](Scalar_t x) {
                  return evaluate_net_weight(net, X, Y, weights, l, k, i, j, x);
               };
               ROOT::Math::Functor1D func(f);
               double dy = deriv.Derivative1(func, W(i, j), dx_eps);

               // Scalar_t dy = finiteDifference(f, dx) / (2.0 * dx);
               Scalar_t dy_ref = dW(i, j);

               std::cout << "dW( " << i << " , " << j << " )  numeric = " << dy << " from BP = " << dy_ref << std::endl;
               // Compute the relative error if dy != 0.
               Scalar_t error;
               if (std::fabs(dy_ref) > 1e-10) {
                  error = std::fabs((dy - dy_ref) / dy_ref);
               } else {
                  error = std::fabs(dy - dy_ref);
               }
               maximum_error = std::max(error, maximum_error);
            }
         }
      }
      std::cout << "max error = " << maximum_error << std::endl;
   }

   std::cout << "\rTesting weight gradients:      ";
   std::cout << "maximum relative error: " << print_error(maximum_error) << std::endl;
   return maximum_error;
}


/*! Generate a random net, perform forward and backward propagation and check
 *  the weight gradients using numerical differentiation. Returns the maximum
 *  relative gradient error and also prints it to stdout. */
//______________________________________________________________________________
template <typename Architecture>
auto testBackpropagationWeights(typename Architecture::Scalar_t dx_eps)
-> typename Architecture::Scalar_t
{
   // TODO pass as function params
   size_t tbatchSize = 3, inputSize = 4, outputSize = 2;

   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using Net_t = TDeepNet<Architecture>;
   // using FCLayer_t  = TDenseLayer<Architecture>;

   // Random net.
   Net_t net(tbatchSize, 1, tbatchSize, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
             EInitialization::kGauss);
   // FCLayer_t* l1 = net.AddDenseLayer(outputSize, EActivationFunction::kIdentity);
   net.AddDenseLayer(outputSize, EActivationFunction::kIdentity);

   net.AddBatchNormLayer();

   //    auto &layers = net.GetLayers();
   // auto bnlayer = new TBatchNormLayer<Architecture>(tbatchSize, 1, 1, outputSize);
   // layers.push_back( bnlayer);

   net.AddDenseLayer(1, EActivationFunction::kIdentity);



   net.Print();

   // Random training data.
   Tensor_t X(1, tbatchSize, inputSize); // T x B x D
   Matrix_t Y(tbatchSize, 1), weights(tbatchSize, 1);
   net.Initialize();
   randomBatch(X);
   // Matrix_t & input = X[0];
   // for (int i = 0; i < tbatchSize; ++i) {
   //    for (int j = 0; j < inputSize; ++j) {
   //       input(i,j) = i*2 + j;
   //    }
   // }
   //input.Print();


   fillMatrix(Y,0.0);
   fillMatrix(weights, 1.0);

   std::cout << "input \n";
   Architecture::PrintTensor(X, "input network X ");
   //X[0].Print();

   std::cout << "weights layer 0" << std::endl;
   using Scalar_t = typename Architecture::Scalar_t;
   TMatrixT<Scalar_t> mw = net.GetLayerAt(0)->GetWeightsAt(0);
   mw.Print();

   // auto & w2 = net.GetLayerAt(2)->GetWeightsAt(0);
   net.GetLayerAt(0)->GetOutput().PrintShape();
   net.GetLayerAt(1)->GetOutput().PrintShape();
   net.GetLayerAt(2)->GetOutput().PrintShape();

   net.Forward(X, true);

   std::cout << "output of layer before BNorm layer \n";
   Tensor_t &inputBN = net.GetLayerAt(0)->GetOutput();
   Architecture::PrintTensor(inputBN, "input BN layer");

   std::cout << "output BN \n";

   std::vector<double> data(tbatchSize);
   std::cout << "Mean and std of BN inputs\n";
   for (size_t k = 0; k < outputSize; ++k) {
      for (size_t i = 0; i < tbatchSize; ++i) data[i] = net.GetLayerAt(0)->GetOutputAt(0)(i,k);
      std::cout << "output DL feature " << k << " mean " <<  TMath::Mean(data.begin(), data.end() ) << "\t";
      std::cout << "output DL std " <<  TMath::RMS(data.begin(), data.end() ) << std::endl;
   }
   std::cout << "output of BN \n";
   Architecture::PrintTensor(net.GetLayerAt(1)->GetOutput(), "output BN layer");

   std::cout << "Mean and std of BN outputs\n";
   for (size_t k = 0; k < outputSize; ++k) {
      for (size_t i = 0; i < tbatchSize; ++i) data[i] = net.GetLayerAt(1)->GetOutputAt(0)(i,k);
      std::cout << "output BN feature " << k << " mean " <<  TMath::Mean(data.begin(), data.end() ) << "\t";
      std::cout << "output BN std " <<  TMath::RMS(data.begin(), data.end() ) << std::endl;
   }

   std::cout << "\nDo backward pass.....\n";

   net.Backward(X, Y, weights);

   // print dx of bnorm layer (activation gradients of previous one)
   Architecture::PrintTensor(net.GetLayerAt(0)->GetActivationGradients(), "Computed BN dx");

   std::cout << "\nTest gradients.....\n";
   double error = check_gradients(net, X, Y, weights, dx_eps);
   return error;
}

/*! Generate a random net, perform forward and backward propagation and check
 *  the weight gradients using numerical differentiation. Returns the maximum
 *  relative gradient error and also prints it to stdout. */
//______________________________________________________________________________
template <typename Architecture>
auto testCNNBackpropagationWeights(typename Architecture::Scalar_t dx_eps) -> typename Architecture::Scalar_t
{
   // TODO pass as function params
   size_t tbatchSize = 10, inputD = 2, inputH = 4, inputW = 4 ;

   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using Net_t = TDeepNet<Architecture>;
   // using FCLayer_t  = TDenseLayer<Architecture>;

   // Random net.
   Net_t net(tbatchSize, inputD, inputH, inputW, tbatchSize, inputD, inputH * inputW, ELossFunction::kMeanSquaredError,
             EInitialization::kGauss);
   // FCLayer_t* l1 = net.AddDenseLayer(outputSize, EActivationFunction::kIdentity);
   net.AddConvLayer(3, 3, 3, 1, 1, 1, 1, EActivationFunction::kRelu);

   //net.AddConvLayer(3, 3, 3, 1, 1, 1, 1, EActivationFunction::kRelu);

   net.AddBatchNormLayer();

   //net.AddMaxPoolLayer(2, 2, 1, 1);

   net.AddReshapeLayer(3, 4, 4, true);

   net.AddDenseLayer(4, EActivationFunction::kSigmoid);

   net.AddBatchNormLayer();

   //    auto &layers = net.GetLayers();
   // auto bnlayer = new TBatchNormLayer<Architecture>(tbatchSize, 1, 1, outputSize);
   // layers.push_back( bnlayer);

   net.AddDenseLayer(1, EActivationFunction::kIdentity);

   net.Print();

   // Random training data.
   Tensor_t X = Architecture::CreateTensor(tbatchSize, inputD, inputH , inputW); // T x B x D
   Matrix_t Y(tbatchSize, 1), weights(tbatchSize, 1);
   net.Initialize();
   randomBatch(X);

   fillMatrix(Y, 0.0);
   fillMatrix(weights, 1.0);

   net.Forward(X, true);

   net.Backward(X, Y, weights);

   double error = check_gradients(net, X, Y, weights, dx_eps);
   return error;
}
