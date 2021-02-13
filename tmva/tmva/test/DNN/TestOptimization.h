// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Various Optimizers for training DeepNet                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ravi Kiran S      <sravikiran0606@gmail.com>  - CERN, Switzerland         *
 *                                                                                *
 * Copyright (c) 2005-2018:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef TMVA_TEST_DNN_TEST_OPTIMIZATION_H
#define TMVA_TEST_DNN_TEST_OPTIMIZATION_H

#include "Utility.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TFormula.h"
#include "TString.h"

#include "TMVA/Configurable.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/IMethod.h"
#include "TMVA/ClassifierFactory.h"

#include "TMVA/DNN/SGD.h"
#include "TMVA/DNN/Adam.h"
#include "TMVA/DNN/Adagrad.h"
#include "TMVA/DNN/RMSProp.h"
#include "TMVA/DNN/Adadelta.h"
#include "TMVA/DNN/TensorDataLoader.h"

#include <limits>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>

using namespace TMVA::DNN;
using TMVA::DNN::EOptimizer;

/** Train a linear neural network on a randomly generated linear mapping
 *  from an 32-dimensional input space to a 1-dimensional output space.
 *  Returns the error of the response of the network to the input containing
 *  only ones to the 1x32 matrix used to generate the training data.
 */
template <typename Architecture_t>
auto testOptimization(typename Architecture_t::Scalar_t momentum, EOptimizer optimizerType, Bool_t debug) ->
   typename Architecture_t::Scalar_t
{
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Tensor_t = typename Architecture_t::Tensor_t;
   using Scalar_t = typename Architecture_t::Scalar_t; 
   using Layer_t = VGeneralLayer<Architecture_t>;
   using DeepNet_t = TDeepNet<Architecture_t, Layer_t>;
   using DataLoader_t = TTensorDataLoader<TensorInput, Architecture_t>;

   size_t nSamples = 256;
   size_t nFeatures = 32;
   size_t batchSize = 32;

   std::chrono::time_point<std::chrono::system_clock> tInitial, tFinal;
   tInitial = std::chrono::system_clock::now();

   // Initialize train and test input
   // XTrain = (1 x nSamples x nFeatures)
   // XTest  = (1 x nSamples x nFeatures)
   std::vector<TMatrixT<Double_t>> XTrain, XTest;

   XTrain.reserve(1);
   XTest.reserve(1);
   XTrain.emplace_back(nSamples, nFeatures);
   XTest.emplace_back(nSamples, nFeatures);

   // Initialize train and test output
   // YTrain = (nSamples x nOutput)
   // YTest  = (nSamples x nOutput)
   size_t nOutput = 1;
   TMatrixT<Double_t> YTrain(nSamples, nOutput), YTest(nSamples, nOutput);

   // Initialize train and test weights
   // WTrain = (nSamples x 1)
   // WTest  = (nSamples x 1)
   TMatrixT<Double_t> WTrain(nSamples, 1), WTest(nSamples, 1);

   // Initialize K
   // K = (nFeatures x nOutput)
   TMatrixT<Double_t> K(nFeatures, nOutput);

   // Use random K to generate linear mapping
   randomMatrix(K);

   randomMatrix(XTrain[0]);
   randomMatrix(XTest[0]);

   // Generate the output
   // YTrain = XTrain[0] * K
   YTrain.Mult(XTrain[0], K);
   YTest.Mult(XTest[0], K);

   // YTrain.Print();
   // YTest.Print();
   // K.Print(); 
   // XTrain[0].Print();
   // XTest[0].Print();

   // Fill-in the batch weights
   fillMatrix(WTrain, 1.0);
   fillMatrix(WTest, 1.0);

   // Construct the deepNet
   size_t inputDepth = 1;
   size_t inputHeight = 1;
   size_t inputWidth = nFeatures;

   size_t batchDepth = 1;
   size_t batchHeight = batchSize;
   size_t batchWidth = nFeatures;

   DeepNet_t deepNet(batchSize, inputDepth, inputHeight, inputWidth, batchDepth, batchHeight, batchWidth,
                     ELossFunction::kMeanSquaredError, EInitialization::kGauss, ERegularization::kNone, 0.0, true);
   deepNet.AddDenseLayer(32, EActivationFunction::kIdentity);
   deepNet.AddDenseLayer(32, EActivationFunction::kIdentity);
   deepNet.AddDenseLayer(1, EActivationFunction::kIdentity);
   deepNet.Initialize();

   if (debug) {
      deepNet.Print();
   }

   // Initialize the tensor inputs
   size_t nThreads = 1;
   TensorInput trainingInput(XTrain, YTrain, WTrain);
   TensorInput testInput(XTest, YTest, WTest);

   std::vector<size_t> shape {inputDepth, inputHeight, inputWidth };
   

   DataLoader_t trainingData(trainingInput, nSamples, batchSize, shape, {batchDepth, batchHeight, batchWidth}, nOutput,
                             nThreads);
   DataLoader_t testingData(testInput,  nSamples, batchSize, shape, {batchDepth, batchHeight, batchWidth}, nOutput, 
                             nThreads);

   // create a pointer to base class VOptimizer
   std::unique_ptr<VOptimizer<Architecture_t, Layer_t, DeepNet_t>> optimizer;

   // Initialize the optimizer
   switch (optimizerType) {
   case EOptimizer::kSGD:
      optimizer = std::unique_ptr<TSGD<Architecture_t, Layer_t, DeepNet_t>>(
         new TSGD<Architecture_t, Layer_t, DeepNet_t>(0.001, deepNet, momentum));
      break;
   case EOptimizer::kAdam:
      optimizer = std::unique_ptr<TAdam<Architecture_t, Layer_t, DeepNet_t>>(
         new TAdam<Architecture_t, Layer_t, DeepNet_t>(deepNet, 0.001));
      break;
   case EOptimizer::kAdagrad:
      optimizer = std::unique_ptr<TAdagrad<Architecture_t, Layer_t, DeepNet_t>>(
         new TAdagrad<Architecture_t, Layer_t, DeepNet_t>(deepNet, 0.01));
      break;
   case EOptimizer::kRMSProp:
      optimizer = std::unique_ptr<TRMSProp<Architecture_t, Layer_t, DeepNet_t>>(
         new TRMSProp<Architecture_t, Layer_t, DeepNet_t>(deepNet, 0.001, momentum));
      break;
   case EOptimizer::kAdadelta:
      optimizer = std::unique_ptr<TAdadelta<Architecture_t, Layer_t, DeepNet_t>>(
         new TAdadelta<Architecture_t, Layer_t, DeepNet_t>(deepNet, 1.0));
      break;
   }

   // Initialize the variables related to training procedure
   bool converged = false;
   size_t testInterval = 1;
   size_t maxEpochs = 200;
   size_t batchesInEpoch = nSamples / deepNet.GetBatchSize();
   size_t convergenceCount = 0;
   size_t convergenceSteps = 100;

   if (debug) {
      std::string separator(62, '-');
      std::cout << separator << std::endl;
      std::cout << std::setw(10) << "Epoch"
                << " | " << std::setw(12) << "Train Err." << std::setw(12) << "Test Err." << std::setw(12)
                << "t(s)/epoch" << std::setw(12) << "Eval t(s)" << std::setw(12) << "nEvents/s" << std::setw(12)
                << "Conv. Steps" << std::endl;
      std::cout << separator << std::endl;
   }

   // start measuring
   std::chrono::time_point<std::chrono::system_clock> tstart, tend;
   tstart = std::chrono::system_clock::now();

   size_t shuffleSeed = 0;
   TMVA::RandomGenerator<TRandom3> rng(shuffleSeed);

   // test the net
   // Logic : Y = X * K
   // Let X = I, Then Y = I * K => Y = K
   // I = (1 x batchSize x nFeatures)

   Matrix_t I( batchSize, nFeatures);
   Architecture_t::InitializeZero(I);
   for (size_t i = 0; i < batchSize; ++i) { 
        I(i,i) =  1.;
   }     
   Tensor_t tI( I );
  
   // do a forward pass to compute initial Mean Error
   deepNet.Forward(tI, false);

   // get the output of the last layer of the deepNet
   TMatrixT<Scalar_t> Ytemp(deepNet.GetLayerAt(deepNet.GetLayers().size() - 1)->GetOutputAt(0));

   std::cout << " Before Training: Mean Absolute Error = " << meanAbsoluteError(Ytemp, K) << ",";

   Double_t minTestError = 0;
   
   while (!converged) {
      optimizer->IncrementGlobalStep();
      trainingData.Shuffle(rng);

      // training process
      for (size_t i = 0; i < batchesInEpoch; i++) {
         auto my_batch = trainingData.GetTensorBatch();

         //Architecture_t::PrintTensor(my_batch.GetInput(), std::string(TString::Format(" input batch %d",i).Data()));

         deepNet.Forward(my_batch.GetInput(), true);
         deepNet.Backward(my_batch.GetInput(), my_batch.GetOutput(), my_batch.GetWeights());
         optimizer->Step();
      }

      // calculating the error
      if ((optimizer->GetGlobalStep() % testInterval) == 0) {
         std::chrono::time_point<std::chrono::system_clock> t1, t2;
         t1 = std::chrono::system_clock::now();

         // compute test error
         Double_t testError = 0.0;
         //int i = 0; 
         for (auto batch : testingData) {
            auto inputTensor = batch.GetInput();
            auto outputMatrix = batch.GetOutput();
            auto weights = batch.GetWeights();
            //Architecture_t::PrintTensor(inputTensor, std::string(TString::Format(" test batch %d",i++).Data()));
            testError += deepNet.Loss(inputTensor, outputMatrix, weights);
         }
         testError /= (Double_t)(nSamples / batchSize);

         t2 = std::chrono::system_clock::now();

         // checking for convergence
         if (testError < minTestError) {
            convergenceCount = 0;
         } else {
            convergenceCount += testInterval;
         }

         // found the minimum test error
         if (testError < minTestError) {
            if (debug) {
               std::cout << std::setw(10) << optimizer->GetGlobalStep() << " Minimum Test error found : " << testError
                         << std::endl;
            }
            minTestError = testError;
         } else if (minTestError <= 0.0)
            minTestError = testError;

         // compute training error
         Double_t trainingError = 0.0;
         for (auto batch : trainingData) {
            auto inputTensor = batch.GetInput();
            auto outputMatrix = batch.GetOutput();
            auto weights = batch.GetWeights();
            trainingError += deepNet.Loss(inputTensor, outputMatrix, weights);
         }
         trainingError /= (Double_t)(nSamples / batchSize);

         // stop measuring
         tend = std::chrono::system_clock::now();

         // compute numerical throughput
         std::chrono::duration<double> elapsed_seconds = tend - tstart;
         std::chrono::duration<double> elapsed1 = t1 - tstart;

         // time to compute training and test errors
         std::chrono::duration<double> elapsed_testing = tend - t1;

         double seconds = elapsed_seconds.count();
         double eventTime = elapsed1.count() / (batchesInEpoch * testInterval * batchSize);

         converged = optimizer->GetGlobalStep() >= maxEpochs || convergenceCount > convergenceSteps;

         if (debug) {
            std::cout << std::setw(10) << optimizer->GetGlobalStep() << " | " << std::setw(12) << trainingError
                      << std::setw(12) << testError << std::setw(12) << seconds / testInterval << std::setw(12)
                      << elapsed_testing.count() << std::setw(12) << 1. / eventTime << std::setw(12) << convergenceCount
                      << std::endl;

            if (converged) {
               std::cout << std::endl;
            }
         }

         tstart = std::chrono::system_clock::now();
      }

      if (converged && debug) {
         std::cout << "Final Deep Net Weights for epoch " << optimizer->GetGlobalStep() << std::endl;
         auto &weights_tensor = deepNet.GetLayerAt(0)->GetWeights();
         auto &bias_tensor = deepNet.GetLayerAt(0)->GetBiases();
         for (size_t l = 0; l < weights_tensor.size(); l++)
            weights_tensor[l].Print();
         bias_tensor[0].Print();
      }
   }


   deepNet.Forward(tI, false);

   // get the output of the last layer of the deepNet
   TMatrixT<Scalar_t> Y(deepNet.GetLayerAt(deepNet.GetLayers().size() - 1)->GetOutputAt(0));

   if (debug) {
      std::cout << "\nY:\n";

      for (auto i = 0; i < Y.GetNrows(); i++) {
         for (auto j = 0; j < Y.GetNcols(); j++) {
            std::cout << Y(i, j) << " ";
         }
         std::cout << std::endl;
      }

      std::cout << "\nK:\n";
      for (auto i = 0; i < K.GetNrows(); i++) {
         for (auto j = 0; j < K.GetNcols(); j++) {
            std::cout << K(i, j) << " ";
         }
         std::cout << std::endl;
      }
   }

   tFinal = std::chrono::system_clock::now();
   std::chrono::duration<double> totalTime = tFinal - tInitial;

   std::cout << " No of Epochs = " << optimizer->GetGlobalStep() << ", total Time(sec) " << totalTime.count() << ", ";

   return meanAbsoluteError(Y, K);
}

#endif
