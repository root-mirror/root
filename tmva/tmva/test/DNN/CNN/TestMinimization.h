// @(#)root/tmva/tmva/cnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Deep Learning Minimizer                                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievski      <ilievski.vladimir@live.com>  - CERN, Switzerland  *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef TMVA_TEST_DNN_TEST_CNN_TEST_DL_MINIMIZER_H
#define TMVA_TEST_DNN_TEST_CNN_TEST_DL_MINIMIZER_H

#include <vector>
#include "TMatrix.h"
#include "TMVA/DNN/DLMinimizers.h"
#include "TMVA/DNN/DeepNet.h"
#include "../Utility.h"
#include "TMVA/DNN/TensorDataLoader.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

template <typename Architecture_t>
auto testMinimization(typename Architecture_t::Scalar_t momentum, bool nestorov) -> void
{
   using DeepNet_t = TDeepNet<Architecture_t>;
   using DataLoader_t = TTensorDataLoader<TensorInput, Architecture_t>;

   size_t nSamples = 1000;
   size_t nChannels = 5;
   size_t imgWidth = 8;
   size_t imgHeight = 8;
   size_t batchSize = 64;

   // Initialize train and test input
   std::vector<TMatrixT<Double_t>> XTrain, XTest;
   XTrain.reserve(nSamples);
   XTest.reserve(nSamples);
   for (size_t i = 0; i < nSamples; i++) {
      XTrain.emplace_back(nChannels, imgWidth * imgHeight);
      XTest.emplace_back(nChannels, imgWidth * imgHeight);
   }

   // Initialize train and test output
   size_t nOutput = 2;
   TMatrixT<Double_t> YTrain(nSamples, nOutput), YTest(nSamples, nOutput);

   // Initialize train and test weights
   TMatrixT<Double_t> WTrain(nSamples, 1), WTest(nSamples, 1);

   // Fill-in the input
   for (size_t i = 0; i < nSamples; i++) {
      randomMatrix(XTrain[i]);
      randomMatrix(XTest[i]);
   }

   // Fill-in the output
   randomMatrix(YTrain);
   randomMatrix(YTest);

   // Fill-in the batch weights
   fillMatrix(WTrain, 1.0);
   fillMatrix(WTest, 1.0);

   // Construct the Linear Net
   size_t batchDepth = batchSize;
   size_t batchHeight = nChannels;
   size_t batchWidth = imgHeight * imgWidth;

   DeepNet_t convNet(batchSize, nChannels, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                     ELossFunction::kMeanSquaredError);

   size_t nThreads = 5;
   std::vector<DeepNet_t> nets{};
   nets.reserve(nThreads);

   for (size_t i = 0; i < nThreads; i++) {
      // create a copies of the master deep net
      nets.push_back(convNet);
   }

   // construct the master and slave conv nets
   constructMasterSlaveConvNets(convNet, nets);

   // Initialize the minimizers
   TDLGradientDescent<Architecture_t> minimizer(0.0001, 5, 1);

   // Initialize the tensor inputs
   TensorInput trainingInput(XTrain, YTrain, WTrain);
   TensorInput testInput(XTest, YTest, WTest);
   std::vector<size_t> shape {batchSize, 1, 1, 1, 1};

   DataLoader_t trainingData(trainingInput, nSamples, batchDepth, batchHeight, batchWidth, nOutput, shape,
                             nThreads);
   DataLoader_t testingData(testInput,  nSamples, batchDepth, batchHeight, batchWidth, nOutput, shape,
                             nThreads);

   // Initialize the vector of batches, one batch for one slave network
   std::vector<TTensorBatch<Architecture_t>> batches{};

   size_t batchesInEpoch = nSamples / convNet.GetBatchSize();

   // execute all epochs
   for (size_t i = 0; i < batchesInEpoch; i += nThreads) {
      batches.clear();
      batches.reserve(nThreads);

      for (size_t j = 0; j < nThreads; j++) {
         batches.push_back(trainingData.GetTensorBatch());
      }

      if (momentum != 0.0) {
         if (nestorov) {
            minimizer.StepNesterov(convNet, nets, batches, momentum);
         } else {
            minimizer.StepMomentum(convNet, nets, batches, momentum);
         }
      } else {
         minimizer.Step(convNet, nets, batches);
      }
   }
   // no returning then putting this void
}

#endif
