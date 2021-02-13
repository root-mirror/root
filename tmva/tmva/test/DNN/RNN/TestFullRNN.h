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

#ifndef TMVA_TEST_DNN_TEST_RNN_TEST_FULL
#define TMVA_TEST_DNN_TEST_RNN_TEST_FULL

#include <iostream>
#include <vector>
#include <string>

#include "../Utility.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"
#include "TMVA/DNN/Net.h"
#include "TMVA/DNN/Architectures/Cpu.h"
#include "TMVA/ROCCurve.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TMath.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

template <typename Architecture>
auto printTensor1(const std::vector<typename Architecture::Matrix_t> &A, const std::string name = "matrix")
-> void
{
  std::cout << name << "\n";
  for (size_t l = 0; l < A.size(); ++l) {
      for (size_t i = 0; i < A[l].GetNrows(); ++i) {
        for (size_t j = 0; j < A[l].GetNcols(); ++j) {
            std::cout << A[l](i, j) << " ";
        }
        std::cout << "\n";
      }
      std::cout << "********\n";
  }
}

template <typename Architecture>
auto printMatrix1(const typename Architecture::Matrix_t &A, const std::string name = "matrix")
-> void
{
  std::cout << name << "\n";
  for (size_t i = 0; i < A.GetNrows(); ++i) {
    for (size_t j = 0; j < A.GetNcols(); ++j) {
        std::cout << A(i, j) << " ";
    }
    std::cout << "\n";
  }
  std::cout << "********\n";
}

/* Generate a full recurrent neural net
 * like a word generative model */
//______________________________________________________________________________
template <typename Architecture>
auto testFullRNN(TString rnnType, size_t batchSize, size_t stateSize,
                 size_t inputSize, size_t outputSize, bool debug = false)
-> bool
{
   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;

   // using RNNLayer_t = TBasicRNNLayer<Architecture>;
   // using FCLayer_t  = TDenseLayer<Architecture>;
   // using Reshape_t  = TReshapeLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;
   using Scalar_t   = typename Architecture::Scalar_t;

   // check, denselayer takes only first one as input,
   // so make sure time = 1, in the current case
   size_t timeSteps = 1;

   // fix random seed
   TRandom3 rndm(12345);
   UInt_t arch_seed = rndm.Integer(TMath::Limits<UInt_t>::Max());
   Architecture::SetRandomSeed(arch_seed);

   Tensor_t XArch = Architecture::CreateTensor(batchSize, timeSteps, inputSize); // B x T x D

   randomBatch(XArch);
   Tensor_t XRef = XArch;

     // B x T x D
   //TMatrixT<Double_t> YRef(batchSize, outputSize);    // B x O  (D = O)

   Matrix_t YArch(batchSize, outputSize);             // B x O  (D = O)

   std::cerr << "Copying output into input\n";
   for (size_t i = 0; i < batchSize; ++i) {
      auto tmp = XArch.At(i);
      for (size_t j = 0; j < outputSize; ++j) {
         YArch(i, j) = tmp(0, j); // time steps is 1
      }
   }

   Net_t rnn(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError, EInitialization::kGauss);

   if (rnnType == "RNN")
      rnn.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, false);  // single output , but not important
   else if (rnnType == "LSTM")
      rnn.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, false);
   else if (rnnType == "GRU")
      rnn.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, false);

   rnn.AddReshapeLayer(1, 1, stateSize, true);
   rnn.AddDenseLayer(outputSize, EActivationFunction::kIdentity);

   Matrix_t W(batchSize, 1);
   for (size_t i = 0; i < batchSize; ++i) W(i, 0) = 1.0;
   rnn.Initialize();

   size_t iter = 0;
   Scalar_t loss = 999;
   while (iter++ < 50) {
      rnn.Forward(XArch);
      loss = rnn.Loss(YArch, W, false);

      //if (iter % 20 == 0) {
         //for (size_t i = 0; i < inputSize; ++i) std::cout << XRef[0](0, i) << " "; std::cout << "\n";
         //for (size_t i = 0; i < inputSize; ++i) std::cout << rnn.GetLayers().back()->GetOutputAt(0)(0, i) << " "; std::cout << "\n";
      //}
      if (debug || iter%10 == 0)
         std::cout << "iteration : " << iter << "  loss: " << loss << std::endl;

      rnn.Backward(XArch, YArch, W);

      rnn.Update(0.1);
   }
   // at the end loos should be < 0.05
   if (loss > 0.05 ) {
      Error("testFullRNN", "%s simple training test failed",rnnType.Data());
      return false;
   }
   Info("testFullRNN", "%s simple training test passed",rnnType.Data());
   return true;
}

/* Generate a full recurrent neural net
   with several time steps and using a dense layer afterwards
   The time steps is fixed to 5
*/
//______________________________________________________________________________
template <typename Architecture>
auto testFullRNN2(TString rnnType, size_t batchSize, size_t stateSize,
                  size_t inputSize, int seed, bool debug = false)
-> bool
{
   using Matrix_t   = typename Architecture::Matrix_t;
   using Tensor_t   = typename Architecture::Tensor_t;
   // using RNNLayer_t = TBasicRNNLayer<Architecture>;
   // using FCLayer_t  = TDenseLayer<Architecture>;
   // using Reshape_t  = TReshapeLayer<Architecture>;
   using Net_t      = TDeepNet<Architecture>;
   using Scalar_t   = typename Architecture::Scalar_t;

   bool saveResult = debug;
   bool useCPUWeights = false;
   if (debug)  useCPUWeights = true;

   TRandom3 rndm(seed);
   UInt_t arch_seed = rndm.Integer(TMath::Limits<UInt_t>::Max());
   Architecture::SetRandomSeed(arch_seed);

   // check, denselayer takes only first one as input,
   // so make sure time = 1, in the current case
   size_t timeSteps = 5;

   Tensor_t XRef = Architecture::CreateTensor(batchSize, timeSteps, inputSize); // B x T x D

  ///std::vector<TMatrixT<Double_t>> XRef(batchSize, TMatrixT<Double_t>(timeSteps, inputSize));    // B x T x D
   //TMatrixT<Double_t> YRef(batchSize, outputSize);    // B x O  (D = O)
   //Tensor_t XArch;
   Matrix_t YArch(batchSize, 1);             // B x O  (D = O)
   for (size_t i = 0; i < batchSize; ++i) {
      // provide input data and labels Yarch
      // define if events is class  0 or 1
      int label = rndm.Integer(2);
      auto xtmp = XRef.At(i);
      YArch(i, 0) = label;
      for (size_t l = 0; l < timeSteps; ++l) {
         double mu = (label == 0) ? 4 : 2*l;
         for (size_t m = 0; m < inputSize; ++m) {
            mu += m;   // shift the varouous inputs
            xtmp(l,m) = rndm.Gaus( mu, 1);
         }
      }
      //std::cerr << "Copying output into input\n";
   }
   Tensor_t XArch = XRef;

   bool useRegularization = false;
   double weightDecay = (useRegularization) ? 1. : 0;


   Net_t rnn(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kCrossEntropy, EInitialization::kGauss,
      (useRegularization) ? ERegularization::kL2 : ERegularization::kNone, weightDecay);

   if (rnnType == "RNN")
      rnn.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, true); // output full sequence
   else if (rnnType == "LSTM")
         rnn.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, true);
   else if (rnnType == "GRU")
      rnn.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, true);

   rnn.AddReshapeLayer(1, 1, timeSteps*stateSize, true);
   rnn.AddDenseLayer(10, EActivationFunction::kTanh);
   rnn.AddDenseLayer(1, EActivationFunction::kIdentity);

   Matrix_t W(batchSize, 1);
   for (size_t i = 0; i < batchSize; ++i) W(i, 0) = 1.0;
   if (useCPUWeights && Architecture::IsCudnn() ) {
      std::cout << "Use Cuddn architecture - create a CPU network to get same weights as CPU run" << std::endl;
      // create a CPU network and initialize and copy weights to have same CPU-GPU network for debugging
      TDeepNet<TCpu<Scalar_t>> rnn2(batchSize, batchSize, timeSteps, inputSize, 0, 0, 0, ELossFunction::kCrossEntropy,
                                 EInitialization::kGauss,
                                 (useRegularization) ? ERegularization::kL2 : ERegularization::kNone, weightDecay);



      if (rnnType == "RNN")
         rnn2.AddBasicRNNLayer(stateSize, inputSize, timeSteps, false, true); // output full sequence
      else if (rnnType == "LSTM")
         rnn2.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, false, true);
      else if (rnnType == "GRU")
         rnn2.AddBasicGRULayer(stateSize, inputSize, timeSteps, false, true);

      rnn2.AddReshapeLayer(1, 1, timeSteps * stateSize, true);
      rnn2.AddDenseLayer(10, EActivationFunction::kTanh);
      rnn2.AddDenseLayer(1, EActivationFunction::kIdentity);
      TCpu<Scalar_t>::SetRandomSeed(arch_seed);
      rnn2.Initialize();
      // print weights
      TCpu<Scalar_t>::PrintTensor(rnn2.GetLayers().front()->GetWeightsAt(0), "RNN2 weight inputs ");

      rnn.Initialize();
      for (size_t i = 0; i < rnn.GetLayers().size(); ++i)
         rnn.GetLayerAt(i)->CopyParameters(*rnn2.GetLayerAt(i));
   } else
      rnn.Initialize();



   if (debug) {
      Architecture::PrintTensor(XArch, " input tensor");
      Architecture::PrintTensor(YArch, " Label tensor");
   }

   size_t iter = 0;
   while (iter++ < 50) {
      rnn.Forward(XArch, true);


      if (debug) {
         Architecture::PrintTensor(rnn.GetLayers().front()->GetOutput(), "Recurrent output");
         Architecture::PrintTensor(rnn.GetLayers().back()->GetOutput(), "Network output");
         // print weights
         Architecture::PrintTensor(rnn.GetLayers().front()->GetWeightsAt(0), "RNN weight inputs ");
         Architecture::PrintTensor(rnn.GetLayers().front()->GetWeightsAt(1), "RNN weight state ");
         Architecture::PrintTensor(rnn.GetLayers()[2]->GetWeightsAt(0), "First Dense weights ");
         Architecture::PrintTensor(rnn.GetLayers().back()->GetWeightsAt(0), "Last Dense weights ");
      }

      Scalar_t loss = rnn.Loss(YArch, W, useRegularization);

      if (debug || iter % 5 == 0)
         std::cout << "iter = " << iter << " loss: " << loss << std::endl;

      rnn.Backward(XArch, YArch, W);

      rnn.Update(0.1);

   }
   // print output
   rnn.Forward(XArch);

   //Matrix_t & out = rnn.GetLayers().back()->GetOutputAt(0);
   //out.Print();
   //YArch.Print();

   // predictions
   Matrix_t yout(batchSize, 1);
   rnn.Prediction( yout, XArch, EOutputFunction::kSigmoid);

   // predicted output and labels
   TMatrixD result(2, batchSize);
   for (size_t i = 0; i < batchSize; ++i) {
      result(0,i) = YArch(i,0);
      result(1,i) = yout(i,0);
   }
   if (debug) result.Print();

   auto h0 = new TH1D("h0","h0",50,0,1);
   auto h1 = new TH1D("h1","h1",50,0,1);
   // build a roc curve
   std::vector<float> values(batchSize);
   // targets must be a bool vector
   std::vector<bool> targets(batchSize);
   for (size_t i = 0; i < batchSize; ++i) {
      values[i] = yout(i,0);
      targets[i] =  YArch(i,0);
      if (targets[i] == 0) h0->Fill(values[i] );
      else h1->Fill(values[i] );
   }
   TMVA::ROCCurve roc(values, targets);
   double auc = roc.GetROCIntegral();


   if (saveResult)  {
      TFile * fout = TFile::Open(TString::Format("test%sResult.root",rnnType.Data()),"RECREATE");
      h0->Write();
      h1->Write();
      auto g = roc.GetROCCurve();
      g->Write("roc");
      fout->Close();
   }

   // compute efficiency at 0.5

   double eff1 = h0->Integral(25,50)/h0->Integral(1,50);
   double eff2 = h1->Integral(25,50)/h1->Integral(1,50);
   std::cout << "Efficiencies are " << eff1 << " and " << eff2 << std::endl;
   std::cout << "ROC integral is " << auc << std::endl;

   bool ok = (eff1 > 0.9 && eff2 < 0.1) || (eff1 < 0.1 && eff2 > 0.9);
   ok &= (auc > 0.95);
   if (ok)
      Info("testFullRNN2","Test trainig full %s passed ",rnnType.Data());
   else
      Error("testFullRNN2","Test training full %s failed !! ",rnnType.Data());
   return ok;
}


#endif
