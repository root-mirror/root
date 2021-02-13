// @(#)root/tmva $Id$
// Author: Saurav Shekhar 16/02/17

/*************************************************************************
 * Copyright (C) 2017, Saurav Shekhar                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Testing RNNLayer backpropagation                               //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Architectures/TCudnn.h"
#include "TestRecurrentBackpropagation.h"
#include "TROOT.h"
#include "TMath.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {

   bool debug = false;
   std::cout << "Testing RNN backward pass on GPU using Cudnn\n";

   //ROOT::EnableImplicitMT(1);

   using Scalar_t = Double_t;
   using Architecture_t = TCudnn<Scalar_t>;

   int seed = 12345;
   gRandom->SetSeed(seed);
   Architecture_t::SetRandomSeed(gRandom->Integer(TMath::Limits<UInt_t>::Max()));

   bool fail = false;

   // timesteps, batchsize, statesize, inputsize  { fixed input, with dense layer, with extra RNN, return full sequence, modify weights }

   fail |= testRecurrentBackpropagation<Architecture_t>(3, 1, 5, 4, 1e-5, {true, true, false}, debug);
   if (fail) {
      Error("testRecurrentBackPropagationCpu", "Fixed test failed");
      if (!debug) {
         Info("testRecurrentBackPropagationCpu", "Rerun test in Debug mode");
         testRecurrentBackpropagation<Architecture_t>(3, 1, 5, 4, 1e-5, {true, true, false}, true);
      }
      return fail;
   }

   fail |= testRecurrentBackpropagation<Architecture_t>(1, 2, 1, 2, 1e-5);

   fail |= testRecurrentBackpropagation<Architecture_t>(1, 2, 3, 2, 1e-5, {}, debug);

   fail |= testRecurrentBackpropagation<Architecture_t>(2, 3, 4, 5, 1e-5, {}, debug);

   // test returning the full sequence
   fail |= testRecurrentBackpropagation<Architecture_t>(2, 3, 4, 5, 1e-5, {false, false, false, true}, debug);


   fail |= testRecurrentBackpropagation<Architecture_t>(4, 2, 10, 5, 1e-5);

   // use batch size <= 16 or get numerical error in numrical gradients
   fail |= testRecurrentBackpropagation<Architecture_t>(5, 16, 10, 5, 1e-5,{});

   // using a fixed input
   fail |= testRecurrentBackpropagation<Architecture_t>(3, 1, 10, 5, 1e-5, {true});

   // with a dense layer (also use not too large batch size or input size to avoid numerical errors)
   fail |= testRecurrentBackpropagation<Architecture_t>(4, 16, 10, 8, 1e-5, {false, true});

   // test returning the full sequence and dense layer
   fail |= testRecurrentBackpropagation<Architecture_t>(3, 8, 5, 4, 1e-5, {false, true, false, true}, debug);

   // with an additional RNN layer and dense layer
   fail |= testRecurrentBackpropagation<Architecture_t>(2, 1, 1, 2, 1e-5, {false, true, true, false}, debug);

   if (fail)
      Error("testRecurrentBackPropagationCudnn", "Test failed !!!");
   else
      Info("testRecurrentPropagationCudnn", "All tests passed !!!");

   return fail;
}
