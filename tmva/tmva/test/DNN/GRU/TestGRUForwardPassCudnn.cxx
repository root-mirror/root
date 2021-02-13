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
// Testing GRU-Layer forward pass for CPU implementation   //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <climits>
#include "TMVA/DNN/Architectures/Cpu.h"
#include "TMVA/DNN/Architectures/TCudnn.h"
#include "TestGRUForwardPass.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {

   using Arch1 = TCpu<Double_t>;
   using Arch2 = TCudnn<Double_t>;

   int seed = 12345;
   gRandom->SetSeed(seed);
   Arch1::SetRandomSeed(gRandom->Integer(INT_MAX));

   bool debug = false;

   std::cout << "Testing GRU Forward pass in both CPU and GPU\n";

   // note need to have afterGate = true in GRU to have same implementation cudnn and cpu

   // timesteps, batchsize, statesize, inputsize , debug , output full seq , use fized input
   bool ok = true;
   ok &= CompareForwardPass<Arch1,Arch2>(1, 2, 3, 2, true, false, debug);
   // using a fixed input with input dim = 1
   ok &= CompareForwardPass<Arch1, Arch2>(2, 1, 2, 1, true, true, debug);
   //return (ok) ? 0 : -1;
   // using a fixed input (input dim = 3)
   ok &= CompareForwardPass<Arch1, Arch2>(2, 1, 2, 3, true, true, debug);

   ok &= CompareForwardPass<Arch1, Arch2>(2, 1, 2, 3, false, false, debug);
   // output full seq
   ok &= CompareForwardPass<Arch1,Arch2>(1, 4, 10, 5, true, false, debug);
   // output only last time
   ok &= CompareForwardPass<Arch1, Arch2>(5, 8, 5, 10, false);
   if (ok) {
      Info("testLSTMForwardPassCudnn", "All GRU Forward tests passed");
   } else {
      Error("testLSTMForwardPassCudnn", "GRU Forward pass test failed !");
   }

   return (ok) ? 0 : -1;
}
