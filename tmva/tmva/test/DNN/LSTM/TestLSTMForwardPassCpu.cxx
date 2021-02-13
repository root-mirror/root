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
// Testing LSTM-Layer forward pass for CPU implementation   //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <climits>
#include "TMVA/DNN/Architectures/Cpu.h"
#include "TestLSTMForwardPass.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {

   using Architecture_t = TCpu<Double_t>;

   int seed = 12345;
   gRandom->SetSeed(seed);
   Architecture_t::SetRandomSeed(gRandom->Integer(INT_MAX));

   std::cout << "Testing LSTM Forward pass\n";
   bool debug = true;
   // timesteps, batchsize, statesize, inputsize
   bool ok = true;
   ok &= testForwardPass<Architecture_t>(1, 2, 3, 2, debug);
   ok &= testForwardPass<Architecture_t>(2, 4, 10, 5);
   ok &= testForwardPass<Architecture_t>(5, 8, 5, 10);
   if (ok) {
      Info("testLSTMForwardPassCpu", "All LSTM Forward tests passed");
   } else {
      Error("testLSTMForwardPassCpu", "LSTM Forward pass test failed !");
   }

   return (ok) ? 0 : -1;
}