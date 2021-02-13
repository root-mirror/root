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
#include "TMVA/DNN/Architectures/Cpu.h"
#include "TestGRUForwardPass.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {

   using Architecture_t = TCpu<Double_t>;

   std::cout << "Testing GRU Forward pass\n";

   // timesteps, batchsize, statesize, inputsize
   bool ok = true;
   ok &= testForwardPass<Architecture_t>(1, 2, 3, 2);
   ok &= testForwardPass<Architecture_t>(2, 4, 10, 5);
   ok &= testForwardPass<Architecture_t>(5, 8, 5, 10);
   if (ok) {
      Info("testGRUForwardPassCpu", "All GRU Forward tests passed");
   } else {
      Error("testGRUForwardPassCpu", "GRU Forward pass test failed !");
   }


   return (ok) ? 0 : -1;
}
