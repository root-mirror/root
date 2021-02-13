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
// Testing GRU-Layer forward pass for Reference implementation   //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Architectures/Reference.h"
#include "TestGRUForwardPass.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {

   std::cout << "Testing GRU Forward pass:\n";
   
   // timesteps, batchsize, statesize, inputsize
   std::cout << testForwardPass<TReference<double>>(1, 2, 3, 2)  << "\n";
   std::cout << testForwardPass<TReference<double>>(1, 8, 100, 50)  << "\n";
   std::cout << testForwardPass<TReference<double>>(5, 9, 128, 64)  << "\n";

   return 0;
}
