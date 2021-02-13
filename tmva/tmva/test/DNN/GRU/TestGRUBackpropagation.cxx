// @(#)root/tmva $Id$
// Author: Surya S Dwivedi 26/06219

/*************************************************************************
 * Copyright (C) 2019, Surya S Dwivedi                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Testing GRULayer backpropagation                               //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Architectures/Reference.h"
#include "TestGRUBackpropagation.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::GRU;

int main() {
   std::cout << "Testing GRU backward pass\n";

   // timesteps, batchsize, statesize, inputsize
   testGRUBackpropagation<TReference<double>>(1, 2, 1, 10, 1e-5);
   
   testGRUBackpropagation<TReference<double>>(2, 1, 1, 2, 1e-5);
  
   testGRUBackpropagation<TReference<double>>(4, 2, 3, 10, 1e-5);
   // using a fixed input 
   testGRUBackpropagation<TReference<double>>(3, 1, 4, 5, 1e-5, {true});
   // with a dense layer 
   testGRUBackpropagation<TReference<double>>(4, 32, 10, 5, 1e-5, {false, true});
   // with an additional GRU layer 
   testGRUBackpropagation<TReference<double>>(4, 32, 10, 5, 1e-5, {false, true, true});

   return 0;
}
