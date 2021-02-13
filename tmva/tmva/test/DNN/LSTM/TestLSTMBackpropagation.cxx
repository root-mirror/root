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
// Testing LSTMLayer backpropagation                               //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <climits>
#include "TMVA/DNN/Architectures/Reference.h"
#include "TMath.h"
#include "TestLSTMBackpropagation.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {
   std::cout << "Testing LSTM backward pass\n";

   // timesteps, batchsize, statesize, inputsize
//   testLSTMBackpropagation<TReference<double>>(1, 2, 1, 10, 1e-5);
   gRandom->SetSeed(12345);
   TReference<double>::SetRandomSeed(gRandom->Integer(TMath::Limits<UInt_t>::Max()));


   testLSTMBackpropagation<TReference<double>>(2, 1, 1, 2, 1e-5);

   return 0;


   testLSTMBackpropagation<TReference<double>>(1, 2, 2, 10, 1e-5);

   testLSTMBackpropagation<TReference<double>>(2, 1, 2, 5, 1e-5);

   testLSTMBackpropagation<TReference<double>>(4, 2, 3, 10, 1e-5);
   // using a fixed input
   testLSTMBackpropagation<TReference<double>>(3, 1, 4, 5, 1e-5, {true});
   // with a dense layer
   testLSTMBackpropagation<TReference<double>>(4, 32, 10, 5, 1e-5, {false, true});
   // with an additional LSTM layer
   testLSTMBackpropagation<TReference<double>>(4, 32, 10, 5, 1e-5, {false, true, true});

   return 0;
}
