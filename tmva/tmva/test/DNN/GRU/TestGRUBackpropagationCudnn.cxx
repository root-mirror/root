// @(#)root/tmva $Id$
// Author: Surya S Dwivedi 26/06/2019

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
#include <climits>
#include "TMVA/DNN/Architectures/TCudnn.h"
#include "TestGRUBackpropagation.h"
#include "TROOT.h"
#include "TMath.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;

int main() {
   std::cout << "Testing GRU backward pass on GPU with CuDNN\n";


   using Scalar_t = Double_t;
   using Architecture_t = TCudnn<Scalar_t>;

   int seed = 12345;
   gRandom->SetSeed(seed);
   Architecture_t::SetRandomSeed(gRandom->Integer(INT_MAX));

   using Scalar_t = Double_t;

   bool debug = false;

   bool iret = false;

   // timesteps, batchsize, statesize, inputsize  { fixed input, with dense layer, with extra GRU , output full sequence}
   iret |= testGRUBackpropagation<Architecture_t>(2, 1, 2, 3, 1e-5,{true,false,false}, debug);

   iret |= testGRUBackpropagation<Architecture_t>(1, 2, 1, 10, 1e-5, {}, debug);


   iret |= testGRUBackpropagation<Architecture_t>(1, 2, 3, 2, 1e-5);

   iret |= testGRUBackpropagation<Architecture_t>(2, 3, 4, 5, 1e-5);

   iret |= testGRUBackpropagation<Architecture_t>(4, 2, 10, 5, 1e-5);

   iret |= testGRUBackpropagation<Architecture_t>(5, 16, 10, 20, 1e-5);

   // using a fixed input (input size must be <=6, time steps <=3 and batch size <=1  )
   iret |= testGRUBackpropagation<Architecture_t>(3, 1, 10, 5, 1e-5, {true, true}, debug);

   // test with a dense layer
   iret |= testGRUBackpropagation<Architecture_t>(4, 8, 20, 10, 1e-5, {false, true, false}, debug);

   // test returning the full sequence and dense layer
   iret |= testGRUBackpropagation<Architecture_t>(3, 8, 5, 4, 1e-5, {false, true, false, true}, debug);

   // with an additional GRU layer
   iret |= testGRUBackpropagation<Architecture_t>(4, 8, 10, 5, 1e-5, {false, true, true});

   if (iret)
      Error("testGRUBackPropagationCudnn", "Test failed !!!");
   else
      Info( "testGRUBackPropagationCudnn", "All tests passed !!!");

   return iret;
}
