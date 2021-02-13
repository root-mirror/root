// @(#)root/tmva $Id$
// Author: Surya S Dwivedi 02/07/19

/*************************************************************************
 * Copyright (C) 2019, Surya S Dwivedi                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Testing full GRU network (for Reference)                      //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Architectures/Reference.h"
#include "TestFullGRU.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::RNN;


int main() {
   std::cout << "Training GRU to identity first";

   //testFullGRU(size_t batchSize, size_t stateSize, size_t inputSize, size_t outputSize)
   // reconstruct 8 bit vector
   // batchsize, statesize, inputsize, outputsize
   testFullGRU<TReference<double>>(2, 3, 2, 2) ;
   //testFullGRU<TReference<double>>(64, 10, 8, 8) ;
   //testFullGRU<TReference<double>>(3, 8, 100, 50) ;

   // test a full GRU with 5 time steps and different signal/backgrund time dependent shapes
   // batchsize, statesize , inputsize, seed
   int seed = 111; 
   testFullGRU2<TReference<double>>(64, 10, 5, seed) ;

   return 0;
}
