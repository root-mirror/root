// @(#)root/tmva/tmva/cnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Downsample method on a CPU architecture                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ashish Kshirsagar       <ashishkshirsagar10@gmail.com>                    *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

////////////////////////////////////////////////////////////////////
// Testing the Transpose Convolutional Layer                      //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "TMVA/DNN/Architectures/Cpu.h"
#include "TestTransConvLayer.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

int main()
{
   using Scalar_t = Double_t;

   std::cout << "Testing Forward Propagation of the Transpose Convolutional Layer on the CPU:" << std::endl;

   bool status = true;

   std::cout << "Test Forward-Propagation 1: " << std::endl;
   status &= testForward1<TCpu<Scalar_t>>();
   if (!status) {
      std::cerr << "ERROR - Forward-Propagation 1 failed " << std::endl;
      return -1;
   }
   
   std::cout << "Test Backward-Propagation 1: " << std::endl;
   status &= testBackward1<TCpu<Scalar_t>>();
   if (!status) {
      std::cerr << "ERROR - Backward-Propagation 1 failed " << std::endl;
      return -1;
   }
   std::cout << "All tests passed!" << std::endl;
}
