// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing various optimizers for Cpu Backend                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ravi Kiran S      <sravikiran0606@gmail.com>  - CERN, Switzerland         *
 *                                                                                *
 * Copyright (c) 2005-2018:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TestOptimization.h"
#include "TMVA/DNN/Architectures/Cuda.h"

#include <iostream>

using namespace TMVA::DNN;
using TMVA::DNN::EOptimizer;

int main()
{
   using Architecture = TCuda<Float_t>;
   gRandom->SetSeed(111);
   Architecture::SetRandomSeed(112);

   std::cout << "Testing optimization: (single precision) for Cuda with cudnn" << std::endl;

   Real_t momentumSinglePrecision = 0.0;
   std::cout << "Stochastic Gradient Descent: ";
   Double_t error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kSGD, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: SGD test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.9;
   std::cout << "Stochastic Gradient Descent with momentum: ";
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kSGD, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: SGD momentum test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.0;
   std::cout << "Adam Optimizer: ";
   // Adam doesn't use momentum. Passing this as a parameter just to match the function prototype.
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kAdam, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: ADAM test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.0;
   std::cout << "Adagrad Optimizer: ";
   // Adagrad doesn't use momentum. Passing this as a parameter just to match the function prototype.
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kAdagrad, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: Adagrad test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.0;
   std::cout << "RMSProp Optimizer without momentum: ";
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kRMSProp, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: RMSPROP test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.9;
   std::cout << "RMSProp Optimizer with momentum: ";
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kRMSProp, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: RMSPROP momentum test failed ! " << std::endl;
      return 1;
   }

   momentumSinglePrecision = 0.0;
   std::cout << "Adadelta Optimizer: ";
   // Adadelta doesn't use momentum. Passing this as a parameter just to match the function prototype.
   error = testOptimization<Architecture>(momentumSinglePrecision, EOptimizer::kAdadelta, false);
   std::cout << "After Training: Mean Absolute error = " << error << std::endl;
   if (error > 1e-1) {
      std::cout << "---->   Error: Adadelta test failed ! " << std::endl;
      return 1;
   }

   return 0;
}
