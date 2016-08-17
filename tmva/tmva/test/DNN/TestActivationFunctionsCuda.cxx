// @(#)root/tmva $Id$
// Author: Simon Pfreundschuh

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////
//  Concrete instantiation of the generic activation function test  //
//  for the TCuda implementation.                                    //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Architectures/Cuda.h"
#include "Utility.h"
#include "TestActivationFunctions.h"

using namespace TMVA::DNN;

int main()
{
    std::cout << "Testing Activation Functions:" << std::endl;

    double error;

    // Identity.

    error = testIdentity<TCuda>(10);
    std::cout << "Testing identity activation:            ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;

    error = testIdentityDerivative<TCuda>(10);
    std::cout << "Testing identity activation derivative: ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;

    // ReLU.

    error = testRelu<TCuda>(10);
    std::cout << "Testing ReLU activation:                ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;

    error = testReluDerivative<TCuda>(10);
    std::cout << "Testing ReLU activation derivative:     ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;

    // Sigmoid.

    error = testSigmoid<TCuda>(10);
    std::cout << "Testing Sigmoid activation:             ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;

    error = testSigmoidDerivative<TCuda>(10);
    std::cout << "Testing Sigmoid activation derivative:  ";
    std::cout << "maximum relative error = " << error << std::endl;
    if (error > 1e-10)
        return 1;
    return 0;
}
