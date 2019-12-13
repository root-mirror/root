// @(#)root/tmva/tmva/cnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing the Upsampling Layer on a CPU architecture                        *
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
// Testing the Upsample function                                //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "TMVA/DNN/Architectures/Cpu.h"
#include "TestUpsampleLayer.h"

int main()
{
    using Scalar_t = Double_t;

    std::cout << "Testing Upsample on the CPU:" << std::endl;

    std::cout << "Test Forward-Propagation 1: " << std::endl;
    if (!testUpsample1<TCpu<Scalar_t>>()) {
        std::cerr << "ERROR - Forward-Propagation 1 failed " << std::endl;
        return -1;
    }

    std::cout << "Test Back-propagation 1: " << std::endl;
    if (!testBackward1<TCpu<Scalar_t>>()) {
        std::cerr << "ERROR - Back-propagation failed " << std::endl;
        return -1;
    }

}
