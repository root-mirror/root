// @(#)root/tmva/tmva/cnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing the Pooling Layer on a CPU architecture                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ashish Kshirsagar      <ashishkshirsagar10@gmail.com>                     *
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
// Testing the Downsample function                                //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "TMVA/DNN/Architectures/Reference.h"
#include "TestUpsampleLayer.h"

int main()
{
    using Scalar_t = Double_t;

    std::cout << "Testing Upsample on the GPU:" << std::endl;

    std::cout << "Test Forward-Propagation 1: " << std::endl;
    if (!testUpsample1<TReference<Scalar_t>>()) {
        std::cerr << "ERROR - Forward-Propagation 1 failed " << std::endl;
        return -1;
    }

    std::cout << "Test Back-propagation: " << std::endl;
    if (!testBackward1<TReference<Scalar_t>>()) {
        std::cerr << "ERROR - Back-propagation failed " << std::endl;
        return -1;
    }

}
