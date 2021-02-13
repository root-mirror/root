// @(#)root/tmva/tmva/cnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Method DL for Conv Net for the Reference backend                  *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievski      <ilievski.vladimir@live.com>  - CERN, Switzerland  *
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

#include "TestMethodDLCNN.h"
#include "TString.h"

int main()
{

   std::cout << "Testing Method DL for CPU backend: " << std::endl;

   TString archCPU = "CPU";
   testMethodDL_CNN(archCPU);

#ifdef R__HAS_TMVAGPU
   std::cout << "Testing Method DL for GPU backend: " << std::endl;
#ifdef R__HAS_CUDNN
   std::cout << "Using cuDNN for the implementation of the convolution operators " << std::endl;
#endif

   TString archGPU  = "GPU";
   testMethodDL_CNN(archGPU);
#endif

   return 0;
}
