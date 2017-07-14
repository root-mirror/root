// @(#)root/tmva $Id$
// Author: Akshay Vashistha (ajatgd)

/*************************************************************************
 * Copyright (C) 2017, ajatgd
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TMVA/DNN/Architectures/Reference.h"
#include "TestStackedLayers.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::DAE;

int main()

{
  std::cout<<"Testing started"<<std::endl;

  testStacked<TReference<double>>();
   //std::cout<<"Transform Layer"<<std::endl;
  // testTransform<TReference<double>>();

 return 0;
}
