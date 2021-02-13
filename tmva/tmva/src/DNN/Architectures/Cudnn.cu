// @(#)root/tmva/tmva/dnn:$Id$
// Author: Joana Niermann 23/07/19

/*************************************************************************
 * Copyright (C) 2019 Joana Niermann                                     *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////
// Explicit instantiation of the TCudnn architecture class with  //
// for Double_t and Real_t floating point types.                 //
///////////////////////////////////////////////////////////////////

// in case we compile C++ code with std-17 and cuda with lower standard
// use experimental string_view, otherwise keep as is
#include "RConfigure.h"
#ifdef R__HAS_STD_STRING_VIEW
#ifndef R__CUDA_HAS_STD_STRING_VIEW
#undef R__HAS_STD_STRING_VIEW
#define R__HAS_STD_EXPERIMENTAL_STRING_VIEW
#endif
#endif

#include "TMVA/DNN/Architectures/TCudnn.h"
#include "Cudnn/Propagate.cu"
#include "Cudnn/Dropout.cu"
#include "Cudnn/Arithmetic.cu"
#include "Cudnn/ActivationFunctions.cu"
#include "Cudnn/Initialization.cu"
#include "Cudnn/LossFunctions.cu"
//#include "Cudnn/Regularization.cu"
#include "Cudnn/OutputFunctions.cu"

//#include "Cudnn/Dropout.cu"
#include "Cudnn/RecurrentPropagation.cu"

namespace TMVA {
namespace DNN  {

template class TCudnn<Float_t>;
template class TCudnn<Double_t>;

#ifndef R__HAS_TMVAGPU
   // if R__HAS_TMVAGPU is not defined this file should not be compiled
   static_assert(false,"GPU/CUDA architecture is not enabled");
#endif

   // CuDNN options


// int TCudnn<Double_t>::CNNOptions::ConvFwdAlgorithm = -1;
// int TCudnn<Double_t>::CNNOptions::ConvBwdDataAlgorithm = -1;
// int TCudnn<Double_t>::CNNOptions::ConvBwdFilterAlgorithm = -1;
// size_t TCudnn<Double_t>::CNNOptions::ConvMaxWorkspaceSize = 0;

} // end namespace DNN
} // end namespace TMVA
