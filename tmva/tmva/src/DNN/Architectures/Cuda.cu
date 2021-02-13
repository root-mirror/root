// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 10/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////////////////////////
// Explicit instantiation of the TCuda architecture class with //
// for Double_t and Float_t floating point types.               //
/////////////////////////////////////////////////////////////////

// in case we compile C++ code with std-17 and cuda with lower standard
// use experimental string_view, otherwise keep as is
#include "RConfigure.h"
#ifdef R__HAS_STD_STRING_VIEW
#ifndef R__CUDA_HAS_STD_STRING_VIEW
#undef R__HAS_STD_STRING_VIEW
#define R__HAS_STD_EXPERIMENTAL_STRING_VIEW
#endif
#endif

#include "TMVA/DNN/Architectures/Cuda.h"
#include "Cuda/Propagation.cu"
#include "Cuda/Arithmetic.cu"
#include "Cuda/ActivationFunctions.cu"
#include "Cuda/OutputFunctions.cu"
#include "Cuda/LossFunctions.cu"
#include "Cuda/Regularization.cu"
#include "Cuda/Initialization.cu"
#include "Cuda/Dropout.cu"
#include "Cuda/RecurrentPropagation.cu"

namespace TMVA {
namespace DNN  {

template class TCuda<Float_t>;
template class TCuda<Double_t>;


#ifndef R__HAS_TMVAGPU
   // if R__HAS_TMVAGPU is not defined this file should not be compiled
   static_assert(false,"GPU/CUDA architecture is not enabled");
#endif


} // namespace tmva
} // namespace dnn
