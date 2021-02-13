// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VOptimizer                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      General Optimizer Class                                                   *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ravi Kiran S      <sravikiran0606@gmail.com>  - CERN, Switzerland         *
 *                                                                                *
 * Copyright (c) 2005-2018 :                                                      *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef TMVA_DNN_OPTIMIZER
#define TMVA_DNN_OPTIMIZER

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/DeepNet.h"
#include <vector>

namespace TMVA {
namespace DNN {

/** \class VOptimizer
    Generic Optimizer class

    This class represents the general class for all optimizers in the Deep Learning
    Module.
 */
template <typename Architecture_t, typename Layer_t = VGeneralLayer<Architecture_t>,
          typename DeepNet_t = TDeepNet<Architecture_t, Layer_t>>
class VOptimizer {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

protected:
   Scalar_t fLearningRate; ///< The learning rate used for training.
   size_t fGlobalStep;     ///< The current global step count during training.
   DeepNet_t &fDeepNet;    ///< The reference to the deep net.

   /*! Update the weights, given the current weight gradients. */
   virtual void
   UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights, const std::vector<Matrix_t> &weightGradients) = 0;

   /*! Update the biases, given the current bias gradients. */
   virtual void
   UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases, const std::vector<Matrix_t> &biasGradients) = 0;

public:
   /*! Constructor. */
   VOptimizer(Scalar_t learningRate, DeepNet_t &deepNet);

   /*! Performs one step of optimization. */
   void Step();

   /*! Virtual Destructor. */
   virtual ~VOptimizer() = default;

   /*! Increments the global step. */
   void IncrementGlobalStep() { this->fGlobalStep++; }

   /*! Getters */
   Scalar_t GetLearningRate() const
   {
      return fLearningRate;
   }
   size_t GetGlobalStep() const { return fGlobalStep; }
   std::vector<Layer_t *> &GetLayers() { return fDeepNet.GetLayers(); }
   Layer_t *GetLayerAt(size_t i) { return fDeepNet.GetLayerAt(i); }

   /*! Setters */
   void SetLearningRate(size_t learningRate) { fLearningRate = learningRate; }
};

//
//
//  The General Optimizer Class - Implementation
//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
VOptimizer<Architecture_t, Layer_t, DeepNet_t>::VOptimizer(Scalar_t learningRate, DeepNet_t &deepNet)
   : fLearningRate(learningRate), fGlobalStep(0), fDeepNet(deepNet)
{
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto VOptimizer<Architecture_t, Layer_t, DeepNet_t>::Step() -> void
{
   for (size_t i = 0; i < this->GetLayers().size(); i++) {
      this->UpdateWeights(i, this->GetLayerAt(i)->GetWeights(), this->GetLayerAt(i)->GetWeightGradients());
      this->UpdateBiases(i, this->GetLayerAt(i)->GetBiases(), this->GetLayerAt(i)->GetBiasGradients());
   }
}

} // namespace DNN
} // namespace TMVA

#endif