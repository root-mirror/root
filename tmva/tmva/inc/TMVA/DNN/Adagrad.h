// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TAdagrad                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Adagrad Optimizer Class                                                      *
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

#ifndef TMVA_DNN_ADAGRAD
#define TMVA_DNN_ADAGRAD

#include "TMatrix.h"
#include "TMVA/DNN/Optimizer.h"
#include "TMVA/DNN/Functions.h"
#include <vector>

namespace TMVA {
namespace DNN {

/** \class TAdagrad
 *  Adagrad Optimizer class
 *
 *  This class represents the Adagrad Optimizer.
 */
template <typename Architecture_t, typename Layer_t = VGeneralLayer<Architecture_t>,
          typename DeepNet_t = TDeepNet<Architecture_t, Layer_t>>
class TAdagrad : public VOptimizer<Architecture_t, Layer_t, DeepNet_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

protected:
   Scalar_t fEpsilon; ///< The Smoothing term used to avoid division by zero.

   std::vector<std::vector<Matrix_t>>
      fPastSquaredWeightGradients; ///< The sum of the square of the past weight gradients associated with the deep net.
   std::vector<std::vector<Matrix_t>>
      fPastSquaredBiasGradients; ///< The sum of the square of the past bias gradients associated with the deep net.
   std::vector<std::vector<Matrix_t>>
      fWorkWeightTensor; ///< working tensor used to keep a temporary copy of weights or weight gradients
   std::vector<std::vector<Matrix_t>>
      fWorkBiasTensor; ///< working tensor used to keep a temporary copy of bias or bias gradients

   /*! Update the weights, given the current weight gradients. */
   void UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights, const std::vector<Matrix_t> &weightGradients);

   /*! Update the biases, given the current bias gradients. */
   void UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases, const std::vector<Matrix_t> &biasGradients);

public:
   /*! Constructor. */
   TAdagrad(DeepNet_t &deepNet, Scalar_t learningRate = 0.01, Scalar_t epsilon = 1e-8);

   /*! Destructor. */
   ~TAdagrad() = default;

   /*! Getters */
   Scalar_t GetEpsilon() const { return fEpsilon; }

   std::vector<std::vector<Matrix_t>> &GetPastSquaredWeightGradients() { return fPastSquaredWeightGradients; }
   std::vector<Matrix_t> &GetPastSquaredWeightGradientsAt(size_t i) { return fPastSquaredWeightGradients[i]; }

   std::vector<std::vector<Matrix_t>> &GetPastSquaredBiasGradients() { return fPastSquaredBiasGradients; }
   std::vector<Matrix_t> &GetPastSquaredBiasGradientsAt(size_t i) { return fPastSquaredBiasGradients[i]; }
};

//
//
//  The Adagrad Optimizer Class - Implementation
//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
TAdagrad<Architecture_t, Layer_t, DeepNet_t>::TAdagrad(DeepNet_t &deepNet, Scalar_t learningRate, Scalar_t epsilon)
   : VOptimizer<Architecture_t, Layer_t, DeepNet_t>(learningRate, deepNet), fEpsilon(epsilon)
{
   std::vector<Layer_t *> &layers = deepNet.GetLayers();
   const size_t layersNSlices = layers.size();
   fPastSquaredWeightGradients.resize(layersNSlices);
   fPastSquaredBiasGradients.resize(layersNSlices);
   fWorkWeightTensor.resize(layersNSlices);
   fWorkBiasTensor.resize(layersNSlices);

   for (size_t i = 0; i < layersNSlices; i++) {
      const size_t weightsNSlices = (layers[i]->GetWeights()).size();

      // weight and weight gradients  tensors should have same
      Architecture_t::CreateWeightTensors( fPastSquaredWeightGradients[i], layers[i]->GetWeights()); 

      for (size_t j = 0; j < weightsNSlices; j++) {
         initialize<Architecture_t>(fPastSquaredWeightGradients[i][j], EInitialization::kZero);
      }

      const size_t biasesNSlices = (layers[i]->GetBiases()).size();

      Architecture_t::CreateWeightTensors( fPastSquaredBiasGradients[i], layers[i]->GetBiases()); 

      for (size_t j = 0; j < biasesNSlices; j++) {
         initialize<Architecture_t>(fPastSquaredBiasGradients[i][j], EInitialization::kZero);
      }

      Architecture_t::CreateWeightTensors(fWorkWeightTensor[i], layers[i]->GetWeights());
      Architecture_t::CreateWeightTensors(fWorkBiasTensor[i], layers[i]->GetBiases());

   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TAdagrad<Architecture_t, Layer_t, DeepNet_t>::UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights,
                                                                 const std::vector<Matrix_t> &weightGradients) -> void
{
   auto &currentLayerPastSquaredWeightGradients = this->GetPastSquaredWeightGradientsAt(layerIndex);
  

   const size_t weightsNSlices = weights.size();
   assert(currentLayerPastSquaredWeightGradients.size() == weightsNSlices);

   for (size_t i = 0; i < weightsNSlices; i++) {

      auto &currentSquaredWeightGradients = fWorkWeightTensor[layerIndex][i];
      // Vt = Vt-1 + currentSquaredWeightGradients
      Architecture_t::Copy(currentSquaredWeightGradients, weightGradients[i]);
      Architecture_t::SquareElementWise(currentSquaredWeightGradients);
      Architecture_t::ScaleAdd(currentLayerPastSquaredWeightGradients[i], currentSquaredWeightGradients, 1.0);

      // updating the weights.
      // theta = theta - learningRate * currentWeightGradients / (sqrt(Vt + epsilon))

      auto &currentWeightUpdates = fWorkWeightTensor[layerIndex][i]; // reuse the work tensor for the weight updates now
      Architecture_t::Copy(currentWeightUpdates, currentLayerPastSquaredWeightGradients[i]);
      Architecture_t::ConstAdd(currentWeightUpdates, this->GetEpsilon());
      Architecture_t::SqrtElementWise(currentWeightUpdates);
      Architecture_t::ReciprocalElementWise(currentWeightUpdates);
      Architecture_t::Hadamard(currentWeightUpdates, weightGradients[i]);
      Architecture_t::ScaleAdd(weights[i], currentWeightUpdates, -this->GetLearningRate());
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TAdagrad<Architecture_t, Layer_t, DeepNet_t>::UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases,
                                                                const std::vector<Matrix_t> &biasGradients) -> void
{
   std::vector<Matrix_t> &currentLayerPastSquaredBiasGradients = this->GetPastSquaredBiasGradientsAt(layerIndex);

   const size_t biasesNSlices = biases.size();
   assert(currentLayerPastSquaredBiasGradients.size() == biasesNSlices);
   for (size_t i = 0; i < biasesNSlices; i++) {

      // Vt = Vt-1 + currentSquaredBiasGradients
      auto &currentSquaredBiasGradients = fWorkBiasTensor[layerIndex][i];
      Architecture_t::Copy(currentSquaredBiasGradients, biasGradients[i]);
      Architecture_t::SquareElementWise(currentSquaredBiasGradients);
      Architecture_t::ScaleAdd(currentLayerPastSquaredBiasGradients[i], currentSquaredBiasGradients, 1.0);

      // updating the biases.
      // theta = theta - learningRate * currentBiasGradients / (sqrt(Vt + epsilon))

      auto &currentBiasUpdates = fWorkBiasTensor[layerIndex][i];  
      Architecture_t::Copy(currentBiasUpdates, currentLayerPastSquaredBiasGradients[i]);
      Architecture_t::ConstAdd(currentBiasUpdates, this->GetEpsilon());
      Architecture_t::SqrtElementWise(currentBiasUpdates);
      Architecture_t::ReciprocalElementWise(currentBiasUpdates);
      Architecture_t::Hadamard(currentBiasUpdates, biasGradients[i]);
      Architecture_t::ScaleAdd(biases[i], currentBiasUpdates, -this->GetLearningRate());
   }
}

} // namespace DNN
} // namespace TMVA

#endif
