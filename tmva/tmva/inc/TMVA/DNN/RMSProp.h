// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TRMSProp                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      RMSProp Optimizer Class                                                   *
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

#ifndef TMVA_DNN_RMSPROP
#define TMVA_DNN_RMSPROP

#include "TMatrix.h"
#include "TMVA/DNN/Optimizer.h"
#include "TMVA/DNN/Functions.h"
#include <vector>

namespace TMVA {
namespace DNN {

/** \class TRMSProp
 *  RMSProp Optimizer class
 *
 *  This class represents the RMSProp Optimizer with options for applying momentum.
 */
template <typename Architecture_t, typename Layer_t = VGeneralLayer<Architecture_t>,
          typename DeepNet_t = TDeepNet<Architecture_t, Layer_t>>
class TRMSProp : public VOptimizer<Architecture_t, Layer_t, DeepNet_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

protected:
   Scalar_t fMomentum; ///< The momentum used for training.
   Scalar_t fRho;      ///< The Rho constant used by the optimizer.
   Scalar_t fEpsilon;  ///< The Smoothing term used to avoid division by zero.
   std::vector<std::vector<Matrix_t>>
      fPastSquaredWeightGradients; ///< The sum of the square of the past weight gradients associated with the deep net.
   std::vector<std::vector<Matrix_t>>
      fPastSquaredBiasGradients; ///< The sum of the square of the past bias gradients associated with the deep net.

   std::vector<std::vector<Matrix_t>> fWeightUpdates; ///< The accumulation of the past Weights for performing updates.
   std::vector<std::vector<Matrix_t>> fBiasUpdates;   ///< The accumulation of the past Biases for performing updates.
   std::vector<std::vector<Matrix_t>>
      fWorkWeightTensor1; ///< working tensor used to keep a temporary copy of weights or weight gradients
   std::vector<std::vector<Matrix_t>>
      fWorkBiasTensor1; ///< working tensor used to keep a temporary copy of bias or bias gradients
   std::vector<std::vector<Matrix_t>>
      fWorkWeightTensor2; ///< working tensor used to keep a temporary copy of weights or weight gradients
   std::vector<std::vector<Matrix_t>>
      fWorkBiasTensor2; ///< working tensor used to keep a temporary copy of bias or bias gradients

   /*! Update the weights, given the current weight gradients. */
   void UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights, const std::vector<Matrix_t> &weightGradients);

   /*! Update the biases, given the current bias gradients. */
   void UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases, const std::vector<Matrix_t> &biasGradients);

public:
   /*! Constructor. */
   TRMSProp(DeepNet_t &deepNet, Scalar_t learningRate = 0.001, Scalar_t momentum = 0.0, Scalar_t rho = 0.9,
            Scalar_t epsilon = 1e-7);

   /*! Destructor. */
   ~TRMSProp() = default;

   /*! Getters */
   Scalar_t GetMomentum() const { return fMomentum; }
   Scalar_t GetRho() const { return fRho; }
   Scalar_t GetEpsilon() const { return fEpsilon; }

   std::vector<std::vector<Matrix_t>> &GetPastSquaredWeightGradients() { return fPastSquaredWeightGradients; }
   std::vector<Matrix_t> &GetPastSquaredWeightGradientsAt(size_t i) { return fPastSquaredWeightGradients[i]; }

   std::vector<std::vector<Matrix_t>> &GetPastSquaredBiasGradients() { return fPastSquaredBiasGradients; }
   std::vector<Matrix_t> &GetPastSquaredBiasGradientsAt(size_t i) { return fPastSquaredBiasGradients[i]; }

   std::vector<std::vector<Matrix_t>> &GetWeightUpdates() { return fWeightUpdates; }
   std::vector<Matrix_t> &GetWeightUpdatesAt(size_t i) { return fWeightUpdates[i]; }

   std::vector<std::vector<Matrix_t>> &GetBiasUpdates() { return fBiasUpdates; }
   std::vector<Matrix_t> &GetBiasUpdatesAt(size_t i) { return fBiasUpdates[i]; }
};

//
//
//  The RMSProp Optimizer Class - Implementation
//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
TRMSProp<Architecture_t, Layer_t, DeepNet_t>::TRMSProp(DeepNet_t &deepNet, Scalar_t learningRate, Scalar_t momentum,
                                                       Scalar_t rho, Scalar_t epsilon)
   : VOptimizer<Architecture_t, Layer_t, DeepNet_t>(learningRate, deepNet), fMomentum(momentum), fRho(rho),
     fEpsilon(epsilon)
{
   std::vector<Layer_t *> &layers = deepNet.GetLayers();
   const size_t layersNSlices = layers.size();
   fPastSquaredWeightGradients.resize(layersNSlices);
   fPastSquaredBiasGradients.resize(layersNSlices);
   fWeightUpdates.resize(layersNSlices);
   fBiasUpdates.resize(layersNSlices);
   fWorkWeightTensor1.resize(layersNSlices);
   fWorkBiasTensor1.resize(layersNSlices);
   fWorkWeightTensor2.resize(layersNSlices);
   fWorkBiasTensor2.resize(layersNSlices);

   for (size_t i = 0; i < layersNSlices; i++) {
      const size_t weightsNSlices = (layers[i]->GetWeights()).size();

      Architecture_t::CreateWeightTensors(fPastSquaredWeightGradients[i], layers[i]->GetWeights());
      Architecture_t::CreateWeightTensors(fWeightUpdates[i], layers[i]->GetWeights());

      for (size_t j = 0; j < weightsNSlices; j++) {
         initialize<Architecture_t>(fPastSquaredWeightGradients[i][j], EInitialization::kZero);
         initialize<Architecture_t>(fWeightUpdates[i][j], EInitialization::kZero);
      }

      const size_t biasesNSlices = (layers[i]->GetBiases()).size();

      Architecture_t::CreateWeightTensors( fPastSquaredBiasGradients[i], layers[i]->GetBiases()); 
      Architecture_t::CreateWeightTensors( fBiasUpdates[i], layers[i]->GetBiases()); 

      for (size_t j = 0; j < biasesNSlices; j++) {
         initialize<Architecture_t>(fPastSquaredBiasGradients[i][j], EInitialization::kZero);
         initialize<Architecture_t>(fBiasUpdates[i][j], EInitialization::kZero);
      }
      Architecture_t::CreateWeightTensors(fWorkWeightTensor1[i], layers[i]->GetWeights());
      Architecture_t::CreateWeightTensors(fWorkBiasTensor1[i], layers[i]->GetBiases());
      Architecture_t::CreateWeightTensors(fWorkWeightTensor2[i], layers[i]->GetWeights());
      Architecture_t::CreateWeightTensors(fWorkBiasTensor2[i], layers[i]->GetBiases());
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TRMSProp<Architecture_t, Layer_t, DeepNet_t>::UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights,
                                                                 const std::vector<Matrix_t> &weightGradients) -> void
{
   std::vector<Matrix_t> &currentLayerPastSquaredWeightGradients = this->GetPastSquaredWeightGradientsAt(layerIndex);
   std::vector<Matrix_t> &currentLayerWeightUpdates = this->GetWeightUpdatesAt(layerIndex);

   for (size_t k = 0; k < currentLayerPastSquaredWeightGradients.size(); k++) {

      // accumulation matrix used for temporary storing of the current accumulation
      auto &accumulation = fWorkWeightTensor1[layerIndex][k];
      auto &currentSquaredWeightGradients = fWorkWeightTensor2[layerIndex][k];

      // Vt = rho * Vt-1 + (1-rho) * currentSquaredWeightGradients
      initialize<Architecture_t>(accumulation, EInitialization::kZero);

      Architecture_t::Copy(currentSquaredWeightGradients, weightGradients[k]);
      Architecture_t::SquareElementWise(currentSquaredWeightGradients);
      Architecture_t::ScaleAdd(accumulation, currentLayerPastSquaredWeightGradients[k], this->GetRho());
      Architecture_t::ScaleAdd(accumulation, currentSquaredWeightGradients, 1 - (this->GetRho()));
      Architecture_t::Copy(currentLayerPastSquaredWeightGradients[k], accumulation);

      // Wt = momentum * Wt-1 + (learningRate * currentWeightGradients) / (sqrt(Vt + epsilon))
      initialize<Architecture_t>(accumulation, EInitialization::kZero);
      auto &dummy = fWorkWeightTensor2[layerIndex][k]; // reuse working tensor
      Architecture_t::Copy(dummy, currentLayerPastSquaredWeightGradients[k]);
      Architecture_t::ConstAdd(dummy, this->GetEpsilon());
      Architecture_t::SqrtElementWise(dummy);
      Architecture_t::ReciprocalElementWise(dummy);
      Architecture_t::Hadamard(dummy, weightGradients[k]);

      Architecture_t::ScaleAdd(accumulation, currentLayerWeightUpdates[k], this->GetMomentum());
      Architecture_t::ScaleAdd(accumulation, dummy, this->GetLearningRate());
      Architecture_t::Copy(currentLayerWeightUpdates[k], accumulation);
   }

   // updating the weights.
   // theta = theta - Wt
   for (size_t i = 0; i < weights.size(); i++) {
      Architecture_t::ScaleAdd(weights[i], currentLayerWeightUpdates[i], -1.0);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TRMSProp<Architecture_t, Layer_t, DeepNet_t>::UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases,
                                                                const std::vector<Matrix_t> &biasGradients) -> void
{
   std::vector<Matrix_t> &currentLayerPastSquaredBiasGradients = this->GetPastSquaredBiasGradientsAt(layerIndex);
   std::vector<Matrix_t> &currentLayerBiasUpdates = this->GetBiasUpdatesAt(layerIndex);

   for (size_t k = 0; k < currentLayerPastSquaredBiasGradients.size(); k++) {

      // accumulation matrix used for temporary storing of the current accumulation
      auto &accumulation = fWorkBiasTensor1[layerIndex][k];
      auto &currentSquaredBiasGradients = fWorkBiasTensor2[layerIndex][k];
      
      // Vt = rho * Vt-1 + (1-rho) * currentSquaredBiasGradients
      initialize<Architecture_t>(accumulation, EInitialization::kZero);
      Architecture_t::Copy(currentSquaredBiasGradients, biasGradients[k]);
      Architecture_t::SquareElementWise(currentSquaredBiasGradients);
      Architecture_t::ScaleAdd(accumulation, currentLayerPastSquaredBiasGradients[k], this->GetRho());
      Architecture_t::ScaleAdd(accumulation, currentSquaredBiasGradients, 1 - (this->GetRho()));
      Architecture_t::Copy(currentLayerPastSquaredBiasGradients[k], accumulation);

      // Wt = momentum * Wt-1 + (learningRate * currentBiasGradients) / (sqrt(Vt + epsilon))
      initialize<Architecture_t>(accumulation, EInitialization::kZero);
      auto &dummy = fWorkBiasTensor2[layerIndex][k]; // reuse working tensor

      Architecture_t::Copy(dummy, currentLayerPastSquaredBiasGradients[k]);
      Architecture_t::ConstAdd(dummy, this->GetEpsilon());
      Architecture_t::SqrtElementWise(dummy);
      Architecture_t::ReciprocalElementWise(dummy);
      Architecture_t::Hadamard(dummy, biasGradients[k]);

      Architecture_t::ScaleAdd(accumulation, currentLayerBiasUpdates[k], this->GetMomentum());
      Architecture_t::ScaleAdd(accumulation, dummy, this->GetLearningRate());
      Architecture_t::Copy(currentLayerBiasUpdates[k], accumulation);
   }

   // updating the Biases.
   // theta = theta - Wt
   for (size_t i = 0; i < biases.size(); i++) {
      Architecture_t::ScaleAdd(biases[i], currentLayerBiasUpdates[i], -1.0);
   }
}

} // namespace DNN
} // namespace TMVA

#endif