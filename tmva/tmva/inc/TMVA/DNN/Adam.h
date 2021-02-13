// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ravi Kiran S

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TAdam                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Adam Optimizer Class                                                      *
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

#ifndef TMVA_DNN_ADAM
#define TMVA_DNN_ADAM

#include "TMatrix.h"
#include "TMVA/DNN/Optimizer.h"
#include "TMVA/DNN/Functions.h"
#include <vector>

namespace TMVA {
namespace DNN {

/** \class TAdam
 *  Adam Optimizer class
 *
 *  This class represents the Adam Optimizer.
 */
template <typename Architecture_t, typename Layer_t = VGeneralLayer<Architecture_t>,
          typename DeepNet_t = TDeepNet<Architecture_t, Layer_t>>
class TAdam : public VOptimizer<Architecture_t, Layer_t, DeepNet_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

protected:
   Scalar_t fBeta1;   ///< The Beta1 constant used by the optimizer.
   Scalar_t fBeta2;   ///< The Beta2 constant used by the optimizer.
   Scalar_t fEpsilon; ///< The Smoothing term used to avoid division by zero.

   std::vector<std::vector<Matrix_t>> fFirstMomentWeights; ///< The decaying average of the first moment of the past
                                                           /// weight gradients associated with the deep net.
   std::vector<std::vector<Matrix_t>> fFirstMomentBiases; ///< The decaying average of the first moment of the past bias
                                                          /// gradients associated with the deep net.

   std::vector<std::vector<Matrix_t>> fSecondMomentWeights; ///< The decaying average of the second moment of the past
                                                            /// weight gradients associated with the deep net.
   std::vector<std::vector<Matrix_t>> fSecondMomentBiases;  ///< The decaying average of the second moment of the past
                                                            /// bias gradients associated with the deep net.

   /*! Update the weights, given the current weight gradients. */
   void UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights, const std::vector<Matrix_t> &weightGradients);

   /*! Update the biases, given the current bias gradients. */
   void UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases, const std::vector<Matrix_t> &biasGradients);

public:
   /*! Constructor. */
   TAdam(DeepNet_t &deepNet, Scalar_t learningRate = 0.001, Scalar_t beta1 = 0.9, Scalar_t beta2 = 0.999,
         Scalar_t epsilon = 1e-7);

   /*! Destructor. */
   ~TAdam() = default;

   /*! Getters */
   Scalar_t GetBeta1() const { return fBeta1; }
   Scalar_t GetBeta2() const { return fBeta2; }
   Scalar_t GetEpsilon() const { return fEpsilon; }

   std::vector<std::vector<Matrix_t>> &GetFirstMomentWeights() { return fFirstMomentWeights; }
   std::vector<Matrix_t> &GetFirstMomentWeightsAt(size_t i) { return fFirstMomentWeights[i]; }

   std::vector<std::vector<Matrix_t>> &GetFirstMomentBiases() { return fFirstMomentBiases; }
   std::vector<Matrix_t> &GetFirstMomentBiasesAt(size_t i) { return fFirstMomentBiases[i]; }

   std::vector<std::vector<Matrix_t>> &GetSecondMomentWeights() { return fSecondMomentWeights; }
   std::vector<Matrix_t> &GetSecondMomentWeightsAt(size_t i) { return fSecondMomentWeights[i]; }

   std::vector<std::vector<Matrix_t>> &GetSecondMomentBiases() { return fSecondMomentBiases; }
   std::vector<Matrix_t> &GetSecondMomentBiasesAt(size_t i) { return fSecondMomentBiases[i]; }
};

//
//
//  The Adam Optimizer Class - Implementation
//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
TAdam<Architecture_t, Layer_t, DeepNet_t>::TAdam(DeepNet_t &deepNet, Scalar_t learningRate, Scalar_t beta1,
                                                 Scalar_t beta2, Scalar_t epsilon)
   : VOptimizer<Architecture_t, Layer_t, DeepNet_t>(learningRate, deepNet), fBeta1(beta1), fBeta2(beta2),
     fEpsilon(epsilon)
{
   std::vector<Layer_t *> &layers = deepNet.GetLayers();
   const size_t layersNSlices = layers.size();
   fFirstMomentWeights.resize(layersNSlices);
   fFirstMomentBiases.resize(layersNSlices);
   fSecondMomentWeights.resize(layersNSlices);
   fSecondMomentBiases.resize(layersNSlices);


   for (size_t i = 0; i < layersNSlices; i++) {

      Architecture_t::CreateWeightTensors( fFirstMomentWeights[i], layers[i]->GetWeights());
      Architecture_t::CreateWeightTensors( fSecondMomentWeights[i], layers[i]->GetWeights());

      const size_t weightsNSlices = (layers[i]->GetWeights()).size();

      for (size_t j = 0; j < weightsNSlices; j++) {
         initialize<Architecture_t>(fFirstMomentWeights[i][j], EInitialization::kZero);
         initialize<Architecture_t>(fSecondMomentWeights[i][j], EInitialization::kZero);
      }

      const size_t biasesNSlices = (layers[i]->GetBiases()).size();

      Architecture_t::CreateWeightTensors( fFirstMomentBiases[i], layers[i]->GetBiases());
      Architecture_t::CreateWeightTensors( fSecondMomentBiases[i], layers[i]->GetBiases());

      for (size_t j = 0; j < biasesNSlices; j++) {
         initialize<Architecture_t>(fFirstMomentBiases[i][j], EInitialization::kZero);
         initialize<Architecture_t>(fSecondMomentBiases[i][j], EInitialization::kZero);
      }
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TAdam<Architecture_t, Layer_t, DeepNet_t>::UpdateWeights(size_t layerIndex, std::vector<Matrix_t> &weights,
                                                              const std::vector<Matrix_t> &weightGradients) -> void
{
   // update of weights using Adam algorithm
   // we use the formulation defined before section 2.1 in the original paper
   // 'Adam: A method for stochastic optimization, D. Kingma, J. Ba, see https://arxiv.org/abs/1412.6980

   std::vector<Matrix_t> &currentLayerFirstMomentWeights = this->GetFirstMomentWeightsAt(layerIndex);
   std::vector<Matrix_t> &currentLayerSecondMomentWeights = this->GetSecondMomentWeightsAt(layerIndex);

   // alpha = learningRate * sqrt(1 - beta2^t) / (1-beta1^t)
   Scalar_t alpha = (this->GetLearningRate()) * (sqrt(1 - pow(this->GetBeta2(), this->GetGlobalStep()))) /
                    (1 - pow(this->GetBeta1(), this->GetGlobalStep()));

   /// Adam update of first and second momentum of the weights
   for (size_t i = 0; i < weights.size(); i++) {
      // Mt = beta1 * Mt-1 + (1-beta1) * WeightGradients
      Architecture_t::AdamUpdateFirstMom(currentLayerFirstMomentWeights[i], weightGradients[i], this->GetBeta1() );
      // Vt = beta2 * Vt-1 + (1-beta2) * WeightGradients^2
      Architecture_t::AdamUpdateSecondMom(currentLayerSecondMomentWeights[i], weightGradients[i], this->GetBeta2() );
      // Weight = Weight - alpha * Mt / (sqrt(Vt) + epsilon)
      Architecture_t::AdamUpdate(weights[i], currentLayerFirstMomentWeights[i], currentLayerSecondMomentWeights[i],
                                 alpha, this->GetEpsilon() );
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t, typename Layer_t, typename DeepNet_t>
auto TAdam<Architecture_t, Layer_t, DeepNet_t>::UpdateBiases(size_t layerIndex, std::vector<Matrix_t> &biases,
                                                             const std::vector<Matrix_t> &biasGradients) -> void
{
   std::vector<Matrix_t> &currentLayerFirstMomentBiases = this->GetFirstMomentBiasesAt(layerIndex);
   std::vector<Matrix_t> &currentLayerSecondMomentBiases = this->GetSecondMomentBiasesAt(layerIndex);

   // alpha = learningRate * sqrt(1 - beta2^t) / (1-beta1^t)
   Scalar_t alpha = (this->GetLearningRate()) * (sqrt(1 - pow(this->GetBeta2(), this->GetGlobalStep()))) /
                    (1 - pow(this->GetBeta1(), this->GetGlobalStep()));

   // updating of the biases.
   for (size_t i = 0; i < biases.size(); i++) {
      // Mt = beta1 * Mt-1 + (1-beta1) * BiasGradients
      Architecture_t::AdamUpdateFirstMom(currentLayerFirstMomentBiases[i], biasGradients[i], this->GetBeta1() );
      // Vt = beta2 * Vt-1 + (1-beta2) * BiasGradients^2
      Architecture_t::AdamUpdateSecondMom(currentLayerSecondMomentBiases[i], biasGradients[i], this->GetBeta2() );
      // theta = theta - alpha * Mt / (sqrt(Vt) + epsilon)
      Architecture_t::AdamUpdate(biases[i], currentLayerFirstMomentBiases[i], currentLayerSecondMomentBiases[i],
                                 alpha, this->GetEpsilon() );
   }
}

} // namespace DNN
} // namespace TMVA

#endif
