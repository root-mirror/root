// @(#)root/tmva/tmva/cnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TDLGradientDescent                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Deel Learning Minimizers                                                  *
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

#ifndef TMVA_DNN_DLMINIMIZERS
#define TMVA_DNN_DLMINIMIZERS

#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

#include <limits>
#include <vector>

namespace TMVA {
namespace DNN {

/*** \class TDLGradientDescent
 *
 *   Generic implementation of gradient descent minimization for the
 *   deep learning neural nets.
 *
 *   The TDLGradientDescent class implements an architecture, input data and
 *   deep learning neural network type independent implementation of the gradient
 *   descent minimization algorithm.
 *
*    This is provided by the Step(...), StepMomentum(...) and
 *   StepNesterov(...) functions that perform a single minimization step.
 *
 *   The main training characteristics are defined by the provided learning rate,
 *   the test interval, and the convergence steps required for convergence. The
 *   test interval defines how often the error on the validation set is computed,
 *   and the values with which the step counter is increased each time the
 *   HasConverged() member function is called. A convergence step is defined as
 *   a step in which the test error is NOT less than 0.999 times the current
 *   minimal test error that has been reached. If between two subsequent calls
 *   to HasConverged(Double_t) the test error has not been sufficiently reduced
 *   it is assumed that a number of convergence steps equal to the test interval
 *   has been performed.
 */

template <typename Architecture_t>
class TDLGradientDescent {
public:
   using DeepNet_t = TDeepNet<Architecture_t>;
   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

private:
   size_t fBatchSize;        ///< Batch size to use for the training.
   size_t fStepCount;        ///< Number of steps performed in the current training session
   size_t fConvergenceSteps; ///< Number of training epochs without considerable
   ///< decrease in the test error for convergence.
   size_t fConvergenceCount; ///< Current number of training epochs without
   ///< considerable decrease in the test error.
   size_t fTestInterval;    ///< Interval for the computation of the test error.
   Scalar_t fTrainingError; ///< Holds the most recently computed training loss.
   Scalar_t fTestError;     ///< Holds the most recently computed test loss.
   Scalar_t fLearningRate;  ///< Learning rate \f$\alpha\f$
   Scalar_t fMinimumError;  ///< The minimum loss achieved on the training set
   ///< during the current traning session.

public:
   TDLGradientDescent();
   TDLGradientDescent(Scalar_t learningRate, size_t convergenceSteps, size_t testInterval);

   /** Reset minimizer object to default state. */
   void Reset()
   {
      fMinimumError = std::numeric_limits<Scalar_t>::infinity();
      fConvergenceCount = 0;
      fStepCount = 0;
   };

   /** Perform a single optimization step on a given batch. Propagates the input
       matrix foward through the net, evaluates the loss and propagates the gradients
       backward through the net. The computed gradients are scaled by the learning
       rate \f$\alpha\f$ and subtracted from the weights and bias values of each
       layer. */
   void Step(DeepNet_t &deepNet, std::vector<Matrix_t> &input, const Matrix_t &output, const Matrix_t &weights);

   /** Does not evaluate the loss and therefore not trigger a possible synchronization
    *  with the device. Trains the weights of each layer, but only the bias terms of
    *  the first layer for compatibility with the previous implementation. */
   void StepReducedWeights(DeepNet_t &deepNet, std::vector<Matrix_t> &input, const Matrix_t &output,
                           const Matrix_t &weights);

   /** Same as Step(...) but also evaluate the loss on the given training data.
    *  Note that this requires synchronization between host and device. */
   Scalar_t StepLoss(DeepNet_t &deepNet, std::vector<Matrix_t> &input, const Matrix_t &output, const Matrix_t &weights);

   /** Similar to StepReducedWeights(...) but also evaluates the loss. May trigger
     * synchronization with the device. */
   Scalar_t StepReducedWeightsLoss(DeepNet_t &deepNet, std::vector<Matrix_t> &input, const Matrix_t &output,
                                   const Matrix_t &weights);

   /** Perform multiple optimization steps simultaneously. Performs the
    *  backprop algorithm on the input batches given in \p batches on
    *  the neural networks given in \p nets. The forward and backward propagation
    *  steps are executed in an interleaving manner in order to exploit potential
    *  batch-level parallelism for asynchronous device calls.
    */
   void Step(DeepNet_t &master, std::vector<DeepNet_t> &nets, std::vector<TTensorBatch<Architecture_t>> &batches);

   /** Same as the Step(...) method for multiple batches but uses momentum. */
   void StepMomentum(DeepNet_t &master, std::vector<DeepNet_t> &nets,
                     std::vector<TTensorBatch<Architecture_t>> &batches, Scalar_t momentum);

   /** Same as the Step(...) method for multiple batches but uses Nesterov
    *  momentum. */
   void StepNesterov(DeepNet_t &master, std::vector<DeepNet_t> &nets,
                     std::vector<TTensorBatch<Architecture_t>> &batches, Scalar_t momentum);

   /** Increases the minimization step counter by the test error evaluation
    *  period and uses the current internal value of the test error to
    *  determine if the minimization has converged. */
   bool HasConverged();

   /** Increases the minimization step counter by the test error evaluation
    *  period and uses the provided test error value to determine if the
    *  minimization has converged. */
   bool HasConverged(Scalar_t testError);

   /** Getters */
   size_t GetConvergenceCount() const { return fConvergenceCount; }
   size_t GetConvergenceSteps() const { return fConvergenceSteps; }
   Scalar_t GetTrainingError() const { return fTrainingError; }
   Scalar_t GetTestError() const { return fTestError; }
   size_t GetTestInterval() const { return fTestInterval; }

   /** Setters */
   void SetConvergenceSteps(size_t steps) { fConvergenceSteps = steps; }
   void SetTestInterval(size_t interval) { fTestInterval = interval; }
   void SetLearningRate(Scalar_t rate) { fLearningRate = rate; }
   void SetBatchSize(Scalar_t rate) { fBatchSize = rate; }
};

//
// Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TDLGradientDescent<Architecture_t>::TDLGradientDescent()
   : fBatchSize(0), fStepCount(0), fConvergenceSteps(0), fConvergenceCount(0), fTestInterval(0), fLearningRate(0),
     fMinimumError(std::numeric_limits<Scalar_t>::infinity())
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t>
TDLGradientDescent<Architecture_t>::TDLGradientDescent(Scalar_t learningRate, size_t convergenceSteps,
                                                       size_t testInterval)
   : fBatchSize(0), fStepCount(0), fConvergenceSteps(convergenceSteps), fConvergenceCount(0),
     fTestInterval(testInterval), fLearningRate(learningRate), fMinimumError(std::numeric_limits<Scalar_t>::infinity())
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDLGradientDescent<Architecture_t>::Step(DeepNet_t &deepNet, std::vector<Matrix_t> &input, const Matrix_t &output,
                                              const Matrix_t &weights)
{
   // Make forward and backward pass and update the net afterwards
   deepNet.Forward(input, true);
   deepNet.Backward(input, output, weights);
   deepNet.Update(fLearningRate);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDLGradientDescent<Architecture_t>::StepReducedWeights(DeepNet_t &deepNet, std::vector<Matrix_t> &input,
                                                            const Matrix_t &output, const Matrix_t &weights)
{
   // Make forward and backward pass and update the net afterwards
   deepNet.Forward(input, true);
   deepNet.Backward(input, output, weights);

   for (size_t i = 0; i < deepNet.GetDepth(); i++) {
      auto *layer = deepNet.GetLayerAt(i);

      layer->UpdateWeights(layer->GetWeightGradients(), fLearningRate);
      if (i == 0) {
         layer->UpdateBiases(layer->GetBiasGradients(), fLearningRate);
      }
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TDLGradientDescent<Architecture_t>::StepLoss(DeepNet_t &deepNet, std::vector<Matrix_t> &input,
                                                  const Matrix_t &output, const Matrix_t &weights) -> Scalar_t
{
   Scalar_t loss = deepNet.Loss(input, output);
   deepNet.Backward(input, output, weights);
   deepNet.Update(fLearningRate);

   return loss;
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TDLGradientDescent<Architecture_t>::StepReducedWeightsLoss(DeepNet_t &deepNet, std::vector<Matrix_t> &input,
                                                                const Matrix_t &output, const Matrix_t &weights)
   -> Scalar_t
{
   Scalar_t loss = deepNet.Loss(input, output);
   fTrainingError = loss;
   deepNet.Backward(input, output, weights);

   for (size_t i = 0; i < deepNet.GetDepth(); i++) {
      auto *layer = deepNet.GetLayerAt(i);

      layer->UpdateWeights(layer->GetWeightGradients(), fLearningRate);
      if (i == 0) {
         layer->UpdateBiases(layer->GetBiasGradients(), fLearningRate);
      }
   }

   return loss;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDLGradientDescent<Architecture_t>::Step(DeepNet_t &master, std::vector<DeepNet_t> &nets,
                                              std::vector<TTensorBatch<Architecture_t>> &batches)
{

   master.ParallelForward(nets, batches);
   master.ParallelBackward(nets, batches, fLearningRate);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDLGradientDescent<Architecture_t>::StepMomentum(DeepNet_t &master, std::vector<DeepNet_t> &nets,
                                                      std::vector<TTensorBatch<Architecture_t>> &batches,
                                                      Scalar_t momentum)
{
   master.ParallelForward(nets, batches);
   master.ParallelBackwardMomentum(nets, batches, fLearningRate, momentum);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDLGradientDescent<Architecture_t>::StepNesterov(DeepNet_t &master, std::vector<DeepNet_t> &nets,
                                                      std::vector<TTensorBatch<Architecture_t>> &batches,
                                                      Scalar_t momentum)
{
   master.ParallelForward(nets, batches);
   master.ParallelBackwardNestorov(nets, batches, fLearningRate, momentum);
}

//______________________________________________________________________________
template <typename Architecture_t>
bool TDLGradientDescent<Architecture_t>::HasConverged()
{
   if (fTestError < fMinimumError * 0.999) {
      fConvergenceCount = 0;
      fMinimumError = fTestError;
   } else {
      fConvergenceCount++;
   }

   return (fConvergenceCount >= fConvergenceSteps);
}

//______________________________________________________________________________
template <typename Architecture_t>
bool TDLGradientDescent<Architecture_t>::HasConverged(Scalar_t testError)
{
   fTestError = testError;
   if (fTestError < fMinimumError * 0.999) {
      fConvergenceCount = 0;
      fMinimumError = fTestError;
   } else {
      fConvergenceCount += fTestInterval;
   }
   return (fConvergenceCount >= fConvergenceSteps);
}

} // namespace DNN
} // namespace TMVA

#endif
