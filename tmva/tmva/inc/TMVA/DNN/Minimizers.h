// @(#)root/tmva $Id$
// Author: Simon Pfreundschuh 21/06/16

/*************************************************************************
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef TMVA_DNN_MINIMIZERS
#define TMVA_DNN_MINIMIZERS

#include "DataLoader.h"
#include "Functions.h"

#include <limits>
#include <vector>

namespace TMVA {
namespace DNN {

//______________________________________________________________________________
//
// Generic Gradient Descent Class
//______________________________________________________________________________
//

/*** \class TGradientDescent
*
*   Generic implementation of gradient descent minimization.
*
*   The TGradientDescent class implements an architecture and input data
*   independent implementation of the gradient descent minimization algorithm.
*
*   Provides Train(...) and TrainMomentum(...) functions that perform a complete
*   training of a neural network. Those are mainly used for testing since for
*   production a more fine grained control of the training process is desirable.
*   This is provided by the Step(...), StepMomentum(...) and StepNesterov(...)
*   functions that perform a single minimization step.
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
*
*/
template<typename Architecture_t>
class TGradientDescent
{
public:
   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

private:
   size_t   fBatchSize; ///< Batch size to use for the training.
   size_t   fStepCount; ///< Number of steps performed in the current training session
   size_t   fConvergenceSteps; ///< Number of training epochs without considerable
   ///< decrease in the test error for convergence.
   size_t   fConvergenceCount; ///< Current number of training epochs without
   ///< considerable decrease in the test error.
   size_t   fTestInterval; ///< Interval for the computation of the test error.
   Scalar_t fTrainingError;///< Holds the most recently computed training loss.
   Scalar_t fTestError;    ///< Holds the most recently computed test loss.
   Scalar_t fLearningRate; ///< Learning rate \f$\alpha\f$
   Scalar_t fMinimumError; ///< The minimum loss achieved on the training set
   ///< during the current traning session.

public:
   TGradientDescent();
   TGradientDescent(Scalar_t learningRate,
                    size_t   convergenceSteps,
                    size_t   testInterval);

   /** Reset minimizer object to default state. */
   void Reset()
   {
      fMinimumError = std::numeric_limits<Scalar_t>::infinity();
      fConvergenceCount = 0;
      fStepCount = 0;
   };

   /** Train the given net using the given training input data (events), training
       output data (labels), test input data (events), test output data (labels). */
   template <typename Data_t, typename Net_t>
   Scalar_t Train(const Data_t & TrainingDataIn, size_t nTrainingSamples,
                  const Data_t & TestDataIn, size_t nTestSamples,
                  Net_t & net, size_t nThreads = 1);

   /** Same as Train(...) but uses the given momentum.*/
   template <typename Data_t, typename Net_t>
   Scalar_t TrainMomentum(const Data_t & TrainingDataIn, size_t nTrainingSamples,
                          const Data_t & TestDataIn, size_t nTestSamples,
                          Net_t & net, Scalar_t momentum, size_t nThreads = 1);

   /** Perform a single optimization step on a given batch. Propagates the input
       matrix foward through the net, evaluates the loss and propagates the gradients
       backward through the net. The computed gradients are scaled by the learning
       rate \f$\alpha\f$ and subtracted from the weights and bias values of each
       layer. */
   template <typename Net_t>
   void Step(Net_t &net, Matrix_t &input, const Matrix_t &output, const Matrix_t &weights);

   /** Same as Step(...) but also evaluate the loss on the given training data.
    *  Note that this requires synchronization between host and device. */
   template <typename Net_t>
   Scalar_t StepLoss(Net_t &net, Matrix_t &input, const Matrix_t &output, const Matrix_t &weights);

   /** Perform multiple optimization steps simultaneously. Performs the
    *  backprop algorithm on the input batches given in \p batches on
    *  the neural networks given in \p nets. The forward and backward propagation
    *  steps are executed in an interleaving manner in order to exploit potential
    *  batch-level parallelism for asynchronous device calls.
    */
   template <typename Net_t>
   void Step(Net_t &master,
             std::vector<Net_t> &nets,
             std::vector<TBatch<Architecture_t>> &batches);

   /** Same as the Step(...) method for multiple batches but uses momentum. */
   template <typename Net_t>
   void StepMomentum(Net_t &master,
                     std::vector<Net_t> &nets,
                     std::vector<TBatch<Architecture_t>> &batches,
                     Scalar_t momentum);
   template <typename Net_t>

   /** Same as the Step(...) method for multiple batches but uses Nesterov
    *  momentum. */
   void StepNesterov(Net_t &master,
                     std::vector<Net_t> &nets,
                     std::vector<TBatch<Architecture_t>> &batches,
                     Scalar_t momentum);

   /** Does not evaluate the loss and therefore not trigger a possible synchronization
    *  with the device. Trains the weights of each layer, but only the bias terms of
    *  the first layer for compatibility with the previous implementation. */
   template <typename Net_t>
   void StepReducedWeights(Net_t &net, Matrix_t &input, const Matrix_t &output);

   /** Similar to StepReducedWeights(...) but also evaluates the loss. May trigger
    * synchronization with the device. */
   template <typename Net_t>
   Scalar_t StepReducedWeightsLoss(Net_t &net, Matrix_t &input, const Matrix_t &output, const Matrix_t &weights);
   /** Increases the minimization step counter by the test error evaluation
    *  period and uses the current internal value of the test error to
    *  determine if the minimization has converged. */
   bool HasConverged();
   /** Increases the minimization step counter by the test error evaluation
    *  period and uses the provided test error value to determine if the
    *  minimization has converged. */
   bool HasConverged(Scalar_t testError);

   size_t   GetConvergenceCount() const {return fConvergenceCount;}
   size_t   GetConvergenceSteps() const {return fConvergenceSteps;}
   Scalar_t GetTrainingError() const {return fTrainingError;}
   Scalar_t GetTestError() const     {return fTestError;}
   size_t   GetTestInterval() const  {return fTestInterval;}

   void SetConvergenceSteps(size_t steps) {fConvergenceSteps = steps;}
   void SetTestInterval(size_t interval)  {fTestInterval = interval;}
   void SetLearningRate(Scalar_t rate)    {fLearningRate = rate;}
   void SetBatchSize(Scalar_t rate)       {fBatchSize    = rate;}
};

//
// Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TGradientDescent<Architecture_t>::TGradientDescent()
   : fBatchSize(0), fStepCount(0), fConvergenceSteps(0), fConvergenceCount(0), fTestInterval(0), fLearningRate(0),
     fMinimumError(std::numeric_limits<Scalar_t>::infinity())
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t>
TGradientDescent<Architecture_t>::TGradientDescent(Scalar_t learningRate, size_t convergenceSteps, size_t testInterval)
   : fBatchSize(0), fStepCount(0), fConvergenceSteps(convergenceSteps), fConvergenceCount(0),
     fTestInterval(testInterval), fLearningRate(learningRate), fMinimumError(std::numeric_limits<Scalar_t>::infinity())
{
   // Nothing to do here.
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Data_t, typename Net_t>
    auto TGradientDescent<Architecture_t>::Train(const Data_t & trainingData,
                                                 size_t nTrainingSamples,
                                                 const Data_t & testData,
                                                 size_t nTestSamples,
                                                 Net_t & net,
                                                 size_t nThreads)
   -> Scalar_t
{
   Reset();

   // Prepare training data.
   TDataLoader<Data_t, Architecture_t> trainLoader(trainingData, nTrainingSamples,
                                                   net.GetBatchSize(),
                                                   net.GetInputWidth(),
                                                   net.GetOutputWidth(), nThreads);
   auto testNet = net.CreateClone(nTestSamples);
   TDataLoader<Data_t, Architecture_t> testLoader(testData, nTestSamples,
                                                  testNet.GetBatchSize(),
                                                  testNet.GetInputWidth(),
                                                  net.GetOutputWidth());
   std::vector<Net_t> nets{};
   nets.reserve(nThreads);
   for (size_t i = 0; i < nThreads; i++) {
       nets.push_back(net);
       for (size_t j = 0; j < net.GetDepth(); j++)
       {
           auto &masterLayer = net.GetLayer(j);
           auto &layer = nets.back().GetLayer(j);
           Architecture_t::Copy(layer.GetWeights(),
                                masterLayer.GetWeights());
           Architecture_t::Copy(layer.GetBiases(),
                                masterLayer.GetBiases());
       }
   }

   size_t batchesInEpoch = nTrainingSamples / net.GetBatchSize();
   std::vector<TBatch<Architecture_t>> batches{};
   batches.reserve(nThreads);

   do {
      for (fStepCount = 0; fStepCount < fTestInterval; fStepCount++) {
         trainLoader.Shuffle();
         for (size_t i = 0; i < batchesInEpoch; i += nThreads) {
            batches.clear();
            for (size_t j = 0; j < nThreads; j++) batches.push_back(trainLoader.GetBatch());
            Step(net, nets, batches);
         }
      }

      auto b = *testLoader.begin();
      auto inputMatrix = b.GetInput();
      auto outputMatrix = b.GetOutput();
      auto weightMatrix = b.GetWeights();
      fTestError = testNet.Loss(inputMatrix, outputMatrix, weightMatrix);

   } while (!HasConverged());

   return fMinimumError;
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Data_t, typename Net_t>
auto TGradientDescent<Architecture_t>::TrainMomentum(const Data_t & trainingData,
                                                     size_t nTrainingSamples,
                                                     const Data_t & testData,
                                                     size_t nTestSamples,
                                                     Net_t & net,
                                                     Scalar_t momentum,
                                                     size_t nThreads)
   -> Scalar_t
{
   Reset();

   // Prepare training data.
   TDataLoader<Data_t, Architecture_t> trainLoader(trainingData, nTrainingSamples,
                                                   net.GetBatchSize(),
                                                   net.GetInputWidth(),
                                                   net.GetOutputWidth(), nThreads);
   auto testNet = net.CreateClone(net.GetBatchSize());
   TDataLoader<Data_t, Architecture_t> testLoader(testData, nTestSamples,
                                                  testNet.GetBatchSize(),
                                                  testNet.GetInputWidth(),
                                                  net.GetOutputWidth());

   net.InitializeGradients();
   std::vector<Net_t> nets{};
   nets.reserve(nThreads);
   for (size_t i = 0; i < nThreads; i++) {
       nets.push_back(net);
       for (size_t j = 0; j < net.GetDepth(); j++)
       {
           auto &masterLayer = net.GetLayer(j);
           auto &layer = nets.back().GetLayer(j);
           Architecture_t::Copy(layer.GetWeights(),
                                masterLayer.GetWeights());
           Architecture_t::Copy(layer.GetBiases(),
                                masterLayer.GetBiases());
       }
   }

   size_t batchesInEpoch = nTrainingSamples / net.GetBatchSize();
   std::vector<TBatch<Architecture_t>> batches{};
   batches.reserve(nThreads);

   do {
      for (fStepCount = 0; fStepCount < fTestInterval; fStepCount++) {
         trainLoader.Shuffle();
         for (size_t i = 0; i < batchesInEpoch; i += nThreads) {
            batches.clear();
            for (size_t j = 0; j < nThreads; j++) batches.push_back(trainLoader.GetBatch());
            if (momentum != 0.0) {
               StepMomentum(net, nets, batches, momentum);
            } else {
               Step(net, nets, batches);
            }
         }
      }

      fTestError = 0.0;
      for (size_t i = 0; i < batchesInEpoch; i++) {
         auto b = testLoader.GetBatch();
         auto inputMatrix = b.GetInput();
         auto outputMatrix = b.GetOutput();
         auto weightMatrix = b.GetWeights();
         fTestError += testNet.Loss(inputMatrix, outputMatrix, weightMatrix);
      }
      fTestError /= (Double_t)batchesInEpoch;
   } while (!HasConverged());
   return fMinimumError;
}

//______________________________________________________________________________
template <typename Architecture_t>
template <typename Net_t>
void inline TGradientDescent<Architecture_t>::Step(Net_t &net, Matrix_t &input, const Matrix_t &output,
                                                   const Matrix_t &weights)
{
   net.Forward(input, true);
   net.Backward(input, output, weights);

   for (size_t i = 0; i < net.GetDepth(); i++)
   {
      auto &layer = net.GetLayer(i);
      Architecture_t::ScaleAdd(layer.GetWeights(),
                               layer.GetWeightGradients(),
                               -fLearningRate);
      Architecture_t::ScaleAdd(layer.GetBiases(),
                               layer.GetBiasGradients(),
                               -fLearningRate);
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
template <typename Net_t>
auto inline TGradientDescent<Architecture_t>::StepLoss(Net_t &net, Matrix_t &input, const Matrix_t &output,
                                                       const Matrix_t &weights) -> Scalar_t
{
   Scalar_t loss = net.Loss(input, output, weights);
   net.Backward(input, output);

   for (size_t i = 0; i < net.GetDepth(); i++)
   {
      auto &layer = net.GetLayer(i);
      Architecture_t::ScaleAdd(layer.GetWeights(),
                               layer.GetWeightGradients(),
                               -fLearningRate);
      Architecture_t::ScaleAdd(layer.GetBiases(),
                               layer.GetBiasGradients(),
                               -fLearningRate);
   }
   return loss;
}

//______________________________________________________________________________
template<typename Architecture_t>
    template <typename Net_t>
    void inline TGradientDescent<Architecture_t>::Step(
        Net_t & master,
        std::vector<Net_t> & nets,
        std::vector<TBatch<Architecture_t>> & batches)
{
   typename Architecture_t::Matrix_t dummy(0,0);
   size_t depth = master.GetDepth();

   // Forward
   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Forward(batches[j].GetInput(), true);
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput(), true);
      }
   }
   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(nets[j].GetLayer(depth - 1).GetActivationGradients(), nets[j].GetLossFunction(),
                                        batches[j].GetOutput(), nets[j].GetLayer(depth - 1).GetOutput(),
                                        batches[j].GetWeights());
   }
   // Backward
   for (size_t i = depth - 1; i > 0; i--)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Backward(nets[j].GetLayer(i-1).GetActivationGradients(),
                                      nets[j].GetLayer(i-1).GetOutput(),
                                      nets[j].GetRegularization(),
                                      nets[j].GetWeightDecay());
      }
   }
   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Backward(dummy,
                                   batches[j].GetInput(),
                                   nets[j].GetRegularization(),
                                   nets[j].GetWeightDecay());
   }

   for (size_t j = 0; j < nets.size(); j++) {
      for (size_t i = 0; i < depth; i++)
      {
         auto &masterLayer = master.GetLayer(i);
         auto &layer       = nets[j].GetLayer(i);
         Architecture_t::ScaleAdd(masterLayer.GetWeights(),
                                  layer.GetWeightGradients(),
                                  -fLearningRate);
         Architecture_t::Copy(layer.GetWeights(),
                              masterLayer.GetWeights());
         Architecture_t::ScaleAdd(masterLayer.GetBiases(),
                                  layer.GetBiasGradients(),
                                  -fLearningRate);
         Architecture_t::Copy(layer.GetBiases(),
                              masterLayer.GetBiases());
      }
   }
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Net_t>
void inline TGradientDescent<Architecture_t>::StepMomentum(
        Net_t & master,
        std::vector<Net_t> & nets,
        std::vector<TBatch<Architecture_t>> & batches,
        Scalar_t momentum)
{
   typename Architecture_t::Matrix_t dummy(0,0);
   size_t depth = master.GetDepth();

   // Forward
   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Forward(batches[j].GetInput(), true);
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput(), true);
      }
   }
   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(nets[j].GetLayer(depth - 1).GetActivationGradients(), nets[j].GetLossFunction(),
                                        batches[j].GetOutput(), nets[j].GetLayer(depth - 1).GetOutput(),
                                        batches[j].GetWeights());
   }
   // Backward
   for (size_t i = depth - 1; i > 0; i--)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Backward(nets[j].GetLayer(i-1).GetActivationGradients(),
                                      nets[j].GetLayer(i-1).GetOutput(),
                                      nets[j].GetRegularization(),
                                      nets[j].GetWeightDecay());
         Architecture_t::ScaleAdd(master.GetLayer(i).GetWeightGradients(),
                                  nets[j].GetLayer(i).GetWeightGradients(),
                                  - fLearningRate / momentum);
         Architecture_t::ScaleAdd(master.GetLayer(i).GetBiasGradients(),
                                  nets[j].GetLayer(i).GetBiasGradients(),
                                  - fLearningRate / momentum);
      }
      Architecture_t::ScaleAdd(master.GetLayer(i).GetWeightGradients(),
                               master.GetLayer(i).GetWeightGradients(),
                               momentum - 1.0);
      Architecture_t::ScaleAdd(master.GetLayer(i).GetBiasGradients(),
                               master.GetLayer(i).GetBiasGradients(),
                               momentum - 1.0);
   }
   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Backward(dummy,
                                   batches[j].GetInput(),
                                   nets[j].GetRegularization(),
                                   nets[j].GetWeightDecay());
      Architecture_t::ScaleAdd(master.GetLayer(0).GetWeightGradients(),
                               nets[j].GetLayer(0).GetWeightGradients(),
                               - fLearningRate / momentum);
      Architecture_t::ScaleAdd(master.GetLayer(0).GetBiasGradients(),
                               nets[j].GetLayer(0).GetBiasGradients(),
                               - fLearningRate / momentum);
   }

   Architecture_t::ScaleAdd(master.GetLayer(0).GetWeightGradients(),
                            master.GetLayer(0).GetWeightGradients(),
                            momentum - 1.0);
   Architecture_t::ScaleAdd(master.GetLayer(0).GetBiasGradients(),
                            master.GetLayer(0).GetBiasGradients(),
                            momentum - 1.0);

   for (size_t i = 0; i < depth; i++)
   {
       auto &masterLayer = master.GetLayer(i);
       Architecture_t::ScaleAdd(masterLayer.GetWeights(),
                                masterLayer.GetWeightGradients(),
                                1.0);
       Architecture_t::ScaleAdd(masterLayer.GetBiases(),
                                masterLayer.GetBiasGradients(),
                                1.0);
       for (size_t j = 0; j < nets.size(); j++) {
         auto &layer       = nets[j].GetLayer(i);
         Architecture_t::Copy(layer.GetWeights(),
                              masterLayer.GetWeights());
         Architecture_t::Copy(layer.GetBiases(),
                              masterLayer.GetBiases());
       }
   }
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Net_t>
void inline TGradientDescent<Architecture_t>::StepNesterov(
        Net_t & master,
        std::vector<Net_t> & nets,
        std::vector<TBatch<Architecture_t>> & batches,
        Scalar_t momentum)
{
   typename Architecture_t::Matrix_t dummy(0,0);
   size_t depth = master.GetDepth();

   // Forward
   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Forward(batches[j].GetInput(), true);
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput(), true);
      }
   }

   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(nets[j].GetLayer(depth - 1).GetActivationGradients(), nets[j].GetLossFunction(),
                                        batches[j].GetOutput(), nets[j].GetLayer(depth - 1).GetOutput(),
                                        batches[j].GetWeights());
   }

   // Backward
   for (size_t i = depth - 1; i > 0; i--)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Backward(nets[j].GetLayer(i-1).GetActivationGradients(),
                                      nets[j].GetLayer(i-1).GetOutput(),
                                      nets[j].GetRegularization(),
                                      nets[j].GetWeightDecay());
      }
   }

   for (size_t j = 0; j < nets.size(); j++) {
      nets[j].GetLayer(0).Backward(dummy,
                                   batches[j].GetInput(),
                                   nets[j].GetRegularization(),
                                   nets[j].GetWeightDecay());
   }

   for (size_t i = 0; i < depth; i++)
   {
      auto &masterLayer = master.GetLayer(i);
      for (size_t j = 0; j < nets.size(); j++) {
         auto &layer       = nets[j].GetLayer(i);
         Architecture_t::Copy(layer.GetWeights(),
                              masterLayer.GetWeights());
         Architecture_t::Copy(layer.GetBiases(),
                              masterLayer.GetBiases());
         Architecture_t::ScaleAdd(layer.GetWeights(),
                                  masterLayer.GetWeightGradients(),
                                  1.0);
         Architecture_t::ScaleAdd(layer.GetBiases(),
                                  masterLayer.GetBiasGradients(),
                                  1.0);
      }
      for (size_t j = 0; j < nets.size(); j++) {
         auto &layer       = nets[j].GetLayer(i);
         Architecture_t::ScaleAdd(masterLayer.GetWeightGradients(),
                                  layer.GetWeightGradients(),
                                  - fLearningRate / momentum);
         Architecture_t::ScaleAdd(masterLayer.GetBiasGradients(),
                                  layer.GetBiasGradients(),
                                  - fLearningRate / momentum);
      }
      Architecture_t::ScaleAdd(masterLayer.GetWeightGradients(),
                               masterLayer.GetWeightGradients(),
                               momentum - 1.0);
      Architecture_t::ScaleAdd(masterLayer.GetBiasGradients(),
                               masterLayer.GetBiasGradients(),
                               momentum - 1.0);
      Architecture_t::ScaleAdd(masterLayer.GetWeights(),
                               masterLayer.GetWeightGradients(),
                               1.0);
      Architecture_t::ScaleAdd(masterLayer.GetBiases(),
                               masterLayer.GetBiasGradients(),
                               1.0);
   }
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Net_t>
void inline TGradientDescent<Architecture_t>::StepReducedWeights(
    Net_t & net,
    Matrix_t &input,
    const Matrix_t &output)
{
   net.Forward(input, true);
   net.Backward(input, output);

   for (size_t i = 0; i < net.GetDepth(); i++)
   {
      auto &layer = net.GetLayer(i);
      Architecture_t::ScaleAdd(layer.GetWeights(),
                               layer.GetWeightGradients(),
                               -fLearningRate);
      if (i == 0) {
         Architecture_t::ScaleAdd(layer.GetBiases(),
                                  layer.GetBiasGradients(),
                                  -fLearningRate);
      }
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
template <typename Net_t>
auto inline TGradientDescent<Architecture_t>::StepReducedWeightsLoss(Net_t &net, Matrix_t &input,
                                                                     const Matrix_t &output, const Matrix_t &weights)
   -> Scalar_t
{
   Scalar_t loss = net.Loss(input, output);
   fTrainingError = loss;
   net.Backward(input, output, weights);

   for (size_t i = 0; i < net.GetDepth(); i++)
   {
      auto &layer = net.GetLayer(i);
      Architecture_t::ScaleAdd(layer.GetWeights(),
                               layer.GetWeightGradients(),
                               -fLearningRate);
      if (i == 0) {
         Architecture_t::ScaleAdd(layer.GetBiases(),
                                  layer.GetBiasGradients(),
                                  -fLearningRate);
      }
   }
   return loss;
}

//______________________________________________________________________________
template<typename Architecture_t>
bool inline TGradientDescent<Architecture_t>::HasConverged()
{
   if (fTestError < fMinimumError * 0.999) {
      fConvergenceCount = 0;
      fMinimumError     = fTestError;
   } else {
      fConvergenceCount++;
   }

   return (fConvergenceCount >= fConvergenceSteps);
}

//______________________________________________________________________________
template<typename Architecture_t>
bool inline TGradientDescent<Architecture_t>::HasConverged(Scalar_t testError)
{
   fTestError = testError;
   if (fTestError < fMinimumError * 0.999) {
      fConvergenceCount = 0;
      fMinimumError     = fTestError;
   } else {
      fConvergenceCount += fTestInterval;
   }
   return (fConvergenceCount >= fConvergenceSteps);
}
} // namespace DNN
} // namespace TMVA

#endif
