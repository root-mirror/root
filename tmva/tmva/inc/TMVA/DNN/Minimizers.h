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
#include <chrono>
#include "tbb/tbb.h"

namespace TMVA {
namespace DNN {

//______________________________________________________________________________
//
// Generic Gradient Descent Class
//______________________________________________________________________________
//

/** \class TGradientDescent

    Generic implementation of gradient descent minimization.

    The TGradientDescent class implements an architecture and input data
    independent implementation of the gradient descent minimization algorithm.

    The Train(...) method trains a given neural network using the provided training
    and test (validation) data. The interface between input data and the matrix
    representation required by the net is given by the prepare_batches(...) and
    prepare_test_data(...) methods as well as the TBatch class of the
    architecture-specific back end.

    The prepare_batches(...) is expected to generate an iterable of batches. On
    each of these batches, the Step(...) routine of the minimizer is called, which
    for the gradient descent method just adds the corresponding gradients scaled
    by \f$-\alpha\f$ to the weights and biases of each layer. Here \f$\alpha\f$ is
    the learning rate of the gradient descent method.

    The prepare_test_data(...) routine should return a batch representing
    the test data, which is used to evaluate the performance of the net
    every testInterval steps.

    \tparam Architecture_t Type representing which implementation of the low-level
    interface to use.
 */
template<typename Architecture_t>
class TGradientDescent
{
public:
    using Scalar_t = typename Architecture_t::Scalar_t;
    using Matrix_t = typename Architecture_t::Matrix_t;

private:
    size_t   fBatchSize; ///< Batch size to use for the training.
    size_t   fStepCount; ///< Number of steps performed in the current training sessiong.
    size_t   fConvergenceSteps; ///< Number of training epochs without considerable decrease in the test error for convergence.
    size_t   fConvergenceCount; ///< Current number of training epochs without considerable decrease in the test error.
    size_t   fTestInterval; ///< Interval for the computation of the test error.
    Scalar_t fTrainingError;///< Holds the most recently computed training loss.
    Scalar_t fTestError;    ///< Holds the most recently computed test loss.
    Scalar_t fLearningRate; ///< Learning rate \f$\alpha\f$
    Scalar_t fMinimumError; ///< The minimum loss achieved on the training set during the current traning session.

public:
    TGradientDescent();
    TGradientDescent(Scalar_t fBatchSize,
                     Scalar_t learningRate,
                     size_t convergenceSteps,
                     size_t testInterval);
    /*! Reset minimizer object to initial state. Does nothing for this minimizer. */
    void Reset() {};
    /*! Train the given net using the given training input data (events), training
      output data (labels), test input data (events), test output data (labels). */
    template <typename Data_t, typename Net_t>
    Scalar_t Train(const Data_t & TrainingDataIn, size_t nTrainingSamples,
                   const Data_t & TestDataIn, size_t nTestSamples,
                   Net_t & net, size_t nThreads = 1);
    template <typename Data_t, typename Net_t>
    Scalar_t TrainTBB(const Data_t & TrainingDataIn, size_t nTrainingSamples,
                      const Data_t & TestDataIn, size_t nTestSamples,
                      Net_t & net, size_t nThreads = 1);
    template <typename Data_t, typename Net_t>
    Scalar_t TrainMomentum(const Data_t & TrainingDataIn, size_t nTrainingSamples,
                           const Data_t & TestDataIn, size_t nTestSamples,
                           Net_t & net, Scalar_t momentum, size_t nThreads = 1);
    /*! Perform a single optimization step on a given batch. Propagates the input
      matrix foward through the net, evaluates the loss and propagates the gradients
      backward through the net. The computed gradients are scaled by the learning
      rate \f$\alpha\f$ and subtracted from the weights and bias values of each
      layer. */
    template <typename Net_t>
    void Step(Net_t &net, Matrix_t &input, const Matrix_t &output);
    template <typename Net_t>
    void Step(Net_t &master,
              std::vector<Net_t> &nets,
              std::vector<TBatch<Architecture_t>> &batches);
    template <typename Data_t, typename Net_t>
    void StepTBB(Net_t &master,
                 std::vector<Net_t> &nets,
                 std::vector<TDataLoader<Data_t,Architecture_t>> &loaders);
    template <typename Net_t>
    void StepMomentum(Net_t &master,
                      std::vector<Net_t> &nets,
                      std::vector<TBatch<Architecture_t>> &batches,
                      Scalar_t momentum);
    template <typename Net_t>
    void StepNesterov(Net_t &master,
                              std::vector<Net_t> &nets,
                              std::vector<TBatch<Architecture_t>> &batches,
                              Scalar_t momentum);
    /** Does not evaluate the loss and therefore not trigger a possible synchronization
     *  with the device. Trains the weights of each layer, but only the bias terms of
     *  the first layer for compatibility with the previous implementation. */
    template <typename Net_t>
    void StepReducedWeights(Net_t &net,
                            Matrix_t &input,
                            const Matrix_t &output);
    /** Similar to StepReducedWeights(...) but also evaluates the loss. May trigger
     * synchronization with the device. */
    template <typename Net_t>
    Scalar_t StepReducedWeightsLoss(Net_t &net,
                                    Matrix_t &input,
                                    const Matrix_t &output);
    template <typename Net_t>
    inline void TestError(Net_t &net,
                          Matrix_t &input,
                          const Matrix_t &output);
    bool HasConverged();

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

//______________________________________________________________________________
//
// Implementation
//______________________________________________________________________________
template<typename Architecture_t>
    TGradientDescent<Architecture_t>::TGradientDescent()
   : fBatchSize(0), fStepCount(0), fConvergenceSteps(0),
     fConvergenceCount(0), fTestInterval(0), fLearningRate(0),
     fMinimumError(1e100)
{
   // Nothing to do here.
}
//______________________________________________________________________________
template<typename Architecture_t>
TGradientDescent<Architecture_t>::TGradientDescent(Scalar_t batchSize,
                                                   Scalar_t learningRate,
                                                   size_t convergenceSteps,
                                                   size_t testInterval)
   : fBatchSize(batchSize), fStepCount(0), fConvergenceSteps(convergenceSteps),
     fConvergenceCount(0), fTestInterval(testInterval), fLearningRate(learningRate),
     fMinimumError(1e100)
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
   // Reset iteration state.
   fMinimumError = 1e100;
   fConvergenceCount = 0;
   fStepCount = 0;

   // Prepare training data.
   bool converged = false;

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

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   while (!converged)
   {
      fStepCount++;

      size_t netIndex = 0;
      std::vector<TBatch<Architecture_t>> batches{};
      for (size_t i = 0; i < nTrainingSamples / net.GetBatchSize(); i += nThreads) {
         batches.clear();
         for (size_t j = 0; j < nThreads; j++) {
            batches.reserve(nThreads);
            batches.push_back(trainLoader.GetBatch());
         }
         Step(net, nets, batches);
      }

      // Compute test error.
      if ((fStepCount % fTestInterval) == 0) {

         end   = std::chrono::system_clock::now();
         std::chrono::duration<double> elapsed_seconds = end - start;
         start = std::chrono::system_clock::now();
         double seconds = elapsed_seconds.count();
         double batchesInEpoch = (double) (nTrainingSamples / net.GetBatchSize());
         double nFlops  = batchesInEpoch * fTestInterval;
         nFlops *= net.GetNFlops();
         std::cout << "Elapsed time for " << fTestInterval << " Epochs: "
                   << seconds << " [s] => " << nFlops * 1e-9 / seconds
                   << " GFlop/s" << std::endl;

         auto b = *testLoader.begin();
         auto inputMatrix  = b.GetInput();
         auto outputMatrix = b.GetOutput();
         Scalar_t loss = testNet.Loss(inputMatrix, outputMatrix);

         std::cout << "Step " << fStepCount << ": Training Error = "
                   << loss << std::endl;
         converged = HasConverged();
      }

   }
   return fMinimumError;
}

//______________________________________________________________________________
template<typename Architecture_t>
template <typename Data_t, typename Net_t>
auto TGradientDescent<Architecture_t>::TrainTBB(const Data_t & trainingData,
                                                size_t nTrainingSamples,
                                                const Data_t & testData,
                                                size_t nTestSamples,
                                                Net_t & net,
                                                size_t nThreads)
    -> typename Architecture_t::Scalar_t
{
   // Reset iteration state.
   fMinimumError = 1e100;
   fConvergenceCount = 0;
   fStepCount = 0;

   // Prepare training data.
   bool converged = false;

   // Create data loaders.
   std::vector<TDataLoader<Data_t, Architecture_t>> loaders;
   loaders.reserve(nThreads);
   for (size_t i = 0; i < nThreads; i++) {
      loaders.emplace_back(trainingData, nTrainingSamples, net.GetBatchSize(),
                           net.GetInputWidth(), net.GetOutputWidth(), 1);
      loaders.back().Shuffle();
   }
   TDataLoader<Data_t, Architecture_t> trainLoader(trainingData, nTrainingSamples,
                                                   net.GetBatchSize(),
                                                   net.GetInputWidth(),
                                                   net.GetOutputWidth(), nThreads);
   auto testNet = net.CreateClone(nTestSamples);
   TDataLoader<Data_t, Architecture_t> testLoader(testData, nTestSamples,
                                                  testNet.GetBatchSize(),
                                                  testNet.GetInputWidth(),
                                                  net.GetOutputWidth());
   // Create nets for parallel compute streams.
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

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   while (!converged)
   {
      fStepCount++;

      for (size_t i = 0; i < nTrainingSamples / net.GetBatchSize(); i += nThreads) {
         StepTBB(net, nets, loaders);
      }

      // Compute test error.
      if ((fStepCount % fTestInterval) == 0) {

         end   = std::chrono::system_clock::now();
         std::chrono::duration<double> elapsed_seconds = end - start;
         start = std::chrono::system_clock::now();
         double seconds = elapsed_seconds.count();
         double batchesInEpoch = (double) (nTrainingSamples / net.GetBatchSize());
         double nFlops  = batchesInEpoch * fTestInterval;
         nFlops *= net.GetNFlops();
         std::cout << "Elapsed time for " << fTestInterval << " Epochs: "
                   << seconds << " [s] => " << nFlops * 1e-9 / seconds
                   << " GFlop/s" << std::endl;

         auto b = *testLoader.begin();
         auto inputMatrix  = b.GetInput();
         auto outputMatrix = b.GetOutput();
         Scalar_t loss = testNet.Loss(inputMatrix, outputMatrix);

         std::cout << "Step " << fStepCount << ": Training Error = "
                   << loss << std::endl;
         converged = HasConverged();
      }

   }
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
   // Reset iteration state.
   fMinimumError = 1e100;
   fConvergenceCount = 0;
   fStepCount = 0;

   // Prepare training data.
   bool converged = false;

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

   while (!converged)
   {
      fStepCount++;

      // Iterate over epoch.
      size_t netIndex = 0;
      std::vector<TBatch<Architecture_t>> batches{};
      for (size_t i = 0; i < nTrainingSamples / net.GetBatchSize(); i += nThreads) {
         batches.clear();
         batches.reserve(nThreads);
         for (size_t j = 0; j < nThreads; j++) {
            batches.push_back(trainLoader.GetBatch());
         }
         if (momentum != 0.0) {
            StepMomentum(net, nets, batches, momentum);
         } else {
            Step(net, nets, batches);
         }
      }

      // Compute test error.
      if ((fStepCount % fTestInterval) == 0) {
         fTestError = 0.0;
         for (size_t i = 0; i < nTestSamples / net.GetBatchSize(); i += nThreads) {
            auto b = testLoader.GetBatch();
            auto inputMatrix  = b.GetInput();
            auto outputMatrix = b.GetOutput();
            fTestError += testNet.Loss(inputMatrix, outputMatrix);
         }
         fTestError /= (Double_t) nTestSamples / net.GetBatchSize();
         converged = HasConverged();
      }

   }
   return fMinimumError;
}

//______________________________________________________________________________
template<typename Architecture_t>
    template <typename Net_t>
    void inline TGradientDescent<Architecture_t>::Step(Net_t & net,
                                                       Matrix_t &input,
                                                       const Matrix_t &output)
{
    //Scalar_t loss = net.Loss(input, output);
    //fTrainingError = loss;
    net.Forward(input);
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
      nets[j].GetLayer(0).Forward(batches[j].GetInput());
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput());
      }
   }
   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(
          nets[j].GetLayer(depth-1).GetActivationGradients(),
          nets[j].GetLossFunction(),
          batches[j].GetOutput(),
          nets[j].GetLayer(depth-1).GetOutput());
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
template <typename Data_t, typename Net_t>
void inline TGradientDescent<Architecture_t>::StepTBB(
    Net_t & master,
    std::vector<Net_t> & nets,
    std::vector<TDataLoader<Data_t, Architecture_t>> & loaders)
{

   auto backprop = [&](size_t i) {
      auto batch = loaders[i].GetBatch();
      auto input = batch.GetInput();
      auto output = batch.GetOutput();

      nets[i].Forward(input);
      nets[i].Backward(input, output);
   };

   auto update = [&](size_t l) {
      for (size_t i = 0; i < nets.size(); i++) {
         Architecture_t::ScaleAdd(master.GetLayer(l).GetWeights(),
                                  nets[i].GetLayer(l).GetWeightGradients(),
                                  - fLearningRate);
         Architecture_t::Copy(nets[i].GetLayer(l).GetWeights(),
                              master.GetLayer(l).GetWeights());
         Architecture_t::ScaleAdd(master.GetLayer(l).GetBiases(),
                                  nets[i].GetLayer(l).GetBiasGradients(),
                                  - fLearningRate);
         Architecture_t::Copy(nets[i].GetLayer(l).GetBiases(),
                              master.GetLayer(l).GetBiases());
      }
   };

   tbb::parallel_for(size_t(0), nets.size(), backprop);
   tbb::parallel_for(size_t(0), master.GetDepth(), update);
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
      nets[j].GetLayer(0).Forward(batches[j].GetInput());
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput());
      }
   }
   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(
          nets[j].GetLayer(depth-1).GetActivationGradients(),
          nets[j].GetLossFunction(),
          batches[j].GetOutput(),
          nets[j].GetLayer(depth-1).GetOutput());
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
      nets[j].GetLayer(0).Forward(batches[j].GetInput());
   }

   for (size_t i = 1; i < depth; i++)
   {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayer(i).Forward(nets[j].GetLayer(i-1).GetOutput());
      }
   }

   // Gradients
   for (size_t j = 0; j < nets.size(); j++) {
      evaluateGradients<Architecture_t>(
          nets[j].GetLayer(depth-1).GetActivationGradients(),
          nets[j].GetLossFunction(),
          batches[j].GetOutput(),
          nets[j].GetLayer(depth-1).GetOutput());
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
   net.Forward(input);
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
template<typename Architecture_t>
    template <typename Net_t>
    auto inline TGradientDescent<Architecture_t>::StepReducedWeightsLoss(
        Net_t & net,
        Matrix_t &input,
        const Matrix_t &output)
    -> Scalar_t
{
   Scalar_t loss = net.Loss(input, output);
   fTrainingError = loss;
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
   return loss;
}

//______________________________________________________________________________
template<typename Architecture_t>
    template <typename Net_t>
    inline void TGradientDescent<Architecture_t>::TestError(Net_t & net,
                                                            Matrix_t &input,
                                                            const Matrix_t &output)
{
   fTestError = net.Loss(input, output, false);
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

} // namespace DNN
} // namespace TMVA

#endif
