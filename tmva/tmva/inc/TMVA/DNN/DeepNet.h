// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TDeepNet                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Deep Neural Network                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Akshay Vashistha     <akshayvashistha1995@gmail.com> - CERN, Switzerland  *
 *      Vladimir Ilievski    <ilievski.vladimir@live.com>  - CERN, Switzerland    *
 *      Saurav Shekhar       <sauravshekhar01@gmail.com> - CERN, Switzerland      *
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

#ifndef TMVA_DNN_DEEPNET
#define TMVA_DNN_DEEPNET

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/TensorDataLoader.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/DenseLayer.h"
#include "TMVA/DNN/ReshapeLayer.h"
#include "TMVA/DNN/BatchNormLayer.h"

#include "TMVA/DNN/CNN/ConvLayer.h"
#include "TMVA/DNN/CNN/MaxPoolLayer.h"

#include "TMVA/DNN/RNN/RNNLayer.h"
#include "TMVA/DNN/RNN/LSTMLayer.h"
#include "TMVA/DNN/RNN/GRULayer.h"

#ifdef HAVE_DAE
#include "TMVA/DNN/DAE/CompressionLayer.h"
#include "TMVA/DNN/DAE/CorruptionLayer.h"
#include "TMVA/DNN/DAE/ReconstructionLayer.h"
#include "TMVA/DNN/DAE/LogisticRegressionLayer.h"
#endif

#include <vector>
#include <cmath>


namespace TMVA {
namespace DNN {

   using namespace CNN;
   using namespace RNN;

   //using namespace DAE;

/** \class TDeepNet
    Generic Deep Neural Network class.
    This classs encapsulates the information for all types of Deep Neural Networks.
    \tparam Architecture The Architecture type that holds the
    architecture-specific data types.
 */
template <typename Architecture_t, typename Layer_t = VGeneralLayer<Architecture_t>>
class TDeepNet {
public:

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;


private:
   bool inline isInteger(Scalar_t x) const { return x == floor(x); }
   size_t calculateDimension(int imgDim, int fltDim, int padding, int stride);

private:
   std::vector<Layer_t *> fLayers; ///< The layers consisting the DeepNet

   size_t fBatchSize;   ///< Batch size used for training and evaluation.
   size_t fInputDepth;  ///< The depth of the input.
   size_t fInputHeight; ///< The height of the input.
   size_t fInputWidth;  ///< The width of the input.

   size_t fBatchDepth;  ///< The depth of the batch used for training/testing.
   size_t fBatchHeight; ///< The height of the batch used for training/testing.
   size_t fBatchWidth;  ///< The width of the batch used for training/testing.

   bool fIsTraining; ///< Is the network training?

   ELossFunction fJ;      ///< The loss function of the network.
   EInitialization fI;    ///< The initialization method of the network.
   ERegularization fR;    ///< The regularization used for the network.
   Scalar_t fWeightDecay; ///< The weight decay factor.

public:
   /*! Default Constructor */
   TDeepNet();

   /*! Constructor */
   TDeepNet(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t BatchDepth,
            size_t BatchHeight, size_t BatchWidth, ELossFunction fJ, EInitialization fI = EInitialization::kZero,
            ERegularization fR = ERegularization::kNone, Scalar_t fWeightDecay = 0.0, bool isTraining = false);

   /*! Copy-constructor */
   TDeepNet(const TDeepNet &);

   /*! Destructor */
   ~TDeepNet();

   /*! Function for adding Convolution layer in the Deep Neural Network,
    *  with a given depth, filter height and width, striding in rows and columns,
    *  the zero paddings, as well as the activation function and the dropout
    *  probability. Based on these parameters, it calculates the width and height
    *  of the convolutional layer. */
   TConvLayer<Architecture_t> *AddConvLayer(size_t depth, size_t filterHeight, size_t filterWidth, size_t strideRows,
                                            size_t strideCols, size_t paddingHeight, size_t paddingWidth,
                                            EActivationFunction f, Scalar_t dropoutProbability = 1.0);

   /*! Function for adding Convolution Layer in the Deep Neural Network,
    *  when the layer is already created.  */
   void AddConvLayer(TConvLayer<Architecture_t> *convLayer);

   /*! Function for adding Pooling layer in the Deep Neural Network,
    *  with a given filter height and width, striding in rows and columns as
    *  well as the dropout probability. The depth is same as the previous
    *  layer depth. Based on these parameters, it calculates the width and
    *  height of the pooling layer. */
   TMaxPoolLayer<Architecture_t> *AddMaxPoolLayer(size_t frameHeight, size_t frameWidth, size_t strideRows,
                                                  size_t strideCols, Scalar_t dropoutProbability = 1.0);
   /*! Function for adding Max Pooling layer in the Deep Neural Network,
    *  when the layer is already created. */
   void AddMaxPoolLayer(CNN::TMaxPoolLayer<Architecture_t> *maxPoolLayer);


   /*! Function for adding Recurrent Layer in the Deep Neural Network,
    * with given parameters */
   TBasicRNNLayer<Architecture_t> *AddBasicRNNLayer(size_t stateSize, size_t inputSize, size_t timeSteps,
                                                    bool rememberState = false,bool returnSequence = false,
                                                    EActivationFunction f = EActivationFunction::kTanh);

   /*! Function for adding Vanilla RNN when the layer is already created
    */
   void AddBasicRNNLayer(TBasicRNNLayer<Architecture_t> *basicRNNLayer);

   /*! Function for adding LSTM Layer in the Deep Neural Network,
    * with given parameters */
   TBasicLSTMLayer<Architecture_t> *AddBasicLSTMLayer(size_t stateSize, size_t inputSize, size_t timeSteps,
                                                    bool rememberState = false, bool returnSequence = false);

   /*! Function for adding LSTM Layer in the Deep Neural Network,
    * when the layer is already created. */
   void AddBasicLSTMLayer(TBasicLSTMLayer<Architecture_t> *basicLSTMLayer);

   /*! Function for adding GRU Layer in the Deep Neural Network,
    * with given parameters */
   TBasicGRULayer<Architecture_t> *AddBasicGRULayer(size_t stateSize, size_t inputSize, size_t timeSteps,
                                                    bool rememberState = false, bool returnSequence = false,
                                                    bool resetGateAfter = false);

   /*! Function for adding GRU Layer in the Deep Neural Network,
    * when the layer is already created. */
   void AddBasicGRULayer(TBasicGRULayer<Architecture_t> *basicGRULayer);

   /*! Function for adding Dense Connected Layer in the Deep Neural Network,
    *  with a given width, activation function and dropout probability.
    *  Based on the previous layer dimensions, it calculates the input width
    *  of the fully connected layer. */
   TDenseLayer<Architecture_t> *AddDenseLayer(size_t width, EActivationFunction f, Scalar_t dropoutProbability = 1.0);

   /*! Function for adding Dense Layer in the Deep Neural Network, when
    *  the layer is already created. */
   void AddDenseLayer(TDenseLayer<Architecture_t> *denseLayer);

   /*! Function for adding Reshape Layer in the Deep Neural Network, with a given
    *  height and width. It will take every matrix from the previous layer and
    *  reshape it to a matrix with new dimensions. */
   TReshapeLayer<Architecture_t> *AddReshapeLayer(size_t depth, size_t height, size_t width, bool flattening);

   /*! Function for adding a Batch Normalization layer with given parameters */
   TBatchNormLayer<Architecture_t> *AddBatchNormLayer(Scalar_t momentum = -1, Scalar_t epsilon = 0.0001);

   /*! Function for adding Reshape Layer in the Deep Neural Network, when
    *  the layer is already created. */
   void AddReshapeLayer(TReshapeLayer<Architecture_t> *reshapeLayer);

#ifdef HAVE_DAE   /// DAE functions
   /*! Function for adding Corruption layer in the Deep Neural Network,
    *  with given number of visibleUnits and hiddenUnits. It corrupts input
    *  according to given corruptionLevel and dropoutProbability. */
   TCorruptionLayer<Architecture_t> *AddCorruptionLayer(size_t visibleUnits, size_t hiddenUnits,
                                                        Scalar_t dropoutProbability, Scalar_t corruptionLevel);

   /*! Function for adding Corruption Layer in the Deep Neural Network,
     *  when the layer is already created.  */
   void AddCorruptionLayer(TCorruptionLayer<Architecture_t> *corruptionLayer);

   /*! Function for adding Compression layer in the Deep Neural Network,
    *  with given number of visibleUnits and hiddenUnits. It compresses the input units
    *   taking weights and biases from prev layers. */
   TCompressionLayer<Architecture_t> *AddCompressionLayer(size_t visibleUnits, size_t hiddenUnits,
                                                          Scalar_t dropoutProbability, EActivationFunction f,
                                                          std::vector<Matrix_t> weights, std::vector<Matrix_t> biases);

   /*! Function for adding Compression Layer in the Deep Neural Network, when
    *  the layer is already created. */
   void AddCompressionLayer(TCompressionLayer<Architecture_t> *compressionLayer);

   /*! Function for adding Reconstruction layer in the Deep Neural Network,
    *  with given number of visibleUnits and hiddenUnits. It reconstructs the input units
    *  taking weights and biases from prev layers. Same corruptionLevel and dropoutProbability
    *  must be passed as in corruptionLayer. */
   TReconstructionLayer<Architecture_t> *AddReconstructionLayer(size_t visibleUnits, size_t hiddenUnits,
                                                                Scalar_t learningRate, EActivationFunction f,
                                                                std::vector<Matrix_t> weights,
                                                                std::vector<Matrix_t> biases, Scalar_t corruptionLevel,
                                                                Scalar_t dropoutProbability);

   /*! Function for adding Reconstruction Layer in the Deep Neural Network, when
    *  the layer is already created. */
   void AddReconstructionLayer(TReconstructionLayer<Architecture_t> *reconstructionLayer);

   /*! Function for adding logisticRegressionLayer in the Deep Neural Network,
    *  with given number of inputUnits and outputUnits. It classifies the outputUnits. */
   TLogisticRegressionLayer<Architecture_t> *AddLogisticRegressionLayer(size_t inputUnits, size_t outputUnits,
                                                                        size_t testDataBatchSize,
                                                                        Scalar_t learningRate);

   /*! Function for adding logisticRegressionLayer in the Deep Neural Network, when
    *  the layer is already created. */
   void AddLogisticRegressionLayer(TLogisticRegressionLayer<Architecture_t> *logisticRegressionLayer);

   /* To train the Deep AutoEncoder network with required number of Corruption, Compression and Reconstruction
    * layers. */
   void PreTrain(std::vector<Matrix_t> &input, std::vector<size_t> numHiddenUnitsPerLayer, Scalar_t learningRate,
                 Scalar_t corruptionLevel, Scalar_t dropoutProbability, size_t epochs, EActivationFunction f,
                 bool applyDropout = false);

   /* To classify outputLabel in Deep AutoEncoder. Should be used after PreTrain if required.
    * Currently, it used Logistic Regression Layer. Otherwise we can use any other classification layer also.
   */
   void FineTune(std::vector<Matrix_t> &input, std::vector<Matrix_t> &testInput, std::vector<Matrix_t> &outputLabel,
                 size_t outputUnits, size_t testDataBatchSize, Scalar_t learningRate, size_t epochs);
#endif

   /*! Function for initialization of the Neural Net. */
   void Initialize();

   /*! Function that executes the entire forward pass in the network. */
   void Forward(Tensor_t &input, bool applyDropout = false);

    /*! Function that reset some training flags after looping all the events but not the weights*/
   void ResetTraining();



   /*! Function that executes the entire backward pass in the network. */
   void Backward(const Tensor_t &input, const Matrix_t &groundTruth, const Matrix_t &weights);


#ifdef USE_PARALLEL_DEEPNET
   /*! Function for parallel forward in the vector of deep nets, where the master
    *  net is the net calling this function. There is one batch for one deep net.*/
   void ParallelForward(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                        std::vector<TTensorBatch<Architecture_t>> &batches, bool applyDropout = false);

   /*! Function for parallel backward in the vector of deep nets, where the master
    *  net is the net calling this function and getting the updates from the other nets.
    * There is one batch for one deep net.*/
   void ParallelBackward(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                         std::vector<TTensorBatch<Architecture_t>> &batches, Scalar_t learningRate);

   /*! Function for parallel backward in the vector of deep nets, where the master
    *  net is the net calling this function and getting the updates from the other nets,
    *  following the momentum strategy. There is one batch for one deep net.*/
   void ParallelBackwardMomentum(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                 std::vector<TTensorBatch<Architecture_t>> &batches, Scalar_t learningRate,
                                 Scalar_t momentum);

   /*! Function for parallel backward in the vector of deep nets, where the master
    *  net is the net calling this function and getting the updates from the other nets,
    *  following the Nestorov momentum strategy. There is one batch for one deep net.*/
   void ParallelBackwardNestorov(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                 std::vector<TTensorBatch<Architecture_t>> &batches, Scalar_t learningRate,
                                 Scalar_t momentum);

#endif // endif use parallel deepnet

   /*! Function that will update the weights and biases in the layers that
    *  contain weights and biases.  */
   void Update(Scalar_t learningRate);

   /*! Function for evaluating the loss, based on the activations stored
    *  in the last layer. */
   Scalar_t Loss(const Matrix_t &groundTruth, const Matrix_t &weights, bool includeRegularization = true) const;

   /*! Function for evaluating the loss, based on the propagation of the given input. */
   Scalar_t Loss(Tensor_t &input, const Matrix_t &groundTruth, const Matrix_t &weights,
                 bool inTraining = false, bool includeRegularization = true);

   /*! Function for computing the regularizaton term to be added to the loss function  */
   Scalar_t RegularizationTerm() const;

   /*! Prediction based on activations stored in the last layer. */
   void Prediction(Matrix_t &predictions, EOutputFunction f) const;

   /*! Prediction for the given inputs, based on what network learned. */
   void Prediction(Matrix_t &predictions, Tensor_t & input, EOutputFunction f);

   /*! Print the Deep Net Info */
   void Print() const;

   /*! Get the layer in the vector of layers at poistion i */
   inline Layer_t *GetLayerAt(size_t i) { return fLayers[i]; }
   inline const Layer_t *GetLayerAt(size_t i) const { return fLayers[i]; }

   /* Depth and the output width of the network. */
   inline size_t GetDepth() const { return fLayers.size(); }
   inline size_t GetOutputWidth() const { return fLayers.back()->GetWidth(); }

   /* Return a reference to the layers. */
   inline std::vector<Layer_t *> &GetLayers() { return fLayers; }
   inline const std::vector<Layer_t *> &GetLayers() const { return fLayers; }

   /*! Remove all layers from the network. */
   inline void Clear() { fLayers.clear(); }

   /*! Getters */
   inline size_t GetBatchSize() const { return fBatchSize; }
   inline size_t GetInputDepth() const { return fInputDepth; }
   inline size_t GetInputHeight() const { return fInputHeight; }
   inline size_t GetInputWidth() const { return fInputWidth; }

   inline size_t GetBatchDepth() const { return fBatchDepth; }
   inline size_t GetBatchHeight() const { return fBatchHeight; }
   inline size_t GetBatchWidth() const { return fBatchWidth; }

   inline bool IsTraining() const { return fIsTraining; }

   inline ELossFunction GetLossFunction() const { return fJ; }
   inline EInitialization GetInitialization() const { return fI; }
   inline ERegularization GetRegularization() const { return fR; }
   inline Scalar_t GetWeightDecay() const { return fWeightDecay; }

   /*! Setters */
   // FIXME many of these won't work as the data structure storing activations
   // and gradients have not changed in all the layers, also params in layers
   // have not changed either
   inline void SetBatchSize(size_t batchSize) { fBatchSize = batchSize; }
   inline void SetInputDepth(size_t inputDepth) { fInputDepth = inputDepth; }
   inline void SetInputHeight(size_t inputHeight) { fInputHeight = inputHeight; }
   inline void SetInputWidth(size_t inputWidth) { fInputWidth = inputWidth; }
   inline void SetBatchDepth(size_t batchDepth) { fBatchDepth = batchDepth; }
   inline void SetBatchHeight(size_t batchHeight) { fBatchHeight = batchHeight; }
   inline void SetBatchWidth(size_t batchWidth) { fBatchWidth = batchWidth; }
   inline void SetLossFunction(ELossFunction J) { fJ = J; }
   inline void SetInitialization(EInitialization I) { fI = I; }
   inline void SetRegularization(ERegularization R) { fR = R; }
   inline void SetWeightDecay(Scalar_t weightDecay) { fWeightDecay = weightDecay; }

   void SetDropoutProbabilities(const std::vector<Double_t> & probabilities);

};

//
//  Deep Net Class - Implementation
//
//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TDeepNet<Architecture_t, Layer_t>::TDeepNet()
   : fLayers(), fBatchSize(0), fInputDepth(0), fInputHeight(0), fInputWidth(0), fBatchDepth(0), fBatchHeight(0),
     fBatchWidth(0), fJ(ELossFunction::kMeanSquaredError), fI(EInitialization::kZero), fR(ERegularization::kNone),
     fIsTraining(true), fWeightDecay(0.0)
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TDeepNet<Architecture_t, Layer_t>::TDeepNet(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                            size_t batchDepth, size_t batchHeight, size_t batchWidth, ELossFunction J,
                                            EInitialization I, ERegularization R, Scalar_t weightDecay, bool isTraining)
   : fLayers(), fBatchSize(batchSize), fInputDepth(inputDepth), fInputHeight(inputHeight), fInputWidth(inputWidth),
     fBatchDepth(batchDepth), fBatchHeight(batchHeight), fBatchWidth(batchWidth), fIsTraining(isTraining), fJ(J), fI(I),
     fR(R), fWeightDecay(weightDecay)
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TDeepNet<Architecture_t, Layer_t>::TDeepNet(const TDeepNet &deepNet)
   : fLayers(), fBatchSize(deepNet.fBatchSize), fInputDepth(deepNet.fInputDepth), fInputHeight(deepNet.fInputHeight),
     fInputWidth(deepNet.fInputWidth), fBatchDepth(deepNet.fBatchDepth), fBatchHeight(deepNet.fBatchHeight),
     fBatchWidth(deepNet.fBatchWidth), fIsTraining(deepNet.fIsTraining), fJ(deepNet.fJ), fI(deepNet.fI), fR(deepNet.fR),
     fWeightDecay(deepNet.fWeightDecay)
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TDeepNet<Architecture_t, Layer_t>::~TDeepNet()
{
   // Relese the layers memory
   for (auto  layer : fLayers)
      delete layer;
   fLayers.clear();
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::calculateDimension(int imgDim, int fltDim, int padding, int stride) -> size_t
{
   Scalar_t dimension = ((imgDim - fltDim + 2 * padding) / stride) + 1;
   if (!isInteger(dimension) || dimension <= 0) {
      this->Print();
      int iLayer = fLayers.size();
      Fatal("calculateDimension","Not compatible hyper parameters for layer %d - (imageDim, filterDim, padding, stride) %d , %d , %d , %d",
            iLayer, imgDim, fltDim, padding, stride);
      // std::cout << " calculateDimension - Not compatible hyper parameters (imgDim, fltDim, padding, stride)"
      //           << imgDim << " , " << fltDim << " , " <<  padding << " , " << stride<< " resulting dim is " << dimension << std::endl;
      // std::exit(EXIT_FAILURE);
   }

   return (size_t)dimension;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TConvLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddConvLayer(size_t depth, size_t filterHeight,
                                                                            size_t filterWidth, size_t strideRows,
                                                                            size_t strideCols, size_t paddingHeight,
                                                                            size_t paddingWidth, EActivationFunction f,
                                                                            Scalar_t dropoutProbability)
{
   // All variables defining a convolutional layer
   size_t batchSize = this->GetBatchSize();
   size_t inputDepth;
   size_t inputHeight;
   size_t inputWidth;
   EInitialization init = this->GetInitialization();
   ERegularization reg = this->GetRegularization();
   Scalar_t decay = this->GetWeightDecay();

   if (fLayers.size() == 0) {
      inputDepth = this->GetInputDepth();
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputDepth = lastLayer->GetDepth();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
   }



   // Create the conv layer
   TConvLayer<Architecture_t> *convLayer = new TConvLayer<Architecture_t>(
           batchSize, inputDepth, inputHeight, inputWidth, depth, init, filterHeight, filterWidth, strideRows,
           strideCols, paddingHeight, paddingWidth, dropoutProbability, f, reg, decay);

   fLayers.push_back(convLayer);
   return convLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddConvLayer(TConvLayer<Architecture_t> *convLayer)
{
   fLayers.push_back(convLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TMaxPoolLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddMaxPoolLayer(size_t frameHeight, size_t frameWidth,
                                                                                  size_t strideRows, size_t strideCols,
                                                                                  Scalar_t dropoutProbability)
{
   size_t batchSize = this->GetBatchSize();
   size_t inputDepth;
   size_t inputHeight;
   size_t inputWidth;

   if (fLayers.size() == 0) {
      inputDepth = this->GetInputDepth();
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputDepth = lastLayer->GetDepth();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
   }

   TMaxPoolLayer<Architecture_t> *maxPoolLayer = new TMaxPoolLayer<Architecture_t>(
      batchSize, inputDepth, inputHeight, inputWidth, frameHeight, frameWidth,
      strideRows, strideCols, dropoutProbability);

   // But this creates a copy or what?
   fLayers.push_back(maxPoolLayer);

   return maxPoolLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddMaxPoolLayer(TMaxPoolLayer<Architecture_t> *maxPoolLayer)
{
   fLayers.push_back(maxPoolLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TBasicRNNLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddBasicRNNLayer(size_t stateSize, size_t inputSize,
                                                                                    size_t timeSteps,
                                                                                    bool rememberState, bool returnSequence,
                                                                                    EActivationFunction f)
{

   // should check if input and time size are consistent

   //std::cout << "Create RNN " << fLayers.size() << "  " << this->GetInputHeight() << "  " << this->GetInputWidth() << std::endl;
   size_t inputHeight, inputWidth, inputDepth;
   if (fLayers.size() == 0) {
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
      inputDepth = this->GetInputDepth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
      inputDepth = lastLayer->GetDepth();
   }
   if (inputSize != inputWidth) {
      Error("AddBasicRNNLayer","Inconsistent input size with input layout  - it should be %zu instead of %zu",inputSize, inputWidth);
   }
   if (timeSteps != inputHeight && timeSteps != inputDepth) {
      Error("AddBasicRNNLayer","Inconsistent time steps with input layout - it should be %zu instead of %zu or %zu",timeSteps, inputHeight,inputDepth);
   }

   TBasicRNNLayer<Architecture_t> *basicRNNLayer =
      new TBasicRNNLayer<Architecture_t>(this->GetBatchSize(), stateSize, inputSize, timeSteps, rememberState, returnSequence,
                                         f, fIsTraining, this->GetInitialization());
   fLayers.push_back(basicRNNLayer);
   return basicRNNLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddBasicRNNLayer(TBasicRNNLayer<Architecture_t> *basicRNNLayer)
{
   fLayers.push_back(basicRNNLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TBasicLSTMLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddBasicLSTMLayer(size_t stateSize, size_t inputSize,
                                                                                      size_t timeSteps, bool rememberState, bool returnSequence)
{
   // should check if input and time size are consistent
   size_t inputHeight, inputWidth, inputDepth;
   if (fLayers.size() == 0) {
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
      inputDepth = this->GetInputDepth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
      inputDepth = lastLayer->GetDepth();
   }
   if (inputSize != inputWidth) {
      Error("AddBasicLSTMLayer", "Inconsistent input size with input layout  - it should be %zu instead of %zu", inputSize, inputWidth);
   }
   if (timeSteps != inputHeight && timeSteps != inputDepth) {
      Error("AddBasicLSTMLayer", "Inconsistent time steps with input layout - it should be %zu instead of %zu", timeSteps, inputHeight);
   }

   TBasicLSTMLayer<Architecture_t> *basicLSTMLayer =
      new TBasicLSTMLayer<Architecture_t>(this->GetBatchSize(), stateSize, inputSize, timeSteps, rememberState, returnSequence,
                                         DNN::EActivationFunction::kSigmoid,
                                         DNN::EActivationFunction::kTanh,
                                         fIsTraining, this->GetInitialization());
   fLayers.push_back(basicLSTMLayer);
   return basicLSTMLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddBasicLSTMLayer(TBasicLSTMLayer<Architecture_t> *basicLSTMLayer)
{
   fLayers.push_back(basicLSTMLayer);
}


//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TBasicGRULayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddBasicGRULayer(size_t stateSize, size_t inputSize,
                                                                                      size_t timeSteps, bool rememberState, bool returnSequence, bool resetGateAfter)
{
   // should check if input and time size are consistent
   size_t inputHeight, inputWidth, inputDepth;
   if (fLayers.size() == 0) {
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
      inputDepth = this->GetInputDepth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
      inputDepth = lastLayer->GetDepth();
   }
   if (inputSize != inputWidth) {
      Error("AddBasicGRULayer", "Inconsistent input size with input layout  - it should be %zu instead of %zu", inputSize, inputWidth);
   }
   if (timeSteps != inputHeight && timeSteps != inputDepth) {
      Error("AddBasicGRULayer", "Inconsistent time steps with input layout - it should be %zu instead of %zu", timeSteps, inputHeight);
   }

   TBasicGRULayer<Architecture_t> *basicGRULayer =
      new TBasicGRULayer<Architecture_t>(this->GetBatchSize(), stateSize, inputSize, timeSteps, rememberState, returnSequence, resetGateAfter,
                                         DNN::EActivationFunction::kSigmoid,
                                         DNN::EActivationFunction::kTanh,
                                         fIsTraining, this->GetInitialization());
   fLayers.push_back(basicGRULayer);
   return basicGRULayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddBasicGRULayer(TBasicGRULayer<Architecture_t> *basicGRULayer)
{
   fLayers.push_back(basicGRULayer);
}



//DAE
#ifdef HAVE_DAE

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TCorruptionLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddCorruptionLayer(size_t visibleUnits,
                                                                                        size_t hiddenUnits,
                                                                                        Scalar_t dropoutProbability,
                                                                                        Scalar_t corruptionLevel)
{
   size_t batchSize = this->GetBatchSize();

   TCorruptionLayer<Architecture_t> *corruptionLayer =
      new TCorruptionLayer<Architecture_t>(batchSize, visibleUnits, hiddenUnits, dropoutProbability, corruptionLevel);
   fLayers.push_back(corruptionLayer);
   return corruptionLayer;
}
//______________________________________________________________________________

template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddCorruptionLayer(TCorruptionLayer<Architecture_t> *corruptionLayer)
{
   fLayers.push_back(corruptionLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TCompressionLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddCompressionLayer(
   size_t visibleUnits, size_t hiddenUnits, Scalar_t dropoutProbability, EActivationFunction f,
   std::vector<Matrix_t> weights, std::vector<Matrix_t> biases)
{
   size_t batchSize = this->GetBatchSize();

   TCompressionLayer<Architecture_t> *compressionLayer = new TCompressionLayer<Architecture_t>(
      batchSize, visibleUnits, hiddenUnits, dropoutProbability, f, weights, biases);
   fLayers.push_back(compressionLayer);
   return compressionLayer;
}
//______________________________________________________________________________

template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddCompressionLayer(TCompressionLayer<Architecture_t> *compressionLayer)
{
   fLayers.push_back(compressionLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TReconstructionLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddReconstructionLayer(
   size_t visibleUnits, size_t hiddenUnits, Scalar_t learningRate, EActivationFunction f, std::vector<Matrix_t> weights,
   std::vector<Matrix_t> biases, Scalar_t corruptionLevel, Scalar_t dropoutProbability)
{
   size_t batchSize = this->GetBatchSize();

   TReconstructionLayer<Architecture_t> *reconstructionLayer = new TReconstructionLayer<Architecture_t>(
      batchSize, visibleUnits, hiddenUnits, learningRate, f, weights, biases, corruptionLevel, dropoutProbability);
   fLayers.push_back(reconstructionLayer);
   return reconstructionLayer;
}
//______________________________________________________________________________

template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddReconstructionLayer(
   TReconstructionLayer<Architecture_t> *reconstructionLayer)
{
   fLayers.push_back(reconstructionLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TLogisticRegressionLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddLogisticRegressionLayer(
   size_t inputUnits, size_t outputUnits, size_t testDataBatchSize, Scalar_t learningRate)
{
   size_t batchSize = this->GetBatchSize();

   TLogisticRegressionLayer<Architecture_t> *logisticRegressionLayer =
      new TLogisticRegressionLayer<Architecture_t>(batchSize, inputUnits, outputUnits, testDataBatchSize, learningRate);
   fLayers.push_back(logisticRegressionLayer);
   return logisticRegressionLayer;
}
//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddLogisticRegressionLayer(
   TLogisticRegressionLayer<Architecture_t> *logisticRegressionLayer)
{
   fLayers.push_back(logisticRegressionLayer);
}
#endif


//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TDenseLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddDenseLayer(size_t width, EActivationFunction f,
                                                                              Scalar_t dropoutProbability)
{
   size_t batchSize = this->GetBatchSize();
   size_t inputWidth;
   EInitialization init = this->GetInitialization();
   ERegularization reg = this->GetRegularization();
   Scalar_t decay = this->GetWeightDecay();

   if (fLayers.size() == 0) {
      inputWidth = this->GetInputWidth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputWidth = lastLayer->GetWidth();
   }

   TDenseLayer<Architecture_t> *denseLayer =
      new TDenseLayer<Architecture_t>(batchSize, inputWidth, width, init, dropoutProbability, f, reg, decay);

   fLayers.push_back(denseLayer);

   return denseLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddDenseLayer(TDenseLayer<Architecture_t> *denseLayer)
{
   fLayers.push_back(denseLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TReshapeLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddReshapeLayer(size_t depth, size_t height,
                                                                                  size_t width, bool flattening)
{
   size_t batchSize = this->GetBatchSize();
   size_t inputDepth;
   size_t inputHeight;
   size_t inputWidth;
   size_t outputNSlices;
   size_t outputNRows;
   size_t outputNCols;

   if (fLayers.size() == 0) {
      inputDepth = this->GetInputDepth();
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputDepth = lastLayer->GetDepth();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
   }

   if (flattening) {
      outputNSlices = 1;
      outputNRows = this->GetBatchSize();
      outputNCols = depth * height * width;
      size_t inputNCols =  inputDepth * inputHeight *  inputWidth;
      if (outputNCols != 0 && outputNCols != inputNCols ) {
         Info("AddReshapeLayer","Dimensions not compatibles - product of input %zu x %zu x %zu should be equal to output %zu x %zu x %zu - Force flattening output to be %zu",
              inputDepth, inputHeight, inputWidth, depth, height, width,inputNCols);
      }
      outputNCols = inputNCols;
      depth = 1;
      height = 1;
      width = outputNCols;
   } else {
      outputNSlices = this->GetBatchSize();
      outputNRows = depth;
      outputNCols = height * width;
   }

   TReshapeLayer<Architecture_t> *reshapeLayer =
      new TReshapeLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, depth, height, width,
                                        outputNSlices, outputNRows, outputNCols, flattening);

   fLayers.push_back(reshapeLayer);

   return reshapeLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
TBatchNormLayer<Architecture_t> *TDeepNet<Architecture_t, Layer_t>::AddBatchNormLayer(Scalar_t momentum, Scalar_t epsilon)
{
   int axis = -1;
   size_t batchSize = this->GetBatchSize();
   size_t inputDepth = 0;
   size_t inputHeight = 0;
   size_t inputWidth = 0;
   // this is the shape of the output tensor (it is columnmajor by default)
   // and it is normally (depth, hw, bsize)  and for dense layers  (bsize, w, 1)
   std::vector<size_t>  shape = {1, 1, 1};
   if (fLayers.size() == 0) {
      inputDepth = this->GetInputDepth();
      inputHeight = this->GetInputHeight();
      inputWidth = this->GetInputWidth();
      // assume that is like for a dense layer
      shape[0] = batchSize;
      shape[1] = inputWidth;
      shape[2] = 1;
   } else {
      Layer_t *lastLayer = fLayers.back();
      inputDepth = lastLayer->GetDepth();
      inputHeight = lastLayer->GetHeight();
      inputWidth = lastLayer->GetWidth();
      shape = lastLayer->GetOutput().GetShape();
      if (dynamic_cast<TConvLayer<Architecture_t> *>(lastLayer) != nullptr ||
          dynamic_cast<TMaxPoolLayer<Architecture_t> *>(lastLayer) != nullptr)
         axis = 1; // use axis = channel axis for convolutional layer
      if (shape.size() > 3) {
         for (size_t i = 3; i < shape.size(); ++i)
            shape[2] *= shape[i];
      }
   }
   // std::cout << "addBNormLayer " << inputDepth << " , " << inputHeight << " , " << inputWidth << " , " << shape[0]
   //           << "  " << shape[1] << "  " << shape[2] << std::endl;

   auto bnormLayer =
      new TBatchNormLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, shape, axis, momentum, epsilon);

   fLayers.push_back(bnormLayer);

   return bnormLayer;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::AddReshapeLayer(TReshapeLayer<Architecture_t> *reshapeLayer)
{
   fLayers.push_back(reshapeLayer);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Initialize() -> void
{
   for (size_t i = 0; i < fLayers.size(); i++) {
      fLayers[i]->Initialize();
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::ResetTraining() -> void
{
   for (size_t i = 0; i < fLayers.size(); i++) {
      fLayers[i]->ResetTraining();
   }
}


//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Forward( Tensor_t &input, bool applyDropout) -> void
{
   fLayers.front()->Forward(input, applyDropout);

   for (size_t i = 1; i < fLayers.size(); i++) {
      fLayers[i]->Forward(fLayers[i - 1]->GetOutput(), applyDropout);
      //std::cout << "forward for layer " << i << std::endl;
      // fLayers[i]->GetOutput()[0].Print();
   }
}


#ifdef HAVE_DAE
//_____________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::PreTrain(std::vector<Matrix_t> &input,
                                                 std::vector<size_t> numHiddenUnitsPerLayer, Scalar_t learningRate,
                                                 Scalar_t corruptionLevel, Scalar_t dropoutProbability, size_t epochs,
                                                 EActivationFunction f, bool applyDropout) -> void
{
   std::vector<Matrix_t> inp1;
   std::vector<Matrix_t> inp2;
   size_t numOfHiddenLayers = sizeof(numHiddenUnitsPerLayer) / sizeof(numHiddenUnitsPerLayer[0]);
   // size_t batchSize = this->GetBatchSize();
   size_t visibleUnits = (size_t)input[0].GetNrows();

   AddCorruptionLayer(visibleUnits, numHiddenUnitsPerLayer[0], dropoutProbability, corruptionLevel);
   fLayers.back()->Initialize();
   fLayers.back()->Forward(input, applyDropout);
   // fLayers.back()->Print();

   AddCompressionLayer(visibleUnits, numHiddenUnitsPerLayer[0], dropoutProbability, f, fLayers.back()->GetWeights(),
                       fLayers.back()->GetBiases());
   fLayers.back()->Initialize();
   fLayers.back()->Forward(fLayers[fLayers.size() - 2]->GetOutput(), applyDropout); // as we have to pass corrupt input

   AddReconstructionLayer(visibleUnits, numHiddenUnitsPerLayer[0], learningRate, f, fLayers.back()->GetWeights(),
                          fLayers.back()->GetBiases(), corruptionLevel, dropoutProbability);
   fLayers.back()->Initialize();
   fLayers.back()->Forward(fLayers[fLayers.size() - 2]->GetOutput(),
                           applyDropout); // as we have to pass compressed Input
   fLayers.back()->Backward(fLayers[fLayers.size() - 2]->GetOutput(), inp1, fLayers[fLayers.size() - 3]->GetOutput(),
                            input);
   // three layers are added, now pointer is on third layer
   size_t weightsSize = fLayers.back()->GetWeights().size();
   size_t biasesSize = fLayers.back()->GetBiases().size();
   for (size_t epoch = 0; epoch < epochs - 1; epoch++) {
      // fLayers[fLayers.size() - 3]->Forward(input,applyDropout);
      for (size_t j = 0; j < weightsSize; j++) {
         Architecture_t::Copy(fLayers[fLayers.size() - 2]->GetWeightsAt(j), fLayers.back()->GetWeightsAt(j));
      }
      for (size_t j = 0; j < biasesSize; j++) {
         Architecture_t::Copy(fLayers[fLayers.size() - 2]->GetBiasesAt(j), fLayers.back()->GetBiasesAt(j));
      }
      fLayers[fLayers.size() - 2]->Forward(fLayers[fLayers.size() - 3]->GetOutput(), applyDropout);
      fLayers[fLayers.size() - 1]->Forward(fLayers[fLayers.size() - 2]->GetOutput(), applyDropout);
      fLayers[fLayers.size() - 1]->Backward(fLayers[fLayers.size() - 2]->GetOutput(), inp1,
                                            fLayers[fLayers.size() - 3]->GetOutput(), input);
   }
   fLayers.back()->Print();

   for (size_t i = 1; i < numOfHiddenLayers; i++) {

      AddCorruptionLayer(numHiddenUnitsPerLayer[i - 1], numHiddenUnitsPerLayer[i], dropoutProbability, corruptionLevel);
      fLayers.back()->Initialize();
      fLayers.back()->Forward(fLayers[fLayers.size() - 3]->GetOutput(),
                              applyDropout); // as we have to pass compressed Input

      AddCompressionLayer(numHiddenUnitsPerLayer[i - 1], numHiddenUnitsPerLayer[i], dropoutProbability, f,
                          fLayers.back()->GetWeights(), fLayers.back()->GetBiases());
      fLayers.back()->Initialize();
      fLayers.back()->Forward(fLayers[fLayers.size() - 2]->GetOutput(), applyDropout);

      AddReconstructionLayer(numHiddenUnitsPerLayer[i - 1], numHiddenUnitsPerLayer[i], learningRate, f,
                             fLayers.back()->GetWeights(), fLayers.back()->GetBiases(), corruptionLevel,
                             dropoutProbability);
      fLayers.back()->Initialize();
      fLayers.back()->Forward(fLayers[fLayers.size() - 2]->GetOutput(),
                              applyDropout); // as we have to pass compressed Input
      fLayers.back()->Backward(fLayers[fLayers.size() - 2]->GetOutput(), inp1, fLayers[fLayers.size() - 3]->GetOutput(),
                               fLayers[fLayers.size() - 5]->GetOutput());

      // three layers are added, now pointer is on third layer
      size_t _weightsSize = fLayers.back()->GetWeights().size();
      size_t _biasesSize = fLayers.back()->GetBiases().size();
      for (size_t epoch = 0; epoch < epochs - 1; epoch++) {
         // fLayers[fLayers.size() - 3]->Forward(input,applyDropout);
         for (size_t j = 0; j < _weightsSize; j++) {
            Architecture_t::Copy(fLayers[fLayers.size() - 2]->GetWeightsAt(j), fLayers.back()->GetWeightsAt(j));
         }
         for (size_t j = 0; j < _biasesSize; j++) {
            Architecture_t::Copy(fLayers[fLayers.size() - 2]->GetBiasesAt(j), fLayers.back()->GetBiasesAt(j));
         }
         fLayers[fLayers.size() - 2]->Forward(fLayers[fLayers.size() - 3]->GetOutput(), applyDropout);
         fLayers[fLayers.size() - 1]->Forward(fLayers[fLayers.size() - 2]->GetOutput(), applyDropout);
         fLayers[fLayers.size() - 1]->Backward(fLayers[fLayers.size() - 2]->GetOutput(), inp1,
                                               fLayers[fLayers.size() - 3]->GetOutput(),
                                               fLayers[fLayers.size() - 5]->GetOutput());
      }
      fLayers.back()->Print();
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::FineTune(std::vector<Matrix_t> &input, std::vector<Matrix_t> &testInput,
                                                 std::vector<Matrix_t> &inputLabel, size_t outputUnits,
                                                 size_t testDataBatchSize, Scalar_t learningRate, size_t epochs) -> void
{
   std::vector<Matrix_t> inp1;
   std::vector<Matrix_t> inp2;
   if (fLayers.size() == 0) // only Logistic Regression Layer
   {
      size_t inputUnits = input[0].GetNrows();

      AddLogisticRegressionLayer(inputUnits, outputUnits, testDataBatchSize, learningRate);
      fLayers.back()->Initialize();
      for (size_t i = 0; i < epochs; i++) {
         fLayers.back()->Backward(inputLabel, inp1, input, inp2);
      }
      fLayers.back()->Forward(input, false);
      fLayers.back()->Print();
   } else { // if used after any other layer
      size_t inputUnits = fLayers.back()->GetOutputAt(0).GetNrows();
      AddLogisticRegressionLayer(inputUnits, outputUnits, testDataBatchSize, learningRate);
      fLayers.back()->Initialize();
      for (size_t i = 0; i < epochs; i++) {
         fLayers.back()->Backward(inputLabel, inp1, fLayers[fLayers.size() - 2]->GetOutput(), inp2);
      }
      fLayers.back()->Forward(testInput, false);
      fLayers.back()->Print();
   }
}
#endif

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Backward(const Tensor_t &input, const Matrix_t &groundTruth,
                                                 const Matrix_t &weights) -> void
{
   //Tensor_t inp1;
   //Tensor_t inp2;
   // Last layer should be dense layer
   Matrix_t last_actgrad = fLayers.back()->GetActivationGradientsAt(0);
   Matrix_t last_output = fLayers.back()->GetOutputAt(0);
   evaluateGradients<Architecture_t>(last_actgrad, this->GetLossFunction(), groundTruth,
                                     last_output, weights);

   for (size_t i = fLayers.size() - 1; i > 0; i--) {
      auto &activation_gradient_backward = fLayers[i - 1]->GetActivationGradients();
      auto &activations_backward = fLayers[i - 1]->GetOutput();
      fLayers[i]->Backward(activation_gradient_backward, activations_backward);
   }

   // need to have a dummy tensor (size=0) to pass for activation gradient backward which
   // are not computed for the first layer
   Tensor_t dummy;
   fLayers[0]->Backward(dummy, input);
}

#ifdef USE_PARALLEL_DEEPNET

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::ParallelForward(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                                        std::vector<TTensorBatch<Architecture_t>> &batches,
                                                        bool applyDropout) -> void
{
   size_t depth = this->GetDepth();

   // The first layer of each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].GetLayerAt(0)->Forward(batches[i].GetInput(), applyDropout);
   }

   // The i'th layer of each deep net
   for (size_t i = 1; i < depth; i++) {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayerAt(i)->Forward(nets[j].GetLayerAt(i - 1)->GetOutput(), applyDropout);
      }
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::ParallelBackward(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                                         std::vector<TTensorBatch<Architecture_t>> &batches,
                                                         Scalar_t learningRate) -> void
{
   std::vector<Matrix_t> inp1;
   std::vector<Matrix_t> inp2;
   size_t depth = this->GetDepth();

   // Evaluate the gradients of the last layers in each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      evaluateGradients<Architecture_t>(nets[i].GetLayerAt(depth - 1)->GetActivationGradientsAt(0),
                                        nets[i].GetLossFunction(), batches[i].GetOutput(),
                                        nets[i].GetLayerAt(depth - 1)->GetOutputAt(0), batches[i].GetWeights());
   }

   // Backpropagate the error in i'th layer of each deep net
   for (size_t i = depth - 1; i > 0; i--) {
      for (size_t j = 0; j < nets.size(); j++) {
         nets[j].GetLayerAt(i)->Backward(nets[j].GetLayerAt(i - 1)->GetActivationGradients(),
                                         nets[j].GetLayerAt(i - 1)->GetOutput(), inp1, inp2);
      }
   }

   std::vector<Matrix_t> dummy;

   // First layer of each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      nets[i].GetLayerAt(0)->Backward(dummy, batches[i].GetInput(), inp1, inp2);
   }

   // Update and copy
   for (size_t i = 0; i < nets.size(); i++) {
      for (size_t j = 0; j < depth; j++) {
         Layer_t *masterLayer = this->GetLayerAt(j);
         Layer_t *layer = nets[i].GetLayerAt(j);

         masterLayer->UpdateWeights(layer->GetWeightGradients(), learningRate);
         layer->CopyWeights(masterLayer->GetWeights());

         masterLayer->UpdateBiases(layer->GetBiasGradients(), learningRate);
         layer->CopyBiases(masterLayer->GetBiases());
      }
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::ParallelBackwardMomentum(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                                                 std::vector<TTensorBatch<Architecture_t>> &batches,
                                                                 Scalar_t learningRate, Scalar_t momentum) -> void
{
   std::vector<Matrix_t> inp1;
   std::vector<Matrix_t> inp2;
   size_t depth = this->GetDepth();

   // Evaluate the gradients of the last layers in each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      evaluateGradients<Architecture_t>(nets[i].GetLayerAt(depth - 1)->GetActivationGradientsAt(0),
                                        nets[i].GetLossFunction(), batches[i].GetOutput(),
                                        nets[i].GetLayerAt(depth - 1)->GetOutputAt(0), batches[i].GetWeights());
   }

   // Backpropagate the error in i'th layer of each deep net
   for (size_t i = depth - 1; i > 0; i--) {
      Layer_t *masterLayer = this->GetLayerAt(i);

      for (size_t j = 0; j < nets.size(); j++) {
         Layer_t *layer = nets[j].GetLayerAt(i);

         layer->Backward(nets[j].GetLayerAt(i - 1)->GetActivationGradients(), nets[j].GetLayerAt(i - 1)->GetOutput(),
                         inp1, inp2);
         masterLayer->UpdateWeightGradients(layer->GetWeightGradients(), learningRate / momentum);
         masterLayer->UpdateBiasGradients(layer->GetBiasGradients(), learningRate / momentum);
      }

      masterLayer->UpdateWeightGradients(masterLayer->GetWeightGradients(), 1.0 - momentum);
      masterLayer->UpdateBiasGradients(masterLayer->GetBiasGradients(), 1.0 - momentum);
   }

   std::vector<Matrix_t> dummy;

   // First layer of each deep net
   Layer_t *masterFirstLayer = this->GetLayerAt(0);
   for (size_t i = 0; i < nets.size(); i++) {
      Layer_t *layer = nets[i].GetLayerAt(0);

      layer->Backward(dummy, batches[i].GetInput(), inp1, inp2);

      masterFirstLayer->UpdateWeightGradients(layer->GetWeightGradients(), learningRate / momentum);
      masterFirstLayer->UpdateBiasGradients(layer->GetBiasGradients(), learningRate / momentum);
   }

   masterFirstLayer->UpdateWeightGradients(masterFirstLayer->GetWeightGradients(), 1.0 - momentum);
   masterFirstLayer->UpdateBiasGradients(masterFirstLayer->GetBiasGradients(), 1.0 - momentum);

   for (size_t i = 0; i < depth; i++) {
      Layer_t *masterLayer = this->GetLayerAt(i);
      masterLayer->Update(1.0);

      for (size_t j = 0; j < nets.size(); j++) {
         Layer_t *layer = nets[j].GetLayerAt(i);

         layer->CopyWeights(masterLayer->GetWeights());
         layer->CopyBiases(masterLayer->GetBiases());
      }
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::ParallelBackwardNestorov(std::vector<TDeepNet<Architecture_t, Layer_t>> &nets,
                                                                 std::vector<TTensorBatch<Architecture_t>> &batches,
                                                                 Scalar_t learningRate, Scalar_t momentum) -> void
{
   std::cout << "Parallel Backward Nestorov" << std::endl;
   std::vector<Matrix_t> inp1;
   std::vector<Matrix_t> inp2;
   size_t depth = this->GetDepth();

   // Evaluate the gradients of the last layers in each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      evaluateGradients<Architecture_t>(nets[i].GetLayerAt(depth - 1)->GetActivationGradientsAt(0),
                                        nets[i].GetLossFunction(), batches[i].GetOutput(),
                                        nets[i].GetLayerAt(depth - 1)->GetOutputAt(0), batches[i].GetWeights());
   }

   // Backpropagate the error in i'th layer of each deep net
   for (size_t i = depth - 1; i > 0; i--) {
      for (size_t j = 0; j < nets.size(); j++) {
         Layer_t *layer = nets[j].GetLayerAt(i);

         layer->Backward(nets[j].GetLayerAt(i - 1)->GetActivationGradients(), nets[j].GetLayerAt(i - 1)->GetOutput(),
                         inp1, inp2);
      }
   }

   std::vector<Matrix_t> dummy;

   // First layer of each deep net
   for (size_t i = 0; i < nets.size(); i++) {
      Layer_t *layer = nets[i].GetLayerAt(0);
      layer->Backward(dummy, batches[i].GetInput(), inp1, inp2);
   }

   for (size_t i = 0; i < depth; i++) {
      Layer_t *masterLayer = this->GetLayerAt(i);
      for (size_t j = 0; j < nets.size(); j++) {
         Layer_t *layer = nets[j].GetLayerAt(i);

         layer->CopyWeights(masterLayer->GetWeights());
         layer->CopyBiases(masterLayer->GetBiases());

         layer->UpdateWeights(masterLayer->GetWeightGradients(), 1.0);
         layer->UpdateBiases(masterLayer->GetBiasGradients(), 1.0);
      }

      for (size_t j = 0; j < nets.size(); j++) {
         Layer_t *layer = nets[j].GetLayerAt(i);

         masterLayer->UpdateWeightGradients(layer->GetWeightGradients(), learningRate / momentum);
         masterLayer->UpdateBiasGradients(layer->GetBiasGradients(), learningRate / momentum);
      }

      masterLayer->UpdateWeightGradients(masterLayer->GetWeightGradients(), 1.0 - momentum);
      masterLayer->UpdateBiasGradients(masterLayer->GetBiasGradients(), 1.0 - momentum);

      masterLayer->Update(1.0);
   }
}
#endif   // use parallel deep net

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Update(Scalar_t learningRate) -> void
{
   for (size_t i = 0; i < fLayers.size(); i++) {
      fLayers[i]->Update(learningRate);
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Loss(const Matrix_t &groundTruth, const Matrix_t &weights,
                                             bool includeRegularization) const -> Scalar_t
{
   // Last layer should not be deep
   auto loss = evaluate<Architecture_t>(this->GetLossFunction(), groundTruth, fLayers.back()->GetOutputAt(0), weights);

   includeRegularization &= (this->GetRegularization() != ERegularization::kNone);
   if (includeRegularization) {
      loss += RegularizationTerm();
   }

   return loss;
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Loss(Tensor_t &input, const Matrix_t &groundTruth,
                                             const Matrix_t &weights, bool inTraining, bool includeRegularization)
   -> Scalar_t
{
   Forward(input, inTraining);
   return Loss(groundTruth, weights, includeRegularization);
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::RegularizationTerm() const -> Scalar_t
{
   Scalar_t reg = 0.0;
   for (size_t i = 0; i < fLayers.size(); i++) {
      for (size_t j = 0; j < (fLayers[i]->GetWeights()).size(); j++) {
         reg += regularization<Architecture_t>(fLayers[i]->GetWeightsAt(j), this->GetRegularization());
      }
   }
   return this->GetWeightDecay() * reg;
}


//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Prediction(Matrix_t &predictions, EOutputFunction f) const -> void
{
   // Last layer should not be deep (assume output is a matrix)
   evaluate<Architecture_t>(predictions, f, fLayers.back()->GetOutputAt(0));
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Prediction(Matrix_t &predictions, Tensor_t & input,
                                                   EOutputFunction f) -> void
{
   Forward(input, false);
   // Last layer should not be deep
   evaluate<Architecture_t>(predictions, f, fLayers.back()->GetOutputAt(0));
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
auto TDeepNet<Architecture_t, Layer_t>::Print() const -> void
{
   std::cout << "DEEP NEURAL NETWORK:   Depth = " << this->GetDepth();
   std::cout << "  Input = ( " << this->GetInputDepth();
   std::cout << ", " << this->GetInputHeight();
   std::cout << ", " << this->GetInputWidth() << " )";
   std::cout << "  Batch size = " << this->GetBatchSize();
   std::cout << "  Loss function = " << static_cast<char>(this->GetLossFunction()) << std::endl;

   //std::cout << "\t Layers: " << std::endl;

   for (size_t i = 0; i < fLayers.size(); i++) {
      std::cout << "\tLayer " << i << "\t";
      fLayers[i]->Print();
   }
}

//______________________________________________________________________________
template <typename Architecture_t, typename Layer_t>
void TDeepNet<Architecture_t, Layer_t>::SetDropoutProbabilities(
    const std::vector<Double_t> & probabilities)
{
   for (size_t i = 0; i < fLayers.size(); i++) {
      if (i < probabilities.size()) {
         fLayers[i]->SetDropoutProbability(probabilities[i]);
      } else {
         fLayers[i]->SetDropoutProbability(1.0);
      }
   }
}


} // namespace DNN
} // namespace TMVA

#endif
