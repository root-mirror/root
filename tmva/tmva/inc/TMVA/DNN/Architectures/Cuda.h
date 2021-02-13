// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 05/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////
// Definition of the TCuda architecture class, which provides an //
// implementation of the low-level functionality for neural      //
// networks for the CUDA computing architectures.                //
///////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_CUDA
#define TMVA_DNN_ARCHITECTURES_CUDA

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN/ContextHandles.h"
#include "TMVA/DNN/BatchNormLayer.h"
#include "TMVA/DNN/CNN/ConvLayer.h"
#include "TMVA/DNN/CNN/MaxPoolLayer.h"


#include "cuda.h"
#include "Cuda/CudaBuffers.h"
#include "Cuda/CudaMatrix.h"
#include "Cuda/CudaTensor.h"
#include "TMVA/DNN/DataLoader.h"
#include <utility>
#include <vector>
#include <string>

class TRandom;

namespace TMVA
{
namespace DNN
{
 struct CudaActivationDescriptor {};
 struct CudaFilterDescriptor {};
 struct CudaConvolutionDescriptor {};
 struct CudaDropoutDescriptor {};
 struct CudaPoolingDescriptor {};
 struct CudaConvolutionFwdAlgo {};
 struct CudaConvolutionBwdDataAlgo {};
 struct CudaConvolutionBwdFilterAlgo {};
 struct CudaDataType {};
 struct DummyCudaDataType {};

 struct CudaEmptyDescriptor {};

/** The TCuda architecture class.
 *
 * Low-level interface class for CUDA computing architectures. Contains as
 * public types the declaration of the scalar, matrix and buffer types
 * for this architecture as well as the remaining functions in the low-level
 * interface in the form of static members.
 */
template<typename AReal = Float_t>
class TCuda
{
private:
   static TRandom * fgRandomGen;
public:

   using AFloat         = AReal;
   using Scalar_t       = AFloat;

   using Matrix_t       = TCudaMatrix<AFloat>;
   using Tensor_t       = TCudaTensor<AFloat>;
   using DeviceBuffer_t = TCudaDeviceBuffer<AFloat>;
   using HostBuffer_t   = TCudaHostBuffer<AFloat>;

   using ActivationDescriptor_t  = CudaActivationDescriptor;
   using ConvolutionDescriptor_t = CudaConvolutionDescriptor;
   using FilterDescriptor_t      = CudaFilterDescriptor;
   using DropoutDescriptor_t     = CudaDropoutDescriptor;
   //using OpTensorDescriptor_t    = CudaOpTensorDescriptor;
   using PoolingDescriptor_t     = CudaPoolingDescriptor;
   using TensorDescriptor_t      = DummyCudaDataType;
   //using ReductionDescriptor_t   = CudaReduceTensorDescriptor;
   using AlgorithmForward_t      = CudaConvolutionFwdAlgo;
   using AlgorithmBackward_t     = CudaConvolutionBwdDataAlgo;
   using AlgorithmHelper_t       = CudaConvolutionBwdFilterAlgo;
   using AlgorithmDataType_t     = DummyCudaDataType;
   using ReduceTensorDescriptor_t = DummyCudaDataType;
   using RecurrentDescriptor_t    = DummyCudaDataType;

   using EmptyDescriptor_t       = CudaEmptyDescriptor;        // Used if a descriptor is not needed in a class

   using GenLayer_t              = VGeneralLayer<TCuda<AReal>>;
   using BNormLayer_t            = TBatchNormLayer<TCuda<AReal>>;
   using BNormDescriptors_t      = TDNNGenDescriptors<BNormLayer_t>;

   using ConvLayer_t             = CNN::TConvLayer<TCuda<AReal>>;
   using ConvDescriptors_t       = CNN::TCNNDescriptors<ConvLayer_t>;
   using ConvWorkspace_t         = CNN::TCNNWorkspace<ConvLayer_t>;
   using PoolingLayer_t          = CNN::TMaxPoolLayer<TCuda<AReal>>;
   using PoolingDescriptors_t    = CNN::TCNNDescriptors<PoolingLayer_t>;
   using PoolingWorkspace_t      = CNN::TCNNWorkspace<PoolingLayer_t>;

   using RNNDescriptors_t = RNN::TRNNDescriptors<TCuda<AReal>>;
   using RNNWorkspace_t = RNN::TRNNWorkspace<TCuda<AReal>>;

   static TMVA::Experimental::MemoryLayout GetTensorLayout() { return TMVA::Experimental::MemoryLayout::ColumnMajor; }

   static Tensor_t CreateTensor(size_t n, size_t c, size_t h, size_t w) {
      return Tensor_t( {c,h*w,n}, GetTensorLayout());
   }
   static Tensor_t CreateTensor(size_t b, size_t t, size_t w)
   {
      return Tensor_t( {t, w, b}, GetTensorLayout());
   }
   static Tensor_t CreateTensor(DeviceBuffer_t buffer, size_t n, size_t c, size_t h, size_t w) {
      return Tensor_t( buffer, {c,h*w, n}, GetTensorLayout(), 0, 0);
   }
   static Tensor_t CreateTensor(DeviceBuffer_t buffer, size_t b, size_t t, size_t w)
   {
      return Tensor_t(buffer, {t, w, b}, GetTensorLayout());
   }

   // create a weight tensor/matrix  from another tensor using its shape
   // static Matrix_t CreateWeightTensor( Matrix_t & A) {
   //    return Matrix_t( A.GetNrows(), A.GetNcols());
   // }
   // create a weight tensor/matrix vector   from another tensor/weight  vector using the given tensor shapes
   // this function is used by the optimizers to stgore intermidiate weights representations
   static void  CreateWeightTensors( std::vector<Matrix_t> & newWeights, const std::vector<Matrix_t> & weights) {
      if (!newWeights.empty()) newWeights.clear();
      size_t n =  weights.size();
      for (size_t i = 0; i < n; ++i)
         newWeights.emplace_back( weights[i].GetNrows(), weights[i].GetNcols());
   }

   static bool IsCudnn() { return false; }
   //____________________________________________________________________________
   //
   // Architecture Initialization
   //____________________________________________________________________________

   /** Initialize CNN data/operator descriptors. Not used at the moment.*/

   static void InitializeBNormDescriptors(TDescriptors * & /*descriptors*/,
                                          BNormLayer_t */*L = nullptr*/) {
      Error("InitializeBNormDescriptrs", "Batch normalization on GPU is supported only with Cudnn");
   }

   static void  InitializeConvDescriptors(TDescriptors *& /*descriptors*/, ConvLayer_t * /*L = nullptr*/) {}

   static void InitializePoolDescriptors(TDescriptors *& /*descriptors*/, PoolingLayer_t * /*L = nullptr*/) {}

   static void InitializeActivationDescriptor(ActivationDescriptor_t &/*descriptors*/, EActivationFunction /*activFunc */ , double /*coef*/ = 0.0) {}

   /** Release CNN data/operator descriptors. Not used at the moment.*/
   static void ReleaseConvDescriptors(TDescriptors * & /*descriptors*/) {}
   static void ReleasePoolDescriptors(TDescriptors * & /*descriptors*/) {}
   static void ReleaseBNormDescriptors(TDescriptors *& /*descriptors*/) {}

   static void InitializeConvWorkspace(TWorkspace * & /*workspace*/,
                                       TDescriptors * & /*descriptors*/,
                                       const DNN::CNN::TConvParams & /*params*/,
                                       ConvLayer_t */*L = nullptr*/) {}
   static void InitializePoolDropoutWorkspace(TWorkspace * & /*workspace*/,
                                       TDescriptors * & /*descriptors*/,
                                       const DNN::CNN::TConvParams & /*params*/,
                                       PoolingLayer_t */*L = nullptr*/) {}

   static void ReleaseDescriptor(ActivationDescriptor_t &  /*activationDescr*/) {}

   static void FreeConvWorkspace(TWorkspace * & /*workspace*/) {}     ///< Only used for certain cudnn on-device memory
   static void FreePoolDropoutWorkspace(TWorkspace * & /*workspace*/) {}

   static void InitializeRNNDescriptors(TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/) {}
   static void InitializeLSTMDescriptors(TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/) {}
   static void InitializeGRUDescriptors(TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/) {}

   static void InitializeRNNWorkspace(TWorkspace *& /*workspace*/, TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/){}
   static void InitializeLSTMWorkspace(TWorkspace *& /*workspace*/, TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/){}
   static void InitializeGRUWorkspace(TWorkspace *& /*workspace*/, TDescriptors *& /*descriptors*/, GenLayer_t * /*L*/){}

   static void InitializeRNNTensors(GenLayer_t * /*layer*/) {}
   static void InitializeLSTMTensors(GenLayer_t * /*layer*/) {}
   static void InitializeGRUTensors(GenLayer_t * /*layer*/) {}

   static void ReleaseRNNDescriptors(TDescriptors *& /*descriptors*/) {}
   static void FreeRNNWorkspace(TWorkspace *& /*workspace*/) {}

   //static void InitializeRNNTensors(RNNLayer_t * /*layer*/) {}

   //____________________________________________________________________________
   //
   // Propagation
   //____________________________________________________________________________

   /** @name Forward Propagation
    * Low-level functions required for the forward propagation of activations
    * through the network.
    */
      ///@{
   /** Matrix-multiply \p input with the transpose of \pweights and
    *  write the results into \p output. */
   static void MultiplyTranspose(Matrix_t &output, const Matrix_t &input, const Matrix_t &weights);

   static void MultiplyTranspose(Tensor_t &output, const Tensor_t &input, const Matrix_t &weights) {
      Matrix_t output_matrix = output.GetMatrix();
      MultiplyTranspose( output_matrix, input.GetMatrix(), weights);
         //ensor_t::MatrixToTensor(output_matrix, output); // this maybe is not needed
   }

   /** Add the vectors biases row-wise to the matrix output */
   static void AddRowWise(Matrix_t &output,const Matrix_t &biases);

   static void AddRowWise(Tensor_t &output, const Matrix_t &biases) {
      Matrix_t output_matrix = output.GetMatrix();
      AddRowWise(output_matrix, biases);
         //Tensor_t::MatrixToTensor(output_matrix, output); // this maybe is not needed
   }

   /** @name Backward Propagation (Dense Layers)
    * Low-level functions required for the forward propagation of activations
    * through the network.
    */
      ///@{
   /** Perform the complete backward propagation step. If the provided
    *  \p activationGradientsBackward matrix is not empty, compute the
    *  gradients of the objective function with respect to the activations
    *  of the previous layer (backward direction).
    *  Also compute the weight and the bias gradients. Modifies the values
    *  in \p df and thus produces only a valid result, if it is applied the
    *  first time after the corresponding forward propagation has been per-
    *  formed. */
   static void Backward(Tensor_t & activationGradientsBackward,
                        Matrix_t & weightGradients,
                        Matrix_t & biasGradients,
                        const Tensor_t & df,
                        const Tensor_t & activationGradients,
                        const Matrix_t & weights,
                        const Tensor_t & activationBackward);

   /** Adds a the elements in matrix B scaled by c to the elements in
    *  the matrix A. This is required for the weight update in the gradient
    *  descent step.*/
   static void ScaleAdd(Matrix_t & A,
                        const Matrix_t & B,
                        Scalar_t beta = 1.0);

   static void Copy(Matrix_t & B,
                    const Matrix_t & A);

      // copy from another type of matrix
   template<typename AMatrix_t>
   static void CopyDiffArch(Matrix_t & B, const AMatrix_t & A);


   /** Above functions extended to vectors */
   static void ScaleAdd(Tensor_t & A,
                        const Tensor_t & B,
                        Scalar_t beta = 1.0);

   static void Copy(Tensor_t & A,
                    const Tensor_t & B);

      // copy from another tensor
   template<typename ATensor_t>
   static void CopyDiffArch(Tensor_t & A,
                            const ATensor_t & B);

      // copy from vector of matrices of different types
   template<typename AMatrix_t>
   static void CopyDiffArch(std::vector<Matrix_t>  & A,
                            const std::vector<AMatrix_t> & B);

      ///@}

      //____________________________________________________________________________
      //
      // Activation Functions
      //____________________________________________________________________________

   /** @name Activation Functions
    * For each activation function, the low-level interface contains two routines.
    * One that applies the acitvation function to a matrix and one that evaluate
    * the derivatives of the activation function at the elements of a given matrix
    * and writes the results into the result matrix.
    */
      ///@{
         /*  impl using Matrix */
   /*inline void evaluate(Matrix_t &A, EActivationFunction f)
   {
    Tensor_t tA(A);
    evaluate<TCuda<AReal>>(tA,f);
   }*/
    static void ActivationFunctionForward(Tensor_t & X, EActivationFunction activFunct,
                          const ActivationDescriptor_t activationDescr,
                          const double coef = 0.0, const AFloat alpha = 1,
                          const AFloat beta = 0);

   /** Computes the gradient of the activation function */
   static void ActivationFunctionBackward(Tensor_t & dX, const Tensor_t & Y,
                                          const Tensor_t & dY,  const Tensor_t & X,
                                          EActivationFunction activFunct,
                                          const ActivationDescriptor_t activationDescr,
                                          const AFloat alpha = 1,
                                          const AFloat beta = 0);

   static void IdentityDerivative(Tensor_t & B,
                                  const Tensor_t &A);

   static void Relu(Tensor_t & B);
   static void ReluDerivative(Tensor_t & B,
                              const Tensor_t & A);

   static void Sigmoid(Tensor_t & B);
   static void SigmoidDerivative(Tensor_t & B,
                                 const Tensor_t & A);

   static void Tanh(Tensor_t & B);
   static void TanhDerivative(Tensor_t & B,
                              const Tensor_t & A);

   static void FastTanh(Tensor_t &B) { return Tanh(B);  }
   static void FastTanhDerivative(Tensor_t &B, const Tensor_t &A) { return TanhDerivative(B, A); }

   static void SymmetricRelu(Tensor_t & B);
   static void SymmetricReluDerivative(Tensor_t & B,
                                       const Tensor_t & A);

   static void SoftSign(Tensor_t & B);
   static void SoftSignDerivative(Tensor_t & B,
                                  const Tensor_t & A);

   static void Gauss(Tensor_t & B);
   static void GaussDerivative(Tensor_t & B,
                               const Tensor_t & A);
      ///@}

      //____________________________________________________________________________
      //
      // Loss Functions
      //____________________________________________________________________________

   /** @name Loss Functions
    * Loss functions compute a scalar value given the \p output of the network
    * for a given training input and the expected network prediction \p Y that
    * quantifies the quality of the prediction. For each function also a routing
    * that computes the gradients (suffixed by Gradients) must be provided for
    * the starting of the backpropagation algorithm.
    */
      ///@{

   static Scalar_t MeanSquaredError(const Matrix_t &Y, const Matrix_t &output,
                                    const Matrix_t &weights);
   static void MeanSquaredErrorGradients(Matrix_t &dY, const Matrix_t &Y,
                                         const Matrix_t &output, const Matrix_t &weights);

   /** Sigmoid transformation is implicitly applied, thus \p output should
    *  hold the linear activations of the last layer in the net. */
   static Scalar_t CrossEntropy(const Matrix_t &Y, const Matrix_t &output,
                                const Matrix_t &weights);

   static void CrossEntropyGradients(Matrix_t &dY, const Matrix_t &Y,
                                     const Matrix_t &output, const Matrix_t &weights);

   /** Softmax transformation is implicitly applied, thus \p output should
    *  hold the linear activations of the last layer in the net. */
   static Scalar_t SoftmaxCrossEntropy(const Matrix_t &Y, const Matrix_t &output,
                                       const Matrix_t &weights);
   static void SoftmaxCrossEntropyGradients(Matrix_t &dY, const Matrix_t &Y,
                                            const Matrix_t &output, const Matrix_t &weights);
      ///@}

      //____________________________________________________________________________
      //
      // Output Functions
      //____________________________________________________________________________

   /** @name Output Functions
    * Output functions transform the activations \p output of the
    * output layer in the network to a valid prediction \p YHat for
    * the desired usage of the network, e.g.  the identity function
    * for regression or the sigmoid transformation for two-class
    * classification.
    */
      ///@{
   static void Sigmoid(Matrix_t &YHat,
                       const Matrix_t & );
   static void Softmax(Matrix_t &YHat,
                       const Matrix_t & );
      ///@}

      //____________________________________________________________________________
      //
      // Regularization
      //____________________________________________________________________________

   /** @name Regularization
    * For each regularization type two functions are required, one named
    * <tt><Type>Regularization</tt> that evaluates the corresponding
    * regularization functional for a given weight matrix and the
    * <tt>Add<Type>RegularizationGradients</tt>, that adds the regularization
    * component in the gradients to the provided matrix.
    */
      ///@{

   static Scalar_t L1Regularization(const Matrix_t & W);
   static void AddL1RegularizationGradients(Matrix_t & A,
                                            const Matrix_t & W,
                                            Scalar_t weightDecay);

   static Scalar_t L2Regularization(const Matrix_t & W);
   static void AddL2RegularizationGradients(Matrix_t & A,
                                            const Matrix_t & W,
                                            Scalar_t weightDecay);
      ///@}

      //____________________________________________________________________________
      //
      // Initialization
      //____________________________________________________________________________

   /** @name Initialization
    * For each initialization method, one function in the low-level interface
    * is provided. The naming scheme is <p>Initialize<Type></p> for a given
    * initialization method Type.
    */
      ///@{

   static void InitializeGauss(Matrix_t & A);
   static void InitializeUniform(Matrix_t & A);
   static void InitializeIdentity(Matrix_t & A);
   static void InitializeZero(Matrix_t & A);
   static void InitializeZero(Tensor_t &A);
   static void InitializeGlorotNormal(Matrix_t & A);
   static void InitializeGlorotUniform(Matrix_t & A);

      // return static instance of random generator used for initialization
      // if generator does not exist it is created the first time with a random seed (e.g. seed = 0)
   static TRandom & GetRandomGenerator();
      // set random seed for the static geenrator
      // if the static geneerator does not exists it is created
   static void SetRandomSeed(size_t seed);
      ///@}

      //____________________________________________________________________________
      //
      // Dropout
      //____________________________________________________________________________

   /** @name Dropout
    */
      ///@{

   /** Apply dropout with activation probability \p p to the given
    *  tensor \p A and scale the result by reciprocal of \p p. */
   static void DropoutForward(Tensor_t & A,
                              TDescriptors * descriptors,
                              TWorkspace   * workspace,
                              Scalar_t p);

   static void DropoutForward(Matrix_t & A, Scalar_t p) {
      Tensor_t tA(A);
      DropoutForward( tA, static_cast<TDescriptors *> (nullptr), static_cast<TWorkspace *> (nullptr), p );
   }

   static void DropoutBackward(Tensor_t & /* A */,
                               TDescriptors * /*descriptors */,
                               TWorkspace   * /* workspace */ ) {}
      ///@}

   //____________________________________________________________________________
   //
   // Batch Normalization
   //____________________________________________________________________________

   /** @name Batch Normalization Layer Propagation
    */
   ///@{

   /** The input from each batch are normalized during training to have zero mean and unit variance
     * and they are then scaled by two parameter, different for each input variable:
     *  - a scale factor \gamma gamma
     *  - an offset \beta beta */

   static void BatchNormLayerForwardTraining(int axis, const Tensor_t &x, Tensor_t &y, Matrix_t &gamma, Matrix_t &beta,
                                             Matrix_t &mean, Matrix_t &, Matrix_t &iVariance, Matrix_t &runningMeans,
                                             Matrix_t &runningVars, Scalar_t nTrainedBatches, Scalar_t momentum,
                                             Scalar_t epsilon, const TensorDescriptor_t &bnParDescriptor);

   /** During inference the inputs are not normalized using the batch mean but the previously computed
    * at  running mean and variance */

   static void BatchNormLayerForwardInference(int axis, const Tensor_t &x, Matrix_t &gamma, Matrix_t &beta, Tensor_t &y,
                                              const Matrix_t &runningMeans, const Matrix_t &runningVars,
                                              Scalar_t epsilon, const TensorDescriptor_t &);

   static void BatchNormLayerBackward(int axis, const Tensor_t &x, const Tensor_t &dy, Tensor_t &dx,
                                      Matrix_t &gamma, //  Matrix_t &beta, (not needed)
                                      Matrix_t &dgamma, Matrix_t &dbeta, const Matrix_t &mean, const Matrix_t &variance,
                                      const Matrix_t &iVariance, Scalar_t epsilon, const TensorDescriptor_t &);

   //____________________________________________________________________________
   //
   //  Convolutional Layer Propagation
   //____________________________________________________________________________

   /** @name Forward Propagation in Convolutional Layer
    */
      ///@{

   /** Calculate how many neurons "fit" in the output layer, given the input as well as the layer's hyperparameters. */
   static size_t calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride);

   /** Transform the matrix B in local view format, suitable for
    *  convolution, and store it in matrix A */
   static void Im2col(Matrix_t &A,
                      const Matrix_t &B,
                      size_t imgHeight,
                      size_t imgWidth,
                      size_t fltHeight,
                      size_t fltWidth,
                      size_t strideRows,
                      size_t strideCols,
                      size_t zeroPaddingHeight,
                      size_t zeroPaddingWidth);

   static void Im2colIndices(std::vector<int> &V, const Matrix_t &B, size_t nLocalViews, size_t imgHeight, size_t imgWidth, size_t fltHeight,
                             size_t fltWidth, size_t strideRows, size_t strideCols, size_t zeroPaddingHeight,
                             size_t zeroPaddingWidth);
   static void Im2colFast(Matrix_t &A, const Matrix_t &B, const std::vector<int> & V);

   /** Rotates the matrix \p B, which is representing a weights,
    *  and stores them in the matrix \p A. */
   static void RotateWeights(Matrix_t &A, const Matrix_t &B, size_t filterDepth, size_t filterHeight,
                             size_t filterWidth, size_t numFilters);

   /** Add the biases in the Convolutional Layer.  */
   static void AddConvBiases(Matrix_t &output, const Matrix_t &biases);
      ///@}

   /** Dummy placeholder - preparation is currently only required for the CUDA architecture. */
   static void PrepareInternals(Tensor_t &) {}

   /** Forward propagation in the Convolutional layer */
   static void ConvLayerForward(Tensor_t & output,
                                Tensor_t & inputActivationFunc,
                                const Tensor_t &input,
                                const Matrix_t &weights, const Matrix_t & biases,
                                const DNN::CNN::TConvParams & params, EActivationFunction activFunc,
                                Tensor_t & /* inputPrime */,
                                const ConvDescriptors_t & /*descriptors*/,   // Empty struct for cuda architecture
                                ConvWorkspace_t & /*workspace*/);      // Empty struct for cuda architecture
                                //void * cudnnWorkspace = nullptr);          // Remains nullptr for cuda architecture
   /** @name Backward Propagation in Convolutional Layer
    */
      ///@{

   /** Perform the complete backward propagation step in a Convolutional Layer.
    *  If the provided \p activationGradientsBackward matrix is not empty, compute the
    *  gradients of the objective function with respect to the activations
    *  of the previous layer (backward direction).
    *  Also compute the weight and the bias gradients. Modifies the values
    *  in \p df and thus produces only a valid result, if it is applied the
    *  first time after the corresponding forward propagation has been per-
    *  formed. */
   static void ConvLayerBackward(Tensor_t &activationGradientsBackward,
                                 Matrix_t &weightGradients, Matrix_t &biasGradients,
                                 Tensor_t &df,
                                 Tensor_t &activationGradients,
                                 const Matrix_t &weights,
                                 const Tensor_t &activationBackward,
                                 const Tensor_t & outputTensor,
                                 EActivationFunction activFunc,
                                 const ConvDescriptors_t & /*descriptors*/,
                                 ConvWorkspace_t & /*workspace*/,
                                 size_t batchSize,   size_t inputHeight,
                                 size_t inputWidth,  size_t depth,
                                 size_t height,      size_t width,
                                 size_t filterDepth, size_t filterHeight,
                                 size_t filterWidth, size_t nLocalViews );

   /** Utility function for calculating the activation gradients of the layer
    *  before the convolutional layer. */
   static void CalculateConvActivationGradients(Tensor_t &activationGradientsBackward,
                                                const Tensor_t &df,
                                                const Matrix_t &weights, size_t batchSize,
                                                size_t inputHeight, size_t inputWidth, size_t depth, size_t height,
                                                size_t width, size_t filterDepth, size_t filterHeight,
                                                size_t filterWidth);

   /** Utility function for calculating the weight gradients of the convolutional
    * layer. */
   static void CalculateConvWeightGradients(Matrix_t &weightGradients,
                                            const Tensor_t &df,
                                            const Tensor_t &activations_backward,
                                            size_t batchSize, size_t inputHeight, size_t inputWidth, size_t depth,
                                            size_t height, size_t width, size_t filterDepth, size_t filterHeight,
                                            size_t filterWidth, size_t nLocalViews);

   /** Utility function for calculating the bias gradients of the convolutional
    *  layer */
   static void CalculateConvBiasGradients(Matrix_t &biasGradients, const Tensor_t &df,
                                          size_t batchSize, size_t depth, size_t nLocalViews);
      ///@}

      //____________________________________________________________________________
      //
      //  Max Pooling Layer Propagation
      //____________________________________________________________________________
   /** @name Forward Propagation in Max Pooling Layer
    */
      ///@{

   /** Downsample the matrix \p C to the matrix \p A, using max
    * operation, such that the winning indices are stored in matrix
    * \p B. */
   static void Downsample(Tensor_t &A, Tensor_t &B, const Tensor_t &C,
                          const PoolingDescriptors_t & /*descriptors*/,
                          PoolingWorkspace_t & /*workspace*/,
                          size_t imgHeight, size_t imgWidth, size_t fltHeight,
                          size_t fltWidth, size_t strideRows, size_t strideCols);

      ///@}

   /** @name Backward Propagation in Max Pooling Layer
    */
      ///@{
   /** Perform the complete backward propagation step in a Pooling Layer. Based on the
    *  winning idices stored in the index matrix, it just forwards the actiovation
    *  gradients to the previous layer. */
   static void MaxPoolLayerBackward(Tensor_t &activationGradientsBackward,
                                    const Tensor_t &activationGradients,
                                    const Tensor_t &indexMatrix,
                                    const Tensor_t & /*inputActivation*/,
                                    const Tensor_t & /*outputTensor*/,
                                    const PoolingDescriptors_t & /*descriptors*/,
                                    PoolingWorkspace_t & /*workspace*/,
                                    size_t imgHeight,
                                    size_t imgWidth,
                                    size_t fltHeight,
                                    size_t fltWidth,
                                    size_t strideRows,
                                    size_t strideCols,
                                    size_t nLocalViews);

      ///@}

      //____________________________________________________________________________
      //
      //  Reshape Layer Propagation
      //____________________________________________________________________________
   /** @name Forward and Backward Propagation in Reshape Layer
    */
      ///@{

   /** Transform the matrix \p B to a matrix with different dimensions \p A */
   static void Reshape(Matrix_t &A, const Matrix_t &B);

   /** Flattens the tensor \p B, such that each matrix, is stretched in
    *  one row, resulting with a matrix \p A. */
   static void Flatten(Tensor_t &A, const Tensor_t &B); // size_t size, size_t nRows, size_t nCols);

   /** Transforms each row of \p B to a matrix and stores it in the
    *  tensor \p B. */
   static void Deflatten(Tensor_t &A, const Tensor_t &B); // size_t index, size_t nRows,size_t nCols);

   /** Rearrage data accoring to time fill B x T x D out with T x B x D matrix in*/
   static void Rearrange(Tensor_t &out, const Tensor_t &in);


   /** Backward pass for Recurrent Networks */
   static Matrix_t & RecurrentLayerBackward(Matrix_t & state_gradients_backward, // BxH
                                            Matrix_t & input_weight_gradients,
                                            Matrix_t & state_weight_gradients,
                                            Matrix_t & bias_gradients,
                                            Matrix_t & df, //DxH
                                            const Matrix_t & state, // BxH
                                            const Matrix_t & weights_input, // HxD
                                            const Matrix_t & weights_state, // HxH
                                            const Matrix_t & input,  // BxD
                                            Matrix_t & input_gradient);

   // dummy RNN functions
   static void RNNForward(const Tensor_t & /* x */, const Matrix_t & /* hx */, const Matrix_t & /* cx */,
                          const Matrix_t & /* weights */, Tensor_t & /* y */, Matrix_t & /* hy */, Matrix_t & /* cy */,
                          const RNNDescriptors_t & /* descr */, RNNWorkspace_t & /* workspace */, bool /* isTraining */)
   {
   }

   static void RNNBackward(const Tensor_t & /* x */, const Matrix_t & /* hx */, const Matrix_t & /* cx */,
                           const Tensor_t & /* y */, const Tensor_t & /* dy */, const Matrix_t & /* dhy */,
                           const Matrix_t & /* dcy */, const Tensor_t & /* weights */, Tensor_t & /* dx */,
                           Matrix_t & /* dhx */, Matrix_t & /* dcx */, Tensor_t & /* dw */,
                           const RNNDescriptors_t & /* desc */, RNNWorkspace_t & /* workspace */)
   {
   }
   static Matrix_t &
   LSTMLayerBackward(Matrix_t &state_gradients_backward, Matrix_t & /*cell_gradients_backward*/,
                     Matrix_t & /*input_weight_gradients*/, Matrix_t & /*forget_weight_gradients*/,
                     Matrix_t & /*candidate_weight_gradients*/, Matrix_t & /*output_weight_gradients*/,
                     Matrix_t & /*input_state_weight_gradients*/, Matrix_t & /*forget_state_weight_gradients*/,
                     Matrix_t & /*candidate_state_weight_gradients*/, Matrix_t & /*output_state_weight_gradients*/,
                     Matrix_t & /*input_bias_gradients*/, Matrix_t & /*forget_bias_gradients*/,
                     Matrix_t & /*candidate_bias_gradients*/, Matrix_t & /*output_bias_gradients*/, Matrix_t & /*di*/,
                     Matrix_t & /*df*/, Matrix_t & /*dc*/, Matrix_t & /*dout*/,
                     const Matrix_t & /*precStateActivations*/, const Matrix_t & /*precCellActivations*/,
                     const Matrix_t & /*fInput*/, const Matrix_t & /*fForget*/, const Matrix_t & /*fCandidate*/,
                     const Matrix_t & /*fOutput*/, const Matrix_t & /*weights_input*/,
                     const Matrix_t & /*weights_forget*/, const Matrix_t & /*weights_candidate*/,
                     const Matrix_t & /*weights_output*/, const Matrix_t & /*weights_input_state*/,
                     const Matrix_t & /*weights_forget_state*/, const Matrix_t & /*weights_candidate_state*/,
                     const Matrix_t & /*weights_output_state*/, const Matrix_t & /*input*/,
                     Matrix_t & /*input_gradient*/, Matrix_t & /*cell_gradient*/, Matrix_t & /*cell_tanh*/)
   {
      Fatal("TCuda::LSTMLayerBackward", "Recurrent layers are not supported in the native Cuda architecture!!!");
      return state_gradients_backward;
   }

   /** Backward pass for GRU Network */
   static Matrix_t &
   GRULayerBackward(Matrix_t &state_gradients_backward, Matrix_t & /*reset_weight_gradients*/,
                    Matrix_t & /*update_weight_gradients*/, Matrix_t & /*candidate_weight_gradients*/,
                    Matrix_t & /*reset_state_weight_gradients*/, Matrix_t & /*update_state_weight_gradients*/,
                    Matrix_t & /*candidate_state_weight_gradients*/, Matrix_t & /*reset_bias_gradients*/,
                    Matrix_t & /*update_bias_gradients*/, Matrix_t & /*candidate_bias_gradients*/, Matrix_t & /*dr*/,
                    Matrix_t & /*du*/, Matrix_t & /*dc*/, const Matrix_t & /*precStateActivations*/,
                    const Matrix_t & /*fReset*/, const Matrix_t & /*fUpdate*/, const Matrix_t & /*fCandidate*/,
                    const Matrix_t & /*weights_reset*/, const Matrix_t & /*weights_update*/,
                    const Matrix_t & /*weights_candidate*/, const Matrix_t & /*weights_reset_state*/,
                    const Matrix_t & /*weights_update_state*/, const Matrix_t & /*weights_candidate_state*/,
                    const Matrix_t & /*input*/, Matrix_t & /*input_gradient*/, bool)
   {
      Fatal("TCuda::GRULayerBackward", "Recurrent layers are not supported in the native Cuda architecture!!!");
      return state_gradients_backward;
   }
   ///@}

   //____________________________________________________________________________
   //
   // Additional Arithmetic Functions
   //____________________________________________________________________________

   /** @name Additional Arithmetic Functions
    *
    * Additional arithmetic on CUDA matrices  used to implement the low-level
    * interface.
    */
      ///@{

   /** Standard multiplication of two matrices \p A and \p B with the result being
    *  written into C.
    */
   static void Multiply(Matrix_t &C,
                        const Matrix_t &A,
                        const Matrix_t &B);
   /** Matrix multiplication of two matrices \p A and \p B^T (transposed) with the
    *  result being written into C.
    */
   static void TransposeMultiply(Matrix_t &output,
                                 const Matrix_t &input,
                                 const Matrix_t &Weights,
                                 Scalar_t alpha = 1.0, Scalar_t beta = 0.);
   /** In-place Hadamard (element-wise) product of matrices \p A and \p B
    *  with the result being written into \p A.
    */
   static void Hadamard(Tensor_t &A,
                        const Tensor_t &B);
   static void Hadamard(Matrix_t &A,
                        const Matrix_t &B);
      // {
      //    Tensor_t tA(A);
      //    Hadamard( tA, Tensor_t(B));
      // }

   /** Sum columns of (m x n) matrixx \p A and write the results into the first
    * m elements in \p A.
    */
   static void SumColumns(Matrix_t &B,
                          const Matrix_t &A,
                          Scalar_t alpha = 1.0, Scalar_t beta = 0.);

   /** Compute the sum of all elements in \p A */
   static Scalar_t Sum(const Matrix_t &A);

   /** Check two matrices for equality, taking floating point arithmetic errors into account. */
   static bool AlmostEquals(const Matrix_t &A, const Matrix_t &B, double epsilon = 0.1);

   /** Add the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstAdd(Matrix_t &A, Scalar_t beta);

   /** Multiply the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstMult(Matrix_t &A, Scalar_t beta);

   /** Reciprocal each element of the matrix \p A and write the result into
    * \p A
    */
   static void ReciprocalElementWise(Matrix_t &A);

   /** Square each element of the matrix \p A and write the result into
    * \p A
    */
   static void SquareElementWise(Matrix_t &A);

   /** Square root each element of the matrix \p A and write the result into
    * \p A
    */
   static void SqrtElementWise(Matrix_t &A);

      // optimizer functions
   static void AdamUpdate(Matrix_t & A, const Matrix_t & M, const Matrix_t & V, Scalar_t alpha, Scalar_t eps);
   static void AdamUpdateFirstMom(Matrix_t & A, const Matrix_t & B, Scalar_t beta);
   static void AdamUpdateSecondMom(Matrix_t & A, const Matrix_t & B, Scalar_t beta);

      // printing of tensor
   static void PrintTensor( const Tensor_t & A, const std::string name = "Cuda-tensor", bool = false);

   ///////////////////////////////////////////////////////////////////////////////
   /// extra functions defined only for CPU architecture !!!
   //////////////////////////////////////////////////////////////////////////////

   /** Sum rows of (m x n) matrix \p A and write the results into the first
   * m elements in \p B.
   */
   static void SumRows(Matrix_t & B, const Matrix_t & A);


};

//____________________________________________________________________________
template <typename AFloat>
template <typename AMatrix_t>
void TCuda<AFloat>::CopyDiffArch(TCudaMatrix<AFloat> &B,
                        const AMatrix_t &A)
{
   // copy from another architecture using the reference one
   // this is not very efficient since creates temporary objects
   TMatrixT<AFloat> tmp = A;
   Copy(B, TCudaMatrix<AFloat>(tmp) );
}

//____________________________________________________________________________
template <typename AFloat>
template <typename AMatrix_t>
void TCuda<AFloat>::CopyDiffArch(std::vector<TCudaMatrix<AFloat>> &B,
                            const std::vector<AMatrix_t> &A)
{
   for (size_t i = 0; i < B.size(); ++i) {
      CopyDiffArch(B[i], A[i]);
   }
}

template <typename AFloat>
void TCuda<AFloat>::PrintTensor(const typename TCuda<AFloat>::Tensor_t & A, const std::string name, bool  )
{
   std::cout << name << "  size = " << A.GetSize() << " shape = { ";
   auto shape = A.GetShape();
   for (size_t k = 0; k < shape.size()-1; ++k)
      std::cout << shape[k] << " , ";
   std::cout << shape.back() << " } ";
   std::cout << " strides = { ";
   auto strides = A.GetStrides();
   for (size_t k = 0; k < strides.size()-1; ++k)
      std::cout << strides[k] << " , ";
   std::cout << strides.back() << " }\n ";

   if (A.GetShape().size() == 2 ) {
      for (size_t i = 0; i < A.GetShape()[0]; ++i) {
         std::cout << "{ ";
         for (size_t j = 0; j < A.GetShape()[1]; ++j) {
            std::cout << A(i,j) << " ";
         }
         std::cout << " } " << std::endl;
      }
   } else if  (A.GetShape().size() == 3 ) {
      for (size_t i = 0; i < A.GetFirstSize(); ++i) {
         std::cout << "{ ";
         for (size_t j = 0; j < A.GetHSize(); ++j) {
            std::cout << "{ ";
            for (size_t k = 0; k < A.GetWSize(); ++k) {
               std::cout << A(i,j,k) << " ";
            }
            std::cout << " } " << std::endl;
         }
         std::cout << " } " << std::endl;
      }
   }
   else {
      for (size_t l = 0; l < A.GetSize(); ++l) {
         std::cout << A.GetData()[l] << " ";
      }
      std::cout << "\n";
   }
}


} // namespace DNN
} // namespace TMVA

#endif
