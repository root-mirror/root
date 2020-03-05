// @(#)root/tmva/tmva/dnn:$Id$
// Author: Joana Niermann 23/07/19

/*************************************************************************
 * Copyright (C) 2019, Joana Niermann                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////
// Definition of the TCudnn architecture class, which provides   //
// a wrapping of the low-level functionality for neural networks //
// in the cuDNN library.                                         //
///////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_CUDNN
#define TMVA_DNN_ARCHITECTURES_CUDNN

#include "RConfigure.h"   // for definition of R__HAS_CUDNN

#ifndef R__HAS_CUDNN
#error This file can be compiled only when cudnn is available in ROOT
#else 

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN/ContextHandles.h"
//#include "TMVA/DNN/CNN/Descriptors.h"
#include "TMVA/DNN/BatchNormLayer.h"
#include "TMVA/DNN/CNN/ConvLayer.h"
#include "TMVA/DNN/CNN/MaxPoolLayer.h"

#include "cudnn.h"
#include "Cuda/CudaBuffers.h"
#include "Cuda/CudaTensor.h"
#include "TMVA/DNN/TensorDataLoader.h"
#include <utility>
#include <vector>

#include "TMVA/DNN/Architectures/Cuda.h"

class TRandom;

namespace TMVA
{
namespace DNN
{

struct TCudnnEmptyDescriptor {};


/** The TCudnn architecture class.
 *
 * Low-level interface class for CUDA computing architectures using the cuDNN
 * library as backend. Contains as public types the declaration of the scalar,
 * matrix and buffer types for this architecture, as well as the remaining
 * functions in the low-level interface in the form of static members.
 */
template<typename AFloat = Float_t>
class TCudnn
{
private:
   static TRandom * fgRandomGen;
public:

   using Scalar_t       = AFloat;
   using Matrix_t       = TCudaTensor<AFloat>;
   using Tensor_t       = TCudaTensor<AFloat>;
   using DeviceBuffer_t = TCudaDeviceBuffer<AFloat>;
   using HostBuffer_t   = TCudaHostBuffer<AFloat>;

   // The descriptors for the (tensor) data are held by the data classes (CudaTensor)
   using ActivationDescriptor_t  = cudnnActivationDescriptor_t;
   using ConvolutionDescriptor_t = cudnnConvolutionDescriptor_t;
   using DropoutDescriptor_t     = cudnnDropoutDescriptor_t;
   using FilterDescriptor_t      = cudnnFilterDescriptor_t;
   //using OpTensorDescriptor_t    = cudnnOpTensorDescriptor_t;
   using PoolingDescriptor_t     = cudnnPoolingDescriptor_t;
   //using ReductionDescriptor_t   = cudnnReduceTensorDescriptor_t;
   using AlgorithmForward_t      = cudnnConvolutionFwdAlgo_t;
   using AlgorithmBackward_t     = cudnnConvolutionBwdDataAlgo_t;
   using AlgorithmHelper_t       = cudnnConvolutionBwdFilterAlgo_t;
   using AlgorithmDataType_t     = cudnnDataType_t;
   using ReduceTensorDescriptor_t = cudnnReduceTensorDescriptor_t;
   using TensorDescriptor_t       = cudnnTensorDescriptor_t;

   using EmptyDescriptor_t       = TCudnnEmptyDescriptor;        // Used if a descriptor is not needed in a class

   using BNormLayer_t            = TBatchNormLayer<TCudnn<AFloat>>;
   using BNormDescriptors_t      = TDNNGenDescriptors<BNormLayer_t>;
   //using BNormWorkspace_t        = CNN::TCNNWorkspace<BNormLayer_t>;*/
   using ConvLayer_t             = CNN::TConvLayer<TCudnn<AFloat>>;
   using ConvDescriptors_t       = CNN::TCNNDescriptors<ConvLayer_t>;
   using ConvWorkspace_t         = CNN::TCNNWorkspace<ConvLayer_t>;
   using PoolingLayer_t          = CNN::TMaxPoolLayer<TCudnn<AFloat>>;
   using PoolingDescriptors_t    = CNN::TCNNDescriptors<PoolingLayer_t>;
   using PoolingWorkspace_t      = CNN::TCNNWorkspace<PoolingLayer_t>;

   // template <typename AFloat>
   // using ConvDescriptors_t = CNN::TCNNDescriptors<CNN::TConvLayer<TCudnn<AFloat>>>;

   // convolution options
   // default is -1 (left to cudnn)
   struct CNNOptions  {

      static int ConvFwdAlgorithm;
      static int ConvBwdDataAlgorithm;
      static int ConvBwdFilterAlgorithm;
      // default is 0 (left to cudnn : a value -1 will indicate to not use any space)
      static Long_t ConvMaxWorkspaceSize;
   }; // namespace DNN

   static TMVA::Experimental::MemoryLayout GetTensorLayout() { return TMVA::Experimental::MemoryLayout::RowMajor; }


   static Tensor_t CreateTensor(size_t n, size_t c, size_t h, size_t w) {
      return Tensor_t( {n,c,h,w}, GetTensorLayout(), 0, 0);
   }

   static Tensor_t CreateTensor(DeviceBuffer_t buffer, size_t n, size_t c, size_t h, size_t w) {
      return Tensor_t( buffer, {n,c,h,w}, GetTensorLayout(), 0, 0);
   }

   // create a weight tensor/matrix vector   from another tensor/weight  vector using the given tensor shapes
   // this function is used by the optimizers to stgore intermidiate weights representations
   static void  CreateWeightTensors( std::vector<Matrix_t> & newWeights, const std::vector<Matrix_t> & weights) {
      if (!newWeights.empty()) newWeights.clear();
      size_t n =  weights.size();
      for (size_t i = 0; i < n; ++i)
         newWeights.emplace_back( weights[i].GetShape(), weights[i].GetLayout(), 0, 0);
   }
   //____________________________________________________________________________
   //
   // Architecture Initialization
   //____________________________________________________________________________

   static void InitializeBNormDescriptors(TDescriptors * & descriptors,
                                          BNormLayer_t *L = nullptr);

   static void InitializeConvDescriptors(TDescriptors * & descriptors,
                                         ConvLayer_t *L = nullptr);

   static void InitializePoolDescriptors(TDescriptors * & descriptors,
                                        PoolingLayer_t *L = nullptr);

   static void InitializeActivationDescriptor(ActivationDescriptor_t & descriptors, EActivationFunction activFunc, double coef = 0.0);

   static void ReleaseConvDescriptors(TDescriptors    * descriptors );
   static void ReleasePoolDescriptors(TDescriptors * descriptors );
   static void ReleaseBNormDescriptors(TDescriptors * descriptors );
   static void ReleaseDescriptor(EmptyDescriptor_t       & emptyDescr) {}        // Does nothing
   static void ReleaseDescriptor(ActivationDescriptor_t  & activationDescr);
   static void ReleaseDescriptor(ConvolutionDescriptor_t & convolutionDescr);
   static void ReleaseDescriptor(DropoutDescriptor_t     & dropoutDescr);
   static void ReleaseDescriptor(FilterDescriptor_t      & filterDescr);
   static void ReleaseDescriptor(PoolingDescriptor_t     & poolingDescr);
   static void ReleaseDescriptor(TensorDescriptor_t      & tensorDescr);


   static void InitializeConvWorkspace(TWorkspace * & workspace,
                                       TDescriptors * & descriptors,
                                       const DNN::CNN::TConvParams & params,
                                       ConvLayer_t *L = nullptr);
   static void InitializePoolDropoutWorkspace(TWorkspace * & workspace,
                                       TDescriptors * & descriptors,
                                       const DNN::CNN::TConvParams & params,
                                       PoolingLayer_t *L = nullptr);

   static void FreeConvWorkspace(TWorkspace * workspace, ConvLayer_t *L = nullptr);
   static void FreePoolDropoutWorkspace(TWorkspace * workspace, PoolingLayer_t *L = nullptr);
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
   static void MultiplyTranspose(Tensor_t &output, const Tensor_t &input, const Matrix_t &weights);

   /** Add the vectors biases row-wise to the matrix output */
   static void AddRowWise(Tensor_t &output,const Matrix_t &biases);

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
                        Tensor_t & df,
                        const Tensor_t & activationGradients,
                        const Matrix_t & weights,
                        const Tensor_t & activationBackward);

   /** Above functions extended to vectors */
   static void ScaleAdd(Tensor_t & A, const Tensor_t & B,
                        Scalar_t alpha = 1.0,
                        Scalar_t beta = 1.0);

   /** Deep copy from B to A. */
   static void Copy(Tensor_t & A, const Tensor_t & B);

   // copy from another tensor
   template<typename ATensor_t>
   static void CopyDiffArch(Tensor_t & A,
                            const ATensor_t & B);

   template <typename ATensor_t>
   static void CopyWeightsDiffArch(Tensor_t &A, const ATensor_t &B);

   //template<>
   static void CopyDiffArch(Tensor_t A, const Tensor_t & B ) { Copy(A,B); }

      // copy from vector of matrices of different types
   template<typename AMatrix_t>
   static void CopyDiffArch(std::vector<Tensor_t>  & A,
                            const std::vector<AMatrix_t> & B);


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
   static void Identity(Tensor_t & X) {}
   static void IdentityDerivative(Tensor_t & dX, Tensor_t& X,
                                  Tensor_t & Y,  Tensor_t & dY,
                                  ActivationDescriptor_t activationDescr,
                                  const AFloat alpha = 1,
                                  const AFloat beta = 1) {}

   static void ActivationFunctionForward(Tensor_t & X, EActivationFunction activFunct,
                          const ActivationDescriptor_t activationDescr,
                          const double coef = 0.0, const AFloat alpha = 1,
                          const AFloat beta = 0);

   // same as above but using different input/output tensors
   static void ActivationFunctionForward(Tensor_t &Y, const Tensor_t & X, EActivationFunction activFunct,
                                         const ActivationDescriptor_t activationDescr, const double coef = 0.0,
                                         const AFloat alpha = 1, const AFloat beta = 0);

   /** Computes the gradient of the activation function */
   static void ActivationFunctionBackward(Tensor_t & dX, const Tensor_t & Y,
                                          const Tensor_t & dY,  const Tensor_t & X,
                                          EActivationFunction activFunct,
                                          const ActivationDescriptor_t activationDescr,
                                          const AFloat alpha = 1,
                                          const AFloat beta = 0);

   //
   // No cudnn implementation for the following activation functions
   //
   //static void SymmetricRelu(Tensor_t & B);
   static void SymmetricReluDerivative(Tensor_t & B,
                                       const Tensor_t & A) {}

   //static void SoftSign(Tensor_t & B);
   static void SoftSignDerivative(Tensor_t & B,
                                  const Tensor_t & A) {}

   //static void Gauss(Tensor_t & B);
   static void GaussDerivative(Tensor_t & B,
                               const Tensor_t & A) {}
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
      // Dropout
      //____________________________________________________________________________

   /** @name Dropout
    */
      ///@{

   /** Apply dropout with activation probability \p p to the given
    *  tensor \p A and scale the result by reciprocal of \p p. */
   static void DropoutForward(Tensor_t & A,
                              TDescriptors * descriptors,
                              TWorkspace         * workspace,
                              Scalar_t p);

   static void DropoutBackward(Tensor_t & A,
                               TDescriptors * descriptors,
                               TWorkspace   * workspace);

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

   static void BatchNormLayerForwardInference(int axis, const Tensor_t &x, Matrix_t &gamma, Matrix_t &beta,
                                              Tensor_t &y, const Matrix_t &runningMeans,
                                              const Matrix_t &runningVars, Scalar_t epsilon,
                                              const TensorDescriptor_t &);

   static void BatchNormLayerBackward(int axis, const Tensor_t &x, const Tensor_t &dy, Tensor_t &dx,
                                      Matrix_t &gamma, //  Matrix_t &beta, (not needed)
                                      Matrix_t &dgamma, Matrix_t &dbeta, const Matrix_t &mean, const Matrix_t &variance,
                                      const Matrix_t &iVariance, Scalar_t epsilon, const TensorDescriptor_t &);

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

   static Scalar_t L1Regularization(const Matrix_t &W)
   {
      TCudaMatrix<AFloat> mW(W.GetDeviceBuffer(), W.GetSize(), 1);
      return TCuda<AFloat>::L1Regularization(mW);
   }
   static void AddL1RegularizationGradients(Matrix_t &A, const Matrix_t &W, Scalar_t weightDecay)
   {
      TCudaMatrix<AFloat> mA(A.GetDeviceBuffer(), A.GetSize(), 1);
      TCudaMatrix<AFloat> mW(W.GetDeviceBuffer(), W.GetSize(), 1);
      return TCuda<AFloat>::AddL1RegularizationGradients(mA, mW, weightDecay);
   }

   static Scalar_t L2Regularization(const Matrix_t &W)
   {
      TCudaMatrix<AFloat> mW(W.GetDeviceBuffer(), W.GetSize(), 1);
      return TCuda<AFloat>::L2Regularization(mW);
   }
   static void AddL2RegularizationGradients(Matrix_t &A, const Matrix_t &W, Scalar_t weightDecay)
   {
      TCudaMatrix<AFloat> mA(A.GetDeviceBuffer(), A.GetSize(), 1);
      TCudaMatrix<AFloat> mW(W.GetDeviceBuffer(), W.GetSize(), 1);
      return TCuda<AFloat>::AddL1RegularizationGradients(mA, mW, weightDecay);
   }
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

   static void InitializeGauss(Matrix_t &A);
   static void InitializeUniform(Matrix_t &A);
   static void InitializeIdentity(Matrix_t &A);
   static void InitializeZero(Matrix_t &A);
   static void InitializeGlorotNormal(Matrix_t &A);
   static void InitializeGlorotUniform(Matrix_t &A);

   // return static instance of random generator used for initialization
   // if generator does not exist it is created the first time with a random seed (e.g. seed = 0)
   static TRandom &GetRandomGenerator();
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
   static void Dropout(Tensor_t &A, Scalar_t p) {}

   ///@}

   //____________________________________________________________________________
   //
   //  Convolutional Layer Propagation
   //____________________________________________________________________________

   /** @name Forward Propagation in Convolutional Layer
    */
   ///@{

   /** Add the biases in the Convolutional Layer.  */
   static void AddConvBiases(Matrix_t &output, const Matrix_t &biases);
   ///@}

   /** Dummy placeholder - preparation is currently only required for the CUDA architecture. */
   static void PrepareInternals(Tensor_t &) {}

   /** Forward propagation in the Convolutional layer */
   static void ConvLayerForward(Tensor_t &output,
                                Tensor_t &inputActivationFunc, // this is output conv w/o activ func.
                                const Tensor_t &input, const Matrix_t &weights, const Matrix_t &biases,
                                const DNN::CNN::TConvParams &params, EActivationFunction activFunc,
                                Tensor_t & /* inputPrime */, const ConvDescriptors_t &descriptors,
                                ConvWorkspace_t &workspace);
   // const AFloat alpha = 1,
   // const AFloat beta  = 1);

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
   static void ConvLayerBackward(Tensor_t &activationGradientsBackward, Matrix_t &weightGradients,
                                 Matrix_t &biasGradients, Tensor_t &inputActivation, Tensor_t &activationGradients,
                                 const Matrix_t &weights, const Tensor_t &activationBackward,
                                 const Tensor_t &outputTensor, EActivationFunction activFunc,
                                 const ConvDescriptors_t &descriptors, ConvWorkspace_t &workspace, size_t /*batchSize*/,
                                 size_t /*inputHeight*/, size_t /*inputWidth*/, size_t /*depth*/, size_t /*height*/,
                                 size_t /*width*/, size_t /*filterDepth*/, size_t /*filterHeight*/,
                                 size_t /*filterWidth*/, size_t /*nLocalViews*/);

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
    * \p B. No winning indices needed for cuDNN. */
   static void Downsample(Tensor_t &A, Tensor_t & /*B*/, const Tensor_t &C, const PoolingDescriptors_t &descriptors,
                          PoolingWorkspace_t &workspace, size_t imgHeight, size_t imgWidth, size_t fltHeight,
                          size_t fltWidth, size_t strideRows, size_t strideCols);

   ///@}

   /** @name Backward Propagation in Max Pooling Layer
    */
   ///@{
   /** Perform the complete backward propagation step in a Pooling Layer. Based on the
    *  input to and output from the MaxPoolLayer, the gradients for the winning pixels
    *  are computed. */
   static void MaxPoolLayerBackward(Tensor_t &activationGradientsBackward, const Tensor_t &activationGradients,
                                    const Tensor_t & /*indexMatrix*/, const Tensor_t &inputActivation,
                                    const Tensor_t &outputTensor, const PoolingDescriptors_t &descriptors,
                                    PoolingWorkspace_t &workspace, size_t imgHeight, size_t imgWidth, size_t fltHeight,
                                    size_t fltWidth, size_t strideRows, size_t strideCols, size_t nLocalViews);

   ///@}

   //____________________________________________________________________________
   //
   //  Reshape Layer Propagation
   //____________________________________________________________________________
   /** @name Forward and Backward Propagation in Reshape Layer
    */
   ///@{

   /** Transform the matrix \p B to a matrix with different dimensions \p A */
   // static void Reshape(Matrix_t &A, const Matrix_t &B);

   /** Flattens the tensor \p B, such that each matrix, is stretched in
    *  one row, resulting with a matrix \p A. */
   static void Flatten(Tensor_t &A, const Tensor_t &B);

   /** Transforms each row of \p B to a matrix and stores it in the
    *  tensor \p B. */
   static void Deflatten(Tensor_t &A, const Tensor_t &B); // size_t index, size_t nRows,size_t nCols);

   /** Rearrage data accoring to time fill B x T x D out with T x B x D matrix in*/
   static void Rearrange(Tensor_t &out, const Tensor_t &in) { TCuda<AFloat>::Rearrange(out, in); }

   /** Backward pass for Recurrent Networks */
   static Matrix_t &RecurrentLayerBackward(Matrix_t &state_gradients_backward, // BxH
                                           Matrix_t & /* input_weight_gradients */, Matrix_t &/* state_weight_gradients */,
                                           Matrix_t &/* bias_gradients */,
                                           Matrix_t &/* df */,                  // DxH
                                           const Matrix_t &/* state */,         // BxH
                                           const Matrix_t &/* weights_input */, // HxD
                                           const Matrix_t &/* weights_state */, // HxH
                                           const Matrix_t &/* input */,         // BxD
                                           Matrix_t &/* input_gradient */)
   {
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
   
   /** In-place Hadamard (element-wise) product of matrices \p A and \p B
    *  with the result being written into \p A.
    */
   static void Hadamard(Tensor_t &A, const Tensor_t &B)
   {
      TCudaMatrix<AFloat> tmpA(A.GetDeviceBuffer(), 1, A.GetSize());
      TCudaMatrix<AFloat> tmpB(B.GetDeviceBuffer(), 1, B.GetSize());
      assert(A.GetSize() == B.GetSize());
      TCuda<AFloat>::Hadamard(tmpA, tmpB);
   }
   // static void Hadamard(Matrix_t &A,
   //                      const Matrix_t &B);*/
   // {
   //    Tensor_t tA(A);
   //    Hadamard( tA, Tensor_t(B));
   // }

   
   /** Compute the sum of all elements in \p A */
   static Scalar_t Sum(const Matrix_t &A, Scalar_t alpha = 1.0, Scalar_t beta = 0.0);

   /** Check two matrices for equality, taking floating point arithmetic errors into account. */
   //static bool AlmostEquals(const Matrix_t &A, const Matrix_t &B, double epsilon = 0.1);

   /** Add the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstAdd(Matrix_t &A, Scalar_t beta) {
      TCudaMatrix<AFloat> tmp(A.GetDeviceBuffer(), 1, A.GetSize());
      TCuda<AFloat>::ConstAdd(tmp,beta);
   }

   /** Multiply the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstMult(Matrix_t &A, Scalar_t beta) {
      TCudaMatrix<AFloat> tmp(A.GetDeviceBuffer(), 1, A.GetSize());
      TCuda<AFloat>::ConstMult(tmp,beta);
   }

   /** Reciprocal each element of the matrix \p A and write the result into
    * \p A
    */
   static void ReciprocalElementWise(Matrix_t &A) {
      TCudaMatrix<AFloat> tmp(A.GetDeviceBuffer(), 1, A.GetSize());
      TCuda<AFloat>::ReciprocalElementWise(tmp);
   }

   /** Square each element of the matrix \p A and write the result into
    * \p A
    */
   static void SquareElementWise(Matrix_t &A) {
      TCudaMatrix<AFloat> tmp(A.GetDeviceBuffer(), 1, A.GetSize());
      TCuda<AFloat>::SquareElementWise(tmp);
   }

   /** Square root each element of the matrix \p A and write the result into
    * \p A
    */
   //static void SqrtElementWise(Matrix_t &A, Scalar_t alpha = 1, Scalar_t beta = 0, Scalar_t gamma = 0) {
   static void SqrtElementWise(Matrix_t &A) {
      TCudaMatrix<AFloat> tmp(A.GetDeviceBuffer(), 1, A.GetSize());
      TCuda<AFloat>::SqrtElementWise(tmp);
   }

      // optimizer functions
   static void AdamUpdate(Matrix_t & A, const Matrix_t & M, const Matrix_t & V, Scalar_t alpha, Scalar_t eps) {
      TCudaMatrix<AFloat> tmpA(A.GetDeviceBuffer(), A.GetSize(),1);
      TCudaMatrix<AFloat> tmpM(M.GetDeviceBuffer(), M.GetSize(),1);
      TCudaMatrix<AFloat> tmpV(V.GetDeviceBuffer(), V.GetSize(),1);
      TCuda<AFloat>::AdamUpdate(tmpA, tmpM, tmpV,alpha, eps);
   }
   static void AdamUpdateFirstMom(Matrix_t & A, const Matrix_t & B, Scalar_t beta) {
      TCudaMatrix<AFloat> tmpA(A.GetDeviceBuffer(), A.GetSize(),1);
      TCudaMatrix<AFloat> tmpB(B.GetDeviceBuffer(), B.GetSize(),1);
      TCuda<AFloat>::AdamUpdateFirstMom(tmpA, tmpB,  beta);
   }
   static void AdamUpdateSecondMom(Matrix_t & A, const Matrix_t & B, Scalar_t beta) {
      TCudaMatrix<AFloat> tmpA(A.GetDeviceBuffer(), A.GetSize(),1);
      TCudaMatrix<AFloat> tmpB(B.GetDeviceBuffer(), B.GetSize(),1);
      TCuda<AFloat>::AdamUpdateSecondMom(tmpA, tmpB,  beta);
   }

      // printing of tensor
   static void PrintTensor( const Tensor_t & A, const std::string name = "tensor", bool = false);



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
template <typename ATensor>
void TCudnn<AFloat>::CopyDiffArch(TCudaTensor<AFloat> &B,
                        const ATensor &A)
{

   // should add static assert that A has not to be same type as B

   // this copying tensors from different architectures
   if (B.GetLayout() == GetTensorLayout() ) {
      assert(B.GetShape().size() == 4);
      for (size_t i = 0; i < A.GetFirstSize(); ++i) {
         TMatrixT<AFloat> matIn = A.At(i).GetMatrix(); // this convert tensor (B,D,HW) in  (D,HW)i -> (D,HW)i
         // TMAtrix has the correct layout (row-wise) no need to traspose in this case
         TCudaTensor<AFloat> tmpOut = B.At(i); // matrix (D,HW)
         // copy will copy the buffer
         TCudaTensor<AFloat> tmpIn(matIn.GetMatrixArray(), tmpOut.GetShape(), tmpOut.GetLayout());
         Copy(tmpOut, tmpIn);
      }
   } else {
      // case of same layout (column major)
      TMatrixT<AFloat> tmp = A;
      TCudaMatrix<AFloat> tmp2(tmp);
      TCudaTensor<AFloat> tA(tmp2);
      Copy(B, tA);
   }
}

//____________________________________________________________________________
template <typename AFloat>
template <typename AMatrix>
void TCudnn<AFloat>::CopyWeightsDiffArch(TCudaTensor<AFloat> &B, const  AMatrix &A)
{
   // copy from another architecture using the reference one
   // this is not very efficient since creates temporary objects
   TMatrixT<AFloat> tmp = A; // .GetMatrix();
   // we need to traspose for different layout
   if (B.GetLayout() == GetTensorLayout()  ) {
      // this is for CNN weights that are in row-major formats
      assert(B.GetShape().size() == 4);  // weights shape should be 4
      tmp.T();
   }
   TCudaMatrix<AFloat> tmp2(tmp);
   TCudaTensor<AFloat> tA(tmp2);
   Copy(B, tA);
}

//____________________________________________________________________________
template <typename AFloat>
template <typename AMatrix_t>
void TCudnn<AFloat>::CopyDiffArch(std::vector<Tensor_t> &B,
                            const std::vector<AMatrix_t> &A)
{
   for (size_t i = 0; i < B.size(); ++i) {
      CopyWeightsDiffArch(B[i], A[i]);
   }
}

template <typename AFloat>
void TCudnn<AFloat>::PrintTensor(const typename TCudnn<AFloat>::Tensor_t & A, const std::string name, bool truncate )
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
         size_t n =  A.GetShape()[1];
         if (truncate) n = std::min(n,size_t(10));
         for (size_t j = 0; j < n; ++j) {
            std::cout << A(i,j) << " ";

         }
         if (truncate && n < A.GetShape()[1]) std::cout << " ...... ";
         std::cout << " } " << std::endl;
      }
   } else if  (A.GetShape().size() == 3 ) {
      for (size_t i = 0; i < A.GetFirstSize(); ++i) {
         std::cout << "{ ";
         for (size_t j = 0; j < A.GetHSize(); ++j) {
            std::cout << "{ ";
            size_t n =  A.GetWSize();
            if (truncate)  n = std::min(n,size_t(10));
            for (size_t k = 0; k < n; ++k) {
               std::cout << A(i,j,k) << " ";
            }
            if (truncate && n < A.GetWSize()) std::cout << " ...... ";
            std::cout << " } " << std::endl;
         }
         std::cout << " } " << std::endl;
      }
   } else if  (A.GetShape().size() == 4 ) {
      for (size_t i = 0; i < A.GetShape()[0]; ++i) {
         std::cout << "{ ";
         for (size_t j = 0; j < A.GetShape()[1]; ++j) {
            std::cout << "{ ";
            for (size_t k = 0; k < A.GetShape()[2]; ++k) {
               size_t n =  A.GetShape()[3];
               if (truncate)  n = std::min(n,size_t(10));
               for (size_t l = 0; l < n; ++l) {
                  std::cout << A(i,j,k,l) << " ";
               }
               if (truncate && n < A.GetShape()[3]) std::cout << " ...... ";
               std::cout << " } " << std::endl;
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

// initialize the CNN options
// possible options for forward (from 0 to 7)
//
//  0 : CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_GEMM;
//  1 : CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_PRECOMP_GEMM;
//  6  : CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD;
//  7 : CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD_NONFUSED;  (lots of memory)

// for backward data (from 0 to 5)
//  1 : CUDNN_CONVOLUTION_BWD_DATA_ALGO_1;
//  5  CUDNN_CONVOLUTION_BWD_DATA_ALGO_WINOGRAD_NONFUSED;

template <typename AFloat>
int TCudnn<AFloat>::CNNOptions::ConvFwdAlgorithm = -1;
template <typename AFloat>
int TCudnn<AFloat>::CNNOptions::ConvBwdDataAlgorithm = -1;
template <typename AFloat>
int TCudnn<AFloat>::CNNOptions::ConvBwdFilterAlgorithm = -1;
template <typename AFloat>
Long_t TCudnn<AFloat>::CNNOptions::ConvMaxWorkspaceSize = -1;  // -1 let use Cudnn defaults

} // namespace DNN
} // namespace TMVA

#endif
#endif
