// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 13/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 //////////////////////////////////////////////////////////////////
 // Implementation of the functions required for the forward and //
 // backward propagation of activations through a neural network //
 // for CUDA architectures using cuDNN library.                  //
 //////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/TCudnn.h"

#include "TMVA/DNN/CNN/ConvLayer.h"

#include "TMVA/DNN/Architectures/Cuda.h"

#include "TRandom.h"

// #include "TMVA/DNN/Architectures/Cuda/Device.h"
// #include "Kernels.cuh"*/
// #include <math.h>

// for std::numeric_limits<T>::max()
#include <limits>

namespace TMVA {
namespace DNN  {


//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::MultiplyTranspose(TCudaTensor<AFloat> &output,
                                       const TCudaTensor<AFloat> &input,
                                       const TCudaTensor<AFloat> &weights)
{
   //PrintTensor(input,"input to MultTrans");
   //PrintTensor(weights,"dense layer  weights");
   TCuda<AFloat>::MultiplyTranspose(output, input, weights.GetMatrix());
   //PrintTensor(input,"output of  MultTrans");
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::AddRowWise(TCudaTensor<AFloat> &output,
                                const TCudaTensor<AFloat> &biases)
{
   TCuda<AFloat>::AddRowWise( output, biases.GetMatrix());
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::Backward(TCudaTensor<AFloat> & activation_gradients_backward,
                              TCudaTensor<AFloat> & weight_gradients,
                              TCudaTensor<AFloat> & bias_gradients,
                              TCudaTensor<AFloat> & df,
                              const TCudaTensor<AFloat> & activation_gradients,
                              const TCudaTensor<AFloat> & weights,
                              const TCudaTensor<AFloat> & activation_backward)
{
   // use implentation from TCuda

   //std::cout << "\n\n ------ Backward--------\n";
   //PrintTensor(activation_backward,"input to backward");
   //PrintTensor(weights,"dense layer  weights");
   //PrintTensor(activation_gradients,"input dy");
   //PrintTensor(activation_gradients,"df");

   TCudaMatrix<AFloat> weightGradMatrix = weight_gradients.GetMatrix();
   TCudaMatrix<AFloat> biasGradMatrix = bias_gradients.GetMatrix();

   TCuda<AFloat>::Backward(activation_gradients_backward,
                              weightGradMatrix,
                              biasGradMatrix,
                              df,
                              activation_gradients,
                              weights.GetMatrix(),
                              activation_backward);

   //PrintTensor(activation_gradients_backward,"computed dx");
   //PrintTensor(weight_gradients,"computed dw");
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::Copy(Tensor_t & B, const Tensor_t & A)
{
   size_t nElements = A.GetSize();
   R__ASSERT(nElements == B.GetSize());

   cudaMemcpyAsync(B.GetDataPointer(), A.GetDataPointer(),
                   nElements * sizeof(AFloat), cudaMemcpyDeviceToDevice, 0);
}



///////////////////////////////////////////////////////////////////////////////
// Initialization of the cuDNN objects for the different layers
// ...
///////////////////////////////////////////////////////////////////////////////
template<typename AFloat>
void TCudnn<AFloat>::InitializeBNormDescriptors(TDescriptors * & descriptors, typename TCudnn<AFloat>::BNormLayer_t *L)
{
   auto bnormDescriptors = new BNormDescriptors_t ();

   // reshaped tensors  of bnorm  layer - look if output tensor has right shape
   Tensor_t &outputTensor = L->GetOutput();
   Tensor_t &data = L->GetReshapedData();
   if (L->GetNormAxis() == -1 && L->GetBatchSize() == outputTensor.GetShape()[0] && L->GetDepth() == 1 && L->GetHeight() == 1 ) {
      // case of dense layer before - need to reshape the data
      R__ASSERT(outputTensor.GetLayout() != GetTensorLayout());  // has to be output column major
      // case of convolutions layer before
      Tensor_t &data = L->GetReshapedData();
      // it is not important which buffer we use. Important is the shape and layout of tensor
      data = Tensor_t(outputTensor.GetDeviceBuffer(), {1, L->GetWidth(), 1, L->GetBatchSize()}, GetTensorLayout(), 0, 0);
   } else if (L->GetNormAxis() == 1 ) {
      // case of convolutional layer  before
      outputTensor.PrintShape("output");
      Tensor_t tmp( {L->GetBatchSize() , L->GetDepth(), L->GetHeight(), L->GetWidth()}, GetTensorLayout(), 0, 0);
      tmp.PrintShape("tmp");
      data = Tensor_t(outputTensor.GetDeviceBuffer(), {L->GetBatchSize() , L->GetDepth(), L->GetHeight(), L->GetWidth() }, GetTensorLayout(), 0, 0);


      // reshape output tensor and activation gradient tensor of pool layer
      outputTensor = Tensor_t(outputTensor.GetDeviceBuffer(),
        {L->GetBatchSize(), L->GetDepth(), L->GetHeight(), L->GetWidth()},
                              GetTensorLayout(), 0, 0 );

      Tensor_t &activationGradients = L->GetActivationGradients();
      activationGradients = Tensor_t(activationGradients.GetDeviceBuffer(),
                                     outputTensor.GetShape(), GetTensorLayout(), 0, 0);
      outputTensor.PrintShape("output2");

   }

   outputTensor.PrintShape("output bnorm");
   data.PrintShape("reshaped data");

   // Tensor_t &activationGradients = L->GetActivationGradients();
   // activationGradients = Tensor_t(activationGradients.GetDeviceBuffer(), outputTensor.GetShape(), GetTensorLayout(), 0, 0);

   CUDNNCHECK(cudnnCreateTensorDescriptor(&bnormDescriptors->HelperDescriptor));

   cudnnBatchNormMode_t bnMode = CUDNN_BATCHNORM_SPATIAL;

   // CUDNN_BATCHNORM_PER_ACTIVATION;
   // if (L->GetNormAxis() == 1)
   //    bnMode = CUDNN_BATCHNORM_SPATIAL;

   CUDNNCHECK(cudnnDeriveBNTensorDescriptor(bnormDescriptors->HelperDescriptor,
                                            data.GetTensorDescriptor(),
                                            bnMode) );

   descriptors = bnormDescriptors;
}

template <typename AFloat>
void TCudnn<AFloat>::InitializeActivationDescriptor(TCudnn<AFloat>::ActivationDescriptor_t &descriptor,
                                                            EActivationFunction activFunc, double coef)
{
   cudnnActivationMode_t activationMode;
   bool isIdentity = false;
   switch (activFunc) {
      case EActivationFunction::kIdentity:
         isIdentity = true;
         break; // Identity activation only works for cudnnConvolutionBiasActivationForward()
      case EActivationFunction::kRelu:
         activationMode = CUDNN_ACTIVATION_RELU;
         break;
      case EActivationFunction::kSigmoid:
         activationMode = CUDNN_ACTIVATION_SIGMOID;
         break;
      case EActivationFunction::kTanh:
         activationMode = CUDNN_ACTIVATION_TANH;
         break;
      case EActivationFunction::kFastTanh:
         activationMode = CUDNN_ACTIVATION_TANH;
         break;
         // The activations otherwise used are not supported by cuDNN
      default:
         activationMode = CUDNN_ACTIVATION_RELU;
   };

   CUDNNCHECK(cudnnCreateActivationDescriptor(&descriptor));

    // Dont set activation function descriptor for identity function
    if (!isIdentity) CUDNNCHECK(cudnnSetActivationDescriptor(descriptor, activationMode, CUDNN_PROPAGATE_NAN, coef));
}

template<typename AFloat>
void TCudnn<AFloat>::InitializeConvDescriptors(TDescriptors * & descriptors, ConvLayer_t *L) {

   auto convDescriptors = new CNN::TCNNDescriptors<typename TCudnn<AFloat>::ConvLayer_t> ();

   //FIXME: Move this to constructor
   cudnnDataType_t   cudnnDataType;
   if      (std::is_same<AFloat, double>::value) { cudnnDataType = CUDNN_DATA_DOUBLE;}
   else if (std::is_same<AFloat, float>::value)  { cudnnDataType = CUDNN_DATA_FLOAT;}

   double coef = 0.0;   // this is a coefficient which can be used for activation (e.g. relu) . it is not yet supported
   InitializeActivationDescriptor(convDescriptors->HelperDescriptor, L->GetActivationFunction(), coef);

   CUDNNCHECK(cudnnCreateConvolutionDescriptor(&convDescriptors->LayerDescriptor));
   CUDNNCHECK(cudnnCreateFilterDescriptor(&convDescriptors->WeightsDescriptor));

   // Set the convolution parameters
   CUDNNCHECK(cudnnSetConvolution2dDescriptor(convDescriptors->LayerDescriptor,
                                              L->GetPaddingHeight(),
                                              L->GetPaddingWidth(),
                                              L->GetStrideRows(),
                                              L->GetStrideCols(),
                                              1,                 //Dilation height
                                              1,                 //Dilation width
                                              CUDNN_CROSS_CORRELATION,
                                              cudnnDataType));

   // Set the  filter parameters
   CUDNNCHECK(cudnnSetFilter4dDescriptor(convDescriptors->WeightsDescriptor,
                                         cudnnDataType,
                                         CUDNN_TENSOR_NCHW,
                                         L->GetDepth(),
                                         L->GetInputDepth(),
                                         L->GetFilterHeight(),
                                         L->GetFilterWidth()));

   descriptors = convDescriptors;
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::InitializePoolDescriptors(TDescriptors * & descriptors,
                                               PoolingLayer_t *L) {
   auto poolDescriptors = new CNN::TCNNDescriptors<typename TCudnn<AFloat>::PoolingLayer_t> ();
   CUDNNCHECK(cudnnCreatePoolingDescriptor(&poolDescriptors->LayerDescriptor));

   CUDNNCHECK(cudnnCreateDropoutDescriptor(&poolDescriptors->HelperDescriptor));

   CUDNNCHECK(cudnnSetPooling2dDescriptor(poolDescriptors->LayerDescriptor,
                                          CUDNN_POOLING_MAX,
                                          CUDNN_PROPAGATE_NAN,
                                          L->GetFilterHeight(),
                                          L->GetFilterWidth(),
                                          0,//L->GetPaddingHeight()
                                          0,//L->GetPaddingWidth()
                                          L->GetStrideRows(),
                                          L->GetStrideCols()));



   descriptors = poolDescriptors;

   // reshape output tensor and activation gradient tensor of pool layer
   Tensor_t &outputTensor = L->GetOutput();
   outputTensor = Tensor_t(outputTensor.GetDeviceBuffer(),
                           {L->GetBatchSize(), L->GetDepth(), L->GetHeight(), L->GetWidth()},
                           GetTensorLayout(), 0, 0);

   Tensor_t &activationGradients = L->GetActivationGradients();
   activationGradients = Tensor_t(activationGradients.GetDeviceBuffer(),
                                  outputTensor.GetShape(), GetTensorLayout(), 0, 0);
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::ReleaseConvDescriptors(TDescriptors * descriptors) {
   auto convDescriptors = static_cast<ConvDescriptors_t *>(descriptors);
   ReleaseDescriptor(convDescriptors->LayerDescriptor);
   ReleaseDescriptor(convDescriptors->HelperDescriptor);
   ReleaseDescriptor(convDescriptors->WeightsDescriptor);
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleasePoolDescriptors(TDescriptors * descriptors) {
   auto poolDescriptors = static_cast<PoolingDescriptors_t *>(descriptors);
   ReleaseDescriptor(poolDescriptors->LayerDescriptor);
   ReleaseDescriptor(poolDescriptors->HelperDescriptor);
   ReleaseDescriptor(poolDescriptors->WeightsDescriptor);
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseBNormDescriptors(TDescriptors * descriptors) {
   auto bnormDescriptors = static_cast<BNormDescriptors_t *>(descriptors);
   ReleaseDescriptor(bnormDescriptors->HelperDescriptor);  // it is a tensor descriptor
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(ActivationDescriptor_t & activationDescr) {
   CUDNNCHECK(cudnnDestroyActivationDescriptor(activationDescr));
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(ConvolutionDescriptor_t & convolutionDescr) {
   CUDNNCHECK(cudnnDestroyConvolutionDescriptor(convolutionDescr));
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(DropoutDescriptor_t & dropoutDescr) {
   CUDNNCHECK(cudnnDestroyDropoutDescriptor(dropoutDescr));
}
//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(TensorDescriptor_t & tensorDescr) {
   CUDNNCHECK(cudnnDestroyTensorDescriptor(tensorDescr));
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(FilterDescriptor_t & filterDescr) {
   CUDNNCHECK(cudnnDestroyFilterDescriptor(filterDescr));
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::ReleaseDescriptor(PoolingDescriptor_t & poolingDescr) {
   CUDNNCHECK(cudnnDestroyPoolingDescriptor(poolingDescr));
}


//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::InitializeConvWorkspace(TWorkspace * & workspace,
                                             TDescriptors * & descriptors,
                                             const DNN::CNN::TConvParams & /*params*/,
                                             ConvLayer_t *L) {
   auto convWorkspace = new ConvWorkspace_t();
   size_t memLimit = (CNNOptions::ConvMaxWorkspaceSize > 0) ? static_cast<size_t>(CNNOptions::ConvMaxWorkspaceSize) : 0;
   auto convDescriptors = static_cast<ConvDescriptors_t *>(descriptors);
   // can we do the following and substitute below???
   // auto weightsDescriptor{convDescriptors->WeightsDescriptor};
   // auto convDescriptor{convDescriptors->LayerDescriptor};

#if (CUDNN_VERSION >= 8000)
   enum algoPreference { no_workspace, fastest, workspace_limit };
   algoPreference algoChoice;
   // C++11 lambdas cannot be templated, so we have to do this HORRIBLE stuff...
   class LocalPerf {
    public:
      LocalPerf(cudnnConvolutionFwdAlgoPerf_t * fwd) {m_fwd = fwd;}
      LocalPerf(cudnnConvolutionBwdFilterAlgoPerf_t * bwdFilter) {m_bwdFilter = bwdFilter;}
      LocalPerf(cudnnConvolutionBwdDataAlgoPerf_t * bwdData) {m_bwdData = bwdData;}
      size_t getMemory(int i) {return m_fwd != nullptr ? m_fwd[i].memory : m_bwdFilter != nullptr ? m_bwdFilter[i].memory : m_bwdData != nullptr ? m_bwdData[i].memory : 0;}
      float getTime(int i) {return m_fwd != nullptr ? m_fwd[i].time : m_bwdFilter != nullptr ? m_bwdFilter[i].time : m_bwdData != nullptr ? m_bwdData[i].time : 0;}
      cudnnStatus_t getStatus(int i) {return m_fwd != nullptr ? m_fwd[i].status : m_bwdFilter != nullptr ? m_bwdFilter[i].status : m_bwdData != nullptr ? m_bwdData[i].status : CUDNN_STATUS_BAD_PARAM;}
      int getIdx(algoPreference const & algoPref, int const algoCount, size_t memLim = std::numeric_limits<size_t>::max()) {
         int algoIdx{0};
         if (algoPref == algoPreference::fastest) {  // prefer fastest
            float temp_runtime{std::numeric_limits<float>::max()};
            for (int i = 0; i < algoCount; ++i) {
               if (getStatus(i) == CUDNN_STATUS_SUCCESS && getTime(i) < temp_runtime) {
                  temp_runtime = getTime(i);
                  algoIdx = i;
               }
            }
         } else if (algoPref == algoPreference::workspace_limit) {  // constrain to workspace size
            float temp_runtime{std::numeric_limits<float>::max()};
            for (int i = 0; i < algoCount; ++i) {
               if (getStatus(i) == CUDNN_STATUS_SUCCESS && getTime(i) < temp_runtime && getMemory(i) <= memLim) {
                  temp_runtime = getTime(i);
                  algoIdx = i;
               }
            }
         } else {  // prefer smallest workspace size
            size_t temp_memsize{std::numeric_limits<size_t>::max()};
            for (int i = 0; i < algoCount; ++i) {
               if (getStatus(i) == CUDNN_STATUS_SUCCESS && getMemory(i) < temp_memsize) {
                  temp_memsize = getMemory(i);
                  algoIdx = i;
               }
            }
         }
         return algoIdx;
      };
    private:
      LocalPerf();
      // these three type are absolutely equivalent
      // and one can access them as they wish to get info
      cudnnConvolutionFwdAlgoPerf_t *m_fwd = nullptr;
      cudnnConvolutionBwdFilterAlgoPerf_t *m_bwdFilter = nullptr;
      cudnnConvolutionBwdDataAlgoPerf_t * m_bwdData = nullptr;
   };
#else
   // More detailed alternative: cudnnFindConvolutionForwardAlgorithm (only option in newer cuDNN versions)
   cudnnConvolutionFwdPreference_t preferenceFwd;
   cudnnConvolutionBwdDataPreference_t preferenceBwdData;
   cudnnConvolutionBwdFilterPreference_t preferenceBwdFilter;
#endif
   // decide on algorithm preference early
   if (CNNOptions::ConvMaxWorkspaceSize < 0) {
#if (CUDNN_VERSION >= 8000)
      // fastest overall
      algoChoice = fastest;
#else
      preferenceFwd = CUDNN_CONVOLUTION_FWD_PREFER_FASTEST;
      preferenceBwdData = CUDNN_CONVOLUTION_BWD_DATA_PREFER_FASTEST;
      preferenceBwdFilter = CUDNN_CONVOLUTION_BWD_FILTER_PREFER_FASTEST;
#endif

   } else if (CNNOptions::ConvMaxWorkspaceSize == 0) {
      // no workspace case
#if (CUDNN_VERSION >= 8000)
      algoChoice = no_workspace;
#else
      preferenceFwd = CUDNN_CONVOLUTION_FWD_NO_WORKSPACE;
      preferenceBwdData = CUDNN_CONVOLUTION_BWD_DATA_NO_WORKSPACE;
      preferenceBwdFilter = CUDNN_CONVOLUTION_BWD_FILTER_NO_WORKSPACE;
#endif

   } else {
      // fastest in memory limit
#if (CUDNN_VERSION >= 8000)
      algoChoice = workspace_limit;
#else
      preferenceFwd = CUDNN_CONVOLUTION_FWD_SPECIFY_WORKSPACE_LIMIT;
      preferenceBwdData = CUDNN_CONVOLUTION_BWD_DATA_SPECIFY_WORKSPACE_LIMIT;
      preferenceBwdFilter = CUDNN_CONVOLUTION_BWD_FILTER_SPECIFY_WORKSPACE_LIMIT;
#endif
   }
   // fix the weight tensor shapes
   // by default the weights are columnmajor, set them to be row major . At this points
   // they are not yet initialized
   Tensor_t & filters = L->GetWeightsAt(0);
   filters = Tensor_t(filters.GetDeviceBuffer(), {L->GetDepth(), L->GetInputDepth(), L->GetFilterHeight(), L->GetFilterWidth()}, MemoryLayout::RowMajor, 0, 0);
   // PrintTensor(L->GetWeightsAt(0));
   Tensor_t & biases = L->GetBiasesAt(0);
   biases = Tensor_t(biases.GetDeviceBuffer(), {1, L->GetDepth(), 1, 1}, GetTensorLayout(), 0, 0);

   Tensor_t & outputTensor = L->GetOutput();
   outputTensor = Tensor_t(outputTensor.GetDeviceBuffer(), {L->GetBatchSize(), L->GetDepth(), L->GetHeight(), L->GetWidth()}, GetTensorLayout(), 0, 0);
   Tensor_t & inputActivation = L->GetInputActivation();
   inputActivation = Tensor_t(inputActivation.GetDeviceBuffer(),outputTensor.GetShape() ,GetTensorLayout(), 0, 0);

   Tensor_t &  activationGradients = L->GetActivationGradients();
   activationGradients = Tensor_t(activationGradients.GetDeviceBuffer(),outputTensor.GetShape(), GetTensorLayout(), 0, 0);

   Tensor_t & weightGradients = L->GetWeightGradientsAt(0);
   weightGradients = Tensor_t(weightGradients.GetDeviceBuffer(), filters.GetShape(), GetTensorLayout(), 0, 0);

   Tensor_t & biasGradients = L->GetBiasGradientsAt(0);
   biasGradients = Tensor_t(biasGradients.GetDeviceBuffer(), biases.GetShape(), GetTensorLayout(), 0, 0);


   // FIXME: Use descriptors instead (Tensor device memory is otherwise allocated during initialization)
   // Tensor_t inputTensor  ({L->GetBatchSize(), L->GetInputDepth(), L->GetInputHeight(), L->GetInputWidth()}, MemoryLayout::RowMajor, 0, 0);
   cudnnTensorDescriptor_t  inputTensorDescriptor;
   CUDNNCHECK(cudnnCreateTensorDescriptor(&inputTensorDescriptor) );
   CUDNNCHECK(cudnnSetTensor4dDescriptor(inputTensorDescriptor,
                                             CUDNN_TENSOR_NCHW,// Layout of the tensor in memory
                                             Tensor_t::GetDataType(),
                                             (int)L->GetBatchSize(),
                                             (int)L->GetInputDepth(),
                                             (int)L->GetInputHeight(),
                                             (int)L->GetInputWidth() ) );


   // size_t outputHeight = ConvLayer_t::calculateDimension(L->GetInputHeight(), L->GetFilterHeight(), L->GetPaddingHeight(), L->GetStrideRows());
   // size_t outputWidth  = ConvLayer_t::calculateDimension(L->GetInputWidth(), L->GetFilterWidth(), L->GetPaddingWidth(),  L->GetStrideCols());
   //Tensor_t outputTensor ({L->GetBatchSize(), L->GetDepth(), outputHeight, outputWidth}, MemoryLayout::RowMajor, 0, 0);

   // Get access to cudnn library handle, which is static for the CudaTensor class
   cudnnHandle_t cudnnHandle = outputTensor.GetCudnnHandle();

   // cuDNN decides which algorithm to use
#if (CUDNN_VERSION >= 8000)
   /**
    * I'm sure there may be a faster way, but this works
    */
   int convRequestedAlgoCount{0};  // requestedAlgoCount is setting how many algorithms to try
   CUDNNCHECK(cudnnGetConvolutionForwardAlgorithmMaxCount(cudnnHandle, &convRequestedAlgoCount))  // ask cuDNN how much it can try

   int algoCount;
   cudnnConvolutionFwdAlgoPerf_t convFwdPerfResults[convRequestedAlgoCount];  // this will store metrics to choose convolution algorithm
   CUDNNCHECK(
      cudnnFindConvolutionForwardAlgorithm(
         cudnnHandle,
         inputTensorDescriptor,
         convDescriptors->WeightsDescriptor,
         convDescriptors->LayerDescriptor,
         outputTensor.GetTensorDescriptor(),
         convRequestedAlgoCount,
         &algoCount,
         convFwdPerfResults
      )
   );
   // we could also do it with the expert mode (cudnnFindConvolutionForwardAlgorithmEx),
   // i.e.
   // create an input tensor before the inputTensorDescriptor
   // and get the descriptor from there
   // Tensor_t inputTensor({L->GetBatchSize(), L->GetInputDepth(), L->GetInputHeight(), L->GetInputWidth()}, MemoryLayout::RowMajor, 0, 0);
   // CUDNNCHECK(cudnnFindConvolutionForwardAlgorithmEx(
   //    cudnnHandle,
   //    inputTensor.GetTensorDescriptor(),
   //    &inputTensor,
   //    convDescriptors->WeightsDescriptor,
   //    &filters,
   //    convDescriptors->LayerDescriptor,
   //    outputTensor.GetTensorDescriptor(),
   //    &outputTensor,
   //    convRequestedAlgoCount,
   //    &algoCount,
   //    &convFwdPerfResults,
   //    &convWorkspace,
   //    memLimit));  // use memLimit for workspace size
   // instead choose either fastest or lowest memory algo as per preference
   LocalPerf fwdPerfResults{convFwdPerfResults};
   convWorkspace->AlgorithmForward = convFwdPerfResults[fwdPerfResults.getIdx(algoChoice, algoCount, memLimit)].algo;
#else
   CUDNNCHECK(cudnnGetConvolutionForwardAlgorithm(
      cudnnHandle, inputTensorDescriptor, convDescriptors->WeightsDescriptor, convDescriptors->LayerDescriptor,
      outputTensor.GetTensorDescriptor(), preferenceFwd,
      memLimit, // Memory limit in bytes for mode CUDNN_CONVOLUTION_FWD_SPECIFY_WORKSPACE_LIMIT
      &convWorkspace->AlgorithmForward));
#endif

   // Allocate memory for the convolution
   //size_t workSpaceSizeInBytes = 0;

   std::cout << "CONV FWD Algo used for convolution of input shape { " << L->GetBatchSize() << " , " <<  L->GetInputDepth() << " , "
                                                <<L->GetInputHeight() << " , " << L->GetInputWidth() << " } is "
             << convWorkspace->AlgorithmForward << std::endl;

   // convWorkspace->AlgorithmForward = CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD;
   // convWorkspace->AlgorithmForward = CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_GEMM;
   // convWorkspace->AlgorithmForward = CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_PRECOMP_GEMM;
   //     CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD_NONFUSED;
   if (CNNOptions::ConvFwdAlgorithm > 0) {
      convWorkspace->AlgorithmForward = (cudnnConvolutionFwdAlgo_t) CNNOptions::ConvFwdAlgorithm;
      std::cout << " but force using " << convWorkspace->AlgorithmForward << std::endl;
   }


   cudnnMathType_t math_type = CUDNN_TENSOR_OP_MATH; // : CUDNN_DEFAULT_MATH);
   // if using tensor math (cudnn version > 7)
   CUDNNCHECK(cudnnSetConvolutionMathType(convDescriptors->LayerDescriptor, math_type));

   CUDNNCHECK(cudnnGetConvolutionForwardWorkspaceSize(cudnnHandle,
                                                      inputTensorDescriptor,
                                                      convDescriptors->WeightsDescriptor,
                                                      convDescriptors->LayerDescriptor,
                                                      outputTensor.GetTensorDescriptor(),
                                                      convWorkspace->AlgorithmForward,
                                                      &convWorkspace->ForwardWorkspaceSize));

   if (convWorkspace->ForwardWorkspaceSize) cudaMalloc(&convWorkspace->ForwardWorkspace, convWorkspace->ForwardWorkspaceSize*sizeof(AFloat));
   if (convWorkspace->ForwardWorkspaceSize > 0 && convWorkspace->ForwardWorkspace == nullptr  ) {
      std::cerr << "Error allocating FWD CONV workspace of size " << convWorkspace->ForwardWorkspaceSize << " - probably running out of memory on the GPU"
      << std::endl;
      std::cout << " layer input shape is  { " << L->GetBatchSize() << " , " <<  L->GetInputDepth() << " , "
                                                <<L->GetInputHeight() << " , " << L->GetInputWidth() << " } " << std::endl;
      //inputTensor.PrintShape("inputTensor");
      R__ASSERT(false);
   }
   //
   // Backward Algorithm
   //

   //Tensor_t activationGradients ({L->GetBatchSize(), L->GetDepth(), outputHeight, outputWidth}, MemoryLayout::RowMajor, 0, 0);
   //Tensor_t activationGradientsBackward ({L->GetBatchSize(), L->GetInputDepth(), L->GetInputHeight(), L->GetInputWidth()}, MemoryLayout::RowMajor, 0, 0);
   cudnnTensorDescriptor_t activationGradientsBackwardDescriptor = inputTensorDescriptor;

   cudnnHandle = activationGradients.GetCudnnHandle();
   // dx : Activation gradient to be computed                               -> activationGradients [in place op]
   // dy : Gradient of activation from the following layer (backpropagation)-> activationGradients

#if (CUDNN_VERSION >= 8000)
   /**
    * I'm sure there may be a faster way, but this works
    */
   CUDNNCHECK(cudnnGetConvolutionBackwardDataAlgorithmMaxCount(cudnnHandle, &convRequestedAlgoCount))  // ask cuDNN how much it can try
   cudnnConvolutionBwdDataAlgoPerf_t convBwdDataPerfResults[convRequestedAlgoCount];  // this will store metrics to choose convolution algorithm
   CUDNNCHECK(cudnnFindConvolutionBackwardDataAlgorithm(
      cudnnHandle,
      convDescriptors->WeightsDescriptor,
      activationGradients.GetTensorDescriptor(),
      convDescriptors->LayerDescriptor,
      activationGradientsBackwardDescriptor,
      convRequestedAlgoCount,
      &algoCount,
      convBwdDataPerfResults));
   // we could also do it with the expert mode (cudnnFindConvolutionForwardAlgorithmEx),
   // i.e.
   // CUDNNCHECK(cudnnFindConvolutionBackwardDataAlgorithmEx(
   //    cudnnHandle,
   //    convDescriptors->WeightsDescriptor,
   //    &filters,
   //    activationGradients.GetTensorDescriptor(),
   //    &activationGradients,
   //    convDescriptors->LayerDescriptor,
   //    activationGradientsBackwardDescriptor,
   //    &inputTensor,
   //    convRequestedAlgoCount,
   //    &algoCount,
   //    &convBwdDataPerfResults,
   //    &convWorkspace,
   //    memLimit));  // use memLimit for workspace size
   // instead choose either fastest or lowest memory algo as per preference
   LocalPerf bwdDataPerfResults{convBwdDataPerfResults};
   convWorkspace->AlgorithmBackward = convBwdDataPerfResults[bwdDataPerfResults.getIdx(algoChoice, algoCount, memLimit)].algo;
#else
   CUDNNCHECK(cudnnGetConvolutionBackwardDataAlgorithm(cudnnHandle,
                                                      convDescriptors->WeightsDescriptor,
                                                      activationGradients.GetTensorDescriptor(),
                                                      convDescriptors->LayerDescriptor,
                                                      activationGradientsBackwardDescriptor,
                                                      preferenceBwdData, memLimit,
                                                      &convWorkspace->AlgorithmBackward));
#endif

   std::cout << "CONV BWD Data Algo used  is "  << convWorkspace->AlgorithmBackward << std::endl;
   //CUDNNCHECK(cudnnSetConvolutionMathType(convDescriptors->LayerDescriptor, CUDNN_TENSOR_OP_MATH));


   if (CNNOptions::ConvBwdDataAlgorithm > 0) {
      convWorkspace->AlgorithmBackward = (cudnnConvolutionBwdDataAlgo_t)CNNOptions::ConvBwdDataAlgorithm;
      std::cout << " but force using " << convWorkspace->AlgorithmBackward << std::endl;
   }

   CUDNNCHECK(cudnnGetConvolutionBackwardDataWorkspaceSize(cudnnHandle,
                                                           convDescriptors->WeightsDescriptor,
                                                           activationGradients.GetTensorDescriptor(),
                                                           convDescriptors->LayerDescriptor,
                                                           activationGradientsBackwardDescriptor,
                                                           convWorkspace->AlgorithmBackward,
                                                           &convWorkspace->BackwardWorkspaceSize));

   if (convWorkspace->BackwardWorkspaceSize) cudaMalloc(&convWorkspace->BackwardWorkspace, convWorkspace->BackwardWorkspaceSize*sizeof(AFloat));
   if (convWorkspace->BackwardWorkspaceSize > 0 && convWorkspace->BackwardWorkspace == nullptr  ) {
      std::cerr << "Error allocating BACKW DATA CONV workspace of size " << convWorkspace->BackwardWorkspaceSize << " - probably running out of memory on the GPU"
      << std::endl;
      std::cout << " layer input shape is  { " << L->GetBatchSize() << " , " <<  L->GetInputDepth() << " , "
                                                <<L->GetInputHeight() << " , " << L->GetInputWidth() << " } " << std::endl;
      //inputTensor.PrintShape("inputTensor");
      R__ASSERT(false);
   }
   // Filter gradient
   //Tensor_t activationBackward ({L->GetBatchSize(), L->GetInputDepth(), L->GetInputHeight(), L->GetInputWidth()}, MemoryLayout::RowMajor, 0, 0);
   // here should be able to use inputTensorDescriptor
   cudnnTensorDescriptor_t activationBackwardDescriptor = inputTensorDescriptor;

#if (CUDNN_VERSION >= 8000)
   /**
    * I'm sure there may be a faster way, but this works
    */
   CUDNNCHECK(cudnnGetConvolutionBackwardFilterAlgorithmMaxCount(cudnnHandle, &convRequestedAlgoCount))  // ask cuDNN how much it can try
   cudnnConvolutionBwdFilterAlgoPerf_t convBwdFilterPerfResults[convRequestedAlgoCount];  // this will store metrics to choose convolution algorithm
   CUDNNCHECK(cudnnFindConvolutionBackwardFilterAlgorithm(
      cudnnHandle,
      activationBackwardDescriptor,
      activationGradients.GetTensorDescriptor(),
      convDescriptors->LayerDescriptor,
      convDescriptors->WeightsDescriptor,
      convRequestedAlgoCount,
      &algoCount,
      convBwdFilterPerfResults));
   // we could also do it with the expert mode (cudnnFindConvolutionForwardAlgorithmEx),
   // i.e.
   // CUDNNCHECK(cudnnFindConvolutionBackwardFilterAlgorithmEx(
   //    cudnnHandle,
   //    activationBackwardDescriptor,
   //    &inputTensor,
   //    activationGradients.GetTensorDescriptor(),
   //    &activationGradients,
   //    convDescriptors->LayerDescriptor,
   //    convDescriptors->WeightsDescriptor,
   //    &filters,
   //    convRequestedAlgoCount,
   //    &algoCount,
   //    &convBwdFilterPerfResults,
   //    &convWorkspace,
   //    memLimit));  // use memLimit for workspace size
   // instead choose either fastest or lowest memory algo as per preference
   LocalPerf bwdFilterPerfResults{convBwdFilterPerfResults};
   convWorkspace->HelperAlgorithm = convBwdFilterPerfResults[bwdFilterPerfResults.getIdx(algoChoice, algoCount, memLimit)].algo;
#else
   CUDNNCHECK(cudnnGetConvolutionBackwardFilterAlgorithm(cudnnHandle,
                                                         activationBackwardDescriptor,
                                                         activationGradients.GetTensorDescriptor(),
                                                         convDescriptors->LayerDescriptor,
                                                         convDescriptors->WeightsDescriptor,
                                                         preferenceBwdFilter,
                                                         memLimit,
                                                         &convWorkspace->HelperAlgorithm));
#endif

   std::cout << "CONV BWD Filter Algo used  is " << convWorkspace->HelperAlgorithm << std::endl;

   if (CNNOptions::ConvBwdFilterAlgorithm > 0) {
      convWorkspace->HelperAlgorithm = (cudnnConvolutionBwdFilterAlgo_t)CNNOptions::ConvBwdFilterAlgorithm;
      std::cout << " but force using " << convWorkspace->HelperAlgorithm << std::endl;
   }

   //CUDNNCHECK(cudnnSetConvolutionMathType(convDescriptors->LayerDescriptor, CUDNN_TENSOR_OP_MATH));

   CUDNNCHECK(cudnnGetConvolutionBackwardFilterWorkspaceSize(
         cudnnHandle, activationBackwardDescriptor, activationGradients.GetTensorDescriptor(),
         convDescriptors->LayerDescriptor, convDescriptors->WeightsDescriptor, convWorkspace->HelperAlgorithm,
         &convWorkspace->HelperWorkspaceSize));

   if (convWorkspace->HelperWorkspaceSize)
         cudaMalloc(&convWorkspace->HelperWorkspace, convWorkspace->HelperWorkspaceSize * sizeof(AFloat));

   if (convWorkspace->HelperWorkspaceSize > 0 && convWorkspace->HelperWorkspace == nullptr) {
      std::cerr << "Error allocating BACKW FILTER CONV workspace of size " << convWorkspace->BackwardWorkspaceSize
                << " - probably running out of memory on the GPU" << std::endl;
      std::cout << " layer input shape is  { " << L->GetBatchSize() << " , " << L->GetInputDepth() << " , "
                << L->GetInputHeight() << " , " << L->GetInputWidth() << " } " << std::endl;
      filters.PrintShape("filterTensor");
      R__ASSERT(false);
   }

   /// allocate workspace and descriptor for reduction operation
   // used to compiute bias gradients
   // try reducing the tensor

   CUDNNCHECK(cudnnCreateReduceTensorDescriptor(&convWorkspace->fReduceTensorDesc));

   auto reduceTensorDesc = convWorkspace->fReduceTensorDesc;
   CUDNNCHECK(cudnnSetReduceTensorDescriptor(reduceTensorDesc, CUDNN_REDUCE_TENSOR_ADD, Tensor_t::GetDataType(),
                                             CUDNN_PROPAGATE_NAN, CUDNN_REDUCE_TENSOR_NO_INDICES, CUDNN_32BIT_INDICES));

   CUDNNCHECK(cudnnGetReductionWorkspaceSize(cudnnHandle, reduceTensorDesc, activationGradients.GetTensorDescriptor(),
                                             biasGradients.GetTensorDescriptor(),
                                             &convWorkspace->fReductionWorkspaceSize));
   if (convWorkspace->fReductionWorkspaceSize > 0)
      cudaMalloc(&convWorkspace->fReductionWorkspace, convWorkspace->fReductionWorkspaceSize);


   // size_t isizeInBytes;
   // void *iSpace = nullptr;
   // CUDNNCHECK(cudnnGetReductionIndicesSize(cudnnHandle, reduceTensorDesc, activationGradients.GetTensorDescriptor(),
   //                                         biasGradients.GetTensorDescriptor(), &isizeInBytes));

   // if (isizeInBytes > 0)
   //    cudaMalloc(&convWorkspace->fIndiceWorkspace, isizeInBytes);

   workspace = convWorkspace;

   CUDNNCHECK(cudnnDestroyTensorDescriptor(inputTensorDescriptor));
}


//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::InitializePoolDropoutWorkspace(TWorkspace * & workspace,
                                             TDescriptors * & descriptors,
                                             const DNN::CNN::TConvParams & /*params*/,
                                             PoolingLayer_t *L) {

   auto poolWorkspace = new PoolingWorkspace_t ();
   auto poolDescriptors = static_cast<PoolingDescriptors_t *>(descriptors);

   //Tensor_t inputTensor ({L->GetBatchSize(), L->GetInputDepth(), L->GetInputHeight(), L->GetInputWidth()}, MemoryLayout::RowMajor, 0, 0);
   cudnnHandle_t cudnnHandle = L->GetOutput().GetCudnnHandle();

   // create tensor descriptors
   cudnnTensorDescriptor_t  inputTensorDescriptor;
   CUDNNCHECK(cudnnCreateTensorDescriptor(&inputTensorDescriptor) );
   CUDNNCHECK(cudnnSetTensor4dDescriptor(inputTensorDescriptor,
                                             CUDNN_TENSOR_NCHW,// Layout of the tensor in memory
                                             Tensor_t::GetDataType(),
                                             (int)L->GetBatchSize(),
                                             (int)L->GetInputDepth(),
                                             (int)L->GetInputHeight(),
                                             (int)L->GetInputWidth() ) );


   // Space needed to execute forward and backward dropout pass
   CUDNNCHECK(cudnnDropoutGetReserveSpaceSize(inputTensorDescriptor,
                                              &poolWorkspace->HelperWorkspaceSize));

   if (poolWorkspace->HelperWorkspaceSize) {
      cudaMalloc(&poolWorkspace->HelperWorkspace, poolWorkspace->HelperWorkspaceSize * sizeof(AFloat));
      if (poolWorkspace->HelperWorkspace == nullptr) {
         std::cerr << "Error allocating POOL reserved droput workspace of size " <<  poolWorkspace->HelperWorkspaceSize
                  << " probably running out of memory on the GPU"
                   << std::endl;
         std::cout << " layer input shape is  { " << L->GetBatchSize() << " , " << L->GetInputDepth() << " , "
                   << L->GetInputHeight() << " , " << L->GetInputWidth() << " } " << std::endl;
         R__ASSERT(false);
      }
   }

   // Space that contain random pass
   CUDNNCHECK(cudnnDropoutGetStatesSize(cudnnHandle,
                                        &poolWorkspace->ForwardWorkspaceSize));

   if (poolWorkspace->ForwardWorkspaceSize) {
      cudaMalloc(&poolWorkspace->ForwardWorkspace, poolWorkspace->ForwardWorkspaceSize * sizeof(AFloat));
      if (poolWorkspace->ForwardWorkspace == nullptr) {
         std::cerr << "Error allocating POOL droput state of size " <<  poolWorkspace->ForwardWorkspaceSize <<
         " probably running out of memory on the GPU"  << std::endl;
         std::cout << " layer input shape is  { " << L->GetBatchSize() << " , " << L->GetInputDepth() << " , "
                   << L->GetInputHeight() << " , " << L->GetInputWidth() << " } " << std::endl;
         R__ASSERT(false);
      }
   }

   // Fill the dropout workspace with random numbers and copy to device
   TRandom &  rand = TCudnn<AFloat>::GetRandomGenerator();
   // create a 64 bit seed using 2 32 bits integers
   unsigned long long seed = (unsigned long long) rand.Integer(UINT_MAX) << 32 + rand.Integer(UINT_MAX);
   // Reset the descriptor at every forward pass, so that random states get newly initialized?
   CUDNNCHECK(cudnnSetDropoutDescriptor(poolDescriptors->HelperDescriptor,
                                        cudnnHandle,
                                        L->GetDropoutProbability(),
                                        poolWorkspace->ForwardWorkspace,
                                        poolWorkspace->ForwardWorkspaceSize,
                                        seed));

   workspace = poolWorkspace;

   CUDNNCHECK(cudnnDestroyTensorDescriptor(inputTensorDescriptor));
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::FreeConvWorkspace(TWorkspace * workspace) {
   if (!workspace) return;
   auto convWorkspace = static_cast<ConvWorkspace_t *>(workspace);

   if(convWorkspace->ForwardWorkspace)  cudaFree(convWorkspace->ForwardWorkspace);
   if(convWorkspace->BackwardWorkspace) cudaFree(convWorkspace->BackwardWorkspace);
   if(convWorkspace->HelperWorkspace)   cudaFree(convWorkspace->HelperWorkspace);

   CUDNNCHECK(cudnnDestroyReduceTensorDescriptor(convWorkspace->fReduceTensorDesc));

   if (convWorkspace->fReductionWorkspace)
      cudaFree(convWorkspace->fReductionWorkspace);

}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::FreePoolDropoutWorkspace(TWorkspace * workspace) {
   if (!workspace) return;
   auto poolWorkspace = static_cast<PoolingWorkspace_t *>(workspace);

   if(poolWorkspace->ForwardWorkspace)  cudaFree(poolWorkspace->ForwardWorkspace);
   if(poolWorkspace->BackwardWorkspace) cudaFree(poolWorkspace->BackwardWorkspace);
   if(poolWorkspace->HelperWorkspace)   cudaFree(poolWorkspace->HelperWorkspace);
}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::BatchNormLayerForwardTraining(int axis, const Tensor_t &x,
                                                    Tensor_t & y,
                                                    Matrix_t &gamma, Matrix_t &beta,
                                                    Matrix_t & mean, Matrix_t &, Matrix_t & iVariance,
                                                    Matrix_t & runningMeans, Matrix_t & runningVars,
                                                    Scalar_t nTrainedBatches, Scalar_t momentum, Scalar_t epsilon,
                                                    const TensorDescriptor_t & bnParDescriptor )

{
   AFloat a = 1.0;
   AFloat b = 0.0;
   //cudnnBatchNormMode_t    bnMode = CUDNN_BATCHNORM_PER_ACTIVATION;
   //if (axis == 1) bnMode = CUDNN_BATCHNORM_SPATIAL;
   cudnnBatchNormMode_t bnMode = CUDNN_BATCHNORM_SPATIAL;

   //x.PrintShape("x");
   //y.PrintShape("y");

   // the factor is defined in Cudnn as 1-momentum
   double exponentialAverageFactor = (momentum < 0.) ? 1. / (1 + nTrainedBatches) :  1. - momentum;
   CUDNNCHECK(cudnnBatchNormalizationForwardTraining(x.GetCudnnHandle(), bnMode,
                                                      &a, &b,
                                                      x.GetTensorDescriptor(), x.GetDataPointer(),
                                                      y.GetTensorDescriptor(), y.GetDataPointer(),
                                                      bnParDescriptor,
                                                      gamma.GetDataPointer(), beta.GetDataPointer(),
                                                      exponentialAverageFactor,
                                                      runningMeans.GetDataPointer(),
                                                      runningVars.GetDataPointer(),
                                                      epsilon, mean.GetDataPointer(), iVariance.GetDataPointer() ) );

}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::BatchNormLayerForwardInference(int axis, const Tensor_t &x, Matrix_t &gamma, Matrix_t &beta,
                                                    Tensor_t &y, const Matrix_t &runningMeans,
                                                    const Matrix_t &runningVars, Scalar_t epsilon,
                                                    const TensorDescriptor_t & bnParDescriptor)

{
   AFloat a = 1.0;
   AFloat b = 0.0;

   cudnnBatchNormMode_t bnMode = CUDNN_BATCHNORM_SPATIAL;


   CUDNNCHECK(cudnnBatchNormalizationForwardInference(x.GetCudnnHandle(), bnMode,
                                                      &a, &b,
                                                      x.GetTensorDescriptor(), x.GetDataPointer(),// pass y as descriptor
                                                      y.GetTensorDescriptor(), y.GetDataPointer(),
                                                      bnParDescriptor,
                                                      gamma.GetDataPointer(), beta.GetDataPointer(),
                                                      runningMeans.GetDataPointer(),
                                                      runningVars.GetDataPointer(),
                                                      epsilon) );

}

//____________________________________________________________________________
template <typename AFloat>
void TCudnn<AFloat>::BatchNormLayerBackward(int axis, const Tensor_t &x, const Tensor_t &dy, Tensor_t &dx,
                                             Matrix_t &gamma, //  Matrix_t &beta, (not needed)
                                             Matrix_t &dgamma, Matrix_t &dbeta, const Matrix_t &mean,
                                             const Matrix_t &variance, const Matrix_t &iVariance,
                                             Scalar_t epsilon, const TensorDescriptor_t & bnParDescriptor)
{
   AFloat a = 1.0;
   AFloat b = 0.0;
   cudnnBatchNormMode_t bnMode = CUDNN_BATCHNORM_SPATIAL;

   CUDNNCHECK(cudnnBatchNormalizationBackward(x.GetCudnnHandle(), bnMode,
                                                      &a, &b, &a, &b,
                                                      x.GetTensorDescriptor(), x.GetDataPointer(),
                                                      dy.GetTensorDescriptor(), dy.GetDataPointer(),
                                                      dx.GetTensorDescriptor(), dx.GetDataPointer(),
                                                      bnParDescriptor, gamma.GetDataPointer(),
                                                      dgamma.GetDataPointer(), dbeta.GetDataPointer(),
                                                      epsilon, mean.GetDataPointer(), iVariance.GetDataPointer() ) );

}


template <typename AFloat>
void TCudnn<AFloat>::ConvLayerForward(Tensor_t & outputTensor,
                                      Tensor_t & inputActivation,
                                      const Tensor_t & input,
                                      const Matrix_t & weights, const Matrix_t & biases,
                                      const DNN::CNN::TConvParams & params,
                                      EActivationFunction activFunc,
                                      Tensor_t & inputPrime,
                                      const ConvDescriptors_t & descriptors,
                                      ConvWorkspace_t & workspace)
//                                    const AFloat alpha,
//                                    const AFloat beta)
{
   //((Tensor_t & )input).Reshape( {params.batchSize, params.inputDepth, params.inputHeight, params.inputWidth});

   assert( input.GetLayout() == GetTensorLayout());

   //size_t outputHeight =  DNN::CNN::TConvLayer<TCudnn<AFloat>>::calculateDimension(params.inputHeight, params.filterHeight, params.paddingHeight, params.strideRows);
   //size_t outputWidth =  DNN::CNN::TConvLayer<TCudnn<AFloat>>::calculateDimension(params.inputWidth, params.filterWidth, params.paddingWidth, params.strideCols);

   // PrintTensor(input,"input");
   // PrintTensor(outputTensor,"output");
   // PrintTensor(weights,"weights");
   // PrintTensor(biases,"biases");
   //((Tensor_t & )weights).Reshape( { params.numberFilters, params.inputDepth, params.filterHeight, params.filterWidth } );
   //((Tensor_t & )biases).Reshape(  { 1,params.numberFilters, 1, 1});
   //biases.Reshape ( { 1,params.numberFilters, 1, 1});

   AFloat alpha = 1.0;
   AFloat beta  = 0.0;
   cudnnHandle_t cudnnHandle = input.GetCudnnHandle();

   // check descriptors
#ifndef NDEBUG
   int n,c,h,w = 0;
   int s1,s2,s3,s4 = 0;
   cudnnDataType_t  dataType;
   cudnnGetTensor4dDescriptor( input.GetTensorDescriptor(), &dataType,&n,&c,&h,&w,&s1,&s2,&s3,&s4 );
   std::vector<size_t>  shape_input = {size_t(n), size_t(c) , size_t(h), size_t(w) };
   assert (shape_input == input.GetShape());

   cudnnGetTensor4dDescriptor( outputTensor.GetTensorDescriptor(), &dataType,&n,&c,&h,&w,&s1,&s2,&s3,&s4 );
   std::vector<size_t>  shape_output = {size_t(n), size_t(c) , size_t(h), size_t(w) };
   assert (shape_output == outputTensor.GetShape());
#endif


   // Perform convolution
   cudnnStatus_t status =  cudnnConvolutionForward(cudnnHandle,
                                      &alpha,
                                      input.GetTensorDescriptor(),
                                      input.GetDataPointer(),
                                      descriptors.WeightsDescriptor,
                                      weights.GetDataPointer(),
                                      descriptors.LayerDescriptor,
                                      workspace.AlgorithmForward,
                                      workspace.ForwardWorkspace,
                                      workspace.ForwardWorkspaceSize,
                                      &beta,
                                      inputActivation.GetTensorDescriptor(),
                                      inputActivation.GetDataPointer());

   // Apply biases
   assert(status == CUDNN_STATUS_SUCCESS);
   CUDNNCHECK(status);

   //PrintTensor(biases,"biases");
   //PrintTensor(outputTensor,"tensor before biases");
   AddConvBiases(inputActivation, biases);

   //PrintTensor(outputTensor,"tensor after biases");

   // Store the conv output before application of activation to use in the backward pass
   //TCudnn<AFloat>::Copy(outputTensor,inputActivation);

   // Apply activation
   TCudnn<AFloat>::ActivationFunctionForward(outputTensor, inputActivation, activFunc, descriptors.HelperDescriptor, 0.0, 1.0, 0.0);


   //perform convolution + biases + activation in one single call
    // could use an extra vector z but here we just pass a dummy tensor
   // AFloat alpha2 = 0;
   // CUDNNCHECK(cudnnConvolutionBiasActivationForward(cudnnHandle,
   //             &alpha,
   //             input.GetTensorDescriptor(),
   //             input.GetDataPointer(),
   //             descriptors.WeightsDescriptor,
   //             weights.GetDataPointer(),
   //             descriptors.LayerDescriptor,
   //             workspace.AlgorithmForward,
   //             workspace.ForwardWorkspace,
   //             workspace.ForwardWorkspaceSize,
   //             &alpha2,
   //             outputTensor.GetTensorDescriptor(),  // z : not used
   //             outputTensor.GetDataPointer(),
   //             biases.GetDescriptor(),
   //             biases.GetDataPointer(),
   //             descriptors.HelperDescriptor,
   //    outputTensor.GetTensorDescriptor(),
   //    outputTensor.GetDataPointer()));

   //TCudnn<AFloat>::PrintTensor(outputTensor, "after activation");

}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::ConvLayerBackward(Tensor_t &activationGradientsBackward,
                                       Matrix_t &weightGradients, Matrix_t &biasGradients,
                                       Tensor_t &inputActivation,
                                       Tensor_t &activationGradients,
                                       const Matrix_t &weights,
                                       const Tensor_t &activationBackward,
                                       const Tensor_t &outputTensor,
                                       EActivationFunction activFunc,
                                       const ConvDescriptors_t & descriptors,
                                       ConvWorkspace_t & workspace,
                                       size_t /*batchSize*/,   size_t /*inputHeight*/,
                                       size_t /*inputWidth*/,  size_t /*depth*/,
                                       size_t /*height*/,      size_t /*width*/,
                                       size_t /*filterDepth*/, size_t /*filterHeight*/,
                                       size_t /*filterWidth*/, size_t /*nLocalViews*/)
{
   // activationGradients.Reshape( outputTensor.GetShape());
   // weightGradients.Reshape( weights.GetShape());
   // biasGradients.Reshape({ 1, outputTensor.GetShape()[1], 1, 1});   // second output dimension is number of filters
   // // activationGradientsBackward.Reshape()
   // activationBackward.Reshape

   //--------------------------------------------------------------------------
   // Activation function gradient
   //--------------------------------------------------------------------------

   // x  : Output of previous layer without activation function             -> inputActivation
   // dx : Activation gradient to be computed                               -> activationGradients [in place op]
   // y  : Ouput of this layer (activation applied)                         -> outputTensor
   // dy : Gradient of activation from the following layer (backpropagation)-> activationGradients

   //if (descriptors.HelperDescriptor)
   ActivationFunctionBackward(activationGradients, outputTensor, activationGradients, inputActivation,
                              activFunc, descriptors.HelperDescriptor);  //y dy x dx

   //--------------------------------------------------------------------------
   // Network Activation gradient
   //--------------------------------------------------------------------------
   const AFloat alpha = 1.0;
   const AFloat beta  = 0.0;

   cudnnHandle_t cudnnHandle = outputTensor.GetCudnnHandle();

   //cudnnMathType_t math_type = CUDNN_TENSOR_OP_MATH; // : CUDNN_DEFAULT_MATH);
   // if using tensor math (cudnn version > 7)
   //CUDNNCHECK(cudnnSetConvolutionMathType(descriptors.LayerDescriptor, math_type));

   // do not compute activation gradients for first layer (i.e. when input activationGradientBackward is a dummy empty tensor)
   if (activationGradientsBackward.GetSize() > 0)
      CUDNNCHECK(cudnnConvolutionBackwardData(cudnnHandle,
                                           &alpha,
                                           descriptors.WeightsDescriptor,
                                           weights.GetDataPointer(),
                                           activationGradients.GetTensorDescriptor(),
                                           activationGradients.GetDataPointer(),
                                           descriptors.LayerDescriptor,
                                           workspace.AlgorithmBackward,
                                           workspace.BackwardWorkspace,
                                           workspace.BackwardWorkspaceSize,
                                           &beta,
                                           activationGradientsBackward.GetTensorDescriptor(),
                                           activationGradientsBackward.GetDataPointer()));

    //--------------------------------------------------------------------------
    // Filter gradient
    //--------------------------------------------------------------------------

   //CUDNNCHECK(cudnnSetConvolutionMathType(descriptors.LayerDescriptor, CUDNN_TENSOR_OP_MATH));

   CUDNNCHECK(cudnnConvolutionBackwardFilter(
      cudnnHandle, &alpha, activationBackward.GetTensorDescriptor(), activationBackward.GetDataPointer(),
      activationGradients.GetTensorDescriptor(), activationGradients.GetDataPointer(), descriptors.LayerDescriptor,
      workspace.HelperAlgorithm, workspace.HelperWorkspace, workspace.HelperWorkspaceSize, &beta,
      descriptors.WeightsDescriptor, weightGradients.GetDataPointer()));

   //--------------------------------------------------------------------------
   // Bias gradient
   //--------------------------------------------------------------------------
#if 0
   CUDNNCHECK(cudnnConvolutionBackwardBias(cudnnHandle, &alpha, activationGradients.GetTensorDescriptor(),
                                           activationGradients.GetDataPointer(), &beta,
                                           biasGradients.GetTensorDescriptor(), biasGradients.GetDataPointer()));

   //PrintTensor(biasGradients,"computed gradient biases");
#endif
#if 0
// try trandforming the activation gradient tensor and reduce it
   auto shape = activationGradients.GetShape();
   Tensor_t actGradTransf({shape[1], shape[0], shape[2], shape[3]}, activationGradients.GetLayout());
   CUDNNCHECK(cudnnTransformTensor(cudnnHandle, &alpha, activationGradients.GetTensorDescriptor(),
                                   activationGradients.GetDataPointer(), &beta, actGradTransf.GetTensorDescriptor(),
                                   actGradTransf.GetDataPointer()));
   // now make the operation
   TCudaMatrix<AFloat> actGradMatrix(actGradTransf.GetDeviceBuffer(), shape[0] * shape[2] * shape[3], shape[1]);
   TCudaMatrix<AFloat> temp(biasGradients.GetDeviceBuffer(), biasGradients.GetShape()[1], 1);
   TCuda<AFloat>::SumColumns(temp, actGradMatrix);
#endif

   // to compute the bias gradients is more efficient to reduce the activation gradient tensor in B,H,W dimensions and
   // adding their elements.
   // we create the descriptor for reduction and necessary workspaces in the initialization workspace function for the
   // convolution layer. Note that the indices workspace is not needed in case of addition operation

   CUDNNCHECK(cudnnReduceTensor(cudnnHandle, workspace.fReduceTensorDesc, nullptr, 0, workspace.fReductionWorkspace,
                                workspace.fReductionWorkspaceSize, &alpha, activationGradients.GetTensorDescriptor(),
                                activationGradients.GetDataPointer(), &beta, biasGradients.GetTensorDescriptor(),
                                biasGradients.GetDataPointer()));


#if 0   // this is very slow for large batches
   biasGradients.Zero();
   TCudaMatrix<AFloat> temp(biasGradients.GetShape()[1], 1);
   TCudaMatrix<AFloat> biasGradMatrix(biasGradients.GetDeviceBuffer(), biasGradients.GetShape()[1], 1);
   size_t batchSize = activationGradients.GetFirstSize();
   for (size_t event = 0; event < batchSize; event++) {
      TCudaTensor<AFloat> actGrad = activationGradients.At(event);  // remember tihs is a rowwise tensor
      TCudaMatrix<AFloat> actGradMatrix(actGrad.GetDeviceBuffer(),
                                        activationGradients.GetShape()[2] * activationGradients.GetShape()[3],
                                        activationGradients.GetShape()[0]);
      TCuda<AFloat>::SumColumns(temp, actGradMatrix);
      TCuda<AFloat>::ScaleAdd(biasGradMatrix, temp);
   }
#endif
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::AddConvBiases(Tensor_t &output,
                                   const Tensor_t &biases)
{
   TCudnn<AFloat>::ScaleAdd(output, biases);
}


//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Downsampling function used as the forward propagation step of a
///        Max-Pooling layer.
///
/// \param[out] A The output matrix. Each row corresponds to a slice and each element
///             is the max within a receptive field.
/// \param[out] B The winning indices matrix. Each element is the index of the max element.
/// \param[in] C The input matrix. Each row is a slice.
/// \param[in] imgHeight The heigh of the input.
/// \param[in] imgWidth The output of the input.
/// \param[in] fltHeight Height of the kernel.
/// \param[in] fltWidth Width of the kernel.
/// \param[in] strideRows stride size in the horizontal dimension.
/// \param[in] strideCols stride size in the vertical dimension.
///
/// Each output element is the maximum of the receptive field. We also save the winning
/// indices to facilitate back-propagation - we need to know which input element influenced
/// the output and only apply the derivative correction to this particular element.
/// The slicing process is the same as in a convolutional layer, however padding is set to 0.
///////////////////////////////////////////////////////////////////////////////////////////////
template<typename AFloat>
void TCudnn<AFloat>::Downsample(Tensor_t &A,
                                Tensor_t &/*B*/,
                                const Tensor_t &C,
                                const PoolingDescriptors_t & descriptors,
                                PoolingWorkspace_t & workspace,
                                size_t imgHeight,
                                size_t imgWidth,
                                size_t fltHeight,
                                size_t fltWidth,
                                size_t strideRows,
                                size_t strideCols)
{
   const AFloat alpha = 1.0;
   const AFloat beta = 0.0;

   CUDNNCHECK(cudnnPoolingForward(C.GetCudnnHandle(),
                                  descriptors.LayerDescriptor,
                                  &alpha,
                                  C.GetTensorDescriptor(),
                                  C.GetDataPointer(),
                                  &beta,
                                  A.GetTensorDescriptor(),
                                  A.GetDataPointer()));
}

//____________________________________________________________________________
template<typename AFloat>
void TCudnn<AFloat>::MaxPoolLayerBackward(Tensor_t & activationGradientsBackward, //dx
                                         const Tensor_t & activationGradients, // dy
                                         const Tensor_t & indexMatrix,
                                         const Tensor_t & activationBackward,  //X
                                         const Tensor_t & outputTensor,        //Y
                                         const PoolingDescriptors_t & descriptors,
                                         PoolingWorkspace_t & workspace,
                                         size_t imgHeight,
                                         size_t imgWidth,
                                         size_t fltHeight,
                                         size_t fltWidth,
                                         size_t strideRows,
                                         size_t strideCols,
                                         size_t /* nLocalViews */)
{
   const AFloat alpha = 1.0;
   const AFloat beta = 0.0;
   // x  : Output of previous layer without                                 -> inputActivation
   // dx : Activation gradient to be computed                               -> activationGradientsBackward
   // y  : Ouput of this layer (activation applied)                         -> outputTensor
   // dy : Gradient of activation from the following layer (backpropagation)-> activationGradients
   CUDNNCHECK(cudnnPoolingBackward(outputTensor.GetCudnnHandle(),
                                   descriptors.LayerDescriptor,
                                   &alpha,
                                   outputTensor.GetTensorDescriptor(),
                                   outputTensor.GetDataPointer(),
                                   activationGradients.GetTensorDescriptor(),
                                   activationGradients.GetDataPointer(),
                                   activationBackward.GetTensorDescriptor(),
                                   activationBackward.GetDataPointer(),
                                   &beta,
                                   activationGradientsBackward.GetTensorDescriptor(),
                                   activationGradientsBackward.GetDataPointer()));
}

//____________________________________________________________________________
/*template<typename AFloat>
void TCudnn<AFloat>::Reshape(TCudaTensor<AFloat> &A, const TCudaTensor<AFloat> &B)
{
    dim3 blockDims = TDevice::BlockDims2D();
    dim3 gridDims  = TDevice::GridDims2D(A);
    cudaStream_t s = A.GetComputeStream();

    ::TMVA::DNN::Cuda::Reshape<<<gridDims, blockDims>>>(A.GetDataPointer(), B.GetDataPointer(),
                                                        A.GetNrows(), A.GetNcols(), B.GetNrows(), B.GetNcols());
}*/

//______________________________________________________________________________
/*template <typename AReal>
void TCudnn<AReal>::Rearrange(std::vector<TCudaTensor<AReal>> &out, const std::vector<TCudaTensor<AReal>> &in)
{
   // B x T x D out --- T x B x D in*/
   /*size_t B = out.size();
   size_t T = out[0].GetNrows();
   size_t D = out[0].GetNcols();
   if ((T != in.size()) || (B != in[0].GetNrows())
       || (D != in[0].GetNcols())) {
      std::cout << "Incompatible Dimensions\n"
         << in.size() << "x" << in[0].GetNrows() << "x" << in[0].GetNcols()
         << " --> " << B << "x" << T << "x" << D << "\n";
      return;
   }
   for (size_t i = 0; i < B; ++i) {
      for (size_t j = 0; j < T; ++j) {
         for (size_t k = 0; k < D; ++k) {
            out[i](j, k) = in[j](i, k);
         }
      }
   }
   return;
}*/

//____________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
/// \brief Flatten a vector of matrices into a single matrix.
///
/// \param[out] A Output matrix.
/// \param[in] B Input vector. Each element is a matrix to be concatenated.
///
//////////////////////////////////////////////////////////////////////////////////
template<typename AFloat>
void TCudnn<AFloat>::Flatten(TCudaTensor<AFloat> &A,
                            const TCudaTensor<AFloat> &B)
{
   TCuda<AFloat>::Flatten(A, B);
}

//_________________________________________TCudaMatrix<AFloat> weightGradMatrix = weight_gradients.GetMatrix(); ____________________________
///////////////////////////////////////////TCudaMatrix<AFloat> weightGradMatrix = weight_gradients.GetMatrix(); //////////////////////////////
/// \brief Deflatten a matrix into a vectorTCudaMatrix<AFloat> weightGradMatrix = weight_gradients.GetMatrix(); rices.
///
/// \param[out] A Output matrices. Each eleTCudaMatrix<AFloat> weightGradMatrix = weight_gradients.GetMatrix(); ll be a part of the input.
/// \param[in] B Input flat matrix.
///
//////////////////////////////////////////////////////////////////////////////////
template<typename AFloat>
void TCudnn<AFloat>::Deflatten(TCudaTensor<AFloat> &A,
                              const TCudaTensor<AFloat> &B)
{
   TCuda<AFloat>::Deflatten(A,B);
}

} // namespace DNN
} // namespace TMVA
