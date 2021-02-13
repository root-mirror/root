// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 13/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////
// Implementation of the TCudaTensor class. //
/////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cuda/CudaTensor.h"
#include "TMVA/DNN/Architectures/Cuda/Device.h"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace TMVA {
namespace DNN  {


// Static members.
//____________________________________________________________________________
#ifdef R__HAS_CUDNN
template<typename AFloat>
std::vector<cudnnHandle_t> TCudaTensor<AFloat>::fCudnnHandle(1);
template<typename AFloat>
cudnnDataType_t          TCudaTensor<AFloat>::fDataType         = CUDNN_DATA_FLOAT;
#endif

template<typename AFloat>
std::vector<int>         TCudaTensor<AFloat>::fInstances(1,0);

/// This information is needed for the multi-dimensional indexing. See here:
/// https://en.wikipedia.org/wiki/Row-_and_column-major_order
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.strides.html
template<typename AFloat>
std::vector<std::size_t> TCudaTensor<AFloat>::ComputeStridesFromShape(const std::vector<std::size_t> &shape,
   bool rowmajorLayout)
{
   const auto size = shape.size();
   std::vector<std::size_t> strides(size);
   if (rowmajorLayout)  {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[size - 1 - i] = 1;
         } else {
            strides[size - 1 - i] = strides[size - 1 - i + 1] * shape[size - 1 - i + 1];
         }
      }
   } else  {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[i] = 1;
         } else {
            strides[i] = strides[i - 1] * shape[i - 1];
         }
      }
   }
   return strides;
}

// Constructors.
//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor()
    : fShape(), fStrides(), fNDim(0), fSize(0), fElementBuffer(), fStreamIndx(0), fTensorDescriptor(nullptr)
{
   //InitializeCuda();
}


//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const std::vector<size_t> & shape,
                                 TCudaTensor::MemoryLayout layout,
                                 int device, int streamIndx)
    : fShape(shape), fStrides(shape.size()), fNDim(shape.size()), fDevice(device), fStreamIndx(streamIndx),
      fTensorDescriptor(nullptr), fMemoryLayout(layout)
{
   fStrides = ComputeStridesFromShape(fShape, layout==MemoryLayout::RowMajor);

   fSize = (layout==MemoryLayout::RowMajor) ? fStrides.front()*fShape.front() :
                                              fStrides.back()*fShape.back();

   // create a new buffer in this case
   fElementBuffer = TCudaDeviceBuffer<AFloat>(fSize, 0);
   // need to initialize Cuda when creating a new Cuda Buffer (e.g. create Tensor descriptor)
   InitializeCuda();
}

//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const AFloat * host_data, const std::vector<size_t> & shape,
                                 TCudaTensor::MemoryLayout layout,
                                 int device, int streamIndx)
   : TCudaTensor(shape, layout, device, streamIndx)
{
   // do I need to allocate this buffer ????
   // is not a mem leak
   // AFloat * buffer = new AFloat[fSize];
   // size_t index = 0;
   // for (size_t j = 0; j < fSize; ++j) {
   //       buffer[j] = static_cast<AFloat>(host_data[j]);
   //    }
   // }

   cudaMemcpy(fElementBuffer, host_data, fSize * sizeof(AFloat),
              cudaMemcpyHostToDevice);

   // no need to initialize cuda. Done in the other constructor that is called before
   //InitializeCuda();
}

//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(TCudaDeviceBuffer<AFloat> buffer,
                                 const std::vector<size_t> & shape,
                                 TMVA::Experimental::MemoryLayout layout,
                                 int device, int streamIndx)
   : fNDim(shape.size()), fElementBuffer(buffer), fShape(shape), fStrides( shape.size()), fDevice(device),
     fStreamIndx(streamIndx), fTensorDescriptor(nullptr), fMemoryLayout(layout)
{
   // constructor from an existing buffer . Buffer size must contain given size
   fStrides = ComputeStridesFromShape(fShape, layout==MemoryLayout::RowMajor);

   fSize = (layout==MemoryLayout::RowMajor) ? fStrides.front()*fShape.front() :
                                              fStrides.back()*fShape.back();
   R__ASSERT(fSize <= buffer.GetSize());

   // need to Initialize Cuda in case device buffer was created separatly
   InitializeCuda();
}

//____________________________________________________________________________
//FIXME: Go to shared_ptr implementation of instance tracking
// template <typename AFloat>
// TCudaTensor<AFloat>::TCudaTensor(const TCudaTensor<AFloat>& oldTensor) :
//    TCudaTensor(oldTensor.fShape, oldTensor.fMemoryLayout, oldTensor.fDevice, oldTensor.fStreamIndx)
// {
//    // No deep copy
//    fStrides       = oldTensor.fStrides;
//    fElementBuffer = oldTensor.fElementBuffer;

//    std::cout << "calling copy constructor of TCuda tensor" << std::endl;

//    InitializeCuda();
// }

//____________________________________________________________________________
template <typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const TCudaMatrix<AFloat>& matrix, size_t dim) :
   TCudaTensor( matrix.GetDeviceBuffer(), {matrix.GetNrows(), matrix.GetNcols()}, MemoryLayout::ColumnMajor)
{
   // No deep copy
   if (dim > 2) {
      // change shape from (nrows,ncols) to (nrows,ncols,1,1)
      // this works onlt for coolum major layout since this is same of TCudaMatrix
      fShape.insert(fShape.end(), dim-2, 1);
      fStrides.insert(fStrides.end(),dim-2,fSize);
      fNDim = dim;
      // need to reset tensor descriptor since we are changing the shape
      SetTensorDescriptor();
   }
}



template<typename AFloat>
TCudaTensor<AFloat>::operator TMatrixT<AFloat>() const
{
   // this should work only for size 2 or 4 tensors
   if (GetLayout() == MemoryLayout::ColumnMajor &&
       (fNDim == 2 || (fNDim == 3 && GetFirstSize() == 1)) ) {
//         return TCudaMatrix<AFloat>(fElementBuffer, GetHSize(), GetWSize());
      TCudaMatrix<AFloat> temp = GetMatrix();
      return temp;
   }
   // we can convert directy to TMatrix
   // assert(fNDim <= 4);
   // size_t nRows = fShape[0]*fShape[1];
   // size_t nCols = fShape[2];
   // if (fNDim == 4) nCols*= fShape[3];

   if (GetLayout() == MemoryLayout::RowMajor) {

     // This assume that tensor D1, D2, D3,D4 is converted in   D1 , D2*D3*D4
     TMatrixT<AFloat> hostMatrix( GetNrows(), GetNcols() );
     cudaMemcpy(hostMatrix.GetMatrixArray(), fElementBuffer, fSize * sizeof(AFloat),
           cudaMemcpyDeviceToHost);
     return hostMatrix;

   }
   // else in case of column major tensor we need to transpose(this is what is done in TCudaMatrix)
   // Here we assume that D1, D2, D3 is converted in   a matrix  (D3, D1*D2)
   TMatrixT<AFloat> hostMatrix( GetNcols(), GetNrows() );
   cudaMemcpy(hostMatrix.GetMatrixArray(), fElementBuffer, fSize * sizeof(AFloat),
              cudaMemcpyDeviceToHost);
   return hostMatrix.T();  // return transpose matrix

}
#ifdef R__HAS_CUDNN
//____________________________________________________________________________
template <typename AFloat>
TCudaTensor<AFloat>::~TCudaTensor()
{
   if (fTensorDescriptor && fTensorDescriptor.use_count() == 1 ) {
      // //std::cout << "Destroy tensor descriptor for shape ";
      // for (int ii = 0; ii < fNDim; ++ii)
      //    std::cout << fShape[ii] << ",";
      // std::cout << std::endl;
      CUDNNCHECK(cudnnDestroyTensorDescriptor(fTensorDescriptor->fCudnnDesc));

      fInstances[fStreamIndx]--;

         // When all tensors in a streamIndx are destroyed, release cudnn resources
      if (fInstances[fStreamIndx] <= 0) {
         std::cout << "All Cuda tensors are -released - destroy cudnn handle " << fInstances[fStreamIndx] << std::endl;
         CUDNNCHECK(cudnnDestroy(fCudnnHandle[fStreamIndx]));
      }

   }
   //std::cout << "Tensor descriptor destroyed - instances are " << fInstances[fStreamIndx] << std::endl;

}
//____________________________________________________________________________
template <typename AFloat>
void TCudaTensor<AFloat>::InitializeCuda()
{
   // descriptor is needed for Cuddn tensor that are rowmajor
   if (!fTensorDescriptor && fSize > 0 && fNDim >= 2) {


      // if ((fInstances[fStreamIndx] < 4 && fInstances[fStreamIndx] > -4) || fInstances[fStreamIndx]%1000 == 0) {
      //    std::cout << " stream index " << fStreamIndx << " instances " << fInstances[fStreamIndx] << std::endl;
      //    PrintShape();
      // }


      // Also check whether a new streamIndx has been opened
      if (fInstances.size() - 1 < fStreamIndx) {
         // If need to resize once, need probably to resize more often
         fInstances.resize(2 * fStreamIndx + 1, 0);
         fCudnnHandle.resize(2 * fStreamIndx + 1, nullptr);
         }
      if (fInstances[fStreamIndx] == 0) {
         std::cout << "TCudaTensor::create cudnn handle ! " << std::endl;
         CUDNNCHECK(cudnnCreate(&fCudnnHandle[fStreamIndx]));
         // CUDNNCHECK(cudnnSetStream(fCudnnHandle[fStreamIndx], fElementBuffer.GetComputeStream()));

         // cublasCreate(&fCublasHandle);
         // CUDACHECK(cudaMalloc(& fDeviceReturn, sizeof(AFloat)));
         // CUDACHECK(cudaMalloc(& fCurandStates, TDevice::NThreads(*this)));
      }
      // if (TDevice::NThreads(*this) > (int) fNCurandStates) {
      //     fNCurandStates = TDevice::NThreads(*this);
      //     if (fCurandStates) {
      //         cudaFree(fCurandStates);
      //     }
      //     cudaMalloc(&fCurandStates, TDevice::NThreads(*this) * sizeof(curandState_t));
      //     InitializeCurandStates();
      // }

      // Prevent template specialization of entire class
      if (std::is_same<AFloat, double>::value) {
         fDataType = CUDNN_DATA_DOUBLE;
      } else if (std::is_same<AFloat, float>::value) {
         fDataType = CUDNN_DATA_FLOAT;
      }

      // create tensor descriptor
      fTensorDescriptor = std::make_shared<TensorDescriptor>();
      // std::cout << "create tensor  descriptor ! " << std::endl;
      CUDNNCHECK(cudnnCreateTensorDescriptor(&(fTensorDescriptor->fCudnnDesc)));

      // we increment instances when we create the descriptor
      fInstances[fStreamIndx]++;
   }

   SetTensorDescriptor();

}
template<typename AFloat>
void TCudaTensor<AFloat>::SetTensorDescriptor() {
      if (!fTensorDescriptor) return;
      if (fSize == 0) return;

      // cuDNN NdTensor format has a minsize of 4 tensor dimensions
      // 4D tensor is more performant on lower dimensions and supports all folowing operations
      // is this really true ???
      if (fNDim == 4 || fNDim > 1 && fMemoryLayout == MemoryLayout::ColumnMajor || fNDim == 2) {
         // pad cudnn tensor column major with extra elements (these are used in the convolutions)
         Shape_t shape = fShape;

         if (fNDim < 4 && fNDim > 1) {
            //    // add 1 to tensor
            if (fMemoryLayout == MemoryLayout::RowMajor)
               shape.insert(shape.end(), 4 - fNDim, 1);
            else
               shape.insert(shape.begin(), 4 - fNDim, 1);
         }

         if (fMemoryLayout == MemoryLayout::RowMajor) {
            auto status = cudnnSetTensor4dDescriptor(fTensorDescriptor->fCudnnDesc,
                                                     CUDNN_TENSOR_NCHW, // Layout of the tensor in memory
                                                     fDataType,
                                                     (int)shape[0],  // batch size
                                                     (int)shape[1],  // no. channels
                                                     (int)shape[2],  // image height
                                                     (int)shape[3]); // image width
            assert(status == CUDNN_STATUS_SUCCESS);
            CUDNNCHECK(status);
         } else {
            CUDNNCHECK(cudnnSetTensor4dDescriptor(fTensorDescriptor->fCudnnDesc,
                                                  CUDNN_TENSOR_NCHW, // Layout of the tensor in memory
                                                  fDataType,
                                                  (int)shape[3],   // batch size
                                                  (int)shape[2],   // no. channels
                                                  (int)shape[1],   // image height
                                                  (int)shape[0])); // image width
         }

         // Some operations in cudnn may not work with this tensor description
         // do not support tensors with dims < 1
      } else if (fNDim >2  || fNDim > 4) {
         // these are used in the RNN layers

         // seems to work for 3d tensor with row major (case of RNN tensors)
         // rnn wnats 3d tensors but it does not work for 2d tensors
         std::vector<int> shape(fShape.begin(), fShape.end());
         std::vector<int> strides(fStrides.begin(), fStrides.end());
         auto status = cudnnSetTensorNdDescriptor(fTensorDescriptor->fCudnnDesc, fDataType, (int)fNDim, shape.data(),
                                                  strides.data());
         assert(status == CUDNN_STATUS_SUCCESS);
         CUDNNCHECK(status);
      }

#ifdef NDEBUG
      size_t tensorSize;
      CUDNNCHECK(cudnnGetTensorSizeInBytes(fTensorDescriptor->fCudnnDesc, &tensorSize));
      assert(fSize == tensorSize/sizeof(AFloat));

        //    int n,c,h,w = 0;
   // int s1,s2,s3,s4 = 0;
   // cudnnDataType_t  dataType;
   // cudnnGetTensor4dDescriptor( fTensorDescriptor, &dataType,&n,&c,&h,&w,&s1,&s2,&s3,&s4 );
   // std::vector<size_t>  shape_input = {n,c,h,w};
   // assert (shape_input == GetShape());

#endif


}
#else // case ROOT has not Cudnn (add dummy implementations)
//____________________________________________________________________________
template <typename AFloat>
TCudaTensor<AFloat>::~TCudaTensor()
{}
//____________________________________________________________________________
template <typename AFloat>
void TCudaTensor<AFloat>::InitializeCuda()
{}
//____________________________________________________________________________
template<typename AFloat>
void TCudaTensor<AFloat>::SetTensorDescriptor()
{}

#endif

//____________________________________________________________________________
template<typename AFloat>
void TCudaTensor<AFloat>::InitializeCurandStates()
{
   // dim3 blockDims = TDevice::BlockDims2D();
   // dim3 gridDims  = TDevice::GridDims2D(*this);
   // CurandInitializationKernel<<<gridDims, blockDims>>>(time(nullptr), fCurandStates);
}

template<typename AFloat>
void TCudaTensor<AFloat>::Print(const char * name, bool truncate) const
{
      //TCudaBuffer<AFloat> hostBuffer (fSize);
      //fElementBuffer.CopyTo(hostBuffer);
    #if 0
      AFloat hostBuffer[fSize];

      cudaMemcpy(hostBuffer, fElementBuffer, fSize * sizeof(AFloat),
                 cudaMemcpyDeviceToHost);

      for (size_t i = 0; i < fSize; i++) std::cout << hostBuffer[i] << "  ";
   #endif
   PrintShape(name);
   size_t n = fSize;
   if (n > 10 && truncate) n = 10;
   std::cout << "Data : { ";
   for (size_t i = 0; i < n; ++i ) {
      AFloat * elementPointer = fElementBuffer + i;
      std::cout << AFloat( TCudaDeviceReference<AFloat>(elementPointer) );
      if (i < n-1) std::cout << " , ";
   }
   if (n < fSize) std::cout << "............   } ";
   std::cout << " } " << std::endl;
}
template<typename AFloat>
void TCudaTensor<AFloat>::PrintShape(const char * name) const
{
      std::string memlayout = (GetLayout() == MemoryLayout::RowMajor) ? "RowMajor" : "ColMajor";
      std::cout << name << " shape : { ";
      for (size_t i = 0; i < fNDim-1; ++i )
         std::cout << fShape[i] << " , ";
      std::cout << fShape.back() << " } " << " Layout : " << memlayout << std::endl;
}
#if 0
// Conversion to RTensor
//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::operator Experimental::RTensor<AFloat>() const
{
   std::vector<size_t> shape(fNDims, fNDims + fDim)

   Experimental::RTensor<AFloat> hostTensor( shape)

   AFloat * buffer = new AFloat[fSize];
   cudaMemcpy(buffer, fElementBuffer, fSize * sizeof(AFloat),
              cudaMemcpyDeviceToHost);

   int index = 0;
   for (int j = 0; j < fSize; j++) {
         hostTensor.GetData()[j] = static_cast<AFloat>(buffer[j]);
      }
   }

   delete[] buffer;
   return hostTensor;
}
#endif
// Explicit Instantiations.

template class TCudaTensor<float>;
template class TCudaTensor<double>;

} // namespace DNN
} // namespace TMVA
