// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 13/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////////
// Contains the TCudaMatrix class for the representation of matrices //
// on CUDA devices as well as the TCudaDeviceReference class which   //
// is a helper class to emulate lvalue references to floating point  //
// values on the device.                                             //
///////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_CUDA_CUDATENSOR
#define TMVA_DNN_ARCHITECTURES_CUDA_CUDATENSOR


#include <vector>
#include <cstring>
#include <cassert>
#include <iostream>

#include "CudaMatrix.h"
#include "TMatrixT.h"
#include "CudaBuffers.h"

//#include "TMVA/RTensor.hxx"

#ifdef R__HAS_CUDNN
#include "cudnn.h"
#define CUDNNCHECK(ans) {cudnnError((ans), __FILE__, __LINE__); }
#endif

namespace TMVA {



#ifndef TMVA_RTENSOR

namespace Experimental {
/// Memory layout type (copy from RTensor.hxx)
enum class MemoryLayout : uint8_t {
   RowMajor = 0x01,
   ColumnMajor = 0x02
};
}
#endif

namespace DNN {

using MemoryLayout = TMVA::Experimental::MemoryLayout;

#ifdef R__HAS_CUDNN
/**
 * Function to handle the status output of cuDNN function calls. See also
 * CUDACHECK in CudaMatrix.h.
 */
inline void cudnnError(cudnnStatus_t status, const char *file, int line, bool abort=true)
{
   if (status != CUDNN_STATUS_SUCCESS) {
      fprintf(stderr, "CUDNN Error: %s %s %d\n", cudnnGetErrorString(status), file, line);
      if (abort)
         exit(status);
   }
}
#endif
//____________________________________________________________________________
//
// Cuda Tensor
//____________________________________________________________________________

/** TCudaTensor Class
 *
 * The TCudaTensor class extends the TCudaMatrix class for dimensions > 2.
 *
 */
template<typename AFloat>
class TCudaTensor
{
public:

   using Shape_t = std::vector<size_t>;
   using MemoryLayout = TMVA::Experimental:: MemoryLayout;
   using Scalar_t = AFloat;


private:

#ifdef R__HAS_CUDNN
   struct TensorDescriptor {
       cudnnTensorDescriptor_t   fCudnnDesc;
   };

   static std::vector<cudnnHandle_t>     fCudnnHandle;      ///< Holds the cuddn library context (one for every CUDA stream)

   static cudnnDataType_t                fDataType;         ///< Cudnn datatype used for the tensor
#else
   struct TensorDescriptor {
   };
#endif

   /** For each GPU device keep the CUDA streams in which tensors are used.
     * Instances belonging to the same stream on the same deviceshare a
     * cudnn library handel to keep cudnn contexts seperated */
   //static std::vector<std::vector<int> > fInstances;
   static std::vector<int> fInstances;

   /** The shape vector (size of dimensions) needs to be ordered as no. channels,
    *  image dimensions.
    */
   Shape_t      fShape;            ///< spatial subdimensions
   Shape_t      fStrides;          ///< Strides between tensor dimensions (always assume dense, non overlapping tensor)
   size_t       fNDim;             ///< Dimension of the tensor (first dimension is the batch size, second is the no. channels)
   size_t       fSize;             ///< No. of elements
   int          fDevice;           ///< Device associated with current tensor instance
   int          fStreamIndx;       ///< Cuda stream associated with current instance

   std::shared_ptr<TensorDescriptor> fTensorDescriptor;
   TCudaDeviceBuffer<AFloat> fElementBuffer;

   MemoryLayout fMemoryLayout;



public:


   //static AFloat * GetOnes() {return fOnes;}

   TCudaTensor();

   TCudaTensor(const AFloat * data,
               const std::vector<size_t> & shape,
               MemoryLayout memlayout = MemoryLayout::ColumnMajor,
               int deviceIndx = 0, int streamIndx = 0);
   TCudaTensor(TCudaDeviceBuffer<AFloat> buffer,
               const std::vector<size_t> & shape,
               MemoryLayout memlayout = MemoryLayout::ColumnMajor,
               int deviceIndx = 0, int streamIndx = 0);
   TCudaTensor(const std::vector<size_t> & shape,
               MemoryLayout memlayout = MemoryLayout::ColumnMajor,
               int deviceIndx = 0, int streamIndx = 0);

   TCudaTensor(size_t bsize, size_t csize, size_t hwsize, MemoryLayout memlayout = MemoryLayout::ColumnMajor,  int deviceIndx = 0, int streamIndx = 0) :
      TCudaTensor( (memlayout == MemoryLayout::ColumnMajor) ? Shape_t({ csize, hwsize, bsize}) : Shape_t({ bsize, csize, hwsize }) , memlayout,
                   deviceIndx, streamIndx)
     {}

   TCudaTensor(size_t bsize, size_t csize, size_t hsize, size_t wsize, MemoryLayout memlayout = MemoryLayout::ColumnMajor,  int deviceIndx = 0, int streamIndx = 0) :

      TCudaTensor( {bsize, csize, hsize, wsize}, memlayout, deviceIndx, streamIndx)
     {
        if (memlayout == MemoryLayout::ColumnMajor)
           *this =  TCudaTensor(fElementBuffer, { csize, hsize, wsize, bsize}, memlayout, deviceIndx, streamIndx);
     }

   TCudaTensor(size_t n, size_t m, MemoryLayout memlayout = MemoryLayout::ColumnMajor,  int deviceIndx = 0, int streamIndx = 0) :
      //   TCudaTensor( {n,m}, memlayout, deviceIndx, streamIndx) :
      TCudaTensor( {n, m}, memlayout, deviceIndx, streamIndx)
     {}

   TCudaTensor(const TCudaMatrix<AFloat> & m, size_t dim = 2);

   TCudaTensor(const TMatrixT<AFloat> & m, size_t dim = 2) :
      TCudaTensor( TCudaMatrix<AFloat>(m), dim)
   {}

   TCudaTensor(TCudaDeviceBuffer<AFloat> buffer, size_t n, size_t m) :
         TCudaTensor( buffer, {n,m}, MemoryLayout::ColumnMajor ,0,0) {}

   TCudaTensor(const TCudaTensor &) = default;
   TCudaTensor(TCudaTensor &&) = default;
   TCudaTensor & operator=(const TCudaTensor  &) = default;
   TCudaTensor & operator=(      TCudaTensor &&) = default;
   ~TCudaTensor();

   /** Convert cuda matrix to Root TMatrix. Performs synchronous data transfer. */
   operator TMatrixT<AFloat>() const;


   MemoryLayout GetLayout() const { return fMemoryLayout; }

   const Shape_t & GetShape() const {return fShape;}
   const Shape_t & GetStrides() const {return fStrides;}
   size_t GetDimAt(size_t i) const {return fShape[i];}
   size_t GetNDim() const {return fNDim;}
   size_t GetSize() const {return fSize;}

   const AFloat * GetDataPointer() const {return fElementBuffer;}
   AFloat       * GetDataPointer()       {return fElementBuffer;}
   const AFloat * GetData() const {return fElementBuffer;}
   AFloat       * GetData()       {return fElementBuffer;}

   const AFloat * GetDataPointerAt(size_t i ) const {
      return (const_cast<TCudaDeviceBuffer<AFloat>&>(fElementBuffer)).GetSubBuffer(i * GetFirstStride(), GetFirstStride() ); }
   AFloat       * GetDataPointerAt(size_t i )       {return fElementBuffer.GetSubBuffer(i * GetFirstStride(), GetFirstStride() ); }


   const TCudaDeviceBuffer<AFloat> & GetDeviceBuffer()     const {return fElementBuffer;}
   TCudaDeviceBuffer<AFloat>       & GetDeviceBuffer()           {return fElementBuffer;}

#ifdef R__HAS_CUDNN
   const cudnnHandle_t             & GetCudnnHandle()      const {return fCudnnHandle[fStreamIndx];}
   const cudnnTensorDescriptor_t   & GetTensorDescriptor() const {return fTensorDescriptor->fCudnnDesc;}
   static cudnnDataType_t   GetDataType() { return fDataType; }
#endif

   cudaStream_t GetComputeStream() const {
      return fElementBuffer.GetComputeStream();
   }
   void         SetComputeStream(cudaStream_t stream) {
       fElementBuffer.SetComputeStream(stream);
   }

   bool isEqual (TCudaTensor<AFloat> & other) {

      if (fSize != other.GetSize()) return false;


      std::unique_ptr<AFloat[]> hostBufferThis(new AFloat[fSize]);
      std::unique_ptr<AFloat[]> hostBufferOther(new AFloat[fSize]);
      cudaMemcpy(hostBufferThis.get(), fElementBuffer, fSize * sizeof(AFloat),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(hostBufferOther.get(), other.GetDeviceBuffer(), fSize * sizeof(AFloat),
                 cudaMemcpyDeviceToHost);

      for (size_t i = 0; i < fSize; i++) {
         if (hostBufferThis[i] != hostBufferOther[i]) return false;
      }
      return true;
   }

   bool isEqual (const AFloat * hostBufferOther, size_t otherSize) {
      if (fSize != otherSize) return false;


      std::unique_ptr<AFloat[]> hostBufferThis(new AFloat[fSize]);
      cudaMemcpy(hostBufferThis.get(), fElementBuffer, fSize * sizeof(AFloat),
                 cudaMemcpyDeviceToHost);

      for (size_t i = 0; i < fSize; i++) {
         if (hostBufferThis[i] != hostBufferOther[i]) return false;
      }

      return true;
   }

   void Print(const char * name = "Tensor", bool truncate = false) const;

   void PrintShape(const char * name="Tensor") const;

   void Zero() {
      cudaMemset(GetDataPointer(), 0, sizeof(AFloat) * GetSize());
   }

   void SetConstVal(const AFloat constVal) {
      TCudaHostBuffer<AFloat> hostBuffer(fSize);
      hostBuffer.SetConstVal(constVal);
      fElementBuffer.CopyFrom(hostBuffer);
   }

   // have this tensor representatrions
   // 2-dimensional tensors  :  NW   where N is batch size W is the feature size . Memory layout should be columnwise in this case
   // 3 -dimensional tensor  : represnetation is NHWC  , tensor should be columnwise storage
   // 4 -dimensional tensor :   representation is NCHW  ande tensor should be row wose
   // a rowmajor tensor with dimension less than trhee should not exist but in case consider as a N, (CHW) for 2d, N, C, (HW) for 3d
   // a columnmajor tensor for dimension >=4 should not exist but in case consider as a N,H,W,C  (i.e. with shape C,W,H,N)

   size_t GetFirstSize() const {
      return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape.back() : fShape.front(); }  // CM order
   size_t GetFirstStride() const {
      return (GetLayout() == MemoryLayout::ColumnMajor ) ?  fStrides.back() : fStrides.front();  } // CM order

   size_t GetCSize() const {
      if  (fNDim == 2) return 1;
      return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape.front() : fShape[1] ; //assume NHWC
   }
   size_t GetHSize() const {
      if  (fNDim == 2) return fShape[0];
      if  (fNDim == 3) return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape[0] : fShape[1] ;// same as C
      if  (fNDim >= 4) return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape[2] : fShape[2] ;
      return 0;
   }
   size_t GetWSize() const {
      if  (fNDim == 2) return fShape[1];
      if  (fNDim == 3) return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape[1] : fShape[2] ;
      if  (fNDim == 4) return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape[3] : fShape[3] ;
      return 0;
   }

   // for backward compatibility (assume column-major
   // for backward compatibility : for CM tensor (n1,n2,n3,n4) -> ( n1*n2*n3, n4)
   //                              for RM tensor (n1,n2,n3,n4) -> ( n2*n3*n4, n1 ) ???
   size_t GetNrows() const { return (GetLayout() == MemoryLayout::ColumnMajor ) ? fStrides.back() : fShape.front();}
   size_t GetNcols() const { return (GetLayout() == MemoryLayout::ColumnMajor ) ? fShape.back() : fStrides.front(); }


   // Matrix conversion for tensors of shape 2
   TCudaMatrix<AFloat> GetMatrix() const  {
      // remember TCudaMatrix is always column-major
      if ( GetLayout() == MemoryLayout::ColumnMajor &&
           (fNDim == 2 || (fNDim == 3 && GetFirstSize() == 1) ) )
         return TCudaMatrix<AFloat>(fElementBuffer, GetHSize(), GetWSize());


      //case of N,M,1,1,..
      bool caseNM11 = true;
      for (size_t i = 2; i < fNDim; ++i)  caseNM11 &= fShape[i] == 1;
      if (caseNM11) {
         return (GetLayout() == MemoryLayout::ColumnMajor ) ?
            TCudaMatrix<AFloat>(fElementBuffer, fShape[0], fShape[1]) :
            TCudaMatrix<AFloat>(fElementBuffer, fShape[1], fShape[0]);
      }
      bool case11NM = true;
      for (size_t i = 0; i < fNDim-2; ++i)  case11NM &= fShape[i] == 1;
      if (case11NM) {
         return  (GetLayout() == MemoryLayout::ColumnMajor ) ?
            TCudaMatrix<AFloat>(fElementBuffer, fShape[fNDim-2], fShape[fNDim-1]) :
            TCudaMatrix<AFloat>(fElementBuffer, fShape[fNDim-1], fShape[fNDim-2]);
      }

      assert(false);
      return TCudaMatrix<AFloat>();
   }

   // for backeard compatibility with old tensor
   TCudaMatrix<AFloat> operator[](size_t i) const {
      //assert(GetLayout() == MemoryLayout::ColumnMajor );
      return At(i).GetMatrix();
   }



   static inline std::vector<std::size_t> ComputeStridesFromShape(const std::vector<std::size_t> &shape,
   bool rowmajorLayout);

   void ReshapeInPlace(const Shape_t & newShape)  {
      fShape   = newShape;
      fStrides = ComputeStridesFromShape(fShape, fMemoryLayout == MemoryLayout::RowMajor);
      fNDim = fShape.size();
      // in principle reshape should not change tensor size
      size_t newSize = (fMemoryLayout == MemoryLayout::RowMajor) ? fStrides.front() * fShape.front() : fStrides.back() * fShape.back();
      R__ASSERT(newSize <= fSize);
      fSize = newSize;
      // reset the descritor for Cudnn
      SetTensorDescriptor();
   }

   TCudaTensor<AFloat> Reshape(const Shape_t & newShape) const {
      TCudaTensor<AFloat> tmp(this->GetDeviceBuffer(), newShape, this->GetLayout(), fDevice, fStreamIndx);
      return tmp;
   }

   void SetTensorDescriptor();

   // return slice of tensor
   // return slices in the first dimension (if row wise) or last dimension if colun wise
   // so single event slides
   TCudaTensor<AFloat> At(size_t i) const {
      Shape_t sliced_shape = (GetLayout() == MemoryLayout::RowMajor)
               ? Shape_t(fShape.begin() + 1, fShape.end()) :
                 Shape_t(fShape.begin(), fShape.end() - 1);


      size_t buffsize = (GetLayout() == MemoryLayout::RowMajor) ?
         fStrides.front() :  fStrides.back();

      size_t offset = i * buffsize;

      return TCudaTensor<AFloat>((const_cast<TCudaDeviceBuffer<AFloat>&>(fElementBuffer)).GetSubBuffer(offset, buffsize), sliced_shape, GetLayout());
   }


   // element access ( for debugging)
   TCudaDeviceReference<AFloat> operator()(size_t i, size_t j) const
   {
      // like this works also for multi-dim tensors
      // and consider the tensor as a multidim one
      size_t nrows = GetNrows();
      size_t ncols = GetNcols();

      size_t offset = (GetLayout() == MemoryLayout::RowMajor) ?
         i * ncols + j  : j * nrows + i;

      AFloat * elementPointer = fElementBuffer + offset;
      return TCudaDeviceReference<AFloat>(elementPointer);
   }
   // element access ( for debugging)
   TCudaDeviceReference<AFloat> operator()(size_t i, size_t j, size_t k) const
   {
      // k is B, i is C, j is HW :
      assert( fNDim >= 3); // || ( k==0 && fNDim == 2 ) );
      //note  for larger dimension k is all other dims collapsed !!!

      size_t offset = (GetLayout() == MemoryLayout::RowMajor) ?
            i * fStrides[0] + j * fStrides[1] + k :
            i * fStrides[2] + k * fStrides[1] + j;

      AFloat * elementPointer = fElementBuffer + offset;

      return TCudaDeviceReference<AFloat>(elementPointer);
   }

    TCudaDeviceReference<AFloat> operator()(size_t i, size_t j, size_t k, size_t l) const
   {
      // for rowsise
      //assert(GetLayout() == MemoryLayout::RowMajor);
      assert( fNDim == 4); // || ( k==0 && fNDim == 2 ) );

      size_t offset = (GetLayout() == MemoryLayout::RowMajor) ?
            i * fStrides[0] + j * fStrides[1] + k * fStrides[2] + l:
            l * fStrides[3] + k * fStrides[2] + j * fStrides[1] + i;

      AFloat * elementPointer = fElementBuffer + offset;

      return TCudaDeviceReference<AFloat>(elementPointer);
   }

private:

   /** Initializes all shared devices resource and makes sure that a sufficient
    *  number of curand states are allocated on the device and initialized as
    *  well as that the one-vector for the summation over columns has the right
    *  size. */
   void InitializeCuda();
   void InitializeCurandStates();

};




} // namespace DNN
} // namespace TMVA

#endif
