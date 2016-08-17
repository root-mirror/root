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
// Implementation of the TCudaMatrix class. //
/////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cuda/CudaMatrix.h"
#include "TMVA/DNN/Architectures/Cuda/Device.h"
#include "TMVA/DNN/Architectures/Cuda/Kernels.h"

namespace TMVA {
namespace DNN  {

// Static members.
//____________________________________________________________________________
size_t          TCudaMatrix::fInstances     = 0;
cublasHandle_t  TCudaMatrix::fCublasHandle  = nullptr;
CudaDouble_t  * TCudaMatrix::fDeviceReturn  = nullptr;
CudaDouble_t  * TCudaMatrix::fOnes          = nullptr;
curandState_t * TCudaMatrix::fCurandStates  = nullptr;
size_t          TCudaMatrix::fNCurandStates = 0;
size_t          TCudaMatrix::fNOnes         = 0;

// Constructors.
//____________________________________________________________________________
TCudaMatrix::TCudaMatrix()
    : fNRows(0), fNCols(0), fElementBuffer()
{
   InitializeCuda();
}

//____________________________________________________________________________
TCudaMatrix::TCudaMatrix(size_t m, size_t n)
    : fNRows(m), fNCols(n), fElementBuffer(m * n, 0)
{
   InitializeCuda();
}

//____________________________________________________________________________
TCudaMatrix::TCudaMatrix(const TMatrixT<CudaDouble_t> & Host)
    : fNRows(Host.GetNrows()), fNCols(Host.GetNcols()),
      fElementBuffer(Host.GetNoElements(), 0)
{
   InitializeCuda();

   CudaDouble_t * buffer = new CudaDouble_t[fNRows * fNCols];
   size_t index = 0;
   for (size_t j = 0; j < fNCols; j++) {
      for (size_t i = 0; i < fNRows; i++) {
         buffer[index] = Host(i, j);
         index++;
      }
   }

   cudaMemcpy(fElementBuffer, buffer, fNRows * fNCols * sizeof(CudaDouble_t),
              cudaMemcpyHostToDevice);
}

//____________________________________________________________________________
TCudaMatrix::TCudaMatrix(TCudaDeviceBuffer buffer,
                         size_t m, size_t n)
    : fNRows(m), fNCols(n), fElementBuffer(buffer)
{
   InitializeCuda();
}

//____________________________________________________________________________
inline void TCudaMatrix::InitializeCuda()
{
   if (fInstances == 0) {
       cublasCreate(&fCublasHandle);
       CUDACHECK(cudaMalloc(& fDeviceReturn, sizeof(CudaDouble_t)));
       CUDACHECK(cudaMalloc(& fCurandStates, TDevice::NThreads(*this)));
   }
   if (TDevice::NThreads(*this) > (int) fNCurandStates) {
       fNCurandStates = TDevice::NThreads(*this);
       if (fCurandStates) {
           cudaFree(fCurandStates);
       }
       cudaMalloc(&fCurandStates, TDevice::NThreads(*this) * sizeof(curandState_t));
       InitializeCurandStates();
   }
   if (fNRows >  fNOnes) {
      fNOnes = fNRows;
      if (fOnes) {
         cudaFree(fOnes);
      }
      cudaMalloc(&fOnes, fNRows * sizeof(CudaDouble_t));
      CudaDouble_t * buffer = new CudaDouble_t[fNRows];
      for (size_t i = 0; i < fNRows; i++) {
         buffer[i] = 1.0;
      }
      cudaMemcpy(fOnes, buffer, fNRows * sizeof(CudaDouble_t),
                 cudaMemcpyHostToDevice);
   }
   fInstances++;
}

//____________________________________________________________________________
void TCudaMatrix::InitializeCurandStates()
{
   dim3 blockDims = TDevice::BlockDims();
   dim3 gridDims  = TDevice::GridDims(*this);
   ::TMVA::DNN::Cuda::InitializeCurandStates<<<gridDims, blockDims>>>(time(nullptr),
                                                                      fCurandStates);

}


// Conversion to TMatrixT.
//____________________________________________________________________________
TCudaMatrix::operator TMatrixT<CudaDouble_t>() const
{
   TMatrixT<CudaDouble_t> hostMatrix(GetNrows(), GetNcols());

   CudaDouble_t * buffer = new CudaDouble_t[fNRows * fNCols];
   cudaMemcpy(buffer, fElementBuffer, fNRows * fNCols * sizeof(CudaDouble_t),
              cudaMemcpyDeviceToHost);

   size_t index = 0;
   for (size_t j = 0; j < fNCols; j++) {
      for (size_t i = 0; i < fNRows; i++) {
         hostMatrix(i, j) = buffer[index];
         index++;
      }
   }

   delete[] buffer;
   return hostMatrix;
}

} // namespace DNN
} // namespace TMVA
