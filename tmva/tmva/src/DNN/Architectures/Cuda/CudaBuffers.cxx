// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 07/08/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////
// Implementation of device and host buffers for CUDA architectures.  //
////////////////////////////////////////////////////////////////////////

#include "TMVA/DataSetInfo.h"
#include "TMVA/DNN/DataLoader.h"

#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Architectures/Cuda.h"
#ifdef R__HAS_CUDNN
#include "TMVA/DNN/Architectures/TCudnn.h"
#endif
#include "TMVA/DNN/Architectures/Cuda/CudaBuffers.h"

#include "cuda_runtime.h"
#include <algorithm>

namespace TMVA {
namespace DNN {

//
// TCudaHostBuffer
//______________________________________________________________________________
template <typename AFloat>
void TCudaHostBuffer<AFloat>::TDestructor::operator()(AFloat **devicePointer)
{
   cudaFreeHost(*devicePointer);
   delete[] devicePointer;
}

//______________________________________________________________________________
template <typename AFloat>
TCudaHostBuffer<AFloat>::TCudaHostBuffer(size_t size) : fOffset(0), fSize(size), fComputeStream(0), fDestructor()
{
   AFloat **pointer = new AFloat *[1];
   cudaMallocHost(pointer, size * sizeof(AFloat));
   fHostPointer = std::shared_ptr<AFloat *>(pointer, fDestructor);
}

//______________________________________________________________________________
template <typename AFloat>
TCudaHostBuffer<AFloat>::operator AFloat *() const
{
   return (fHostPointer) ? *fHostPointer + fOffset : nullptr;
}

//______________________________________________________________________________
template <typename AFloat>
TCudaHostBuffer<AFloat> TCudaHostBuffer<AFloat>::GetSubBuffer(size_t offset, size_t size)
{
   TCudaHostBuffer buffer = *this;
   buffer.fOffset = offset;
   buffer.fSize = size;
   return buffer;
}

//______________________________________________________________________________
template <typename AFloat>
void TCudaHostBuffer<AFloat>::SetConstVal(const AFloat constVal)
{
   std::fill(*fHostPointer, *fHostPointer+fSize, constVal);
}

//
// TCudaDevicePointer
//______________________________________________________________________________
template <typename AFloat>
void TCudaDeviceBuffer<AFloat>::TDestructor::operator()(AFloat **devicePointer)
{
   cudaFree(*devicePointer);
   delete[] devicePointer;
}

//______________________________________________________________________________
template <typename AFloat>
TCudaDeviceBuffer<AFloat>::TCudaDeviceBuffer(size_t size) : fOffset(0), fSize(size), fDestructor()
{
   AFloat **pointer = new AFloat *[1];
   cudaMalloc(pointer, size * sizeof(AFloat));
   fDevicePointer = std::shared_ptr<AFloat *>(pointer, fDestructor);
   cudaStreamCreate(&fComputeStream);
}

//______________________________________________________________________________
template <typename AFloat>
TCudaDeviceBuffer<AFloat>::TCudaDeviceBuffer(size_t size, cudaStream_t stream)
   : fOffset(0), fSize(size), fComputeStream(stream), fDestructor()
{
   AFloat **pointer = new AFloat *[1];
   cudaMalloc(pointer, size * sizeof(AFloat));
   fDevicePointer = std::shared_ptr<AFloat *>(pointer, fDestructor);
}

//______________________________________________________________________________
template <typename AFloat>
TCudaDeviceBuffer<AFloat>::TCudaDeviceBuffer(AFloat *devicePointer, size_t size, cudaStream_t stream)
   : fOffset(0), fSize(size), fComputeStream(stream), fDestructor()
{
   AFloat **pointer = new AFloat *[1];
   *pointer = devicePointer;
   fDevicePointer = std::shared_ptr<AFloat *>(pointer, fDestructor);
}

//______________________________________________________________________________
template <typename AFloat>
TCudaDeviceBuffer<AFloat> TCudaDeviceBuffer<AFloat>::GetSubBuffer(size_t offset, size_t size)
{
   TCudaDeviceBuffer buffer = *this;
   buffer.fOffset = offset;
   buffer.fSize = size;
   return buffer;
}

//______________________________________________________________________________
template <typename AFloat>
TCudaDeviceBuffer<AFloat>::operator AFloat *() const
{
   return (fDevicePointer) ? *fDevicePointer + fOffset : nullptr;
}

//______________________________________________________________________________
template <typename AFloat>
void TCudaDeviceBuffer<AFloat>::CopyFrom(const TCudaHostBuffer<AFloat> &buffer) const
{
   cudaStreamSynchronize(fComputeStream);
   cudaMemcpyAsync(*this, buffer, fSize * sizeof(AFloat), cudaMemcpyHostToDevice, fComputeStream);
}

//______________________________________________________________________________
template <typename AFloat>
void TCudaDeviceBuffer<AFloat>::CopyTo(const TCudaHostBuffer<AFloat> &buffer) const
{
   cudaMemcpyAsync(buffer, *this, fSize * sizeof(AFloat), cudaMemcpyDeviceToHost, fComputeStream);
   buffer.fComputeStream = fComputeStream;
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<float>>::CopyInput(TCudaHostBuffer<float> &buffer, IndexIterator_t sampleIterator,
                                                         size_t batchSize)
{
   const TMatrixT<Double_t> &inputMatrix = std::get<0>(fData);
   size_t n = inputMatrix.GetNcols();

   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = static_cast<float>(inputMatrix(sampleIndex, j));
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<float>>::CopyOutput(TCudaHostBuffer<float> &buffer,
                                                          IndexIterator_t sampleIterator, size_t batchSize)
{
   const TMatrixT<Double_t> &outputMatrix = std::get<1>(fData);
   size_t n = outputMatrix.GetNcols();

   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = static_cast<float>(outputMatrix(sampleIndex, j));
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<float>>::CopyWeights(TCudaHostBuffer<float> &buffer,
                                                           IndexIterator_t sampleIterator, size_t batchSize)
{
   const TMatrixT<Double_t> &weightMatrix = std::get<2>(fData);
   for (size_t i = 0; i < batchSize; i++) {
      buffer[i] = static_cast<float>(weightMatrix(*sampleIterator, 0));
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<float>>::CopyInput(TCudaHostBuffer<float> &buffer, IndexIterator_t sampleIterator,
                                                       size_t batchSize)
{
   Event *event = std::get<0>(fData)[0];
   size_t n  = event->GetNVariables();
   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = * sampleIterator++;
      event = std::get<0>(fData)[sampleIndex];
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = static_cast<float>(event->GetValue(j));
      }
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<float>>::CopyOutput(TCudaHostBuffer<float> &buffer, IndexIterator_t sampleIterator,
                                                        size_t batchSize)
{
  const DataSetInfo &info = std::get<1>(fData);
  size_t n = buffer.GetSize() / batchSize;

  // Copy target(s).

  for (size_t i = 0; i < batchSize; i++) {
    size_t sampleIndex = *sampleIterator++;
    Event *event = std::get<0>(fData)[sampleIndex];
    for (size_t j = 0; j < n; j++) {
      // Copy output matrices.
      size_t bufferIndex = j * batchSize + i;
      // Classification
      if (event->GetNTargets() == 0) {
        if (n == 1) {
          // Binary.
          buffer[bufferIndex] = (info.IsSignal(event)) ? 1.0 : 0.0;
        } else {
          // Multiclass.
          buffer[bufferIndex] = 0.0;
          if (j == event->GetClass()) {
            buffer[bufferIndex] = 1.0;
          }
        }
      } else {
        buffer[bufferIndex] = static_cast<float>(event->GetTarget(j));
      }
    }
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<float>>::CopyWeights(TCudaHostBuffer<float> &buffer, IndexIterator_t sampleIterator,
                                                         size_t batchSize)
{
   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      buffer[i] = static_cast<float>(event->GetWeight());
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<double>>::CopyInput(TCudaHostBuffer<double> &buffer,
                                                          IndexIterator_t sampleIterator, size_t batchSize)
{
   const TMatrixT<Double_t> &inputMatrix = std::get<0>(fData);
   size_t n = inputMatrix.GetNcols();

   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = inputMatrix(sampleIndex, j);
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<double>>::CopyOutput(TCudaHostBuffer<double> &buffer,
                                                           IndexIterator_t sampleIterator, size_t batchSize)
{
   const TMatrixT<Double_t> &outputMatrix = std::get<1>(fData);
   size_t n = outputMatrix.GetNcols();

   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = outputMatrix(sampleIndex, j);
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<MatrixInput_t, TCuda<double>>::CopyWeights(TCudaHostBuffer<double> &buffer,
                                                            IndexIterator_t sampleIterator, size_t batchSize)
{
   const TMatrixT<Double_t> &weightMatrix = std::get<2>(fData);
   for (size_t i = 0; i < batchSize; i++) {
      buffer[i] = static_cast<double>(weightMatrix(*sampleIterator, 0));
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<double>>::CopyInput(TCudaHostBuffer<double> &buffer, IndexIterator_t sampleIterator,
                                                        size_t batchSize)
{
   Event *event = std::get<0>(fData)[0];
   size_t n  = event->GetNVariables();
   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = * sampleIterator++;
      event = std::get<0>(fData)[sampleIndex];
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * batchSize + i;
         buffer[bufferIndex] = event->GetValue(j);
      }
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<double>>::CopyOutput(TCudaHostBuffer<double> &buffer,
                                                         IndexIterator_t sampleIterator, size_t batchSize)
{
  const DataSetInfo &info = std::get<1>(fData);
  size_t n = buffer.GetSize() / batchSize;

  // Copy target(s).

  for (size_t i = 0; i < batchSize; i++) {
    size_t sampleIndex = *sampleIterator++;
    Event *event = std::get<0>(fData)[sampleIndex];
    for (size_t j = 0; j < n; j++) {
      // Copy output matrices.
      size_t bufferIndex = j * batchSize + i;
      // Classification
      if (event->GetNTargets() == 0) {
        // Binary.
        if (n == 1) {
          buffer[bufferIndex] = (info.IsSignal(event)) ? 1.0 : 0.0;
        } else {
          // Multiclass.
          buffer[bufferIndex] = 0.0;
          if (j == event->GetClass()) {
            buffer[bufferIndex] = 1.0;
          }
        }
      } else {
        buffer[bufferIndex] = event->GetTarget(j);
      }
    }
   }
}

//______________________________________________________________________________
template <>
void TDataLoader<TMVAInput_t, TCuda<double>>::CopyWeights(TCudaHostBuffer<double> &buffer,
                                                          IndexIterator_t sampleIterator, size_t batchSize)
{
   for (size_t i = 0; i < batchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      buffer[i] = static_cast<double>(event->GetWeight());
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<float>>::CopyTensorInput(TCudaHostBuffer<float> &buffer,
                                                                   IndexIterator_t sampleIterator)
{
   const std::vector<TMatrixT<Double_t>> &inputTensor = std::get<0>(fData);

   if (fBatchDepth == 1) {
      for (size_t i = 0; i < fBatchHeight; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchWidth; j++) {
            size_t bufferIndex = j * fBatchHeight + i;
            buffer[bufferIndex] = static_cast<float>(inputTensor[0](sampleIndex, j));
         }
         sampleIterator++;
      }
   } else {
      for (size_t i = 0; i < fBatchDepth; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + k * fBatchHeight + j;
               buffer[bufferIndex] = static_cast<float>(inputTensor[sampleIndex](j, k));
            }
         }
         sampleIterator++;
      }
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<float>>::CopyTensorOutput(TCudaHostBuffer<float> &buffer,
                                                                    IndexIterator_t sampleIterator)
{
   const TMatrixT<Double_t> &outputMatrix = std::get<1>(fData);
   size_t n = outputMatrix.GetNcols();

   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * fBatchSize + i;
         buffer[bufferIndex] = static_cast<float>(outputMatrix(sampleIndex, j));
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<float>>::CopyTensorWeights(TCudaHostBuffer<float> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   const TMatrixT<Double_t> &weightMatrix = std::get<2>(fData);
   for (size_t i = 0; i < fBatchSize; i++) {
      buffer[i] = static_cast<float>(weightMatrix(*sampleIterator, 0));
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<float>>::CopyTensorInput(TCudaHostBuffer<float> &buffer,
                                                                   IndexIterator_t sampleIterator)
{
   // one event, one  example in the batch

   if (fBatchDepth == 1 && fBatchHeight == fBatchSize) {
      for (size_t i = 0; i < fBatchHeight; i++) {
         size_t sampleIndex = *sampleIterator;
         Event * event = std::get<0>(fData)[sampleIndex];
         for (size_t j = 0; j < fBatchWidth; j++) {
            size_t bufferIndex = j * fBatchHeight + i;
            buffer[bufferIndex] = event->GetValue(j);
         }
         sampleIterator++;
      }
   } else if (fBatchDepth == fBatchSize) {
      // batchDepth is batch size
      for (size_t i = 0; i < fBatchDepth; i++) {
         size_t sampleIndex = *sampleIterator;
         Event * event = std::get<0>(fData)[sampleIndex];
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               // because of the column-major ordering
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + k * fBatchHeight + j;
               buffer[bufferIndex] = event->GetValue(j * fBatchWidth + k);
            }
         }
         sampleIterator++;
      }
   }
   else {
      std::cout  << fBatchDepth << fBatchSize << fBatchHeight << std::endl;
      Error("TTensorDataLoader","Inconsistency between batch depth and batch size");
      R__ASSERT(0);
   }
}
//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<float>>::CopyTensorOutput(TCudaHostBuffer<float> &buffer,
                                                                    IndexIterator_t sampleIterator)
{
   const DataSetInfo &info = std::get<1>(fData);
   size_t n = buffer.GetSize() / fBatchSize;

   // Copy target(s).

   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      for (size_t j = 0; j < n; j++) {
         // Copy output matrices.
         size_t bufferIndex = j * fBatchSize + i;
         // Classification
         if (event->GetNTargets() == 0) {
            if (n == 1) {
               // Binary.
               buffer[bufferIndex] = (info.IsSignal(event)) ? 1.0 : 0.0;
            } else {
               // Multiclass.
               buffer[bufferIndex] = 0.0;
               if (j == event->GetClass()) {
                  buffer[bufferIndex] = 1.0;
               }
            }
         } else {
            buffer[bufferIndex] = static_cast<Float_t>(event->GetTarget(j));
         }
      }
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<float>>::CopyTensorWeights(TCudaHostBuffer<float> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      buffer[i] = event->GetWeight();
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<Double_t>>::CopyTensorInput(TCudaHostBuffer<double> &buffer,
                                                                    IndexIterator_t sampleIterator)
{
   const std::vector<TMatrixT<Double_t>> &inputTensor = std::get<0>(fData);

   if (fBatchDepth == 1) {
      for (size_t i = 0; i < fBatchHeight; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchWidth; j++) {
            size_t bufferIndex = j * fBatchHeight + i;
            buffer[bufferIndex] = static_cast<float>(inputTensor[0](sampleIndex, j));
         }
         sampleIterator++;
      }
   } else {
      for (size_t i = 0; i < fBatchDepth; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + k * fBatchHeight + j;
               buffer[bufferIndex] = static_cast<float>(inputTensor[sampleIndex](j, k));
            }
         }
         sampleIterator++;
      }
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<Double_t>>::CopyTensorOutput(TCudaHostBuffer<double> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   const TMatrixT<Double_t> &outputMatrix = std::get<1>(fData);
   size_t n = outputMatrix.GetNcols();

   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator;
      for (size_t j = 0; j < n; j++) {
         size_t bufferIndex = j * fBatchSize + i;
         buffer[bufferIndex] = outputMatrix(sampleIndex, j);
      }
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCuda<Double_t>>::CopyTensorWeights(TCudaHostBuffer<double> &buffer,
                                                                      IndexIterator_t sampleIterator)
{
   const TMatrixT<Double_t> &weightMatrix = std::get<2>(fData);

   for (size_t i = 0; i < fBatchSize; i++) {
      buffer[i] = weightMatrix(*sampleIterator, 0);
      sampleIterator++;
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<Double_t>>::CopyTensorInput(TCudaHostBuffer<double> &buffer,
                                                                    IndexIterator_t sampleIterator)
{
   // one event, one  example in the batch

   if (fBatchDepth == 1 && fBatchHeight == fBatchSize) {
      for (size_t i = 0; i < fBatchHeight; i++) {
         size_t sampleIndex = *sampleIterator;
         Event * event = std::get<0>(fData)[sampleIndex];
         for (size_t j = 0; j < fBatchWidth; j++) {
            size_t bufferIndex = j * fBatchHeight + i;
            buffer[bufferIndex] = event->GetValue(j);
         }
         sampleIterator++;
      }
   } else if (fBatchDepth == fBatchSize) {
      // batchDepth is batch size
      for (size_t i = 0; i < fBatchDepth; i++) {
         size_t sampleIndex = *sampleIterator;
         Event * event = std::get<0>(fData)[sampleIndex];
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               // because of the column-major ordering
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + k * fBatchHeight + j;
               buffer[bufferIndex] = event->GetValue(j * fBatchWidth + k);
            }
         }
         sampleIterator++;
      }
   }
   else {
      std::cout  << fBatchDepth << fBatchSize << fBatchHeight << std::endl;
      Error("TTensorDataLoader","Inconsistency between batch depth and batch size");
      R__ASSERT(0);
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<Double_t>>::CopyTensorOutput(TCudaHostBuffer<double> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   const DataSetInfo &info = std::get<1>(fData);
   size_t n = buffer.GetSize() / fBatchSize;

   // Copy target(s).

   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      for (size_t j = 0; j < n; j++) {
         // Copy output matrices.
         size_t bufferIndex = j * fBatchSize + i;
         // Classification
         if (event->GetNTargets() == 0) {
            if (n == 1) {
               // Binary.
               buffer[bufferIndex] = (info.IsSignal(event)) ? 1.0 : 0.0;
            } else {
               // Multiclass.
               buffer[bufferIndex] = 0.0;
               if (j == event->GetClass()) {
                  buffer[bufferIndex] = 1.0;
               }
            }
         } else {
            buffer[bufferIndex] = static_cast<Double_t>(event->GetTarget(j));
         }
      }
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCuda<Double_t>>::CopyTensorWeights(TCudaHostBuffer<double> &buffer,
                                                                      IndexIterator_t sampleIterator)
{
   for (size_t i = 0; i < fBatchSize; i++) {
      size_t sampleIndex = *sampleIterator++;
      Event *event = std::get<0>(fData)[sampleIndex];
      buffer[i] = event->GetWeight();
   }
}

#if 0
//______________________________________________________________________________
template <>
TTensorBatch<TCuda<float> > TTensorDataLoader<TensorInput, TCuda<float> >::GetTensorBatch()
{
   // After copying the data to the device, wrap the device buffer in the respective
   // architectures matrix type
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<Matrix_t> inputTensor(std::get<0>(DeviceBuffers), fBatchSize, )
   size_t jump = fBatchHeight * fBatchWidth;
   for (size_t i = 0; i < fBatchSize; i++) {
      DeviceBuffer_t subInputDeviceBuffer = std::get<0>(DeviceBuffers).GetSubBuffer(i * jump, jump);
      inputTensor.emplace_back(subInputDeviceBuffer, fBatchHeight, fBatchWidth);
   }
   Matrix_t outputMatrix(std::get<1>(DeviceBuffers), fBatchSize, fNOutputFeatures);
   Matrix_t weightMatrix(std::get<2>(DeviceBuffers), fBatchSize, fNOutputFeatures);

   fBatchIndex++;
   return TTensorBatch<TCuda<float>>(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCuda<double> > TTensorDataLoader<TensorInput, TCuda<double> >::GetTensorBatch()
{
   // After copying the data to the device, wrap the device buffer in the respective
   // architectures matrix type
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<Matrix_t> inputTensor;
   size_t jump = fBatchHeight * fBatchWidth;
   for (size_t i = 0; i < fBatchSize; i++) {
      DeviceBuffer_t subInputDeviceBuffer = std::get<0>(DeviceBuffers).GetSubBuffer(i * jump, jump);
      inputTensor.emplace_back(subInputDeviceBuffer, fBatchHeight, fBatchWidth);
   }
   Matrix_t outputMatrix(std::get<1>(DeviceBuffers), fBatchSize, fNOutputFeatures);
   Matrix_t weightMatrix(std::get<2>(DeviceBuffers), fBatchSize, fNOutputFeatures);

   fBatchIndex++;
   return TTensorBatch<TCuda<double>>(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCuda<float> > TTensorDataLoader<TMVAInput_t, TCuda<float> >::GetTensorBatch()
{
   // After copying the data to the device, wrap the device buffer in the respective
   // architectures matrix type
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<Matrix_t> inputTensor;
   size_t jump = fBatchHeight * fBatchWidth;
   for (size_t i = 0; i < fBatchSize; i++) {
      DeviceBuffer_t subInputDeviceBuffer = std::get<0>(DeviceBuffers).GetSubBuffer(i * jump, jump);
      inputTensor.emplace_back(subInputDeviceBuffer, fBatchHeight, fBatchWidth);
   }
   Matrix_t outputMatrix(std::get<1>(DeviceBuffers), fBatchSize, fNOutputFeatures);
   Matrix_t weightMatrix(std::get<2>(DeviceBuffers), fBatchSize, fNOutputFeatures);

   fBatchIndex++;
   return TTensorBatch<TCuda<float>>(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCuda<double> > TTensorDataLoader<TMVAInput_t, TCuda<double> >::GetTensorBatch()
{
   // After copying the data to the device, wrap the device buffer in the respective
   // architectures matrix type
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<Matrix_t> inputTensor;
   size_t jump = fBatchHeight * fBatchWidth;
   for (size_t i = 0; i < fBatchSize; i++) {
      DeviceBuffer_t subInputDeviceBuffer = std::get<0>(DeviceBuffers).GetSubBuffer(i * jump, jump);
      inputTensor.emplace_back(subInputDeviceBuffer, fBatchHeight, fBatchWidth);
   }
   Matrix_t outputMatrix(std::get<1>(DeviceBuffers), fBatchSize, fNOutputFeatures);
   Matrix_t weightMatrix(std::get<2>(DeviceBuffers), fBatchSize, fNOutputFeatures);

   fBatchIndex++;
   return TTensorBatch<TCuda<double>>(inputTensor, outputMatrix, weightMatrix);
}
#endif


// see file Cudnn/TensorDataLoader.cxx for Cudnn definitions

//______________________________________________________________________________
// Explicit Instantiations.

template class TCudaDeviceBuffer<float>;
template class TCudaDeviceBuffer<double>;

template class TCudaHostBuffer<float>;
template class TCudaHostBuffer<double>;

template class TDataLoader<MatrixInput_t, TCuda<float>>;
template class TDataLoader<TMVAInput_t, TCuda<float>>;
template class TDataLoader<MatrixInput_t, TCuda<double>>;
template class TDataLoader<TMVAInput_t, TCuda<double>>;

template class TTensorDataLoader<TensorInput, TCuda<float> >;
template class TTensorDataLoader<TMVAInput_t, TCuda<float> >;
template class TTensorDataLoader<TensorInput, TCuda<double >>;
template class TTensorDataLoader<TMVAInput_t, TCuda<double> >;


} // TMVA
} // DNN
