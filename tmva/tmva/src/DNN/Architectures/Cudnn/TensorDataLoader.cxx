// @(#)root/tmva/tmva/dnn:$Id$
// Author: Lorenzo Moneta,


////////////////////////////////////////////////////////////////////////
// Implementation of TensorDataLoader functions for CUDA with CuDNN architecture.  //
////////////////////////////////////////////////////////////////////////

#include "TMVA/DataSetInfo.h"

#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Architectures/Cuda/CudaBuffers.h"

#include "TMVA/DNN/Architectures/TCudnn.h"



#include "cuda_runtime.h"
#include <algorithm>

namespace TMVA {
namespace DNN {

//______________________________________________________________________________
//
// cuDNN
//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCudnn<float> >::CopyTensorInput(TCudaHostBuffer<float> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   const std::vector<TMatrixT<Double_t> > &inputTensor = std::get<0>(fData);

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
void TTensorDataLoader<TensorInput, TCudnn<float> >::CopyTensorOutput(TCudaHostBuffer<float> &buffer,
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
void TTensorDataLoader<TensorInput, TCudnn<float> >::CopyTensorWeights(TCudaHostBuffer<float> &buffer,
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
void TTensorDataLoader<TMVAInput_t, TCudnn<float> >::CopyTensorInput(TCudaHostBuffer<float> &buffer,
                                                                     IndexIterator_t sampleIterator)
{
   // Image has channel depth 1 -> they are ordered as row-vectors in a matrix (batchHeight = batchSize)
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
   // A batch is made up by a single image with its channels
   } else if (fBatchDepth == fBatchSize) {
      for (size_t i = 0; i < fBatchSize; i++) {
         size_t sampleIndex = *sampleIterator;
         Event * event = std::get<0>(fData)[sampleIndex];
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               // Cudnn order is NCHW
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + j * fBatchWidth + k;
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
void TTensorDataLoader<TMVAInput_t, TCudnn<float> >::CopyTensorOutput(TCudaHostBuffer<float> &buffer,
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
void TTensorDataLoader<TMVAInput_t, TCudnn<float> >::CopyTensorWeights(TCudaHostBuffer<float> &buffer,
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
void TTensorDataLoader<TensorInput, TCudnn<double> >::CopyTensorInput(TCudaHostBuffer<double> &buffer,
                                                                      IndexIterator_t sampleIterator)
{
   const std::vector<TMatrixT<Double_t> > &inputTensor = std::get<0>(fData);

   if (fBatchDepth == 1) {
      for (size_t i = 0; i < fBatchHeight; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchWidth; j++) {
            size_t bufferIndex = j * fBatchHeight + i;
            buffer[bufferIndex] = static_cast<double>(inputTensor[0](sampleIndex, j));
         }
         sampleIterator++;
      }
   } else {
      for (size_t i = 0; i < fBatchDepth; i++) {
         size_t sampleIndex = *sampleIterator;
         for (size_t j = 0; j < fBatchHeight; j++) {
            for (size_t k = 0; k < fBatchWidth; k++) {
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + k * fBatchHeight + j;
               buffer[bufferIndex] = static_cast<double>(inputTensor[sampleIndex](j, k));
            }
         }
         sampleIterator++;
      }
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TensorInput, TCudnn<double> >::CopyTensorOutput(TCudaHostBuffer<double> &buffer,
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
void TTensorDataLoader<TensorInput, TCudnn<double> >::CopyTensorWeights(TCudaHostBuffer<double> &buffer,
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
void TTensorDataLoader<TMVAInput_t, TCudnn<double> >::CopyTensorInput(TCudaHostBuffer<double> &buffer,
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
               size_t bufferIndex = i * fBatchHeight * fBatchWidth + j * fBatchWidth + k;
               buffer[bufferIndex] = event->GetValue(j * fBatchWidth + k);
            }
         }
         sampleIterator++;
      }
   }
   else {
      Error("TTensorDataLoader","Inconsistency between batch depth and batch size");
      R__ASSERT(0);
   }
}

//______________________________________________________________________________
template <>
void TTensorDataLoader<TMVAInput_t, TCudnn<double> >::CopyTensorOutput(TCudaHostBuffer<double> &buffer,
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
void TTensorDataLoader<TMVAInput_t, TCudnn<double> >::CopyTensorWeights(TCudaHostBuffer<double> &buffer,
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
TTensorBatch<TCudnn<float> > TTensorDataLoader<TensorInput, TCudnn<float> >::GetTensorBatch()
{
   // Get buffer tuple on device that contains the data
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<size_t> outputShape  {fBatchSize, 1, fNOutputFeatures, 1};
   std::vector<size_t> wheightShape {fBatchSize, 1, 1, 1};
   std::vector<TCudaTensor<float> > inputTensor(1, TCudaTensor<float>(std::get<0>(DeviceBuffers),
                                                this->GetTensorDim(),  fInputShape));
   TCudaTensor<float> outputMatrix(std::get<1>(DeviceBuffers), this->GetTensorDim(), outputShape);
   TCudaTensor<float> weightMatrix(std::get<2>(DeviceBuffers), this->GetTensorDim(), wheightShape);

   fBatchIndex++;
   return TTensorBatch<TCudnn<float> >(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCudnn<double> > TTensorDataLoader<TensorInput, TCudnn<double> >::GetTensorBatch()
{
   // Get buffer tuple on device that contains the data
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<size_t> outputShape  {fBatchSize, 1, fNOutputFeatures, 1};
   std::vector<size_t> wheightShape {fBatchSize, 1, 1, 1};
   std::vector<TCudaTensor<double> > inputTensor(1, TCudaTensor<double>(std::get<0>(DeviceBuffers),
                                                 this->GetTensorDim(),  fInputShape));
   TCudaTensor<double> outputMatrix(std::get<1>(DeviceBuffers), this->GetTensorDim(), outputShape);
   TCudaTensor<double> weightMatrix(std::get<2>(DeviceBuffers), this->GetTensorDim(), wheightShape);

   fBatchIndex++;
   return TTensorBatch<TCudnn<double> >(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCudnn<float> > TTensorDataLoader<TMVAInput_t, TCudnn<float> >::GetTensorBatch()
{
   // Get buffer tuple on device that contains the data
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<size_t> outputShape  {fBatchSize, 1, fNOutputFeatures, 1};
   std::vector<size_t> wheightShape {fBatchSize, 1, 1, 1};
   std::vector<TCudaTensor<float> > inputTensor(1, TCudaTensor<float>(std::get<0>(DeviceBuffers),
                                                this->GetTensorDim(),  fInputShape));
   TCudaTensor<float> outputMatrix(std::get<1>(DeviceBuffers), this->GetTensorDim(), outputShape);
   TCudaTensor<float> weightMatrix(std::get<2>(DeviceBuffers), this->GetTensorDim(), wheightShape);

   fBatchIndex++;
   return TTensorBatch<TCudnn<float> >(inputTensor, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template <>
TTensorBatch<TCudnn<double> > TTensorDataLoader<TMVAInput_t, TCudnn<double> >::GetTensorBatch()
{
   // Get buffer tuple on device that contains the data
   DeviceBufferTuple DeviceBuffers = CopyTensorBatches();

   std::vector<size_t> outputShape  {fBatchSize, 1, fNOutputFeatures, 1};
   std::vector<size_t> wheightShape {fBatchSize, 1, 1, 1};
   std::vector<TCudaTensor<double> > inputTensor(1, TCudaTensor<double>(std::get<0>(DeviceBuffers),
                                                 this->GetTensorDim(),  fInputShape));
   TCudaTensor<double> outputMatrix(std::get<1>(DeviceBuffers), fNOutputFeatures + 2, outputShape);
   TCudaTensor<double> weightMatrix(std::get<2>(DeviceBuffers), 3, wheightShape);

   fBatchIndex++;
   return TTensorBatch<TCudnn<double> >(inputTensor, outputMatrix, weightMatrix);
}
#endif


//______________________________________________________________________________
// Explicit Instantiations.

template class TTensorDataLoader<TensorInput, TCudnn<float> >;
template class TTensorDataLoader<TMVAInput_t, TCudnn<float> >;
template class TTensorDataLoader<TensorInput, TCudnn<double> >;
template class TTensorDataLoader<TMVAInput_t, TCudnn<double> >;

} // TMVA
} // DNN
