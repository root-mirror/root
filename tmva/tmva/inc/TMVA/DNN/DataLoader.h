// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 08/08/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////
// Generic data loader for neural network input data. Provides a   //
// high level abstraction for the transfer of training data to the //
// device.                                                         //
/////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_DATALOADER
#define TMVA_DNN_DATALOADER

#include "TMatrix.h"
#include "TMVA/Event.h"

#include <algorithm>
#include <random>
#include <vector>
#include <utility>

namespace TMVA {

class DataSetInfo;

namespace DNN  {

//
// Input Data Types
//______________________________________________________________________________
using MatrixInput_t = std::tuple<const TMatrixT<Double_t> &, const TMatrixT<Double_t> &, const TMatrixT<Double_t> &>;
using TMVAInput_t =
    std::tuple<const std::vector<Event *> &, const DataSetInfo &>;

using IndexIterator_t = typename std::vector<size_t>::iterator;

/** TBatch
 *
 * Class representing training batches consisting of a matrix of input data
 * and a matrix of output data. The input and output data can be accessed using
 * the GetInput() and GetOutput() member functions.
 *
 * \tparam AArchitecture The underlying architecture.
 */
//______________________________________________________________________________
template <typename AArchitecture>
class TBatch
{
private:

   using Matrix_t       = typename AArchitecture::Matrix_t;

   Matrix_t fInputMatrix;
   Matrix_t fOutputMatrix;
   Matrix_t fWeightMatrix;

public:
   TBatch(Matrix_t &, Matrix_t &, Matrix_t &);
   TBatch(const TBatch  &) = default;
   TBatch(      TBatch &&) = default;
   TBatch & operator=(const TBatch  &) = default;
   TBatch & operator=(      TBatch &&) = default;

   /** Return the matrix representing the input data. */
   Matrix_t &GetInput() { return fInputMatrix; }
   /** Return the matrix representing the output data. */
   Matrix_t &GetOutput() { return fOutputMatrix; }
   /** Return the matrix holding the event weights. */
   Matrix_t &GetWeights() { return fWeightMatrix; }
};

template<typename Data_t, typename AArchitecture> class TDataLoader;

/** TBatchIterator
 *
 * Simple iterator class for the iterations over the training batches in
 * a given data set represented by a TDataLoader object.
 *
 * \tparam AData         The input data type.
 * \tparam AArchitecture The underlying architecture type.
 */
template<typename Data_t, typename AArchitecture>
class TBatchIterator
{
private:

   TDataLoader<Data_t, AArchitecture> & fDataLoader;
   size_t fBatchIndex;

public:

TBatchIterator(TDataLoader<Data_t, AArchitecture> & dataLoader, size_t index = 0)
: fDataLoader(dataLoader), fBatchIndex(index)
{
   // Nothing to do here.
}

   TBatch<AArchitecture> operator*() {return fDataLoader.GetBatch();}
   TBatchIterator operator++() {fBatchIndex++; return *this;}
   bool operator!=(const TBatchIterator & other) {
      return fBatchIndex != other.fBatchIndex;
   }
};

/** TDataLoader
 *
 * Service class managing the streaming of the training data from the input data
 * type to the accelerator device or the CPU. A TDataLoader object manages a number
 * of host and device buffer pairs that are used in a round-robin manner for the
 * transfer of batches to the device.
 *
 * Each TDataLoader object has an associated batch size and a number of total
 * samples in the dataset. One epoch is the number of buffers required to transfer
 * the complete training set. Using the begin() and end() member functions allows
 * the user to iterate over the batches in one epoch.
 *
 * \tparam AData The input data type.
 * \tparam AArchitecture The achitecture class of the underlying architecture.
 */
template<typename Data_t, typename AArchitecture>
class TDataLoader
{
private:

   using HostBuffer_t    = typename AArchitecture::HostBuffer_t;
   using DeviceBuffer_t  = typename AArchitecture::DeviceBuffer_t;
   using Matrix_t        = typename AArchitecture::Matrix_t;
   using BatchIterator_t = TBatchIterator<Data_t, AArchitecture>;

   const Data_t &fData;

   size_t fNSamples;
   size_t fBatchSize;
   size_t fNInputFeatures;
   size_t fNOutputFeatures;
   size_t fBatchIndex;

   size_t fNStreams;                            ///< Number of buffer pairs.
   std::vector<DeviceBuffer_t> fDeviceBuffers;
   std::vector<HostBuffer_t>   fHostBuffers;

   std::vector<size_t> fSampleIndices; ///< Ordering of the samples in the epoch.

public:

   TDataLoader(const Data_t & data, size_t nSamples, size_t batchSize,
               size_t nInputFeatures, size_t nOutputFeatures, size_t nStreams = 1);
   TDataLoader(const TDataLoader  &) = default;
   TDataLoader(      TDataLoader &&) = default;
   TDataLoader & operator=(const TDataLoader  &) = default;
   TDataLoader & operator=(      TDataLoader &&) = default;

   /** Copy input matrix into the given host buffer. Function to be specialized by
    *  the architecture-specific backend. */
   void  CopyInput(HostBuffer_t &buffer, IndexIterator_t begin, size_t batchSize);
   /** Copy output matrix into the given host buffer. Function to be specialized
    * by the architecture-spcific backend. */
   void CopyOutput(HostBuffer_t &buffer, IndexIterator_t begin, size_t batchSize);
   /** Copy weight matrix into the given host buffer. Function to be specialized
    * by the architecture-spcific backend. */
   void CopyWeights(HostBuffer_t &buffer, IndexIterator_t begin, size_t batchSize);

   BatchIterator_t begin() {return TBatchIterator<Data_t, AArchitecture>(*this);}
   BatchIterator_t end()
   {
      return TBatchIterator<Data_t, AArchitecture>(*this, fNSamples / fBatchSize);
   }

   /** Shuffle the order of the samples in the batch. The shuffling is indirect,
    *  i.e. only the indices are shuffled. No input data is moved by this
    * routine. */
   void Shuffle();

   /** Return the next batch from the training set. The TDataLoader object
    *  keeps an internal counter that cycles over the batches in the training
    *  set. */
   TBatch<AArchitecture> GetBatch();

};

//
// TBatch Class.
//______________________________________________________________________________
template <typename AArchitecture>
TBatch<AArchitecture>::TBatch(Matrix_t &inputMatrix, Matrix_t &outputMatrix, Matrix_t &weightMatrix)
   : fInputMatrix(inputMatrix), fOutputMatrix(outputMatrix), fWeightMatrix(weightMatrix)
{
    // Nothing to do here.
}

//
// TDataLoader Class.
//______________________________________________________________________________
template<typename Data_t, typename AArchitecture>
TDataLoader<Data_t, AArchitecture>::TDataLoader(
    const Data_t & data, size_t nSamples, size_t batchSize,
    size_t nInputFeatures, size_t nOutputFeatures, size_t nStreams)
    : fData(data), fNSamples(nSamples), fBatchSize(batchSize),
      fNInputFeatures(nInputFeatures), fNOutputFeatures(nOutputFeatures),
      fBatchIndex(0), fNStreams(nStreams), fDeviceBuffers(), fHostBuffers(),
      fSampleIndices()
{
   size_t inputMatrixSize  = fBatchSize * fNInputFeatures;
   size_t outputMatrixSize = fBatchSize * fNOutputFeatures;
   size_t weightMatrixSize = fBatchSize;

   for (size_t i = 0; i < fNStreams; i++)
   {
      fHostBuffers.push_back(HostBuffer_t(inputMatrixSize + outputMatrixSize + weightMatrixSize));
      fDeviceBuffers.push_back(DeviceBuffer_t(inputMatrixSize + outputMatrixSize + weightMatrixSize));
   }

   fSampleIndices.reserve(fNSamples);
   for (size_t i = 0; i < fNSamples; i++) {
      fSampleIndices.push_back(i);
   }
}

//______________________________________________________________________________
template<typename Data_t, typename AArchitecture>
TBatch<AArchitecture> TDataLoader<Data_t, AArchitecture>::GetBatch()
{
   fBatchIndex %= (fNSamples / fBatchSize); // Cycle through samples.


   size_t inputMatrixSize  = fBatchSize * fNInputFeatures;
   size_t outputMatrixSize = fBatchSize * fNOutputFeatures;
   size_t weightMatrixSize = fBatchSize;

   size_t streamIndex = fBatchIndex % fNStreams;
   HostBuffer_t   & hostBuffer   = fHostBuffers[streamIndex];
   DeviceBuffer_t & deviceBuffer = fDeviceBuffers[streamIndex];

   HostBuffer_t inputHostBuffer  = hostBuffer.GetSubBuffer(0, inputMatrixSize);
   HostBuffer_t outputHostBuffer = hostBuffer.GetSubBuffer(inputMatrixSize,
                                                           outputMatrixSize);
   HostBuffer_t weightHostBuffer = hostBuffer.GetSubBuffer(inputMatrixSize + outputMatrixSize, weightMatrixSize);

   DeviceBuffer_t inputDeviceBuffer  = deviceBuffer.GetSubBuffer(0, inputMatrixSize);
   DeviceBuffer_t outputDeviceBuffer = deviceBuffer.GetSubBuffer(inputMatrixSize,
                                                                 outputMatrixSize);
   DeviceBuffer_t weightDeviceBuffer = deviceBuffer.GetSubBuffer(inputMatrixSize + outputMatrixSize, weightMatrixSize);

   size_t sampleIndex = fBatchIndex * fBatchSize;
   IndexIterator_t sampleIndexIterator = fSampleIndices.begin() + sampleIndex;

   CopyInput(inputHostBuffer,   sampleIndexIterator, fBatchSize);
   CopyOutput(outputHostBuffer, sampleIndexIterator, fBatchSize);
   CopyWeights(weightHostBuffer, sampleIndexIterator, fBatchSize);

   deviceBuffer.CopyFrom(hostBuffer);
   Matrix_t  inputMatrix(inputDeviceBuffer,  fBatchSize, fNInputFeatures);
   Matrix_t outputMatrix(outputDeviceBuffer, fBatchSize, fNOutputFeatures);
   Matrix_t weightMatrix(weightDeviceBuffer, fBatchSize, fNOutputFeatures);

   fBatchIndex++;
   return TBatch<AArchitecture>(inputMatrix, outputMatrix, weightMatrix);
}

//______________________________________________________________________________
template<typename Data_t, typename AArchitecture>
void TDataLoader<Data_t, AArchitecture>::Shuffle()
{
   std::shuffle(fSampleIndices.begin(), fSampleIndices.end(), std::default_random_engine{});
}

} // namespace DNN
} // namespace TMVA

#endif
