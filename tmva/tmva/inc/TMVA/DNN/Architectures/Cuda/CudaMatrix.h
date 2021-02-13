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

#ifndef TMVA_DNN_ARCHITECTURES_CUDA_CUDAMATRIX
#define TMVA_DNN_ARCHITECTURES_CUDA_CUDAMATRIX

// in case we compile C++ code with std-17 and cuda with lower standard
// use experimental string_view, otherwise keep as is
#include "RConfigure.h"
#ifdef R__HAS_STD_STRING_VIEW
#ifndef R__CUDA_HAS_STD_STRING_VIEW
#undef R__HAS_STD_STRING_VIEW
#define R__HAS_STD_EXPERIMENTAL_STRING_VIEW
#endif
#endif

#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "curand_kernel.h"

#include "TMatrixT.h"
#include "CudaBuffers.h"

#define CUDACHECK(ans) {cudaError((ans), __FILE__, __LINE__); }

namespace TMVA {
namespace DNN {

/** Function to check cuda return code. Taken from
 * http://stackoverflow.com/questions/14038589/
 */
inline void cudaError(cudaError_t code, const char *file, int line, bool abort=true);

//____________________________________________________________________________
//
// Cuda Device Reference
//____________________________________________________________________________

/** TCudaDeviceReference
 *
 * Helper class emulating lvalue references for AFloat values that are
 * physically on the device. Allows for example to assign to matrix elements.
 * Note that device access through CudaDeviceReferences enforces synchronization
 * with all streams and thus qualifies as performance killer. Only used for
 * testing.
 */
template<typename AFloat>
class TCudaDeviceReference
{
private:

    AFloat * fDevicePointer;

public:

    TCudaDeviceReference(AFloat * devicePointer);

    operator AFloat();

    void operator=(const TCudaDeviceReference &other);
    void operator=(AFloat value);
    void operator+=(AFloat value);
    void operator-=(AFloat value);
};

//____________________________________________________________________________
//
// Cuda Matrix
//____________________________________________________________________________

/** TCudaMatrix Class
 *
 * The TCudaMatrix class represents matrices on a CUDA device. The elements
 * of the matrix are stored in a TCudaDeviceBuffer object which takes care of
 * the allocation and freeing of the device memory. TCudaMatrices are lightweight
 * object, that means on assignment and copy creation only a shallow copy is
 * performed and no new element buffer allocated. To perform a deep copy use
 * the static Copy method of the TCuda architecture class.
 *
 * The TCudaDeviceBuffer has an associated cuda stream, on which the data is
 * transferred to the device. This stream can be accessed through the
 * GetComputeStream member function and used to synchronize computations.
 *
 * The TCudaMatrix class also holds static references to CUDA resources.
 * Those are the cublas handle, a buffer of curand states for the generation
 * of random numbers as well as a vector containing ones, which is used for
 * summing column matrices using matrix-vector multiplication. The class also
 * has a static buffer for returning results from the device.
 *
 */
template<typename AFloat>
class TCudaMatrix
{
public:

private:

   static size_t          fInstances;    ///< Current number of matrix instances.
   static cublasHandle_t  fCublasHandle;
   static AFloat        * fDeviceReturn; ///< Buffer for kernel return values.
   static AFloat        * fOnes;         ///< Vector used for summations of columns.
   static size_t          fNOnes;        ///< Current length of the one vector.
   static curandState_t * fCurandStates;
   static size_t          fNCurandStates;


   size_t                    fNRows;
   size_t                    fNCols;
   TCudaDeviceBuffer<AFloat> fElementBuffer;

public:

   static Bool_t gInitializeCurand;

   static AFloat * GetOnes() {return fOnes;}

   TCudaMatrix();
   TCudaMatrix(size_t i, size_t j);
   TCudaMatrix(const TMatrixT<AFloat> &);
   TCudaMatrix(TCudaDeviceBuffer<AFloat> buffer, size_t m, size_t n);

   TCudaMatrix(const TCudaMatrix  &) = default;
   TCudaMatrix(      TCudaMatrix &&) = default;
   TCudaMatrix & operator=(const TCudaMatrix  &) = default;
   TCudaMatrix & operator=(      TCudaMatrix &&) = default;
   ~TCudaMatrix() = default;

   /** Convert cuda matrix to Root TMatrix. Performs synchronous data transfer. */
   operator TMatrixT<AFloat>() const;

   inline cudaStream_t GetComputeStream() const;
   inline void         SetComputeStream(cudaStream_t stream);
   /** Set the return buffer on the device to the specified value. This is
    * required for example for reductions in order to initialize the
    * accumulator. */
   inline static void ResetDeviceReturn(AFloat value = 0.0);
   /** Transfer the value in the device return buffer to the host. This
    *  tranfer is synchronous */
   inline static AFloat GetDeviceReturn();
   /** Return device pointer to the device return buffer */
   inline static AFloat *        GetDeviceReturnPointer() {return fDeviceReturn;}
   inline static curandState_t * GetCurandStatesPointer() {return fCurandStates;}

   /** Blocking synchronization with the associated compute stream, if it's
    * not the default stream. */
   inline void Synchronize(const TCudaMatrix &) const;

   static size_t GetNDim() {return 2;}
   size_t GetNrows() const {return fNRows;}
   size_t GetNcols() const {return fNCols;}
   size_t GetNoElements() const {return fNRows * fNCols;}

   const AFloat * GetDataPointer() const {return fElementBuffer;}
   AFloat *       GetDataPointer()       {return fElementBuffer;}
   const cublasHandle_t & GetCublasHandle() const    {return fCublasHandle;}

   inline  TCudaDeviceBuffer<AFloat> GetDeviceBuffer() const { return fElementBuffer;}

   /** Access to elements of device matrices provided through TCudaDeviceReference
    *  class. Note that access is synchronous end enforces device synchronization
    *  on all streams. Only used for testing. */
   TCudaDeviceReference<AFloat> operator()(size_t i, size_t j) const;

   void Print() const {
      TMatrixT<AFloat> mat(*this);
      mat.Print();
   }

   void Zero() {
      cudaMemset(GetDataPointer(), 0, sizeof(AFloat) * GetNoElements());
   }


private:

   /** Initializes all shared devices resource and makes sure that a sufficient
    *  number of curand states are allocated on the device and initialized as
    *  well as that the one-vector for the summation over columns has the right
    *  size. */
   void InitializeCuda();
   void InitializeCurandStates();

};

//
// Inline Functions.
//______________________________________________________________________________
inline void cudaError(cudaError_t code, const char *file, int line, bool abort)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//______________________________________________________________________________
template<typename AFloat>
TCudaDeviceReference<AFloat>::TCudaDeviceReference(AFloat * devicePointer)
    : fDevicePointer(devicePointer)
{
   // Nothing to do here.
}

//______________________________________________________________________________
template<typename AFloat>
TCudaDeviceReference<AFloat>::operator AFloat()
{
    AFloat buffer;
    cudaMemcpy(& buffer, fDevicePointer, sizeof(AFloat),
               cudaMemcpyDeviceToHost);
    return buffer;
}

//______________________________________________________________________________
template<typename AFloat>
void TCudaDeviceReference<AFloat>::operator=(const TCudaDeviceReference &other)
{
   cudaMemcpy(fDevicePointer, other.fDevicePointer, sizeof(AFloat),
              cudaMemcpyDeviceToDevice);
}

//______________________________________________________________________________
template<typename AFloat>
void TCudaDeviceReference<AFloat>::operator=(AFloat value)
{
   AFloat buffer = value;
   cudaMemcpy(fDevicePointer, & buffer, sizeof(AFloat),
              cudaMemcpyHostToDevice);
}

//______________________________________________________________________________
template<typename AFloat>
void TCudaDeviceReference<AFloat>::operator+=(AFloat value)
{
   AFloat buffer;
   cudaMemcpy(& buffer, fDevicePointer, sizeof(AFloat),
              cudaMemcpyDeviceToHost);
   buffer += value;
   cudaMemcpy(fDevicePointer, & buffer, sizeof(AFloat),
              cudaMemcpyHostToDevice);
}

//______________________________________________________________________________
template<typename AFloat>
void TCudaDeviceReference<AFloat>::operator-=(AFloat value)
{
   AFloat buffer;
   cudaMemcpy(& buffer, fDevicePointer, sizeof(AFloat),
              cudaMemcpyDeviceToHost);
   buffer -= value;
   cudaMemcpy(fDevicePointer, & buffer, sizeof(AFloat),
              cudaMemcpyHostToDevice);
}

//______________________________________________________________________________
template<typename AFloat>
inline cudaStream_t TCudaMatrix<AFloat>::GetComputeStream() const
{
   return fElementBuffer.GetComputeStream();
}

//______________________________________________________________________________
template<typename AFloat>
inline void TCudaMatrix<AFloat>::SetComputeStream(cudaStream_t stream)
{
   return fElementBuffer.SetComputeStream(stream);
}

//______________________________________________________________________________
template<typename AFloat>
inline void TCudaMatrix<AFloat>::Synchronize(const TCudaMatrix &A) const
{
   cudaEvent_t event;
   cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
   cudaEventRecord(event, A.GetComputeStream());
   cudaStreamWaitEvent(fElementBuffer.GetComputeStream(), event, 0);
   cudaEventDestroy(event);
}

//______________________________________________________________________________
template<typename AFloat>
inline void TCudaMatrix<AFloat>::ResetDeviceReturn(AFloat value)
{
   AFloat buffer = value;
   cudaMemcpy(fDeviceReturn, & buffer, sizeof(AFloat), cudaMemcpyHostToDevice);
}

//______________________________________________________________________________
template<typename AFloat>
inline AFloat TCudaMatrix<AFloat>::GetDeviceReturn()
{
   AFloat buffer;
   cudaMemcpy(& buffer, fDeviceReturn, sizeof(AFloat), cudaMemcpyDeviceToHost);
   return buffer;
}

//______________________________________________________________________________
template<typename AFloat>
TCudaDeviceReference<AFloat> TCudaMatrix<AFloat>::operator()(size_t i, size_t j) const
{
    AFloat * elementPointer = fElementBuffer;
    elementPointer += j * fNRows + i;
    return TCudaDeviceReference<AFloat>(elementPointer);
}

} // namespace DNN
} // namespace TMVA

#endif
