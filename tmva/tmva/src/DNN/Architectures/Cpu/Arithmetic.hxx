// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 20/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////
//  Implementation of Helper arithmetic functions for the //
// multi-threaded CPU implementation of DNNs.             //
////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cpu.h"

#ifdef R__HAS_TMVACPU
#include "TMVA/DNN/Architectures/Cpu/Blas.h"
#else
#include "TMVA/DNN/Architectures/Reference.h"
#include "TVectorT.h"
#endif

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"

//#include "tbb/tbb.h"

#pragma GCC diagnostic pop
#endif

namespace TMVA
{
namespace DNN
{

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::Multiply(TCpuMatrix<AReal> &C,
                            const TCpuMatrix<AReal> &A,
                            const TCpuMatrix<AReal> &B)
{
    int m = (int) A.GetNrows();
    int k = (int) A.GetNcols();
    int n = (int) B.GetNcols();

    R__ASSERT((int) C.GetNrows() == m);
    R__ASSERT((int) C.GetNcols() == n);
    R__ASSERT((int) B.GetNrows() == k);

#ifdef R__HAS_TMVACPU

    char transa = 'N';
    char transb = 'N';

    AReal alpha = 1.0;
    AReal beta  = 0.0;

    const AReal * APointer = A.GetRawDataPointer();
    const AReal * BPointer = B.GetRawDataPointer();
          AReal * CPointer = C.GetRawDataPointer();

    ::TMVA::DNN::Blas::Gemm(&transa, &transb, &m, &n, &k, &alpha,
                            APointer, &m, BPointer, &k, &beta, CPointer, &m);
#else
   TMatrixT<AReal> tmp(C.GetNrows(), C.GetNcols());
   tmp.Mult(A,B);
   C = tmp;
#endif
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::TransposeMultiply(TCpuMatrix<AReal> &C,
                                     const TCpuMatrix<AReal> &A,
                                     const TCpuMatrix<AReal> &B,
                                     AReal alpha, AReal beta)
{
#ifdef R__HAS_TMVACPU
    int m = (int) A.GetNcols();
    int k = (int) A.GetNrows();
    int n = (int) B.GetNcols();

    R__ASSERT((int) C.GetNrows() == m);
    R__ASSERT((int) C.GetNcols() == n);
    R__ASSERT((int) B.GetNrows() == k);

    char transa = 'T';
    char transb = 'N';

    //AReal alpha = 1.0;
    //AReal beta  = 0.0;

    const AReal *APointer = A.GetRawDataPointer();
    const AReal *BPointer = B.GetRawDataPointer();
          AReal *CPointer = C.GetRawDataPointer();

    ::TMVA::DNN::Blas::Gemm(&transa, &transb, &m, &n, &k, &alpha,
                            APointer, &k, BPointer, &k, &beta, CPointer, &m);
#else
   TMatrixT<AReal> tmp(C.GetNrows(), C.GetNcols());
   tmp.TMult(A, B);
   tmp = alpha * tmp;
   if (beta != 0.0) {
      TMatrixT<AReal> tmp0(C);
      tmp = tmp + beta * tmp0;
   }
   C = tmp;
#endif
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::Hadamard(TCpuMatrix<AReal> &B,
                            const TCpuMatrix<AReal> &A)
{
   const AReal *dataA      = A.GetRawDataPointer();
   AReal *dataB      = B.GetRawDataPointer();

   size_t nElements =  A.GetNoElements();
   R__ASSERT(B.GetNoElements() == nElements);
   size_t nSteps = TCpuMatrix<AReal>::GetNWorkItems(nElements);

   auto f = [&](UInt_t workerID)
   {
      for (size_t j = 0; j < nSteps; ++j) {
         size_t idx = workerID+j;
         if (idx >= nElements) break;
         dataB[idx] *= dataA[idx];
      }
      return 0;
   };

   if (nSteps < nElements) {
#ifdef DL_USE_MTE
      B.GetThreadExecutor().Foreach(f, ROOT::TSeqI(0,nElements,nSteps));
#else
      for (size_t i = 0;  i < nElements ; i+= nSteps)
         f(i);
#endif
   }
   else {
      f(0);
   }
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::Hadamard(TCpuTensor<AReal> &B,
                            const TCpuTensor<AReal> &A)
{
   const AReal *dataA      = A.GetRawDataPointer();
   AReal *dataB      = B.GetRawDataPointer();

   size_t nElements =  A.GetNoElements();
   R__ASSERT(B.GetNoElements() == nElements);
   size_t nSteps = TCpuMatrix<AReal>::GetNWorkItems(nElements);

   auto f = [&](UInt_t workerID)
   {
      for (size_t j = 0; j < nSteps; ++j) {
         size_t idx = workerID+j;
         if (idx >= nElements) break;
         dataB[idx] *= dataA[idx];
      }
      return 0;
   };

   if (nSteps < nElements) {
#ifdef DL_USE_MTE
      TMVA::Config::Instance().GetThreadExecutor().Foreach(f, ROOT::TSeqI(0,nElements,nSteps));
#else
      for (size_t i = 0;  i < nElements ; i+= nSteps)
         f(i);
#endif
   }
   else {
      f(0);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Checks two matrices for element-wise equality.
/// \tparam AReal An architecture-specific floating point number type.
/// \param A The first matrix.
/// \param B The second matrix.
/// \param epsilon Equality tolerance, needed to address floating point arithmetic.
/// \return Whether the two matrices can be considered equal element-wise
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename AReal>
bool TCpu<AReal>::AlmostEquals(const TCpuMatrix<AReal> &A, const TCpuMatrix<AReal> &B, double epsilon)
{
    if (A.GetNrows() != B.GetNrows() || A.GetNcols() != B.GetNcols()) {
        Fatal("AlmostEquals", "The passed matrices have unequal shapes.");
    }

    const AReal *dataA = A.GetRawDataPointer();
    const AReal *dataB = B.GetRawDataPointer();
    size_t nElements =  A.GetNoElements();

    for(size_t i = 0; i < nElements; i++) {
        if(fabs(dataA[i] - dataB[i]) > epsilon) return false;
    }
    return true;
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::SumColumns(TCpuMatrix<AReal> &B,
                              const TCpuMatrix<AReal> &A,
                              AReal alpha, AReal beta)
{

   int m = (int) A.GetNrows();
   int n = (int) A.GetNcols();

   assert((int) B.GetNoElements() >= n);

#ifdef R__HAS_TMVACPU
   int inc = 1;
   char   trans   = 'T';

   const AReal * APointer = A.GetRawDataPointer();
         AReal * BPointer = B.GetRawDataPointer();

   // compute B = alpha * A * I + beta * B

   ::TMVA::DNN::Blas::Gemv(&trans, &m, &n, &alpha, APointer, &m,
                           TCpuMatrix<AReal>::GetOnePointer(), &inc,
                           &beta, BPointer, &inc);
#else
   TMatrixT<AReal> tA(A);
   tA.T();
   TVectorT<AReal> ones(m, TCpuMatrix<AReal>::GetOnePointer());
   TVectorT<AReal> tmp(n, B.GetRawDataPointer());
   assert(B.GetNrows() == 1 || B.GetNcols() == 1);
   tmp = alpha * tA * ones + beta * tmp;
   // copy result buffer in B matrix
   std::copy(tmp.GetMatrixArray(), tmp.GetMatrixArray() + n, B.GetRawDataPointer());
#endif
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::ScaleAdd(TCpuMatrix<AReal> &B,
                            const TCpuMatrix<AReal> &A,
                            AReal alpha)
{
#ifdef R__HAS_TMVACPU
   int n = (int) (A.GetNcols() * A.GetNrows());
   int inc = 1;

   const AReal *x = A.GetRawDataPointer();
   AReal *y = B.GetRawDataPointer();

   ::TMVA::DNN::Blas::Axpy(&n, &alpha, x, &inc, y, &inc);
#else
   TMatrixT<AReal> tmp(B);
   TReference<AReal>::ScaleAdd(tmp, A, alpha);
   B = tmp;
#endif
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::Copy(TCpuMatrix<AReal> &B,
                        const TCpuMatrix<AReal> &A)
{
   auto f = [](AReal x) {return x;};
   B.MapFrom(f, A);
}


//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::ScaleAdd(TCpuTensor<AReal> &B,
                            const TCpuTensor<AReal> &A,
                            AReal alpha)
{
   // should re-implemented at tensor level
   for (size_t i = 0; i < B.GetFirstSize(); ++i) {
      TCpuMatrix<AReal> B_m = B.At(i).GetMatrix();
      ScaleAdd(B_m, A.At(i).GetMatrix(), alpha);
   }
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::Copy(TCpuTensor<AReal> &B,
                            const TCpuTensor<AReal> &A)
{

   auto f = [](AReal x) {return x;};
   B.MapFrom(f, A);
}

//____________________________________________________________________________
template <typename AReal>
void TCpu<AReal>::ConstAdd(TCpuMatrix<AReal> &A, AReal beta)
{
   auto f = [beta](AReal x) { return x + beta; };
   A.Map(f);
}

//____________________________________________________________________________
template <typename AReal>
void TCpu<AReal>::ConstMult(TCpuMatrix<AReal> &A, AReal beta)
{
   auto f = [beta](AReal x) { return x * beta; };
   A.Map(f);
}

//____________________________________________________________________________
template <typename AReal>
void TCpu<AReal>::ReciprocalElementWise(TCpuMatrix<AReal> &A)
{
   auto f = [](AReal x) { return 1.0 / x; };
   A.Map(f);
}

//____________________________________________________________________________
template <typename AReal>
void TCpu<AReal>::SquareElementWise(TCpuMatrix<AReal> &A)
{
   auto f = [](AReal x) { return x * x; };
   A.Map(f);
}

//____________________________________________________________________________
template <typename AReal>
void TCpu<AReal>::SqrtElementWise(TCpuMatrix<AReal> &A)
{
   auto f = [](AReal x) { return sqrt(x); };
   A.Map(f);
}

/// Adam updates
//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::AdamUpdate(TCpuMatrix<AReal> &A, const TCpuMatrix<AReal> & M, const TCpuMatrix<AReal> & V, AReal alpha, AReal eps)
{
   // ADAM update the weights.
   // Weight = Weight - alpha * M / (sqrt(V) + epsilon)
   AReal * a = A.GetRawDataPointer();
   const AReal * m = M.GetRawDataPointer();
   const AReal * v = V.GetRawDataPointer();
   for (size_t index = 0; index < A.GetNoElements() ; ++index) {
      a[index] = a[index] - alpha * m[index]/( sqrt(v[index]) + eps);
   }
}

//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::AdamUpdateFirstMom(TCpuMatrix<AReal> &A, const TCpuMatrix<AReal> & B, AReal beta)
{
   // First momentum weight gradient update for ADAM
   // Mt = beta1 * Mt-1 + (1-beta1) * WeightGradients
   AReal * a = A.GetRawDataPointer();
   const AReal * b = B.GetRawDataPointer();
   for (size_t index = 0; index < A.GetNoElements() ; ++index) {
      a[index] = beta * a[index] + (1.-beta) * b[index];
   }
}
//____________________________________________________________________________
template<typename AReal>
void TCpu<AReal>::AdamUpdateSecondMom(TCpuMatrix<AReal> &A, const TCpuMatrix<AReal> & B, AReal beta)
{
   // Second momentum weight gradient update for ADAM
   // Vt = beta2 * Vt-1 + (1-beta2) * WeightGradients^2
   AReal * a = A.GetRawDataPointer();
   const AReal * b = B.GetRawDataPointer();
   for (size_t index = 0; index < A.GetNoElements() ; ++index) {
      a[index] = beta * a[index] + (1.-beta) * b[index] * b[index];
   }
}


} // DNN
} // TMVA
