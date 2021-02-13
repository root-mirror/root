// @(#)root/tmva/tmva/dnn:$Id$ // Author: Simon Pfreundschuh 20/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////
// Test arithmetic functions defined on matrices and compare the //
// results to the reference implementation.                      //
///////////////////////////////////////////////////////////////////

#include "TMatrix.h"
#include "Utility.h"
#include "TMVA/DNN/Architectures/Reference.h"

#include <cmath>
#include <vector>

/** Test multiplication (standard, transpose ) operation on
 *  architecture specific matrix types and compare with results
 *  obtained with TMatrixT.
 */
//______________________________________________________________________________
template<typename Architecture_t>
auto testMultiplication(size_t ntests)
    -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;
   Scalar_t mean = 5.0, sigma = 2.0;

   for (size_t t = 0; t < ntests; t++) {
      size_t m, n, k;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;
      k = rand() % 100 + 1;

      TMatrixT<Scalar_t> ARef(m,k), A2Ref(m,k), ATRef(k,m) , BRef(k,n), BTRef(n,k), CRef(m,n);
      TMVA::DNN::randomMatrix(ARef, mean, sigma);
      TMVA::DNN::randomMatrix(A2Ref, mean, sigma);
      TMVA::DNN::randomMatrix(ATRef, mean, sigma);
      TMVA::DNN::randomMatrix(BRef, mean, sigma);
      TMVA::DNN::randomMatrix(BTRef, mean, sigma);
      Matrix_t A(ARef), A2(A2Ref), AT(ATRef), B(BRef), BT(BTRef),  C(CRef);

      // A * B
      CRef.Mult(ARef,BRef);
      Architecture_t::Multiply(C, A, B);
      Scalar_t error = TMVA::DNN::maximumRelativeError(C, CRef);
      maximumError   = std::max(error, maximumError);

      // A^T * B
      CRef.TMult(ATRef,BRef);
      Architecture_t::TransposeMultiply(C, AT, B);
      error = TMVA::DNN::maximumRelativeError(C, CRef);
      maximumError   = std::max(error, maximumError);

      // A * B^T
      CRef.MultT(ARef,BTRef);
      Architecture_t::MultiplyTranspose(C, A, BT);
      error = TMVA::DNN::maximumRelativeError(C, CRef);
      maximumError   = std::max(error, maximumError);

      // A .* B
      for (size_t i = 0; i < (size_t) ARef.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t) ARef.GetNcols(); j++) {
            ARef(i,j) *= A2Ref(i,j);
         }
      }
      Architecture_t::Hadamard(A, A2);
      error = TMVA::DNN::maximumRelativeError(A, ARef);
      maximumError   = std::max(error, maximumError);
   }

   return maximumError;
}

/** Test  hadamard operation on
 *  architecture specific matrix and tensor types and compare with results
 *  obtained with TMatrixT.
 */
//______________________________________________________________________________
template<typename Architecture_t>
auto testHadamrdMultiplication(size_t ntests)
    -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Tensor_t = typename Architecture_t::Tensor_t;

   Scalar_t maximumError = 0.0;
   Scalar_t mean = 5.0, sigma = 2.0;

   for (size_t t = 0; t < ntests; t++) {
      size_t m = rand() % 100 + 1;
      size_t n = rand() % 100 + 1;
      //k = rand() % 100 + 1;

      TMatrixT<Scalar_t> ARef(m,n), A2Ref(m,n);
      TMVA::DNN::randomMatrix(ARef, mean, sigma);
      TMVA::DNN::randomMatrix(A2Ref, mean, sigma);
     
      Matrix_t A(ARef), A2(A2Ref);


      // A .* B
      for (size_t i = 0; i < (size_t) ARef.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t) ARef.GetNcols(); j++) {
            ARef(i,j) *= A2Ref(i,j);
         }
      }
      Architecture_t::Hadamard(A, A2);
      Scalar_t error = TMVA::DNN::maximumRelativeError(A, ARef);
      maximumError   = std::max(error, maximumError);
      if (maximumError > 1.E-3) { 
         std::cout << "Hadamard Matrix test - max error =  " << maximumError << std::endl;
      }   
   }
   // tensor test
   //ntests = 1;
   for (size_t t = 0; t < ntests; t++) {
      size_t m, n, k;
      m = rand() % 10 + 1;
      n = rand() % 100 + 1;
      k = rand() % 100 + 1;
      // m = 1; n= 2;  
      // k = 3; 

      std::vector<TMatrixT<Scalar_t>> ARef, A2Ref; 
      Tensor_t A(m,n,k);
      Tensor_t A2(m,n,k);
      for (size_t i = 0; i < m; ++i) { 
         ARef.emplace_back(n,k);
         A2Ref.emplace_back(n,k);
         TMVA::DNN::randomMatrix(ARef[i], mean, sigma);
         TMVA::DNN::randomMatrix(A2Ref[i], mean, sigma);

         Matrix_t mA = A.At(i).GetMatrix(); 
         Matrix_t mA2 = A2.At(i).GetMatrix(); 
         Architecture_t::CopyDiffArch(mA, ARef[i]);
         Architecture_t::CopyDiffArch(mA2, A2Ref[i]);
      }
      


      // A .* B
      for (size_t b = 0; b < m; ++b) { 
         for (size_t i = 0; i < (size_t) ARef[b].GetNrows(); i++) {
            for (size_t j = 0; j < (size_t) ARef[b].GetNcols(); j++) {
               ARef[b](i,j) *= A2Ref[b](i,j);
            }
         }
      }
      // convert ARef to the specific architecture
      Tensor_t tARef(m,n,k);
      for (size_t b = 0; b < m; ++b) { 
         Matrix_t mARef = tARef.At(b).GetMatrix(); 
         Architecture_t::CopyDiffArch(mARef, ARef[b]);
      }

      Architecture_t::Hadamard(A, A2);

      
      
      Scalar_t error = TMVA::DNN::maximumRelativeErrorTensor(A, tARef);
      maximumError   = std::max(error, maximumError);

      if (maximumError > 1.E-3) { 
         //Architecture_t::PrintTensor(A0, "input1");
         Architecture_t::PrintTensor(A2, "input2");

         Architecture_t::PrintTensor(A,"output");
         Architecture_t::PrintTensor(tARef,"ref output");
      }
   }

   return maximumError;
}

/** Test the summing over columns by summing by the sums obtained
 *  from a matrix filled with column indices as elements.
 */
//______________________________________________________________________________
template<typename Architecture_t>
auto testSumColumns(size_t ntests)
    -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;
   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m,n), BRef(n,1);

      for (size_t i = 0; i < (size_t) ARef.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t) ARef.GetNcols(); j++) {
            ARef(i,j) = j;
            if (i == 0) BRef(j, 0) = m * j;
         }
      }

      Matrix_t A(ARef), B(n, 1);
      Architecture_t::SumColumns(B, A);

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>) B ,BRef);
      maximumError   = std::max(error, maximumError);
   }
   return maximumError;
}

/** Test the addition of a constant to every element of
 *  the given matrix.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testConstAdd(size_t ntests) -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;
   Scalar_t beta;

   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      beta = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m, n);
      TMVA::DNN::randomMatrix(ARef);

      Matrix_t A(ARef);
      Architecture_t::ConstAdd(A, beta);

      for (size_t i = 0; i < (size_t)ARef.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t)ARef.GetNcols(); j++) {
            ARef(i, j) += (Double_t)beta;
         }
      }

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>)A, ARef);
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
}

/** Test the multiplication of a constant to every element of
 *  the given matrix.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testConstMult(size_t ntests) -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;
   Scalar_t beta;

   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      beta = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m, n);
      TMVA::DNN::randomMatrix(ARef);

      Matrix_t A(ARef);
      Architecture_t::ConstMult(A, beta);

      for (size_t i = 0; i < (size_t)ARef.GetNrows(); i++) {
         for (size_t j = 0; j < (size_t)ARef.GetNcols(); j++) {
            ARef(i, j) *= (Double_t)beta;
         }
      }

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>)A, ARef);
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
}

/** Test the application of reciprocating every element of
 *  the given matrix.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testReciprocalElementWise(size_t ntests) -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;

   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m, n);
      TMVA::DNN::randomMatrix(ARef);

      Matrix_t A(ARef);
      Architecture_t::ReciprocalElementWise(A);

      TMVA::DNN::applyMatrix(ARef, [](double x) { return 1.0 / x; });

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>)A, ARef);
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
}

/** Test the application of squaring every element of
 *  the given matrix.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testSquareElementWise(size_t ntests) -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;

   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m, n);
      TMVA::DNN::randomMatrix(ARef);

      Matrix_t A(ARef);
      Architecture_t::SquareElementWise(A);

      TMVA::DNN::applyMatrix(ARef, [](double x) { return x * x; });

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>)A, ARef);
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
}

/** Test the application of taking square root of every element of
 *  the given matrix.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testSqrtElementWise(size_t ntests) -> typename Architecture_t::Scalar_t
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;

   Scalar_t maximumError = 0.0;
   Scalar_t mean = 5.0, sigma = 2.0;

   for (size_t t = 0; t < ntests; t++) {

      Scalar_t error;

      size_t m, n;
      m = rand() % 100 + 1;
      n = rand() % 100 + 1;

      TMatrixT<Double_t> ARef(m, n);
      TMVA::DNN::randomMatrix(ARef, mean, sigma);

      Matrix_t A(ARef);
      Architecture_t::SqrtElementWise(A);

      TMVA::DNN::applyMatrix(ARef, [](double x) { return sqrt(x); });

      error = TMVA::DNN::maximumRelativeError((TMatrixT<Scalar_t>)A, ARef);
      maximumError = std::max(error, maximumError);
   }
   return maximumError;
}
