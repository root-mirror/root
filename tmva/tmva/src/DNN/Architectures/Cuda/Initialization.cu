// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 14/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 /////////////////////////////////////////////////////////////
 // Implementation of the initialization functions for CUDA //
 // Architectures                                           //
 /////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include "TMatrix.h"
#include "TMVA/DNN/Architectures/Cuda.h"
#include "Kernels.cuh"

namespace TMVA
{
namespace DNN
{

template <typename AFloat>
TRandom * TCuda<AFloat>::fgRandomGen = nullptr;
//______________________________________________________________________________
template<typename AFloat>
void TCuda<AFloat>::SetRandomSeed(size_t seed)
{
   if (!fgRandomGen) fgRandomGen = new TRandom3();
   fgRandomGen->SetSeed(seed);
}
template<typename AFloat>
TRandom & TCuda<AFloat>::GetRandomGenerator()
{
   if (!fgRandomGen) fgRandomGen = new TRandom3(0);
   return *fgRandomGen;
}
//______________________________________________________________________________
template<typename AFloat>
void TCuda<AFloat>::InitializeGauss(TCudaMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();
   TMatrixT<AFloat> B(m, n);

   Double_t sigma = sqrt(2.0 / ((Double_t) n));

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         B(i,j) = rand.Gaus(0.0, sigma);
      }
   }
   A = B;
}

//______________________________________________________________________________
template<typename AFloat>
void TCuda<AFloat>::InitializeUniform(TCudaMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();
   TMatrixT<AFloat> B(m, n);

   Double_t range = sqrt(2.0 / ((Double_t) n));

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         B(i,j) = rand.Uniform(-range, range);
      }
   }
   A = B;
}

//______________________________________________________________________________
///  Truncated normal initialization (Glorot, called also Xavier normal)
///  The values are sample with a normal distribution with stddev = sqrt(2/N_input + N_output) and
///   values larger than 2 * stddev are discarded
///  See Glorot & Bengio, AISTATS 2010 - http://jmlr.org/proceedings/papers/v9/glorot10a/glorot10a.pdf
template<typename AFloat>
void TCuda<AFloat>::InitializeGlorotNormal(TCudaMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();
   TMatrixT<AFloat> B(m, n);

   AFloat sigma = sqrt(2.0 /( ((AFloat) n) + ((AFloat) m)) );

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         AFloat value = 0;
         do {
            value = rand.Gaus(0.0, sigma);
         } while ( std::abs(value) > 2*sigma);
         B(i,j) = value;
      }
   }
   A = B;
}

//______________________________________________________________________________
/// Sample from a uniform distribution in range [ -lim,+lim] where
///  lim = sqrt(6/N_in+N_out).
/// This initialization is also called Xavier uniform
/// see Glorot & Bengio, AISTATS 2010 - http://jmlr.org/proceedings/papers/v9/glorot10a/glorot10a.pdf
template<typename AFloat>
void TCuda<AFloat>::InitializeGlorotUniform(TCudaMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();
   TMatrixT<AFloat> B(m, n);

   AFloat range = sqrt(6.0 /( ((AFloat) n) + ((AFloat) m)) );

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         B(i,j) = rand.Uniform(-range, range);
      }
   }
   printf("initialize glorotuniform \n");
   B.Print();
   A = B;
}

//______________________________________________________________________________
template<typename AFloat>
void TCuda<AFloat>::InitializeIdentity(TCudaMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();
   TMatrixT<AFloat> B(m, n);

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n ; j++) {
         B(i,j) = 0.0;
      }

      if (i < n) {
         B(i,i) = 1.0;
      }
   }
   A = B;
}

//______________________________________________________________________________
template<typename AFloat>
void TCuda<AFloat>::InitializeZero(TCudaMatrix<AFloat> & A)
{
   // use fast zero initialization on the device
   A.Zero();
}
//______________________________________________________________________________
template <typename AFloat>
void TCuda<AFloat>::InitializeZero(TCudaTensor<AFloat> &T)
{
   // use fast zero initialization on the device
   T.Zero();
}

} // namespace DNN
} // namespace TMVA
