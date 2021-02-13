// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 21/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

 //////////////////////////////////////////////////////////////
 // Implementation of the DNN initialization methods for the //
 // multi-threaded CPU backend.                              //
 //////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include "TMVA/DNN/Architectures/Cpu.h"

namespace TMVA
{
namespace DNN
{

template <typename AFloat_t>
TRandom * TCpu<AFloat_t>::fgRandomGen = nullptr;
//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::SetRandomSeed(size_t seed)
{
   if (!fgRandomGen) fgRandomGen = new TRandom3();
   fgRandomGen->SetSeed(seed);
}
template<typename AFloat>
TRandom & TCpu<AFloat>::GetRandomGenerator()
{
   if (!fgRandomGen) fgRandomGen = new TRandom3(0);
   return *fgRandomGen;
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::InitializeGauss(TCpuMatrix<AFloat> & A)
{
   size_t n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat sigma = sqrt(2.0 / ((AFloat) n));

   for (size_t i = 0; i < A.GetSize(); ++i) {
      A.GetRawDataPointer()[i] = rand.Gaus(0.0, sigma);
   }
   // for (size_t i = 0; i < A.GetSize(); ++i) {
   //    A.GetRawDataPointer()[i] = rand.Gaus(0.0, sigma);
   // }
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::InitializeUniform(TCpuMatrix<AFloat> & A)
{
   //size_t m = A.GetNrows();
   size_t n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat range = sqrt(2.0 / ((AFloat) n));

   // for debugging
   //range = 1;
   //rand.SetSeed(111);

   for (size_t i = 0; i < A.GetSize(); ++i) {
      A.GetRawDataPointer()[i] = rand.Uniform(-range, range);
   }
}

 //______________________________________________________________________________
///  Truncated normal initialization (Glorot, called also Xavier normal)
///  The values are sample with a normal distribution with stddev = sqrt(2/N_input + N_output) and
///   values larger than 2 * stddev are discarded
///  See Glorot & Bengio, AISTATS 2010 - http://jmlr.org/proceedings/papers/v9/glorot10a/glorot10a.pdf
template<typename AFloat>
void TCpu<AFloat>::InitializeGlorotNormal(TCpuMatrix<AFloat> & A)
{
   size_t m,n;
   // for conv layer weights output m is only output depth. It shouild ne multiplied also by filter sizes
   // e.g. 9 for a 3x3 filter. But this information is lost if we use Tensors of dims 2
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat sigma = sqrt(6.0 /( ((AFloat) n) + ((AFloat) m)) );
   // AFloat sigma = sqrt(2.0 /( ((AFloat) m)) );

   size_t nsize = A.GetSize();
   for (size_t i = 0; i < nsize; i++) {
      AFloat value = 0;
      do {
         value = rand.Gaus(0.0, sigma);
      } while (std::abs(value) > 2 * sigma);
      A.GetRawDataPointer()[i] = value;
   }
}

//______________________________________________________________________________
/// Sample from a uniform distribution in range [ -lim,+lim] where
///  lim = sqrt(6/N_in+N_out).
/// This initialization is also called Xavier uniform
/// see Glorot & Bengio, AISTATS 2010 - http://jmlr.org/proceedings/papers/v9/glorot10a/glorot10a.pdf
template<typename AFloat>
void TCpu<AFloat>::InitializeGlorotUniform(TCpuMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows(); // output size
   n = A.GetNcols();  // input size
   // Note that m and n are inverted with respect to cudnn because tensor is here column-wise

   TRandom &  rand = GetRandomGenerator();

   AFloat range = sqrt(6.0 /( ((AFloat) n) + ((AFloat) m)) );

   size_t nsize = A.GetSize();
   for (size_t i = 0; i < nsize; i++) {
      A.GetRawDataPointer()[i] = rand.Uniform(-range, range);
   }
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::InitializeIdentity(TCpuMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j <  n; j++) {
         A(i,j) = 0.0;
         //A(i,j) = 1.0;
      }

      if (i < n) {
         A(i,i) = 1.0;
      }
   }
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::InitializeZero(TCpuMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n ; j++) {
         A(i,j) = 0.0;
      }
   }
}
//______________________________________________________________________________
template <typename AFloat>
void TCpu<AFloat>::InitializeZero(TCpuTensor<AFloat> &A)
{
   size_t n = A.GetSize();

   for (size_t i = 0; i < n; i++) {
      A.GetRawDataPointer()[i] = 0.0;
   }
}

} // namespace DNN
} // namespace TMVA
