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
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();
 
   AFloat sigma = sqrt(2.0 / ((AFloat) n));

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         A(i,j) = rand.Gaus(0.0, sigma);
      }
   }
}

//______________________________________________________________________________
template<typename AFloat>
void TCpu<AFloat>::InitializeUniform(TCpuMatrix<AFloat> & A)
{
   size_t m,n;
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat range = sqrt(2.0 / ((AFloat) n));

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         A(i,j) = rand.Uniform(-range, range);
      }
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
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat sigma = sqrt(2.0 /( ((AFloat) n) + ((AFloat) m)) );

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         AFloat value = rand.Gaus(0.0, sigma);
         if ( std::abs(value) > 2*sigma) continue; 
         A(i,j) = rand.Gaus(0.0, sigma);
      }
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
   m = A.GetNrows();
   n = A.GetNcols();

   TRandom &  rand = GetRandomGenerator();

   AFloat range = sqrt(6.0 /( ((AFloat) n) + ((AFloat) m)) );

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         A(i,j) = rand.Uniform(-range, range);
      }
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
      for (size_t j = 0; j < n ; j++) {
         A(i,j) = 0.0;
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

} // namespace DNN
} // namespace TMVA
