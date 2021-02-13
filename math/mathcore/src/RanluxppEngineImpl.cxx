// @(#)root/mathcore:$Id$
// Author: Jonas Hahnfeld 11/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class ROOT::Math::RanluxppEngine
Implementation of the RANLUX++ generator

RANLUX++ is an LCG equivalent of RANLUX using 576 bit numbers.

Described in
A. Sibidanov, *A revision of the subtract-with-borrow random numbergenerators*,
*Computer Physics Communications*, 221(2017), 299-303,
preprint https://arxiv.org/pdf/1705.03123.pdf

The code is loosely based on the Assembly implementation by A. Sibidanov
available at https://github.com/sibidanov/ranluxpp/.
*/

#include "Math/RanluxppEngine.h"

#include "mulmod.h"

#include <cassert>
#include <cstdint>

namespace {

// Variable templates are a feature of C++14, use the older technique of having
// a static member in a template class.

template <int p>
struct RanluxppData;

template <>
struct RanluxppData<24> {
   static const uint64_t kA[9];
};
const uint64_t RanluxppData<24>::kA[] = {
   0x0000000000000000, 0x0000000000000000, 0x0000000000010000, 0xfffe000000000000, 0xffffffffffffffff,
   0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffeffffffff, 0xffffffffffffffff,
};

template <>
struct RanluxppData<2048> {
   static const uint64_t kA[9];
};
const uint64_t RanluxppData<2048>::kA[] = {
   0xed7faa90747aaad9, 0x4cec2c78af55c101, 0xe64dcb31c48228ec, 0x6d8a15a13bee7cb0, 0x20b2ca60cb78c509,
   0x256c3d3c662ea36c, 0xff74e54107684ed2, 0x492edfcc0cc8e753, 0xb48c187cf5b22097,
};

} // end anonymous namespace

namespace ROOT {
namespace Math {

template <int w, int p>
class RanluxppEngineImpl {

private:
   uint64_t fState[9]; ///< State of the generator
   int fPosition = 0;  ///< Current position in bits

   static constexpr const uint64_t *kA = RanluxppData<p>::kA;
   static constexpr int kMaxPos = 9 * 64;

   /// Produce next block of random bits
   void Advance()
   {
      mulmod(kA, fState);
      fPosition = 0;
   }

public:
   /// Return the next random bits, generate a new block if necessary
   uint64_t NextRandomBits()
   {
      if (fPosition + w > kMaxPos) {
         Advance();
      }

      int idx = fPosition / 64;
      int offset = fPosition % 64;
      int numBits = 64 - offset;

      uint64_t bits = fState[idx] >> offset;
      if (numBits < w) {
         bits |= fState[idx + 1] << numBits;
      }
      bits &= ((uint64_t(1) << w) - 1);

      fPosition += w;
      assert(fPosition <= kMaxPos && "position out of range!");

      return bits;
   }

   /// Initialize and seed the state of the generator
   void SetSeed(uint64_t s)
   {
      fState[0] = 1;
      for (int i = 1; i < 9; i++) {
         fState[i] = 0;
      }

      uint64_t a_seed[9];
      // Skip 2 ** 96 states.
      powermod(kA, a_seed, uint64_t(1) << 48);
      powermod(a_seed, a_seed, uint64_t(1) << 48);
      // Skip another s states.
      powermod(a_seed, a_seed, s);
      mulmod(a_seed, fState);

      fPosition = 0;
   }

   /// Skip `n` random numbers without generating them
   void Skip(uint64_t n)
   {
      int left = (kMaxPos - fPosition) / w;
      assert(left >= 0 && "position was out of range!");
      if (n < (uint64_t)left) {
         // Just skip the next few entries in the currently available bits.
         fPosition += n * w;
         assert(fPosition <= kMaxPos && "position out of range!");
         return;
      }

      n -= left;
      // Need to advance and possibly skip over blocks.
      int nPerState = kMaxPos / w;
      int skip = (n / nPerState);

      uint64_t a_skip[9];
      powermod(kA, a_skip, skip + 1);
      mulmod(a_skip, fState);

      // Potentially skip numbers in the freshly generated block.
      int remaining = n - skip * nPerState;
      assert(remaining >= 0 && "should not end up at a negative position!");
      fPosition = remaining * w;
      assert(fPosition <= kMaxPos && "position out of range!");
   }
};

template <int p>
RanluxppEngine<p>::RanluxppEngine(uint64_t seed) : fImpl(new RanluxppEngineImpl<52, p>)
{
   fImpl->SetSeed(seed);
}

template <int p>
RanluxppEngine<p>::~RanluxppEngine() = default;

template <int p>
double RanluxppEngine<p>::Rndm()
{
   return (*this)();
}

template <int p>
double RanluxppEngine<p>::operator()()
{
   // Get 52 bits of randomness.
   uint64_t bits = fImpl->NextRandomBits();

   // Construct the double in [1, 2), using the random bits as mantissa.
   static constexpr uint64_t exp = 0x3ff0000000000000;
   union {
      double dRandom;
      uint64_t iRandom;
   };
   iRandom = exp | bits;

   // Shift to the right interval of [0, 1).
   return dRandom - 1;
}

template <int p>
uint64_t RanluxppEngine<p>::IntRndm()
{
   return fImpl->NextRandomBits();
}

template <int p>
void RanluxppEngine<p>::SetSeed(uint64_t seed)
{
   fImpl->SetSeed(seed);
}

template <int p>
void RanluxppEngine<p>::Skip(uint64_t n)
{
   fImpl->Skip(n);
}

template class RanluxppEngine<24>;
template class RanluxppEngine<2048>;

} // end namespace Math
} // end namespace ROOT
