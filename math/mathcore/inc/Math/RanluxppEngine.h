// @(#)root/mathcore:$Id$
// Author: Jonas Hahnfeld 11/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Math_RanluxppEngine
#define ROOT_Math_RanluxppEngine

#include "Math/TRandomEngine.h"

#include <cstdint>
#include <memory>

namespace ROOT {
namespace Math {

template <int w, int p>
class RanluxppEngineImpl;

template <int p>
class RanluxppEngine final : public TRandomEngine {

private:
   std::unique_ptr<RanluxppEngineImpl<52, p>> fImpl;

public:
   RanluxppEngine(uint64_t seed = 314159265);
   virtual ~RanluxppEngine();

   /// Generate a double-precision random number with 52 bits of randomness
   double Rndm() override;
   /// Generate a double-precision random number (non-virtual method)
   double operator()();
   /// Generate a random integer value with 52 bits
   uint64_t IntRndm();

   /// Initialize and seed the state of the generator
   void SetSeed(uint64_t seed);
   /// Skip `n` random numbers without generating them
   void Skip(uint64_t n);

   /// Get name of the generator
   static const char *Name() { return "RANLUX++"; }
};

using RanluxppEngine24 = RanluxppEngine<24>;
using RanluxppEngine2048 = RanluxppEngine<2048>;

extern template class RanluxppEngine<24>;
extern template class RanluxppEngine<2048>;

} // end namespace Math
} // end namespace ROOT

#endif /* ROOT_Math_RanluxppEngine */
