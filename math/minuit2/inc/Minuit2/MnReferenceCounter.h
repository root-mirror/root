// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnReferenceCounter
#define ROOT_Minuit2_MnReferenceCounter

#include <cassert>

#include "StackAllocator.h"

namespace ROOT {

namespace Minuit2 {

// extern StackAllocator gStackAllocator;

class MnReferenceCounter {

public:
   MnReferenceCounter() : fReferences(0) {}

   MnReferenceCounter(const MnReferenceCounter &other) : fReferences(other.fReferences) {}

   MnReferenceCounter &operator=(const MnReferenceCounter &other)
   {
      fReferences = other.fReferences;
      return *this;
   }

   ~MnReferenceCounter() { assert(fReferences == 0); }

   void *operator new(size_t nbytes) { return StackAllocatorHolder::Get().Allocate(nbytes); }

   void operator delete(void *p, size_t /*nbytes */) { StackAllocatorHolder::Get().Deallocate(p); }

   unsigned int References() const { return fReferences; }

   void AddReference() const { fReferences++; }

   void RemoveReference() const { fReferences--; }

private:
   mutable unsigned int fReferences;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnReferenceCounter
