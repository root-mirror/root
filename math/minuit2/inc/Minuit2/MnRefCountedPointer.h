// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnRefCountedPointer
#define ROOT_Minuit2_MnRefCountedPointer

#include "MnReferenceCounter.h"

namespace ROOT {

namespace Minuit2 {

template <class T>
class MnRefCountedPointer {

public:
   // Default constructor needed for use inside array, vector, etc.
   MnRefCountedPointer() : fPtr(0), fCounter(0) {}

   MnRefCountedPointer(T *pt) : fPtr(pt), fCounter(new MnReferenceCounter()) { AddReference(); }

   MnRefCountedPointer(const MnRefCountedPointer<T> &other) : fPtr(other.fPtr), fCounter(other.fCounter)
   {
      AddReference();
   }

   ~MnRefCountedPointer()
   {
      /*
      if(References() == 0) {
        if(fPtr) delete fPtr;
        if(fCounter) delete fCounter;
      }
      else RemoveReference();
      */
      if (References() != 0)
         RemoveReference();
   }

   bool IsValid() const { return fPtr != 0; }

   MnRefCountedPointer &operator=(const MnRefCountedPointer<T> &other)
   {
      if (this != &other && fPtr != other.fPtr) {
         RemoveReference();
         fPtr = other.fPtr;
         fCounter = other.fCounter;
         AddReference();
      }
      return *this;
   }

   MnRefCountedPointer &operator=(T *ptr)
   {
      if (fPtr != ptr) {
         fPtr = ptr;
         fCounter = new MnReferenceCounter();
      }
      return *this;
   }

   T *Get() const { return fPtr; }

   T *operator->() const
   {
      DoCheck();
      return fPtr;
   }

   T &operator*() const
   {
      DoCheck();
      return *fPtr;
   }

   bool operator==(const T *otherP) const { return fPtr == otherP; }

   bool operator<(const T *otherP) const { return fPtr < otherP; }

   unsigned int References() const { return fCounter->References(); }

   void AddReference() const { fCounter->AddReference(); }

   void RemoveReference()
   {
      fCounter->RemoveReference();
      if (References() == 0) {
         delete fPtr;
         fPtr = 0;
         delete fCounter;
         fCounter = 0;
      }
   }

private:
   T *fPtr;
   MnReferenceCounter *fCounter;

private:
   void DoCheck() const { assert(IsValid()); }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnRefCountedPointer
