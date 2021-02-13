// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_StackAllocator
#define ROOT_Minuit2_StackAllocator

#include "Minuit2/MnConfig.h"

// comment out this line and recompile if you want to gain additional
// performance (the gain is mainly for "simple" functions which are easy
// to calculate and vanishes quickly if going to cost-intensive functions)
// the library is no longer thread save however

#ifdef MN_USE_STACK_ALLOC
#define _MN_NO_THREAD_SAVE_
#endif

#include <cstdlib>
#include <new>

namespace ROOT {

namespace Minuit2 {

/// define stack allocator symbol

class StackOverflow {
};
class StackError {
};
//  using namespace std;

/** StackAllocator controls the memory allocation/deallocation of Minuit. If
    _MN_NO_THREAD_SAVE_ is defined, memory is taken from a pre-allocated piece
    of heap memory which is then used like a stack, otherwise via standard
    malloc/free. Note that defining _MN_NO_THREAD_SAVE_ makes the code thread-
    unsave. The gain in performance is mainly for cost-cheap FCN functions.
 */

class StackAllocator {

public:
   //   enum {default_size = 1048576};
   enum { default_size = 524288 };

   StackAllocator() : fStack(0)
   {
#ifdef _MN_NO_THREAD_SAVE_
      // std::cout<<"StackAllocator Allocate "<<default_size<<std::endl;
      fStack = new unsigned char[default_size];
#endif
      fStackOffset = 0;
      fBlockCount = 0;
   }

   ~StackAllocator()
   {
#ifdef _MN_NO_THREAD_SAVE_
      // std::cout<<"StackAllocator destruct "<<fStackOffset<<std::endl;
      if (fStack)
         delete[] fStack;
#endif
   }

   void *Allocate(size_t nBytes)
   {
#ifdef _MN_NO_THREAD_SAVE_
      if (fStack == 0)
         fStack = new unsigned char[default_size];
      int nAlloc = AlignedSize(nBytes);
      CheckOverflow(nAlloc);

      //       std::cout << "Allocating " << nAlloc << " bytes, requested = " << nBytes << std::endl;

      // write the start position of the next block at the start of the block
      WriteInt(fStackOffset, fStackOffset + nAlloc);
      // write the start position of the new block at the end of the block
      WriteInt(fStackOffset + nAlloc - sizeof(int), fStackOffset);

      void *result = fStack + fStackOffset + sizeof(int);
      fStackOffset += nAlloc;
      fBlockCount++;

#ifdef DEBUG_ALLOCATOR
      CheckConsistency();
#endif

#else
      void *result = malloc(nBytes);
      if (!result)
         throw std::bad_alloc();
#endif

      return result;
   }

   void Deallocate(void *p)
   {
#ifdef _MN_NO_THREAD_SAVE_
      // int previousOffset = ReadInt( fStackOffset - sizeof(int));
      int delBlock = ToInt(p);
      int nextBlock = ReadInt(delBlock);
      int previousBlock = ReadInt(nextBlock - sizeof(int));
      if (nextBlock == fStackOffset) {
         // deallocating last allocated
         fStackOffset = previousBlock;
      } else {
         // overwrite previous adr of next block
         int nextNextBlock = ReadInt(nextBlock);
         WriteInt(nextNextBlock - sizeof(int), previousBlock);
         // overwrite head of deleted block
         WriteInt(previousBlock, nextNextBlock);
      }
      fBlockCount--;

#ifdef DEBUG_ALLOCATOR
      CheckConsistency();
#endif
#else
      free(p);
#endif
      // std::cout << "Block at " << delBlock
      //   << " deallocated, fStackOffset = " << fStackOffset << std::endl;
   }

   int ReadInt(int offset)
   {
      int *ip = (int *)(fStack + offset);

      // std::cout << "read " << *ip << " from offset " << offset << std::endl;

      return *ip;
   }

   void WriteInt(int offset, int Value)
   {

      // std::cout << "writing " << Value << " to offset " << offset << std::endl;

      int *ip = reinterpret_cast<int *>(fStack + offset);
      *ip = Value;
   }

   int ToInt(void *p)
   {
      unsigned char *pc = static_cast<unsigned char *>(p);

      // std::cout << "toInt: p = " << p << " fStack = " << (void*) fStack << std::endl;
      // VC 7.1 warning:conversin from __w64 int to int
      int userBlock = pc - fStack;
      return userBlock - sizeof(int); // correct for starting int
   }

   int AlignedSize(int nBytes)
   {
      const int fAlignment = 4;
      int needed = nBytes % fAlignment == 0 ? nBytes : (nBytes / fAlignment + 1) * fAlignment;
      return needed + 2 * sizeof(int);
   }

   void CheckOverflow(int n)
   {
      if (fStackOffset + n >= default_size) {
         // std::cout << " no more space on stack allocator" << std::endl;
         throw StackOverflow();
      }
   }

   bool CheckConsistency()
   {

      // std::cout << "checking consistency for " << fBlockCount << " blocks"<< std::endl;

      // loop over all blocks
      int beg = 0;
      int end = fStackOffset;
      int nblocks = 0;
      while (beg < fStackOffset) {
         end = ReadInt(beg);

         // std::cout << "beg = " << beg << " end = " << end
         //     << " fStackOffset = " << fStackOffset << std::endl;

         int beg2 = ReadInt(end - sizeof(int));
         if (beg != beg2) {
            // std::cout << "  beg != beg2 " << std::endl;
            return false;
         }
         nblocks++;
         beg = end;
      }
      if (end != fStackOffset) {
         // std::cout << " end != fStackOffset" << std::endl;
         return false;
      }
      if (nblocks != fBlockCount) {
         // std::cout << "nblocks != fBlockCount" << std::endl;
         return false;
      }
      // std::cout << "Allocator is in consistent state, nblocks = " << nblocks << std::endl;
      return true;
   }

private:
   unsigned char *fStack;
   //   unsigned char fStack[default_size];
   int fStackOffset;
   int fBlockCount;
};

class StackAllocatorHolder {

   // t.b.d need to use same trick as  Boost singleton.hpp to be sure that
   // StackAllocator is created before main()

public:
   static StackAllocator &Get()
   {
      static StackAllocator gStackAllocator;
      return gStackAllocator;
   }
};

} // namespace Minuit2

} // namespace ROOT

#endif
