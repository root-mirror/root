/// \file ROOT/RPage.hxx
/// \ingroup Forest ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPage
#define ROOT7_RPage

#include <ROOT/RForestUtil.hxx>

#include <cstddef>

namespace ROOT {
namespace Experimental {

namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RPage
\ingroup Forest
\brief A page is a fixed size slice of an column that is mapped into memory

The page provides a fixed-size opaque memory buffer for uncompressed data. It does not know how to interpret
the contents but it does now about the size (and thus the number) of the elements inside as well as the element
number range within the backing column. The memory buffer is not managed by the page but normally by the page pool.
*/
// clang-format on
class RPage {
   void* fBuffer;
   std::size_t fCapacity;
   std::size_t fSize;
   std::size_t fElementSize;
   TreeIndex_t fRangeStart;

public:
   RPage() : fBuffer(nullptr), fCapacity(0), fSize(0), fElementSize(0), fRangeStart(0) {}
   RPage(void* buffer, std::size_t capacity, std::size_t elementSize)
      : fBuffer(buffer), fCapacity(capacity), fSize(0), fElementSize(elementSize), fRangeStart(0) {}
   ~RPage() = default;

   /// The total space available in the page
   std::size_t GetCapacity() const { return fCapacity; }
   /// The space taken by column elements in the buffer
   std::size_t GetSize() const { return fSize; }
   TreeIndex_t GetNElements() const { return fSize / fElementSize; }
   TreeIndex_t GetRangeStart() const { return fRangeStart; }
   void* GetBuffer() const { return fBuffer; }
   /// Return a pointer after the last element that has space for nElements new elements. If there is not enough capacity,
   /// return nullptr
   void* Reserve(std::size_t nElements) {
      size_t offset = fSize;
      size_t nbyte = nElements * fElementSize;
      if (offset + nbyte > fCapacity) {
        return nullptr;
      }
      fSize += nbyte;
      return static_cast<unsigned char *>(fBuffer) + offset;
   }
   /// Forget all currently stored elements (size == 0) and set a new starting index.
   void Reset(TreeIndex_t rangeStart) { fSize = 0; fRangeStart = rangeStart; }

   bool IsNull() const { return fBuffer == nullptr; }
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
