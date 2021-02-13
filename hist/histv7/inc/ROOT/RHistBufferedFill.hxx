/// \file ROOT/RHistBufferedFill.hxx
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-03
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistBufferedFill
#define ROOT7_RHistBufferedFill

#include "ROOT/RSpan.hxx"

namespace ROOT {
namespace Experimental {

namespace Internal {
template <class DERIVED, class HIST, int SIZE>
class RHistBufferedFillBase {
public:
   using CoordArray_t = typename HIST::CoordArray_t;
   using Weight_t = typename HIST::Weight_t;

private:
   size_t fCursor = 0;
   std::array<CoordArray_t, SIZE> fXBuf;
   std::array<Weight_t, SIZE> fWBuf;

public:
   RHistBufferedFillBase() {}
   ~RHistBufferedFillBase() { toDerived().Flush(); }

   DERIVED &toDerived() { return *static_cast<DERIVED *>(this); }
   const DERIVED &toDerived() const { return *static_cast<const DERIVED *>(this); }

   std::span<const CoordArray_t> GetCoords() const
   {
      return std::span<const CoordArray_t>(fXBuf.begin(), fXBuf.begin() + fCursor);
   }
   std::span<const Weight_t> GetWeights() const
   {
      return std::span<const Weight_t>(fWBuf.begin(), fWBuf.begin() + fCursor);
   }

   void Fill(const CoordArray_t &x, Weight_t weight = 1.)
   {
      fXBuf[fCursor] = x;
      fWBuf[fCursor++] = weight;
      if (fCursor == SIZE) {
         Flush();
      }
   }

   void Flush() {
      toDerived().FlushImpl();
      fCursor = 0;
   }
};

} // namespace Internal

/** \class RHistBufferedFill
 Buffers calls to Fill().

 Once the buffer is full, on destruction of when calling Flush(), it sends the
 buffers off as an ideally vectorizable FillN() operation. It also serves as a
 multi-threaded way of filling the same histogram, reducing the locking
 frequency.

 The HIST template can be either a RHist instance, a RHistImpl instance, or
 a RHistLockedFill instance.
 **/

template <class HIST, int SIZE = 1024>
class RHistBufferedFill: public Internal::RHistBufferedFillBase<RHistBufferedFill<HIST, SIZE>, HIST, SIZE> {
public:
   using Hist_t = HIST;
   using CoordArray_t = typename HIST::CoordArray_t;
   using Weight_t = typename HIST::Weight_t;

private:
   HIST &fHist;

   friend class Internal::RHistBufferedFillBase<RHistBufferedFill<HIST, SIZE>, HIST, SIZE>;
   void FlushImpl() { fHist.FillN(this->GetCoords(), this->GetWeights()); }

public:
   RHistBufferedFill(Hist_t &hist): fHist{hist} {}

   void FillN(const std::span<const CoordArray_t> xN, const std::span<const Weight_t> weightN)
   {
      fHist.FillN(xN, weightN);
   }

   void FillN(const std::span<const CoordArray_t> xN) { fHist.FillN(xN); }

   HIST &GetHist()
   {
      this->Flush(); // synchronize!
      return fHist;
   }
   operator HIST &() { return GetHist(); }

   static constexpr int GetNDim() { return HIST::GetNDim(); }
};
} // namespace Experimental
} // namespace ROOT

#endif
