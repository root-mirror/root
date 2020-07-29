/// \file ROOT/RHistData.h
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-06-14
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistData
#define ROOT7_RHistData

#include <cmath>
#include <vector>
#include "ROOT/RSpan.hxx"
#include "ROOT/RHistUtils.hxx"

namespace ROOT {
namespace Experimental {

namespace Internal {

/**
 Empty implementation not keeping track of the Poisson uncertainty per bin.
 */
template <int Dimensions, class WeightType, template <class EL> class Container, bool>
class RHistStatUncertainty {
public:
   RHistStatUncertainty() = default;
   RHistStatUncertainty(size_t /*bin_size*/, size_t /*overflow_size*/) {}

   WeightType GetBinArray(int binidx) const { return  0.; }
   void Fill(int /*binidx*/, WeightType /*weight*/ = 1.) const {}

   /// Merge with other `RHistStatUncertainty` data, assuming same bin configuration.
   void Add(const RHistStatUncertainty& other) const {}
};

/**
 \class RHistStatUncertainty
 Histogram statistics to keep track of the Poisson uncertainty per bin.
 */
template <int Dimensions, class WeightType, template <class EL> class Container>
class RHistStatUncertainty<Dimensions, WeightType, Container, true> {
public:
   /// Type of the bin content array.
   using Content_t = Container<WeightType>;

private:
   /// Uncertainty of the content for each bin excluding under-/overflow.
   Content_t fSumWeightsSquared; ///< Sum of squared weights.
   /// Uncertainty of the under-/overflow content.
   Content_t fOverflowSumWeightsSquared; ///< Sum of squared weights for under-/overflow.

public:
   RHistStatUncertainty() = default;
   RHistStatUncertainty(size_t bin_size, size_t overflow_size): fSumWeightsSquared(bin_size), fOverflowSumWeightsSquared(overflow_size) {}

   /// Get a reference to the bin corresponding to `binidx` of the correct bin
   /// content array 
   /// i.e. depending if `binidx` is a regular bin or an under- / overflow bin.
   WeightType GetBinArray(int binidx) const
   {
      if (binidx < 0){
         return fOverflowSumWeightsSquared[-binidx - 1];
      } else {
         return fSumWeightsSquared[binidx - 1];
      }
   }

   /// Get a reference to the bin corresponding to `binidx` of the correct bin
   /// content array (non-const)
   /// i.e. depending if `binidx` is a regular bin or an under- / overflow bin.
   WeightType& GetBinArray(int binidx)
   {
      if (binidx < 0){
         return fOverflowSumWeightsSquared[-binidx - 1];
      } else {
         return fSumWeightsSquared[binidx - 1];
      }
   }

   /// Add weight to the bin at `binidx`; the coordinate was `x`.
   void Fill(int binidx, WeightType weight = 1.)
   {
      GetBinArray(binidx) += weight * weight;
   }

   /// Calculate a bin's (Poisson) uncertainty of the bin content as the
   /// square-root of the bin's sum of squared weights.
   double GetBinUncertaintyImpl(int binidx) const { return std::sqrt(GetBinArray(binidx)); }

   /// Get a bin's sum of squared weights.
   WeightType GetSumOfSquaredWeights(int binidx) const { return GetBinArray(binidx); }
   /// Get a bin's sum of squared weights.
   WeightType &GetSumOfSquaredWeights(int binidx) { return GetBinArray(binidx); }

   /// Get the structure holding the sum of squares of weights.
   const Content_t &GetSumOfSquaredWeights() const { return fSumWeightsSquared; }
   /// Get the structure holding the sum of squares of weights (non-const).
   Content_t &GetSumOfSquaredWeights() { return fSumWeightsSquared; }

   /// Get the structure holding the under-/overflow sum of squares of weights.
   const Content_t &GetOverflowSumOfSquaredWeights() const { return fOverflowSumWeightsSquared; }
   /// Get the structure holding the under-/overflow sum of squares of weights (non-const).
   Content_t &GetOverflowSumOfSquaredWeights() { return fOverflowSumWeightsSquared; }

   /// Merge with other `RHistStatUncertainty` data, assuming same bin configuration.
   void Add(const RHistStatUncertainty& other) {
      assert(fSumWeightsSquared.size() == other.fSumWeightsSquared.size()
               && "this and other have incompatible bin configuration!");
      assert(fOverflowSumWeightsSquared.size() == other.fOverflowSumWeightsSquared.size()
               && "this and other have incompatible bin configuration!");
      for (size_t b = 0; b < fSumWeightsSquared.size(); ++b)
         fSumWeightsSquared[b] += other.fSumWeightsSquared[b];
      for (size_t b = 0; b < fOverflowSumWeightsSquared.size(); ++b)
         fOverflowSumWeightsSquared[b] += other.fOverflowSumWeightsSquared[b];
   }
};

/**
 Empty implementation not keeping track of moments.
 */
template <int Dimensions, class WeightType, bool>
class RHistDataMomentUncert {
   RHistDataMomentUncert() = default;
   RHistDataMomentUncert(size_t, size_t) {}

   void FillMoment(const std::array<double, Dimensions> &/*x*/, WeightType /*weight*/ = 1.) const {}

   void Add(const RHistDataMomentUncert& other) const {}
};

/** \class RHistDataMomentUncert
  For now do as `RH1`: calculate first (xw) and second (x^2w) moment.
*/
template <int Dimensions, class WeightType>
class RHistDataMomentUncert<Dimensions, WeightType, true> {
public:
   /// Type of the moments array.
   using MomentArr_t = std::array<WeightType, Dimensions>;

private:
   MomentArr_t fMomentXW;
   MomentArr_t fMomentX2W;
   // FIXME: Add sum(w.x.y)-style stats.

public:
   RHistDataMomentUncert() = default;
   RHistDataMomentUncert(size_t, size_t) {}

   /// Add weight to the bin at binidx; the coordinate was x.
   void FillMoment(const std::array<double, Dimensions> &x, WeightType weight = 1.)
   {
      for (int idim = 0; idim < Dimensions; ++idim) {
         const WeightType xw = x[idim] * weight;
         fMomentXW[idim] += xw;
         fMomentX2W[idim] += x[idim] * xw;
      }
   }

   // FIXME: Add a way to query the inner data

   /// Merge with other RHistDataMomentUncert data, assuming same bin configuration.
   void Add(const RHistDataMomentUncert& other) {
      for (size_t d = 0; d < Dimensions; ++d) {
         fMomentXW[d] += other.fMomentXW[d];
         fMomentX2W[d] += other.fMomentX2W[d];
      }
   }
};
} // namespace Internal

namespace Detail {

/** \class RHistData
  A `RHistImplBase`'s data, provides accessors to all its statistics.
  */
template <int Dimensions, class WeightType, int StatConfig, template <class EL> class Container>
class RHistData:
public Internal::RHistStatUncertainty<Dimensions, WeightType, Container,
   (bool)(StatConfig & Hist::Stat::kUncertainty)>,
Internal::RHistDataMomentUncert<Dimensions, WeightType,
   (bool)(StatConfig & (Hist::Stat::k1stMoment | Hist::Stat::k2ndMoment))>
{
public:
   /// Type of the bin content array.
   using Content_t = Container<WeightType>;
   using HistStatUncertainty_t = Internal::RHistStatUncertainty<Dimensions, WeightType, Container,
      (bool)(StatConfig & (int)Hist::Stat::kUncertainty)>;
   using HistDataMomentUncert_t = Internal::RHistDataMomentUncert<Dimensions, WeightType,
      (bool)(StatConfig & ((int)Hist::Stat::k1stMoment | (int)Hist::Stat::k2ndMoment))>;
private:
   /// Number of calls to Fill().
   ssize_t fEntries = 0;

   /// Sum of weights.
   WeightType fSumWeights = 0;

   /// Sum of (weights^2).
   WeightType fSumWeights2 = 0;

   /// Bin content.
   Content_t fBinContent;

   /// Under- and overflow bin content.
   Content_t fOverflowBinContent;

public:
   RHistData() = default;
   RHistData(size_t bin_size, size_t overflow_size):
      HistStatUncertainty_t(bin_size, overflow_size),
      HistDataMomentUncert_t(bin_size, overflow_size),
    fBinContent(bin_size), fOverflowBinContent(overflow_size) {}

   /// Get a reference to the bin corresponding to `binidx` of the correct bin
   /// content array 
   /// i.e. depending if `binidx` is a regular bin or an under- / overflow bin.
   WeightType GetBinArray(int binidx) const
   {
      if (binidx < 0){
         return fOverflowBinContent[-binidx - 1];
      } else {
         return fBinContent[binidx - 1];
      }
   }

   /// Get a reference to the bin corresponding to `binidx` of the correct bin
   /// content array (non-const)
   /// i.e. depending if `binidx` is a regular bin or an under- / overflow bin.
   WeightType& GetBinArray(int binidx)
   {
      if (binidx < 0){
         return fOverflowBinContent[-binidx - 1];
      } else {
         return fBinContent[binidx - 1];
      }
   }

   /// Add weight to the bin content at `binidx`.
   void Fill(int binidx, WeightType weight = 1.)
   {
      GetBinArray(binidx) += weight;
      ++fEntries;
   }

   /// Get the number of entries filled into the histogram - i.e. the number of
   /// calls to Fill().
   int64_t GetEntries() const { return fEntries; }

   /// Get the number of bins exluding under- and overflow.
   size_t sizeNoOver() const noexcept { return fBinContent.size(); }

   /// Get the number of bins including under- and overflow..
   size_t size() const noexcept { return fBinContent.size() + fOverflowBinContent.size(); }

   /// Get the number of bins including under- and overflow..
   size_t sizeUnderOver() const noexcept { return fOverflowBinContent.size(); }

   /// Get the bin content for the given bin.
   WeightType operator[](int binidx) const { return GetBinArray(binidx); }
   /// Get the bin content for the given bin (non-const).
   WeightType &operator[](int binidx) { return GetBinArray(binidx); }

   /// Get the bin content for the given bin.
   WeightType GetBinContent(int binidx) const { return GetBinArray(binidx); }
   /// Get the bin content for the given bin (non-const).
   WeightType &GetBinContent(int binidx) { return GetBinArray(binidx); }

   /// Retrieve the content array.
   const Content_t &GetContentArray() const { return fBinContent; }
   /// Retrieve the content array (non-const).
   Content_t &GetContentArray() { return fBinContent; }

   /// Retrieve the under-/overflow content array.
   const Content_t &GetOverflowContentArray() const { return fOverflowBinContent; }
   /// Retrieve the under-/overflow content array (non-const).
   Content_t &GetOverflowContentArray() { return fOverflowBinContent; }

   /// Merge with other RHistStatContent, assuming same bin configuration.
   void Add(const RHistData& other) {
      assert(fBinContent.size() == other.fBinContent.size()
               && "this and other have incompatible bin configuration!");
      assert(fOverflowBinContent.size() == other.fOverflowBinContent.size()
               && "this and other have incompatible bin configuration!");
      fEntries += other.fEntries;
      for (size_t b = 0; b < fBinContent.size(); ++b)
         fBinContent[b] += other.fBinContent[b];
      for (size_t b = 0; b < fOverflowBinContent.size(); ++b)
         fOverflowBinContent[b] += other.fOverflowBinContent[b];
      HistStatUncertainty_t::Add(other);
      HistDataMomentUncert_t::AdD(other);
   }

};
} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
