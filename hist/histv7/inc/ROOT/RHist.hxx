/// \file ROOT/RHist.h
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-03-23
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHist
#define ROOT7_RHist

#include "ROOT/RSpan.hxx"
#include "ROOT/RAxis.hxx"
#include "ROOT/RAxisLayout.hxx"
#include "ROOT/RHistData.hxx"
#include "ROOT/RHistUtils.hxx"
#include <initializer_list>
#include <stdexcept>

namespace ROOT {
namespace Experimental {

/// Using declaration for the default container for bin content and uncertainties.
template <class El>
using StdVector_t = std::vector<El>;

/**
 \class RHistBase
 Histogram base class for histograms with `Dimensions` dimensions, where each
 bin count is stored by a value of type `WeightType`. STAT stores statistical
 data of the entries filled into the histogram (bin content, uncertainties etc).

 A histogram counts occurrences of values or n-dimensional combinations thereof.
 Contrary to for instance a `RTree`, a histogram combines adjacent values. The
 resolution of this combination is defined by the axis binning, see e.g.
 http://www.wikiwand.com/en/Histogram

 Histogram axes are only known to the derived classes.
 */

template <int Dimensions, class WeightType, int StatConfig, template <class EL> class Container = StdVector_t>
class RHistBase {
public:
   static constexpr int kNDim = Dimensions;
   using Weight_t = WeightType;
   static constexpr int kStatConfig = StatConfig;
   using Data_t = Detail::RHistData<kNDim, WeightType, kStatConfig, Container>;

private:
   Data_t fData;

protected:
   RHistBase() = default;
   RHistBase(size_t numBins);

   Data_t &GetData() { return fData; }

public:
   /// Number of dimensions of the coordinates
   static constexpr int GetNDim() noexcept { return kNDim; }

   std::array<RAxisBase*, kNDim> GetAxes() const = 0;

   WeightType GetBinContent(size_t bin) const { return fData[bin]; }
   WeightType GetStatContent(size_t bin, Hist::Stat::EStat stat) const { return fData.GetStat(stat)[bin]; }
   double_t GetUncertainty(size_t bin) const { return fData.GetCombinedUncertainty(bin); }
   size_t GetNBins() const { return fData.GetNBins(); }
   const Data_t &GetData() const { return fData; }
};

template <int Dimensions, class WeightType, class AxisTuple, int StatConfig = 0, template <class EL> class Container = StdVector_t >
class RHist: public RHistBase<Dimensions, WeightType, StatConfig, Container> {
   static_assert(! (std::tuple_size<AxisTuple>::value < Dimensions), "Too few histogram axes specified for a histogram of this dimension!");
   static_assert(! (std::tuple_size<AxisTuple>::value > Dimensions), "Too many histogram axes specified for a histogram of this dimension!");

public:
   using AxisLayout_t = Internal::RAxisLayout<AxisTuple>;
   using CoordTuple_t = typename AxisLayout_t::CoordTuple_t;
   static constexpr int kNDim = Dimensions;
   static constexpr int kStatConfig = StatConfig;

private:
   using Base_t = RHistBase<Dimensions, WeightType, StatConfig, Container>;

   AxisLayout_t fAxisLayout;

   ///\{
   /// Internal interfaces to convert a coord tuple to an array of doubles.
   template <class...Axis>
   struct TupleElementsConvertibleToDouble:
      std::integral_constant<bool, (std::is_convertible<double, Axis>::value && ...)>
   {};

   std::array<double, Dimensions> CoordTupleToDoubleArr(...) { return {}; } // FIXME: more helpful diags for users!

   template <class...Coord, class T = std::enable_if_t<TupleElementsConvertibleToDouble<Coord...>::value>>
   std::array<double, Dimensions> CoordTupleToDoubleArr(const Coord &...coord) {
      return {((double)coord),...};
   }
   ///\}

public:
   /// Create a histogram from an `array` of axes (`RAxisConfig`s). Example code:
   ///
   /// Construct a 1-dimensional histogram that can be filled with `floats`s.
   /// The axis has 10 bins between 0. and 1. The two outermost sets of curly
   /// braces are to reach the initialization of the `std::array` elements; the
   /// inner one is for the initialization of a `RAxisCoordinate`.
   ///
   ///     RHist<1,float> h1f({{ {10, 0., 1.} }});
   ///
   /// Construct a 2-dimensional histogram, with the first axis as before, and
   /// the second axis having non-uniform ("irregular") binning, where all bin-
   /// edges are specified. As this is itself an array it must be enclosed by
   /// double curlies.
   ///
   ///     RHist<2,int> h2i({{ {10, 0., 1.}, {{-1., 0., 1., 10., 100.}} }});
   explicit RHist(const AxisTuple &axes);

   /// Constructor taking the single axis for a 1D histogram.
   template <int EnableIfNDim = Dimensions, class = typename std::enable_if<EnableIfNDim == 1>::type>
   explicit RHist(const std::tuple_element<0, AxisTuple> &axis): RHist(std::forward_as_tuple(axis)) {}

   /// Constructor taking the two axes for a 2D histogram.
   template <int EnableIfNDim = Dimensions, class = typename std::enable_if<EnableIfNDim == 2>::type>
   explicit RHist(const std::tuple_element<0, AxisTuple> &axis0,
      const std::tuple_element<1, AxisTuple> &axis1):
   RHist(std::forward_as_tuple(axis0, axis1)) {}

   /// Add `weight` to the bin containing coordinate `x`.
   void Fill(const CoordTuple_t &x, WeightType weight = (WeightType)1) noexcept {
      // FIXME: handle growth!
      ssize_t bin = fAxisLayout.FindBin(x);
      if /*constexpr*/ (StatConfig >= Hist::Stat::k1stMoment)
         Base_t::GetData().FillMoment(bin, std::apply(CoordTupleToDoubleArr, x), weight);
      Base_t::GetData().Fill(bin, weight);
   }

   /// Add `weight` to the bin containing coordinate `x`, non-tuple overload for 1D histogram.
   template <int EnableIfNDim = Dimensions, class = typename std::enable_if<EnableIfNDim == 1>::type>
   void Fill(const typename std::tuple_element<0, AxisTuple>::Coord_t &x, WeightType weight = (WeightType)1) {
      Fill(CoordTuple_t{x}, weight);
   }

   /// Add `weight` to the bin containing coordinate `x`, non-tuple overload for 2D histogram.
   template <int EnableIfNDim = Dimensions, class = typename std::enable_if<EnableIfNDim == 2>::type>
   void Fill(const typename std::tuple_element<0, AxisTuple>::Coord_t &x,
             const typename std::tuple_element<1, AxisTuple>::Coord_t &y,
             WeightType weight = (WeightType)1) {
      Fill(CoordTuple_t{x, y}, weight);
   }

   /// For each coordinate in `xN`, add `weightN[i]` to the bin at coordinate
   /// `xN[i]`. The sizes of `xN` and `weightN` must be the same. This is more
   /// efficient than many separate calls to `Fill()`.
   void FillN(const std::span<const CoordTuple_t> xN, const std::span<const WeightType> weightN) noexcept
   {
      std::vector<ssize_t> binN = fAxisLayout.FindBin(xN);
      if /*constexpr*/ (StatConfig >= Hist::Stat::k1stMoment) {
         std::vector<std::array<double, Dimensions>> coordValN;
         coordValN.reserve(xN.size());
         for (auto &&x: xN)
            coordValN.emplace_back(std::apply(CoordTupleToDoubleArr, x));
         Base_t::GetData().FillMoment(binN, coordValN, weightN);
      }
      Base_t::GetData().Fill(binN, weightN);
   }

   /// Convenience overload: `FillN()` with weight 1.
   void FillN(const std::span<const CoordTuple_t> xN) noexcept {
      std::vector<size_t> binN = fAxisLayout.FindBin(xN);
      if /*constexpr*/ (StatConfig >= Hist::Stat::k1stMoment) {
         std::vector<std::array<double, Dimensions>> coordValN;
         coordValN.reserve(xN.size());
         for (auto &&x: xN)
            coordValN.emplace_back(std::apply(CoordTupleToDoubleArr, x));
         Base_t::GetData().FillMoment(binN, coordValN);
      }
      Base_t::GetData().Fill(binN);
   }

   /// Get the content of the bin at `x`.
   WeightType GetBinContent(const CoordTuple_t &x) const {
      return Base_t::GetBinContent(fAxisLayout.FindBin(x));
   }

   /// Get the uncertainty on the content of the bin at `x`.
   double GetBinUncertainty(const CoordTuple_t &x) const {
      return Base_t::GetBinUncertainty(fAxisLayout.FindBin(x));
   }

   /// TODO:
   /*
   const_iterator begin() const { return const_iterator(*fImpl); }

   const_iterator end() const { return const_iterator(*fImpl, fImpl->GetNBinsNoOver()); }

   /// Swap *this and other.
   ///
   /// Very efficient; swaps the `fImpl` pointers.
   void swap(RHist<Dimensions, WeightType, STAT...> &other) noexcept
   {
      std::swap(fImpl, other.fImpl);
      std::swap(fFillFunc, other.fFillFunc);
   }
   */
};

/// \name RHist Typedefs
///\{ Convenience typedefs (ROOT6-compatible type names)

// Keep them as typedefs, to make sure old-style documentation tools can understand them.
/// TODO: Figure out how to provide aliases given the axis tuple - maybe `RH1DBase` and `template <class Axis> RH1D`?
/*
using RH1D = RHist<1, double, Hist::Stat::kUncertainty>;
using RH1F = RHist<1, float, Hist::Stat::kUncertainty>;
using RH1C = RHist<1, char>;
using RH1I = RHist<1, int>;
using RH1LL = RHist<1, int64_t>;

using RH2D = RHist<2, double, Hist::Stat::kUncertainty>;
using RH2F = RHist<2, float, Hist::Stat::kUncertainty>;
using RH2C = RHist<2, char>;
using RH2I = RHist<2, int>;
using RH2LL = RHist<2, int64_t>;

using RH3D = RHist<3, double, Hist::Stat::kUncertainty>;
using RH3F = RHist<3, float, Hist::Stat::kUncertainty>;
using RH3C = RHist<3, char>;
using RH3I = RHist<3, int>;
using RH3LL = RHist<3, int64_t>;
///\}
*/

/// TODO:
/*
/// Add two histograms.
///
/// This operation may currently only be performed if the two histograms have
/// the same axis configuration, use the same WeightType, and if `from` records
/// at least the same statistics as `to` (recording more stats is fine).
///
/// Adding histograms with incompatible axis binning will be reported at runtime
/// with an `std::runtime_error`. Insufficient statistics in the source
/// histogram will be detected at compile-time and result in a compiler error.
///
/// In the future, we may either adopt a more relaxed definition of histogram
/// addition or provide a mechanism to convert from one histogram type to
/// another. We currently favor the latter path.
template <int Dimensions, class WeightType,
          template <int D_, class P_> class... STAT_TO,
          template <int D_, class P_> class... STAT_FROM>
void Add(RHist<Dimensions, WeightType, STAT_TO...> &to, const RHist<Dimensions, WeightType, STAT_FROM...> &from)
{
   // Enforce "same axis configuration" policy.
   auto& toImpl = *to.GetImpl();
   const auto& fromImpl = *from.GetImpl();
   for (int dim = 0; dim < Dimensions; ++dim) {
      if (!toImpl.GetAxis(dim).HasSameBinningAs(fromImpl.GetAxis(dim))) {
         throw std::runtime_error("Attempted to add RHists with incompatible axis binning");
      }
   }

   // Now that we know that the two axes have the same binning, we can just add
   // the statistics directly.
   toImpl.GetStat().Add(fromImpl.GetStat());
}
*/

} // namespace Experimental
} // namespace ROOT

#endif
