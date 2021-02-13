/// \file ROOT/RHistImpl.hxx
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

#ifndef ROOT7_RHistImpl
#define ROOT7_RHistImpl

#include <cassert>
#include <cctype>
#include <functional>
#include "ROOT/RSpan.hxx"
#include "ROOT/RTupleApply.hxx"

#include "ROOT/RAxis.hxx"
#include "ROOT/RHistBinIter.hxx"
#include "ROOT/RHistUtils.hxx"
#include "ROOT/RLogger.hxx"

class TRootIOCtor;

namespace ROOT {
namespace Experimental {

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
class RHist;

namespace Hist {
/// Iterator over n dimensional axes - an array of n axis iterators.
template <int NDIMS>
using AxisIter_t = std::array<RAxisBase::const_iterator, NDIMS>;
/// Range over n dimensional axes - a pair of arrays of n axis iterators.
template <int NDIMS>
using AxisIterRange_t = std::array<AxisIter_t<NDIMS>, 2>;

/// Kinds of under- and overflow handling.
enum class EOverflow {
   kNoOverflow = 0x0, ///< Exclude under- and overflows
   kUnderflow = 0x1,  ///< Include underflows
   kOverflow = 0x2,   ///< Include overflows
   kUnderOver = 0x3,  ///< Include both under- and overflows
};

inline bool operator&(EOverflow a, EOverflow b)
{
   return static_cast<int>(a) & static_cast<int>(b);
}
} // namespace Hist

namespace Detail {

/**
 \class RHistImplPrecisionAgnosticBase
 Base class for `RHistImplBase` that abstracts out the histogram's `PRECISION`.

 For operations such as painting a histogram, the `PRECISION` (type of the bin
 content) is not relevant; painting will cast the underlying bin type to double.
 To facilitate this, `RHistImplBase` itself inherits from the
 `RHistImplPrecisionAgnosticBase` interface.
 */
template <int DIMENSIONS>
class RHistImplPrecisionAgnosticBase {
public:
   /// Type of the coordinates.
   using CoordArray_t = Hist::CoordArray_t<DIMENSIONS>;
   /// Type of the local per-axis bin indices.
   using BinArray_t = std::array<int, DIMENSIONS>;
   /// Range type.
   using AxisIterRange_t = Hist::AxisIterRange_t<DIMENSIONS>;

   RHistImplPrecisionAgnosticBase() = default;
   RHistImplPrecisionAgnosticBase(const RHistImplPrecisionAgnosticBase &) = default;
   RHistImplPrecisionAgnosticBase(RHistImplPrecisionAgnosticBase &&) = default;
   RHistImplPrecisionAgnosticBase(std::string_view title): fTitle(title) {}
   virtual ~RHistImplPrecisionAgnosticBase() {}

   /// Number of dimensions of the coordinates.
   static constexpr int GetNDim() { return DIMENSIONS; }
   /// Number of bins of this histogram, including all overflow and underflow
   /// bins. Simply the product of all axes' total number of bins.
   virtual int GetNBins() const noexcept = 0;
   /// Number of bins of this histogram, excluding all overflow and underflow
   /// bins. Simply the product of all axes' number of regular bins.
   virtual int GetNBinsNoOver() const noexcept = 0;
   /// Number of under- and overflow bins of this histogram, excluding all
   /// regular bins.
   virtual int GetNOverflowBins() const noexcept = 0;

   /// Get the histogram title.
   const std::string &GetTitle() const { return fTitle; }

   /// Given the coordinate `x`, determine the index of the bin.
   virtual int GetBinIndex(const CoordArray_t &x) const = 0;
   /// Given the coordinate `x`, determine the index of the bin, possibly
   /// growing axes for which `x` is out of range.
   virtual int GetBinIndexAndGrow(const CoordArray_t &x) const = 0;

   /// Given the local per-axis bins `x`, determine the index of the bin.
   virtual int GetBinIndexFromLocalBins(const BinArray_t &x) const = 0;
   /// Given the index of the bin, determine the local per-axis bins `x`.
   virtual BinArray_t GetLocalBins(int binidx) const = 0;

   /// Get the center in all dimensions of the bin with index `binidx`.
   virtual CoordArray_t GetBinCenter(int binidx) const = 0;
   /// Get the lower edge in all dimensions of the bin with index `binidx`.
   virtual CoordArray_t GetBinFrom(int binidx) const = 0;
   /// Get the upper edge in all dimensions of the bin with index `binidx`.
   virtual CoordArray_t GetBinTo(int binidx) const = 0;

   /// Get the uncertainty of the bin with index `binidx`.
   virtual double GetBinUncertainty(int binidx) const = 0;

   /// Whether this histogram's statistics provide storage for uncertainties, or
   /// whether uncertainties are determined as poisson uncertainty of the content.
   virtual bool HasBinUncertainty() const = 0;

   /// The bin content, cast to double.
   virtual double GetBinContentAsDouble(int binidx) const = 0;

   /// Get a base-class view on axis with index `iAxis`.
   ///
   /// \param iAxis - index of the axis, must be `0 <= iAxis < DIMENSION`.
   virtual const RAxisBase &GetAxis(int iAxis) const = 0;

   /// Get an `AxisIterRange_t` for the whole histogram,
   /// excluding under- and overflow.
   virtual AxisIterRange_t GetRange() const = 0;

private:
   std::string fTitle; ///< The histogram's title.
};

/**
 \class RHistImplBase
 Interface class for `RHistImpl`.

 `RHistImpl` is templated for a specific configuration of axes. To enable access
 through `RHist`, `RHistImpl` inherits from `RHistImplBase`, exposing only dimension
 (`DIMENSION`) and bin type (`PRECISION`).
 */
template <class DATA>
class RHistImplBase: public RHistImplPrecisionAgnosticBase<DATA::GetNDim()> {
public:
   /// Type of the statistics (bin content, uncertainties etc).
   using Stat_t = DATA;
   /// Type of the coordinates.
   using CoordArray_t = Hist::CoordArray_t<DATA::GetNDim()>;
   /// Type of the local per-axis bin indices.
   using BinArray_t = std::array<int, DATA::GetNDim()>;
   /// Type of the bin content (and thus weights).
   using Weight_t = typename DATA::Weight_t;

   /// Type of the `Fill(x, w)` function
   using FillFunc_t = void (RHistImplBase::*)(const CoordArray_t &x, Weight_t w);

private:
   /// The histogram's bin content, uncertainties etc.
   Stat_t fStatistics;

public:
   RHistImplBase() = default;
   RHistImplBase(size_t numBins, size_t numOverflowBins): fStatistics(numBins, numOverflowBins) {}
   RHistImplBase(std::string_view title, size_t numBins, size_t numOverflowBins)
      : RHistImplPrecisionAgnosticBase<DATA::GetNDim()>(title), fStatistics(numBins, numOverflowBins)
   {}
   RHistImplBase(const RHistImplBase &) = default;
   RHistImplBase(RHistImplBase &&) = default;

   virtual std::unique_ptr<RHistImplBase> Clone() const = 0;

   /// Interface function to fill a vector or array of coordinates with
   /// corresponding weights.
   /// \note the size of `xN` and `weightN` must be the same!
   virtual void FillN(const std::span<const CoordArray_t> xN, const std::span<const Weight_t> weightN) = 0;

   /// Interface function to fill a vector or array of coordinates.
   virtual void FillN(const std::span<const CoordArray_t> xN) = 0;

   /// Retrieve the pointer to the overridden `Fill(x, w)` function.
   virtual FillFunc_t GetFillFunc() const = 0;

   /// Get the bin content (sum of weights) for the bin at coordinate `x`.
   virtual Weight_t GetBinContent(const CoordArray_t &x) const = 0;

   using RHistImplPrecisionAgnosticBase<DATA::GetNDim()>::GetBinUncertainty;

   /// Get the bin uncertainty for the bin at coordinate x.
   virtual double GetBinUncertainty(const CoordArray_t &x) const = 0;

   /// Get the number of bins in this histogram, including possible under- and
   /// overflow bins.
   int GetNBins() const noexcept final { return fStatistics.size(); }

   /// Get the number of bins in this histogram, excluding possible under- and
   /// overflow bins.
   int GetNBinsNoOver() const noexcept final { return fStatistics.sizeNoOver(); }

   /// Get the number of under- and overflow bins of this histogram, excluding all
   /// regular bins.
   int GetNOverflowBins() const noexcept final { return fStatistics.sizeUnderOver(); }

   /// Get the bin content (sum of weights) for bin index `binidx`.
   Weight_t GetBinContent(int binidx) const 
   {
      assert(binidx != 0);
      return fStatistics[binidx]; 
   }

   /// Get the bin content (sum of weights) for bin index `binidx` (non-const).
   Weight_t &GetBinContent(int binidx) 
   {
      assert(binidx != 0);
      return fStatistics[binidx]; 
   }

   /// Const access to statistics.
   const Stat_t &GetStat() const noexcept { return fStatistics; }

   /// Non-const access to statistics.
   Stat_t &GetStat() noexcept { return fStatistics; }

   /// Get the bin content (sum of weights) for bin index `binidx`, cast to
   /// `double`.
   double GetBinContentAsDouble(int binidx) const final { return (double)GetBinContent(binidx); }

   /// Add `w` to the bin at index `bin`.
   void AddBinContent(int binidx, Weight_t w) 
   {
      assert(binidx != 0);
      fStatistics[binidx] += w; 
   }
};
} // namespace Detail

namespace Internal {
/** \name Histogram traits
    Helper traits for histogram operations.
 */
///\{

/// Specifies if the wanted result is the bin's lower edge, center or higher edge.
enum class EBinCoord {
   kBinFrom,   ///< Get the lower bin edge
   kBinCenter, ///< Get the bin center
   kBinTo      ///< Get the bin high edge
};

/// Status of FindBin(x) and FindAdjustedBin(x)
enum class EFindStatus {
   kCanGrow, ///< The coordinate could fit after growing the axis
   kValid    ///< The returned bin index is valid
};

/// \name Axis tuple operations
/// Template operations on axis tuple.
///@{

/// Recursively gets the total number of bins in whole hist, excluding under- and overflow.
/// Each call gets the current axis' number of bins (excluding under- and overflow) multiplied
/// by that of the next axis.
template <int IDX, class AXISTUPLE>
struct RGetNBinsNoOverCount;

template <class AXES>
struct RGetNBinsNoOverCount<0, AXES> {
   int operator()(const AXES &axes) const { return std::get<0>(axes).GetNBinsNoOver(); }
};

template <int I, class AXES>
struct RGetNBinsNoOverCount {
   int operator()(const AXES &axes) const { return std::get<I>(axes).GetNBinsNoOver() * RGetNBinsNoOverCount<I - 1, AXES>()(axes); }
};

/// Get the number of bins in whole hist, excluding under- and overflow.
template <class... AXISCONFIG>
int GetNBinsNoOverFromAxes(AXISCONFIG... axisArgs)
{
   using axesTuple = std::tuple<AXISCONFIG...>;
   return RGetNBinsNoOverCount<sizeof...(AXISCONFIG) - 1, axesTuple>()(axesTuple{axisArgs...});
}

/// Recursively gets the total number of bins in whole hist, including under- and overflow.
/// Each call gets the current axis' number of bins (including under- and overflow) multiplied
/// by that of the next axis.
template <int IDX, class AXISTUPLE>
struct RGetNBinsCount;

template <class AXES>
struct RGetNBinsCount<0, AXES> {
   int operator()(const AXES &axes) const { return std::get<0>(axes).GetNBins(); }
};

template <int I, class AXES>
struct RGetNBinsCount {
   int operator()(const AXES &axes) const { return std::get<I>(axes).GetNBins() * RGetNBinsCount<I - 1, AXES>()(axes); }
};

/// Get the number of bins in whole hist, including under- and overflow.
template <class... AXISCONFIG>
int GetNBinsFromAxes(AXISCONFIG... axisArgs)
{
   using axesTuple = std::tuple<AXISCONFIG...>;
   return RGetNBinsCount<sizeof...(AXISCONFIG) - 1, axesTuple>()(axesTuple{axisArgs...});
}

/// Get the number of under- and overflow bins in whole hist, excluding regular bins.
template <class... AXISCONFIG>
int GetNOverflowBinsFromAxes(AXISCONFIG... axisArgs)
{
   using axesTuple = std::tuple<AXISCONFIG...>;
   return RGetNBinsCount<sizeof...(AXISCONFIG) - 1, axesTuple>()(axesTuple{axisArgs...}) - RGetNBinsNoOverCount<sizeof...(AXISCONFIG) - 1, axesTuple>()(axesTuple{axisArgs...});
}

/// Recursively fills the ranges of all axes, excluding under- and overflow.
/// Each call fills `range` with `begin()` and `end()` of the current axis, excluding
/// under- and overflow.
template <int I, class AXES>
struct RFillIterRange;

template <class AXES>
struct RFillIterRange<-1, AXES> {
   void operator()(Hist::AxisIterRange_t<std::tuple_size<AXES>::value> & /*range*/, const AXES & /*axes*/) const
   {}
};

template <int I, class AXES>
struct RFillIterRange {
   void operator()(Hist::AxisIterRange_t<std::tuple_size<AXES>::value> &range, const AXES &axes) const
   {
      range[0][I] = std::get<I>(axes).begin();
      range[1][I] = std::get<I>(axes).end();
      RFillIterRange<I - 1, AXES>()(range, axes);
   }
};

/// Recursively gets the number of regular bins just before the current dimension.
/// Each call gets the previous axis' number of regular bins multiplied
/// by the number of regular bins before the previous axis.
template <int I, int NDIMS, typename BINS, class AXES>
struct RGetNRegularBinsBefore;

template <int NDIMS, typename BINS, class AXES>
struct RGetNRegularBinsBefore<-1, NDIMS, BINS, AXES> {
   void operator()(BINS &/*binSizes*/, const AXES &/*axes*/) const
   {}
};

template <int I, int NDIMS, typename BINS, class AXES>
struct RGetNRegularBinsBefore {
   void operator()(BINS &binSizes, const AXES &axes) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      binSizes[thisAxis] = binSizes[thisAxis-1] * std::get<thisAxis-1>(axes).GetNBinsNoOver();
      RGetNRegularBinsBefore<I - 1, NDIMS, BINS, AXES>()(binSizes, axes);
   }
};

/// Recursively gets the total number of regular bins before the current dimension,
/// when computing a global bin that is in under- or overflow in at least one
/// dimension. That global bin's local per-axis bin indices are passed through
/// the `localBins` parameter. These `localBins` were translated to 0-based bins,
/// which is more convenient for some operations and which are the `virtualBins`
/// parameter.
/// Each call gets the current axis' number of regular bins before the global_bin
/// in the current dimension multiplied by the number of regular bins before the
/// current axis.
/// If the global_bin is in under- or overflow in the current dimension (local bin),
/// there is no need to process further.

//  - We want to know how many regular bins lie before the current overflow bin in the
// histogram's global binning order (which so far I thought was row-major, but now I'm
// not sure, maybe it's actually column-major... it doesn't matter, we don't need to spell out what is the global binning order anyway).

template <int I, int NDIMS, typename BINS, class AXES>
struct RComputeGlobalBin;

template <int NDIMS, typename BINS, class AXES>
struct RComputeGlobalBin<-1, NDIMS, BINS, AXES> {
   int operator()(int total_regular_bins_before, const AXES &/*axes*/, const BINS &/*virtualBins*/, const BINS &/*binSizes*/, const BINS &/*localBins*/) const
   {
      return total_regular_bins_before;
   }
};

template <int I, int NDIMS, typename BINS, class AXES>
struct RComputeGlobalBin {
   int operator()(int total_regular_bins_before, const AXES &axes, const BINS &virtualBins, const BINS &binSizes, const BINS &localBins) const
   {
      // We can tell how many regular bins lie before us on this axis,
      // accounting for the underflow bin of this axis if it has one.
      const int num_underflow_bins = static_cast<int>(!std::get<I>(axes).CanGrow());
      const int num_regular_bins_before =
         std::max(virtualBins[I] - num_underflow_bins, 0);
      total_regular_bins_before += num_regular_bins_before * binSizes[I];

      // If we are on an overflow or underflow bin on this axis, we know that
      // we don't need to look at the remaining axes. Projecting on those
      // dimensions would only take us into an hyperplane of over/underflow
      // bins for the current axis, and we know that by construction there
      // will be no regular bins in there.
      if (localBins[I] < 1)
         return total_regular_bins_before;
      
      return RComputeGlobalBin<I - 1, NDIMS, BINS, AXES>()(total_regular_bins_before, axes, virtualBins, binSizes, localBins);
   }
};

/// Recursively compute some quantities needed for `ComputeLocalBins`, namely
/// the total number of bins per hyperplane (overflow and regular) and the
/// number of regular bins per hyperplane on the hyperplanes that have them.
template <int I, int NDIMS, class AXES>
struct RComputeLocalBinsInitialisation;

template <int NDIMS, class AXES>
struct RComputeLocalBinsInitialisation<0, NDIMS, AXES> {
   void operator()(std::array<int, NDIMS-1> /* bins_per_hyperplane */, std::array<int, NDIMS-1> /* regular_bins_per_hyperplane */, const AXES & /*axes*/) const
   {}
};

template <int I, int NDIMS, class AXES>
struct RComputeLocalBinsInitialisation {
   void operator()(std::array<int, NDIMS-1>& bins_per_hyperplane, std::array<int, NDIMS-1>& regular_bins_per_hyperplane, const AXES &axes) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      bins_per_hyperplane[thisAxis] = Internal::RGetNBinsCount<thisAxis, AXES>()(axes);
      regular_bins_per_hyperplane[thisAxis] = Internal::RGetNBinsNoOverCount<thisAxis, AXES>()(axes);
      RComputeLocalBinsInitialisation<I - 1, NDIMS, AXES>()(bins_per_hyperplane, regular_bins_per_hyperplane, axes);
   }
};

/// Recursively computes the number of regular bins before the current dimension,
/// as well as the number of under- and overflow bins left to account for, after
/// the current dimension. If the latter is equal to 0, there is no need to process
/// further.
/// It is computing local bins that are in under- or overflow in at least one
/// dimension.
/// Starting at the highest dimension, it examines how many full hyperplanes of
/// regular bins lie before, then projects on the remaining dimensions.

template <int I, int NDIMS, class AXES>
struct RComputeLocalBins;

template <int NDIMS, class AXES>
struct RComputeLocalBins<0, NDIMS, AXES> {
   void operator()(const AXES &/*axes*/, int &/*unprocessed_previous_overflow_bin*/,
         int &/*num_regular_bins_before*/, std::array<int, NDIMS-1> /* bins_per_hyperplane */,
         std::array<int, NDIMS-1> /* regular_bins_per_hyperplane */, int /* curr_bins_per_hyperplane */,
         int /* curr_regular_bins_per_hyperplane */) const
   {}
};

template <int I, int NDIMS, class AXES>
struct RComputeLocalBins {
   void operator()(const AXES &axes, int &unprocessed_previous_overflow_bin,
         int &num_regular_bins_before, std::array<int, NDIMS-1> bins_per_hyperplane,
         std::array<int, NDIMS-1> regular_bins_per_hyperplane, int curr_bins_per_hyperplane,
         int curr_regular_bins_per_hyperplane) const
   {
      // Let's start by computing the contribution of the underflow
      // hyperplane (if any), in which we know there will be no regular bins
      const int num_underflow_hyperplanes =
            static_cast<int>(!std::get<I>(axes).CanGrow());
      const int bins_in_underflow_hyperplane =
            num_underflow_hyperplanes * bins_per_hyperplane[I-1];

      // Next, from the total number of bins per hyperplane and the number of
      // regular bins per hyperplane that has them, we deduce the number of
      // overflow bins per hyperplane that has regular bins.
      const int overflow_bins_per_regular_hyperplane =
            bins_per_hyperplane[I-1] - regular_bins_per_hyperplane[I-1];

      // This allows us to answer a key question: are there any under/overflow
      // bins on the hyperplanes that have regular bins? It may not be the
      // case if all of their axes are growable, and thus don't have overflow bins.
      if (overflow_bins_per_regular_hyperplane != 0) {
            // If so, we start by cutting off the contribution of the underflow
            // and overflow hyperplanes, to focus specifically on regular bins.
            const int overflow_bins_in_regular_hyperplanes =
               std::min(
                  std::max(
                     unprocessed_previous_overflow_bin
                           - bins_in_underflow_hyperplane,
                     0
                  ),
                  overflow_bins_per_regular_hyperplane
                        * std::get<I>(axes).GetNBinsNoOver()
               );

            // We count how many _complete_ "regular" hyperplanes that leaves
            // before us, and account for those in our regular bin count.
            const int num_regular_hyperplanes_before =
               overflow_bins_in_regular_hyperplanes
                  / overflow_bins_per_regular_hyperplane;
            num_regular_bins_before +=
               num_regular_hyperplanes_before
                  * regular_bins_per_hyperplane[I-1];

            // This only leaves the _current_ hyperplane as a possible source of
            // more regular bins that we haven't accounted for yet. We'll take
            // those into account while processing previous dimensions.
            unprocessed_previous_overflow_bin =
               overflow_bins_in_regular_hyperplanes
                  % overflow_bins_per_regular_hyperplane;
      } else {
            // If there are no overflow bins in regular hyperplane, then the
            // rule changes: observing _one_ overflow bin after the underflow
            // hyperplane means that _all_ regular hyperplanes on this axis are
            // already before us.
            if (unprocessed_previous_overflow_bin >= bins_in_underflow_hyperplane) {
               num_regular_bins_before +=
                  std::get<I>(axes).GetNBinsNoOver()
                        * regular_bins_per_hyperplane[I-1];
            }

            // In this case, we're done, because the current bin may only lie
            // in the underflow or underflow hyperplane of this axis. Which
            // means that there are no further regular bins to be accounted for
            // in the current hyperplane.
            unprocessed_previous_overflow_bin = 0;
      }

      // No need to continue this loop if we've taken into account all
      // overflow bins that were associated with regular bins.
      if (unprocessed_previous_overflow_bin == 0)
         return;
      
      return Internal::RComputeLocalBins<I - 1, NDIMS, AXES>()
                                 (axes, unprocessed_previous_overflow_bin, num_regular_bins_before, bins_per_hyperplane,
                                 regular_bins_per_hyperplane, curr_bins_per_hyperplane, curr_regular_bins_per_hyperplane);
   }
};

/// Recursively computes zero-based local bin indices, given...
///
/// - A zero-based global bin index
/// - The number of considered bins on each axis (can be either `GetNBinsNoOver`
///   or `GetNBins` depending on what you are trying to do)
/// - A policy of treating all bins as regular (i.e. no negative indices)
template <int I, int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeLocalBinsRaw;

template <int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeLocalBinsRaw<-1, NDIMS, BINS, AXES, BINTYPE> {
   void operator()(BINS & /*virtualBins*/, const AXES & /*axes*/, int /*zeroBasedGlobalBin*/, BINTYPE /*GetNBinType*/) const
   {}
};

template <int I, int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeLocalBinsRaw {
   void operator()(BINS &virtualBins, const AXES &axes, int zeroBasedGlobalBin, BINTYPE GetNBinType) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      virtualBins[thisAxis] = zeroBasedGlobalBin % (std::get<thisAxis>(axes).*GetNBinType)();
      RComputeLocalBinsRaw<I - 1, NDIMS, BINS, AXES, BINTYPE>()(virtualBins, axes, zeroBasedGlobalBin / (std::get<thisAxis>(axes).*GetNBinType)(), GetNBinType);
   }
};

/// Recursively computes a zero-based global bin index, given...
///
/// - A set of zero-based per-axis bin indices
/// - The number of considered bins on each axis (can be either `GetNBinsNoOver`
///   or `GetNBins` depending on what you are trying to do)
/// - A policy of treating all bins qs regular (i.e. no negative indices)
template <int I, int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeGlobalBinRaw;

template <int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeGlobalBinRaw<-1, NDIMS, BINS, AXES, BINTYPE> {
   int operator()(int globalVirtualBin, const AXES & /*axes*/, const BINS & /*zeroBasedLocalBins*/, int /*binSize*/, BINTYPE /*GetNBinType*/) const
   {
      return globalVirtualBin;
   }
};

template <int I, int NDIMS, typename BINS, class AXES, class BINTYPE>
struct RComputeGlobalBinRaw {
   int operator()(int globalVirtualBin, const AXES &axes, const BINS &zeroBasedLocalBins, int binSize, BINTYPE GetNBinType) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      globalVirtualBin += zeroBasedLocalBins[thisAxis] * binSize;
      binSize *= (std::get<thisAxis>(axes).*GetNBinType)();
      return Internal::RComputeGlobalBinRaw<I - 1, NDIMS, BINS, AXES, BINTYPE>()(globalVirtualBin, axes, zeroBasedLocalBins, binSize, GetNBinType);
   }
};

/// Recursively converts zero-based virtual bins where the underflow bin
/// has index `0` and the overflow bin has index `N+1` where `N` is the axis'
/// number of regular bins, to the standard `kUnderflowBin`/`kOverflowBin` for under/overflow
/// bin indexing convention.
///
/// For growable axes, must add 1 to go back to standard indices as their virtual
/// indexing convention is also 0-based, with zero designating the first regular bin.
template <int I, int NDIMS, typename BINS, class AXES>
struct RVirtualBinsToLocalBins;

template <int NDIMS, typename BINS, class AXES>
struct RVirtualBinsToLocalBins<-1, NDIMS, BINS, AXES> {
   void operator()(BINS & /*localBins*/, const AXES & /*axes*/, const BINS & /*virtualBins*/) const
   {}
};

template <int I, int NDIMS, typename BINS, class AXES>
struct RVirtualBinsToLocalBins {
   void operator()(BINS &localBins, const AXES &axes, const BINS &virtualBins) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      if ((!std::get<thisAxis>(axes).CanGrow()) && (virtualBins[thisAxis] == 0)) {
         localBins[thisAxis] = RAxisBase::kUnderflowBin;
      } else if ((!std::get<thisAxis>(axes).CanGrow()) && (virtualBins[thisAxis] == (std::get<thisAxis>(axes).GetNBins() - 1))) {
         localBins[thisAxis] = RAxisBase::kOverflowBin;
      } else {
         const int regular_bin_offset = -static_cast<int>(std::get<thisAxis>(axes).CanGrow());
         localBins[thisAxis] = virtualBins[thisAxis] - regular_bin_offset;
      }
      RVirtualBinsToLocalBins<I - 1, NDIMS, BINS, AXES>()(localBins, axes, virtualBins);
   }
};

/// Recursively converts local axis bins from the standard `kUnderflowBin`/`kOverflowBin` for under/overflow
/// bin indexing convention, to a "virtual bin" convention where the underflow bin
/// has index `0` and the overflow bin has index `N+1` where `N` is the axis'
/// number of regular bins.
///
/// For growable axes, subtract 1 from regular indices so that the indexing
/// convention remains zero-based (this means that there will be no "holes" in
/// global binning, which matters more than the choice of regular index base)
template <int I, int NDIMS, typename BINS, class AXES>
struct RLocalBinsToVirtualBins;

template <int NDIMS, typename BINS, class AXES>
struct RLocalBinsToVirtualBins<-1, NDIMS, BINS, AXES> {
   void operator()(BINS & /*virtualBins*/, const AXES & /*axes*/, const BINS & /*localBins*/) const
   {}
};

template <int I, int NDIMS, typename BINS, class AXES>
struct RLocalBinsToVirtualBins {
   void operator()(BINS &virtualBins, const AXES &axes, const BINS &localBins) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      switch (localBins[thisAxis]) {
         case RAxisBase::kUnderflowBin:
               virtualBins[thisAxis] = 0; break;
         case RAxisBase::kOverflowBin:
               virtualBins[thisAxis] = std::get<thisAxis>(axes).GetNBins() - 1; break;
         default:
               virtualBins[thisAxis] = localBins[thisAxis] - static_cast<int>(std::get<thisAxis>(axes).CanGrow());
      }
      RLocalBinsToVirtualBins<I - 1, NDIMS, BINS, AXES>()(virtualBins, axes, localBins);
   }
};

/// Find the per-axis local bin indices associated with a certain set of coordinates.
template <int I, int NDIMS, typename BINS, typename COORD, class AXES>
struct RFindLocalBins;

template <int NDIMS, typename BINS, typename COORD, class AXES>
struct RFindLocalBins<-1, NDIMS, BINS, COORD, AXES> {
   void operator()(BINS & /*localBins*/, const AXES & /*axes*/, const COORD & /*coords*/) const
   {}
};

template <int I, int NDIMS, typename BINS, typename COORD, class AXES>
struct RFindLocalBins {
   void operator()(BINS &localBins, const AXES &axes, const COORD &coords) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      localBins[thisAxis] = std::get<thisAxis>(axes).FindBin(coords[thisAxis]);
      RFindLocalBins<I - 1, NDIMS, BINS, COORD, AXES>()(localBins, axes, coords);
   }
};

/// Recursively converts local axis bins from the standard `kUnderflowBin`/`kOverflowBin` for
/// under/overflow bin indexing convention, to the corresponding bin coordinates.
template <int I, int NDIMS, typename BINS, typename COORD, class AXES>
struct RLocalBinsToCoords;

template <int NDIMS, typename BINS, typename COORD, class AXES>
struct RLocalBinsToCoords<-1, NDIMS, BINS, COORD, AXES> {
   void operator()(COORD & /*coords*/, const AXES & /*axes*/, const BINS & /*localBins*/, EBinCoord /*kind*/) const
   {}
};

template <int I, int NDIMS, typename BINS, typename COORD, class AXES>
struct RLocalBinsToCoords {
   void operator()(COORD &coords, const AXES &axes, const BINS &localBins, EBinCoord kind) const
   {
      constexpr const int thisAxis = NDIMS - I - 1;
      int axisbin = localBins[thisAxis];
      switch (kind) {
         case EBinCoord::kBinFrom: coords[thisAxis] = std::get<thisAxis>(axes).GetBinFrom(axisbin); break;
         case EBinCoord::kBinCenter: coords[thisAxis] = std::get<thisAxis>(axes).GetBinCenter(axisbin); break;
         case EBinCoord::kBinTo: coords[thisAxis] = std::get<thisAxis>(axes).GetBinTo(axisbin); break;
      }
      RLocalBinsToCoords<I - 1, NDIMS, BINS, COORD, AXES>()(coords, axes, localBins, kind);
   }
};

template <class... AXISCONFIG>
static std::array<const RAxisBase *, sizeof...(AXISCONFIG)> GetAxisView(const AXISCONFIG &... axes) noexcept
{
   std::array<const RAxisBase *, sizeof...(AXISCONFIG)> axisViews{{&axes...}};
   return axisViews;
}

///\}
} // namespace Internal

namespace Detail {

template <class DATA, class... AXISCONFIG>
class RHistImpl final: public RHistImplBase<DATA> {
   static_assert(sizeof...(AXISCONFIG) == DATA::GetNDim(), "Number of axes must equal histogram dimension");

   friend typename DATA::Hist_t;

public:
   using ImplBase_t = RHistImplBase<DATA>;
   using CoordArray_t = typename ImplBase_t::CoordArray_t;
   using BinArray_t = typename ImplBase_t::BinArray_t;
   using Weight_t = typename ImplBase_t::Weight_t;
   using typename ImplBase_t::FillFunc_t;
   template <int NDIMS = DATA::GetNDim()>
   using AxisIterRange_t = typename Hist::AxisIterRange_t<NDIMS>;

private:
   std::tuple<AXISCONFIG...> fAxes; ///< The histogram's axes

public:
   RHistImpl(TRootIOCtor *);
   RHistImpl(AXISCONFIG... axisArgs);
   RHistImpl(std::string_view title, AXISCONFIG... axisArgs);

   std::unique_ptr<ImplBase_t> Clone() const override {
      return std::unique_ptr<ImplBase_t>(new RHistImpl(*this));
   }

   /// Retrieve the fill function for this histogram implementation, to prevent
   /// the virtual function call for high-frequency fills.
   FillFunc_t GetFillFunc() const final { 
      return (FillFunc_t)&RHistImpl::Fill; 
   }

   /// Get the axes of this histogram.
   const std::tuple<AXISCONFIG...> &GetAxes() const { return fAxes; }

   /// Normalized axes access, converting from actual axis type to base class.
   const RAxisBase &GetAxis(int iAxis) const final { return *std::apply(Internal::GetAxisView<AXISCONFIG...>, fAxes)[iAxis]; }

   /// Computes a zero-based global bin index, given...
   ///
   /// - A set of zero-based per-axis bin indices
   /// - The number of considered bins on each axis (can be either `GetNBinsNoOver`
   ///   or `GetNBins` depending on what you are trying to do)
   /// - A policy of treating all bins qs regular (i.e. no negative indices)
   template <int NDIMS, typename BINTYPE>
   int ComputeGlobalBinRaw(const BinArray_t& zeroBasedLocalBins, BINTYPE GetNBinType) const {
      int result = 0;
      int binSize = 1;
      return Internal::RComputeGlobalBinRaw<NDIMS - 1, NDIMS, BinArray_t, decltype(fAxes), BINTYPE>()(result, fAxes, zeroBasedLocalBins, binSize, GetNBinType);
   }

   /// Computes zero-based local bin indices, given...
   ///
   /// - A zero-based global bin index
   /// - The number of considered bins on each axis (can be either `GetNBinsNoOver`
   ///   or `GetNBins` depending on what you are trying to do)
   /// - A policy of treating all bins as regular (i.e. no negative indices)
   template <int NDIMS, typename BINTYPE>
   BinArray_t ComputeLocalBinsRaw(int zeroBasedGlobalBin, BINTYPE GetNBinType) const {
      BinArray_t result;
      Internal::RComputeLocalBinsRaw<NDIMS - 1, NDIMS, BinArray_t, decltype(fAxes), BINTYPE>()(result, fAxes, zeroBasedGlobalBin, GetNBinType);
      return result;
   }

   /// Converts local axis bins from the standard `kUnderflowBin`/`kOverflowBin` for under/overflow
   /// bin indexing convention, to a "virtual bin" convention where the underflow bin
   /// has index `0` and the overflow bin has index `N+1` where `N` is the axis'
   /// number of regular bins.
   template <int NDIMS>
   BinArray_t LocalBinsToVirtualBins(const BinArray_t& localBins) const {
      BinArray_t virtualBins;
      Internal::RLocalBinsToVirtualBins<NDIMS - 1, NDIMS, BinArray_t, decltype(fAxes)>()(virtualBins, fAxes, localBins);
      return virtualBins;
   }

   /// Converts zero-based virtual bins where the underflow bin has
   /// index `0` and the overflow bin has index `N+1` where `N` is the axis'
   /// number of regular bins, to the standard `kUnderflowBin`/`kOverflowBin` for under/overflow
   /// bin indexing convention.
   template <int NDIMS>
   BinArray_t VirtualBinsToLocalBins(const BinArray_t& virtualBins) const {
      BinArray_t localBins = {};
      Internal::RVirtualBinsToLocalBins<NDIMS - 1, NDIMS, BinArray_t, decltype(fAxes)>()(localBins, fAxes, virtualBins);
      return localBins;
   }

   /// Computes the global index of a certain bin on an `NDIMS`-dimensional histogram,
   /// knowing the local per-axis bin indices as returned by calling `FindBin()` on each axis.
   template <int NDIMS>
   int ComputeGlobalBin(BinArray_t& local_bins) const {
      // Get regular bins out of the way
      if (std::all_of(local_bins.cbegin(), local_bins.cend(),
                     [](int bin) { return bin >= 1; })) {
         for (int bin = 0; bin < NDIMS; bin++)
            local_bins[bin] -= 1;
         return ComputeGlobalBinRaw<NDIMS>(local_bins, &ROOT::Experimental::RAxisBase::GetNBinsNoOver) + 1;
      }

      // Convert bin indices to a zero-based coordinate system where the underflow
      // bin (if any) has coordinate 0 and the overflow bin (if any) has
      // coordinate N-1, where N is the axis' total number of bins.
      BinArray_t virtual_bins = LocalBinsToVirtualBins<NDIMS>(local_bins);

      // Deduce what the global bin index would be in this coordinate system that
      // unifies regular and overflow bins.
      const int global_virtual_bin = ComputeGlobalBinRaw<NDIMS>(virtual_bins, &ROOT::Experimental::RAxisBase::GetNBins);

      // Move to 1-based and negative indexing
      const int neg_1based_virtual_bin = -global_virtual_bin - 1;

      // At this point, we have an index that represents a count of all bins, both
      // regular and overflow, that are located before the current bin when
      // enumerating histogram bins in row-major order.
      //
      // We will next count the number of _regular_ bins which are located before
      // the current bin, and by removing this offset from the above index, we
      // will get a count of overflow bins that are located before the current bin
      // in row-major order. Which is what we want as our overflow bin index.
      //
      int total_regular_bins_before = 0;

      // First, we need to know how many regular bins we leave behind us for each
      // step on each axis, assuming that the bin from which we come was regular.
      //
      // If mathematically inclined, you can also think of this as the size of an
      // hyperplane of regular bins when projecting on lower-numbered dimensions.
      // See the docs of ComputeLocalBins for more on this worldview.
      //
      BinArray_t bin_sizes;
      bin_sizes[0] = 1;
      Internal::RGetNRegularBinsBefore<NDIMS - 2, NDIMS, BinArray_t, decltype(fAxes)>()(bin_sizes, fAxes);

      // With that, we can deduce how many regular bins lie before us.
      total_regular_bins_before = Internal::RComputeGlobalBin<NDIMS - 1, NDIMS, BinArray_t, decltype(fAxes)>()
         (total_regular_bins_before, fAxes, virtual_bins, bin_sizes, local_bins);

      // Now that we know how many bins lie before us, and how many of those are
      // regular bins, we can trivially deduce how many overflow bins lie before
      // us, and emit that as our global overflow bin index.
      return neg_1based_virtual_bin + total_regular_bins_before;
   }

   /// Computes the local per-axis bin indices of a certain bin on an `NDIMS`-dimensional histogram,
   /// knowing the global histogram bin index.
   template <int NDIMS>
   BinArray_t ComputeLocalBins(int global_bin) const {
      // Get regular bins out of the way
      if (global_bin >= 1) {
         BinArray_t computed_bins = ComputeLocalBinsRaw<NDIMS>(global_bin - 1, &ROOT::Experimental::RAxisBase::GetNBinsNoOver);
         for (int bin = 0; bin < NDIMS; ++bin)
            computed_bins[bin] += 1;
         return computed_bins;
      }

      // Convert our negative index to something positive and 0-based, as that is
      // more convenient to work with. Note, however, that this is _not_
      // equivalent to the virtual_bin that we had before, because what we have
      // here is a count of overflow bins, not of all bins...
      const int corrected_virtual_overflow_bin = -global_bin - 1;

      // ...so we need to retrieve and bring back the regular bin count, and this
      // is where the fun begins.
      //
      // The main difficulty is that the number of regular bins is not fixed as
      // one slides along a histogram axis. Using a 2D binning case as a simple
      // motivating example...
      //
      //    -1   -2   -3   -4   <- No regular bins on the underflow line of axis 1
      //    -5    1    2   -6   <- Some of them on middle lines of axis 1
      //    -7    3    4   -8
      //    -9   -10  -11  -12  <- No regular bins on the overflow line of axis 1
      //
      // As we go to higher dimensions, the geometry becomes more complex, but
      // if we replace "line" with "plane", we get a similar picture in 3D when we
      // slide along axis 2:
      //
      //  No regular bins on the    Some of them on the     No regular bins again
      //    UF plane of axis 2    regular planes of ax.2   on the OF plane of ax.2
      //
      //    -1   -2   -3   -4       -17  -18  -19  -20      -29  -30  -31  -32
      //    -5   -6   -7   -8       -21   1    2   -22      -33  -34  -35  -36
      //    -9   -10  -11  -12      -23   3    4   -24      -37  -37  -39  -40
      //    -13  -14  -15  -16      -25  -26  -27  -28      -41  -42  -43  -44
      //
      // We can generalize this to N dimensions by saying that as we slide along
      // the last axis of an N-d histogram, we see an hyperplane full of overflow
      // bins, then some hyperplanes with regular bins in the "middle" surrounded
      // by overflow bins, then a last hyperplane full of overflow bins.
      //
      // From this, we can devise a recursive algorithm to recover the number of
      // regular bins before the overflow bin we're currently looking at:
      //
      // - Start by processing the last histogram axis.
      // - Ignore the first and last hyperplane on this axis, which only contain
      //   underflow and overflow bins respectively.
      // - Count how many complete hyperplanes of regular bins lie before us on
      //   this axis, which we can do indirectly in our overflow bin based
      //   reasoning by computing the perimeter of the regular region and dividing
      //   our "regular" overflow bin count by that amount.
      // - Now we counted previous hyperplanes on this last histogram axis, but
      //   we need to process the hyperplane that our bin is located in, if any.
      //      * For this, we reduce our overflow bin count to a count of
      //        _unaccounted_ overflow bins in the current hyperplane...
      //      * ...which allows us to recursively continue the computation by
      //        processing the next (well, previous) histogram axis in the context
      //        of this hyperplane, in the same manner as above.
      //
      // Alright, now that the general plan is sorted out, let's compute some
      // quantities that we are going to need, namely the total number of bins per
      // hyperplane (overflow and regular) and the number of regular bins per
      // hyperplane on the hyperplanes that have them.
      //
      std::array<int, NDIMS - 1> bins_per_hyperplane;
      std::array<int, NDIMS - 1> regular_bins_per_hyperplane;
      Internal::RComputeLocalBinsInitialisation<NDIMS - 1, NDIMS, decltype(fAxes)>()(bins_per_hyperplane, regular_bins_per_hyperplane, fAxes);
      
      int curr_bins_per_hyperplane = Internal::RGetNBinsCount<NDIMS - 1, decltype(fAxes)>()(fAxes);
      int curr_regular_bins_per_hyperplane = Internal::RGetNBinsNoOverCount<NDIMS - 1, decltype(fAxes)>()(fAxes);

      // Given that, we examine each axis, starting from the last one.
      int unprocessed_previous_overflow_bin = corrected_virtual_overflow_bin;
      int num_regular_bins_before = 0;
      Internal::RComputeLocalBins<NDIMS - 1, NDIMS, decltype(fAxes)>()
                                 (fAxes, unprocessed_previous_overflow_bin, num_regular_bins_before, bins_per_hyperplane,
                                 regular_bins_per_hyperplane, curr_bins_per_hyperplane, curr_regular_bins_per_hyperplane);

      // By the time we reach the first axis, there should only be at most one
      // full row of regular bins before us:
      //
      //    -1  1  2  3  -2
      //     ^            ^
      //     |            |
      //     |        Option 2: one overflow bin before us
      //     |
      // Option 1: no overflow bin before us
      //
      num_regular_bins_before +=
         unprocessed_previous_overflow_bin * std::get<0>(fAxes).GetNBinsNoOver();

      // Now that we know the number of regular bins before us, we can add this to
      // to the zero-based overflow bin index that we started with to get a global
      // zero-based bin index accounting for both under/overflow bins and regular
      // bins, just like what we had in the ComputeGlobalBin<DATA::GetNDim()>() implementation.
      const int global_virtual_bin =
         corrected_virtual_overflow_bin + num_regular_bins_before;

      // We can then easily go back to zero-based "virtual" bin indices...
      const BinArray_t virtual_bins = ComputeLocalBinsRaw<NDIMS>(global_virtual_bin, &ROOT::Experimental::RAxisBase::GetNBins);

      // ...and from that go back to the -1/-2 overflow bin indexing convention.
      return VirtualBinsToLocalBins<NDIMS>(virtual_bins);
   }

   /// Get the bin index for the given coordinates `x`. The use of `RFindLocalBins`
   /// allows to convert the coordinates to local per-axis bin indices before using
   /// `ComputeGlobalBin()`.
   int GetBinIndex(const CoordArray_t &x) const final
   {
      BinArray_t localBins = {};
      Internal::RFindLocalBins<DATA::GetNDim() - 1, DATA::GetNDim(), BinArray_t, CoordArray_t, decltype(fAxes)>()(localBins, fAxes, x);
      int result = ComputeGlobalBin<DATA::GetNDim()>(localBins);
      return result;
   }

   /// Get the bin index for the given coordinates `x`, growing the axes as needed.
   /// The use of `RFindLocalBins` allows to convert the coordinates to local
   /// per-axis bin indices before using `ComputeGlobalBin()`.
   ///
   /// TODO: implement growable behavior
   int GetBinIndexAndGrow(const CoordArray_t &x) const final
   {
      Internal::EFindStatus status = Internal::EFindStatus::kCanGrow;
      int ret = 0;
      BinArray_t localBins = {};
      while (status == Internal::EFindStatus::kCanGrow) {
         Internal::RFindLocalBins<DATA::GetNDim() - 1, DATA::GetNDim(), BinArray_t, CoordArray_t, decltype(fAxes)>()(localBins, fAxes, x);
         ret = ComputeGlobalBin<DATA::GetNDim()>(localBins);
         status = Internal::EFindStatus::kValid;
      }
      return ret;
   }

   /// Get the bin index for the given local per-axis bin indices `x`, using
   /// `ComputeGlobalBin()`.
   int GetBinIndexFromLocalBins(const BinArray_t &x) const final
   {
      BinArray_t localBins = x;
      int result = ComputeGlobalBin<DATA::GetNDim()>(localBins);
      return result;
   }

   /// Get the local per-axis bin indices `x` for the given bin index, using
   /// `ComputeLocalBins()`.
   BinArray_t GetLocalBins(int binidx) const final
   {
      BinArray_t localBins = ComputeLocalBins<DATA::GetNDim()>(binidx);
      return localBins;
   }

   /// Get the center coordinates of the bin with index `binidx`.
   CoordArray_t GetBinCenter(int binidx) const final
   {
      BinArray_t localBins = ComputeLocalBins<DATA::GetNDim()>(binidx);
      CoordArray_t coords;
      Internal::RLocalBinsToCoords<DATA::GetNDim() - 1, DATA::GetNDim(), BinArray_t, CoordArray_t, decltype(fAxes)>()(coords, fAxes, localBins, Internal::EBinCoord::kBinCenter);
      return coords;
   }

   /// Get the coordinates of the low limit of the bin with index `binidx`.
   CoordArray_t GetBinFrom(int binidx) const final
   {
      BinArray_t localBins = ComputeLocalBins<DATA::GetNDim()>(binidx);
      CoordArray_t coords;
      Internal::RLocalBinsToCoords<DATA::GetNDim() - 1, DATA::GetNDim(), BinArray_t, CoordArray_t, decltype(fAxes)>()(coords, fAxes, localBins, Internal::EBinCoord::kBinFrom);
      return coords;
   }

   /// Get the coordinates of the high limit of the bin with index `binidx`.
   CoordArray_t GetBinTo(int binidx) const final
   {
      BinArray_t localBins = ComputeLocalBins<DATA::GetNDim()>(binidx);
      CoordArray_t coords;
      Internal::RLocalBinsToCoords<DATA::GetNDim() - 1, DATA::GetNDim(), BinArray_t, CoordArray_t, decltype(fAxes)>()(coords, fAxes, localBins, Internal::EBinCoord::kBinTo);
      return coords;
   }

   /// Fill an array of `weightN` to the bins specified by coordinates `xN`.
   /// For each element `i`, the weight `weightN[i]` will be added to the bin
   /// at the coordinate `xN[i]`
   /// \note `xN` and `weightN` must have the same size!
   void FillN(const std::span<const CoordArray_t> xN, const std::span<const Weight_t> weightN) final
   {
#ifndef NDEBUG
      if (xN.size() != weightN.size()) {
         R__LOG_ERROR(HistLog()) << "Not the same number of points and weights!";
         return;
      }
#endif

      for (size_t i = 0; i < xN.size(); ++i) {
         Fill(xN[i], weightN[i]);
      }
   }

   /// Fill an array of `weightN` to the bins specified by coordinates `xN`.
   /// For each element `i`, the weight `weightN[i]` will be added to the bin
   /// at the coordinate `xN[i]`
   void FillN(const std::span<const CoordArray_t> xN) final
   {
      for (auto &&x: xN) {
         Fill(x);
      }
   }

   /// Add a single weight `w` to the bin at coordinate `x`.
   void Fill(const CoordArray_t &x, Weight_t w = 1.)
   {
      int bin = GetBinIndexAndGrow(x);
      this->GetStat().Fill(x, bin, w);
   }

   /// Get the content of the bin at position `x`.
   Weight_t GetBinContent(const CoordArray_t &x) const final
   {
      int bin = GetBinIndex(x);
      return ImplBase_t::GetBinContent(bin);
   }

   /// Return the uncertainties for the given bin index.
   double GetBinUncertainty(int binidx) const final { return this->GetStat().GetBinUncertainty(binidx); }

   /// Get the bin uncertainty for the bin at coordinate `x`.
   double GetBinUncertainty(const CoordArray_t &x) const final
   {
      const int bin = GetBinIndex(x);
      return this->GetBinUncertainty(bin);
   }

   /// Whether this histogram's statistics provide storage for uncertainties, or
   /// whether uncertainties are determined as poisson uncertainty of the content.
   bool HasBinUncertainty() const final { return this->GetStat().HasBinUncertainty(); }

   /// Get the begin() and end() for each axis.
   AxisIterRange_t<DATA::GetNDim()>
   GetRange() const final
   {
      std::array<std::array<RAxisBase::const_iterator, DATA::GetNDim()>, 2> ret;
      Internal::RFillIterRange<DATA::GetNDim() - 1, decltype(fAxes)>()(ret, fAxes);
      return ret;
   }

   /// Grow the axis number `iAxis` to fit the coordinate `x`.
   ///
   /// The histogram (conceptually) combines pairs of bins along this axis until
   /// `x` is within the range of the axis.
   /// The axis must support growing for this to work (e.g. a `RAxisGrow`).
   void GrowAxis(int /*iAxis*/, double /*x*/)
   {
      // TODO: Implement GrowAxis()
   }

   /// \{
   /// \name Iterator interface
   using const_iterator = RHistBinIter<const ImplBase_t>;
   using iterator = RHistBinIter<ImplBase_t>;
   iterator begin() noexcept { return iterator(*this); }
   const_iterator begin() const noexcept { return const_iterator(*this); }
   iterator end() noexcept { return iterator(*this, this->GetNBinsNoOver()); }
   const_iterator end() const noexcept { return const_iterator(*this, this->GetNBinsNoOver()); }
   /// \}
};

template <class DATA, class... AXISCONFIG>
RHistImpl<DATA, AXISCONFIG...>::RHistImpl(TRootIOCtor *)
{}

template <class DATA, class... AXISCONFIG>
RHistImpl<DATA, AXISCONFIG...>::RHistImpl(AXISCONFIG... axisArgs)
   : ImplBase_t(Internal::GetNBinsNoOverFromAxes(axisArgs...), Internal::GetNOverflowBinsFromAxes(axisArgs...)), fAxes{axisArgs...}
{}

template <class DATA, class... AXISCONFIG>
RHistImpl<DATA, AXISCONFIG...>::RHistImpl(std::string_view title, AXISCONFIG... axisArgs)
   : ImplBase_t(title, Internal::GetNBinsNoOverFromAxes(axisArgs...), Internal::GetNOverflowBinsFromAxes(axisArgs...)), fAxes{axisArgs...}
{}

#if 0
// In principle we can also have a runtime version of RHistImpl, that does not
// contain a tuple of concrete axis types but a vector of `RAxisConfig`.
template <class DATA>
class RHistImplRuntime: public RHistImplBase<DATA> {
public:
  RHistImplRuntime(std::array<RAxisConfig, DATA::GetNDim()>&& axisCfg);
};
#endif

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
