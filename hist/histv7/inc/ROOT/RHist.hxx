/// \file ROOT/RHist.hxx
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
#include "ROOT/RHistBinIter.hxx"
#include "ROOT/RHistImpl.hxx"
#include "ROOT/RHistData.hxx"
#include "ROOT/RLogger.hxx"
#include <initializer_list>
#include <stdexcept>

namespace ROOT {
namespace Experimental {

// fwd declare for fwd declare for friend declaration in RHist...
template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
class RHist;

// fwd declare for friend declaration in RHist.
template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
class RHist<DIMENSIONS, PRECISION, STAT...>
HistFromImpl(std::unique_ptr<typename RHist<DIMENSIONS, PRECISION, STAT...>::ImplBase_t> pHistImpl);

/**
 \class RHist
 Histogram class for histograms with `DIMENSIONS` dimensions, where each
 bin count is stored by a value of type `PRECISION`. STAT stores statistical
 data of the entries filled into the histogram (bin content, uncertainties etc).

 A histogram counts occurrences of values or n-dimensional combinations thereof.
 Contrary to for instance a `RTree`, a histogram combines adjacent values. The
 resolution of this combination is defined by the axis binning, see e.g.
 http://www.wikiwand.com/en/Histogram
 */

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
class RHist {
public:
   /// The type of the `Detail::RHistImplBase` of this histogram.
   using ImplBase_t =
      Detail::RHistImplBase<Detail::RHistData<DIMENSIONS, PRECISION, std::vector<PRECISION>, STAT...>>;
   /// The coordinates type: a `DIMENSIONS`-dimensional `std::array` of `double`.
   using CoordArray_t = typename ImplBase_t::CoordArray_t;
   /// The type of weights
   using Weight_t = PRECISION;
   /// Pointer type to `HistImpl_t::Fill`, for faster access.
   using FillFunc_t = typename ImplBase_t::FillFunc_t;
   /// Range.
   using AxisRange_t = typename ImplBase_t::AxisIterRange_t;

   using const_iterator = Detail::RHistBinIter<ImplBase_t>;

   /// Number of dimensions of the coordinates
   static constexpr int GetNDim() noexcept { return DIMENSIONS; }

   RHist() = default;
   RHist(RHist &&) = default;
   RHist(const RHist &other): fImpl(other.fImpl->Clone()), fFillFunc(other.fFillFunc)
   {}

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
   explicit RHist(std::array<RAxisConfig, DIMENSIONS> axes);

   /// Constructor overload taking the histogram title.
   RHist(std::string_view histTitle, std::array<RAxisConfig, DIMENSIONS> axes);

   /// Constructor overload that's only available for a 1-dimensional histogram.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 1>::type>
   explicit RHist(const RAxisConfig &xaxis): RHist(std::array<RAxisConfig, 1>{{xaxis}})
   {}

   /// Constructor overload that's only available for a 1-dimensional histogram,
   /// also passing the histogram title.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 1>::type>
   RHist(std::string_view histTitle, const RAxisConfig &xaxis): RHist(histTitle, std::array<RAxisConfig, 1>{{xaxis}})
   {}

   /// Constructor overload that's only available for a 2-dimensional histogram.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 2>::type>
   RHist(const RAxisConfig &xaxis, const RAxisConfig &yaxis): RHist(std::array<RAxisConfig, 2>{{xaxis, yaxis}})
   {}

   /// Constructor overload that's only available for a 2-dimensional histogram,
   /// also passing the histogram title.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 2>::type>
   RHist(std::string_view histTitle, const RAxisConfig &xaxis, const RAxisConfig &yaxis)
      : RHist(histTitle, std::array<RAxisConfig, 2>{{xaxis, yaxis}})
   {}

   /// Constructor overload that's only available for a 3-dimensional histogram.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 3>::type>
   RHist(const RAxisConfig &xaxis, const RAxisConfig &yaxis, const RAxisConfig &zaxis)
      : RHist(std::array<RAxisConfig, 3>{{xaxis, yaxis, zaxis}})
   {}

   /// Constructor overload that's only available for a 3-dimensional histogram,
   /// also passing the histogram title.
   template <int ENABLEIF_NDIM = DIMENSIONS, class = typename std::enable_if<ENABLEIF_NDIM == 3>::type>
   RHist(std::string_view histTitle, const RAxisConfig &xaxis, const RAxisConfig &yaxis, const RAxisConfig &zaxis)
      : RHist(histTitle, std::array<RAxisConfig, 3>{{xaxis, yaxis, zaxis}})
   {}

   /// Access the ImplBase_t this RHist points to.
   ImplBase_t *GetImpl() const noexcept { return fImpl.get(); }

   /// "Steal" the ImplBase_t this RHist points to.
   std::unique_ptr<ImplBase_t> TakeImpl() && noexcept { return std::move(fImpl); }

   /// Add `weight` to the bin containing coordinate `x`.
   void Fill(const CoordArray_t &x, Weight_t weight = (Weight_t)1) noexcept { (fImpl.get()->*fFillFunc)(x, weight); }

   /// For each coordinate in `xN`, add `weightN[i]` to the bin at coordinate
   /// `xN[i]`. The sizes of `xN` and `weightN` must be the same. This is more
   /// efficient than many separate calls to `Fill()`.
   void FillN(const std::span<const CoordArray_t> xN, const std::span<const Weight_t> weightN) noexcept
   {
      fImpl->FillN(xN, weightN);
   }

   /// Convenience overload: `FillN()` with weight 1.
   void FillN(const std::span<const CoordArray_t> xN) noexcept { fImpl->FillN(xN); }

   /// Get the number of entries this histogram was filled with.
   int64_t GetEntries() const noexcept { return fImpl->GetStat().GetEntries(); }

   /// Get the content of the bin at `x`.
   Weight_t GetBinContent(const CoordArray_t &x) const { return fImpl->GetBinContent(x); }

   /// Get the uncertainty on the content of the bin at `x`.
   double GetBinUncertainty(const CoordArray_t &x) const { return fImpl->GetBinUncertainty(x); }

   const_iterator begin() const { return const_iterator(*fImpl, 1); }

   const_iterator end() const { return const_iterator(*fImpl, fImpl->GetNBinsNoOver() + 1); }

   /// Swap *this and other.
   ///
   /// Very efficient; swaps the `fImpl` pointers.
   void swap(RHist<DIMENSIONS, PRECISION, STAT...> &other) noexcept
   {
      std::swap(fImpl, other.fImpl);
      std::swap(fFillFunc, other.fFillFunc);
   }

private:
   /// The actual histogram implementation.
   std::unique_ptr<ImplBase_t> fImpl;

   /// Pointer to RHistImpl::Fill() member function.
   FillFunc_t fFillFunc = nullptr;    //!

   friend RHist HistFromImpl<>(std::unique_ptr<ImplBase_t>);
};

/// RHist with no STAT parameter uses RHistStatContent by default.
template <int DIMENSIONS, class PRECISION>
class RHist<DIMENSIONS, PRECISION>: public RHist<DIMENSIONS, PRECISION, RHistStatContent> {
   using RHist<DIMENSIONS, PRECISION, RHistStatContent>::RHist;
};

/// Swap two histograms.
///
/// Very efficient; swaps the `fImpl` pointers.
template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
void swap(RHist<DIMENSIONS, PRECISION, STAT...> &a, RHist<DIMENSIONS, PRECISION, STAT...> &b) noexcept
{
   a.swap(b);
};

namespace Internal {
/**
 Generate RHist::fImpl from RHist constructor arguments.
 */
template <int NDIM, int IDIM, class DATA, class... PROCESSEDAXISCONFIG>
struct RHistImplGen {
   /// Select the template argument for the next axis type, and "recurse" into
   /// RHistImplGen for the next axis.
   template <RAxisConfig::EKind KIND>
   std::unique_ptr<Detail::RHistImplBase<DATA>>
   MakeNextAxis(std::string_view title, const std::array<RAxisConfig, NDIM> &axes,
                PROCESSEDAXISCONFIG... processedAxisArgs)
   {
      using NextAxis_t = typename AxisConfigToType<KIND>::Axis_t;
      NextAxis_t nextAxis = AxisConfigToType<KIND>()(axes[IDIM]);
      using HistImpl_t = RHistImplGen<NDIM, IDIM + 1, DATA, PROCESSEDAXISCONFIG..., NextAxis_t>;
      return HistImpl_t()(title, axes, processedAxisArgs..., nextAxis);
   }

   /// Make a RHistImpl-derived object reflecting the RAxisConfig array.
   ///
   /// Delegate to the appropriate MakeNextAxis instantiation, depending on the
   /// axis type selected in the RAxisConfig.
   /// \param title - title of the derived object.
   /// \param axes - `RAxisConfig` objects describing the axis of the resulting
   ///   RHistImpl.
   /// \param processedAxisArgs - the RAxisBase-derived axis objects describing the
   ///   axes of the resulting RHistImpl. There are `IDIM` of those; in the end
   /// (`IDIM` == `GetNDim()`), all `axes` have been converted to
   /// `processedAxisArgs` and the RHistImpl constructor can be invoked, passing
   /// the `processedAxisArgs`.
   std::unique_ptr<Detail::RHistImplBase<DATA>> operator()(std::string_view title,
                                                           const std::array<RAxisConfig, NDIM> &axes,
                                                           PROCESSEDAXISCONFIG... processedAxisArgs)
   {
      switch (axes[IDIM].GetKind()) {
      case RAxisConfig::kEquidistant: return MakeNextAxis<RAxisConfig::kEquidistant>(title, axes, processedAxisArgs...);
      case RAxisConfig::kGrow: return MakeNextAxis<RAxisConfig::kGrow>(title, axes, processedAxisArgs...);
      case RAxisConfig::kIrregular: return MakeNextAxis<RAxisConfig::kIrregular>(title, axes, processedAxisArgs...);
      default: R__LOG_ERROR(HistLog()) << "Unhandled axis kind";
      }
      return nullptr;
   }
};

/// Generate RHist::fImpl from constructor arguments; recursion end.
template <int NDIM, class DATA, class... PROCESSEDAXISCONFIG>
/// Create the histogram, now that all axis types and initializer objects are
/// determined.
struct RHistImplGen<NDIM, NDIM, DATA, PROCESSEDAXISCONFIG...> {
   using HistImplBase_t = ROOT::Experimental::Detail::RHistImplBase<DATA>;
   std::unique_ptr<HistImplBase_t>
   operator()(std::string_view title, const std::array<RAxisConfig, DATA::GetNDim()> &, PROCESSEDAXISCONFIG... axisArgs)
   {
      using HistImplt_t = Detail::RHistImpl<DATA, PROCESSEDAXISCONFIG...>;
      return std::make_unique<HistImplt_t>(title, axisArgs...);
   }
};
} // namespace Internal

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
RHist<DIMENSIONS, PRECISION, STAT...>::RHist(std::string_view title, std::array<RAxisConfig, DIMENSIONS> axes)
   : fImpl{std::move(
        Internal::RHistImplGen<RHist::GetNDim(), 0,
                               Detail::RHistData<DIMENSIONS, PRECISION, std::vector<PRECISION>, STAT...>>()(
           title, axes))}
{
   fFillFunc = fImpl->GetFillFunc();
}

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
RHist<DIMENSIONS, PRECISION, STAT...>::RHist(std::array<RAxisConfig, DIMENSIONS> axes): RHist("", axes)
{}

/// Adopt an external, stand-alone RHistImpl. The RHist will take ownership.
template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
RHist<DIMENSIONS, PRECISION, STAT...>
HistFromImpl(std::unique_ptr<typename RHist<DIMENSIONS, PRECISION, STAT...>::ImplBase_t> pHistImpl)
{
   RHist<DIMENSIONS, PRECISION, STAT...> ret;
   ret.fFillFunc = pHistImpl->GetFillFunc();
   std::swap(ret.fImpl, pHistImpl);
   return ret;
};

/// \name RHist Typedefs
///\{ Convenience typedefs (ROOT6-compatible type names)

// Keep them as typedefs, to make sure old-style documentation tools can understand them.
using RH1D = RHist<1, double, RHistStatContent, RHistStatUncertainty>;
using RH1F = RHist<1, float, RHistStatContent, RHistStatUncertainty>;
using RH1C = RHist<1, char, RHistStatContent>;
using RH1I = RHist<1, int, RHistStatContent>;
using RH1LL = RHist<1, int64_t, RHistStatContent>;

using RH2D = RHist<2, double, RHistStatContent, RHistStatUncertainty>;
using RH2F = RHist<2, float, RHistStatContent, RHistStatUncertainty>;
using RH2C = RHist<2, char, RHistStatContent>;
using RH2I = RHist<2, int, RHistStatContent>;
using RH2LL = RHist<2, int64_t, RHistStatContent>;

using RH3D = RHist<3, double, RHistStatContent, RHistStatUncertainty>;
using RH3F = RHist<3, float, RHistStatContent, RHistStatUncertainty>;
using RH3C = RHist<3, char, RHistStatContent>;
using RH3I = RHist<3, int, RHistStatContent>;
using RH3LL = RHist<3, int64_t, RHistStatContent>;
///\}

/// Add two histograms.
///
/// This operation may currently only be performed if the two histograms have
/// the same axis configuration, use the same precision, and if `from` records
/// at least the same statistics as `to` (recording more stats is fine).
///
/// Adding histograms with incompatible axis binning will be reported at runtime
/// with an `std::runtime_error`. Insufficient statistics in the source
/// histogram will be detected at compile-time and result in a compiler error.
///
/// In the future, we may either adopt a more relaxed definition of histogram
/// addition or provide a mechanism to convert from one histogram type to
/// another. We currently favor the latter path.
template <int DIMENSIONS, class PRECISION,
          template <int D_, class P_> class... STAT_TO,
          template <int D_, class P_> class... STAT_FROM>
void Add(RHist<DIMENSIONS, PRECISION, STAT_TO...> &to, const RHist<DIMENSIONS, PRECISION, STAT_FROM...> &from)
{
   // Enforce "same axis configuration" policy.
   auto& toImpl = *to.GetImpl();
   const auto& fromImpl = *from.GetImpl();
   for (int dim = 0; dim < DIMENSIONS; ++dim) {
      if (!toImpl.GetAxis(dim).HasSameBinningAs(fromImpl.GetAxis(dim))) {
         throw std::runtime_error("Attempted to add RHists with incompatible axis binning");
      }
   }

   // Now that we know that the two axes have the same binning, we can just add
   // the statistics directly.
   toImpl.GetStat().Add(fromImpl.GetStat());
}

} // namespace Experimental
} // namespace ROOT

#endif
