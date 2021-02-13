/// \file ROOT/RAxis.hxx
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

#ifndef ROOT7_RAxis
#define ROOT7_RAxis

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "ROOT/RAxisConfig.hxx"
#include "ROOT/RStringView.hxx"
#include "ROOT/RLogger.hxx"

namespace ROOT {
namespace Experimental {

/**
 \class RAxisBase
 Histogram axis base class. Keeps track of the number of bins and overflow
 handling. Offers bin iteration.

 Regular bin indices are starting from 1, up to N + 1 for an axis with N bins.
 Index -1 is for the underflow bin, representing values that are lower than
 the axis range. Index -2 is the overflow bin for values larger than the axis
 range.
 Growable axes do not have underflow or overflow bins, as they don't need them.
 */
class RAxisBase {
protected:
   ///\name Inaccessible copy, assignment
   /// The copy and move constructors and assignment operators are protected to
   /// prevent slicing.
   ///\{
   RAxisBase(const RAxisBase &) = default;
   RAxisBase(RAxisBase &&) = default;
   RAxisBase &operator=(const RAxisBase &) = default;
   RAxisBase &operator=(RAxisBase &&) = default;
   ///\}

   /// Default construct a RAxisBase (for use by derived classes for I/O)
   RAxisBase() noexcept(noexcept(std::string())) = default;

   /// Virtual destructor needed in this inheritance-based design
   virtual ~RAxisBase();

   /// Construct a RAxisBase.
   ///
   ///\param[in] title - axis title used for graphics and text representation.
   RAxisBase(std::string_view title) noexcept: fTitle(title) {}

   /// Given rawbin (`<0` for underflow, `>=GetNBinsNoOver()` for overflow),
   /// determine the bin number taking into account how over/underflow
   /// should be handled.
   ///
   /// \param[in] rawbin for which to determine the bin number.
   /// \return Returns the bin number adjusted for potential over- and underflow
   /// bins. Returns `kInvalidBin` if the axis cannot handle the over- / underflow.
   ///
   int AdjustOverflowBinNumber(double rawbin) const
   {
      ++rawbin;

      // Underflow: Put in underflow bin if any, otherwise ignore
      if (rawbin < GetFirstBin())
         return CanGrow() ? kInvalidBin : GetUnderflowBin();

      // Overflow: Put in overflow bin if any, otherwise ignore
      // `rawbin` is not an integer, cannot compare `rawbin > GetLastBin()`.
      if (rawbin >= GetLastBin() + 1)
         return CanGrow() ? kInvalidBin : GetOverflowBin();

      // Bin index is in range and has been corrected for over/underflow
      return (int)rawbin;
   }

   /// Check if two axis have the same bin borders
   ///
   /// Default implementation should work for any RAxis type, but is quite
   /// inefficient as it does virtual GetBinFrom calls in a loop. RAxis
   /// implementations are encouraged to provide optimized overrides for common
   /// axis binning comparison scenarios.
   virtual bool HasSameBinBordersAs(const RAxisBase& other) const {
      // Axis growability (and thus under/overflow bin existence) must match
      if (CanGrow() != other.CanGrow())
         return false;

      // Number of normal bins must match
      if (GetNBinsNoOver() != other.GetNBinsNoOver())
         return false;

      // Left borders of normal bins must match
      for (int bin: *this)
         if (GetBinFrom(bin) != other.GetBinFrom(bin))
            return false;

      // Right border of the last normal bin (aka maximum) must also match
      if (GetMaximum() != other.GetMaximum())
         return false;

      // If all of these checks passed, the two axes have the same bin borders
      return true;
   }

public:
   /**
    \class const_iterator
    Random const_iterator through bins. Represents the bin index, not a bin
    content: the axis has no notion of any content.
    */
   class const_iterator: public std::iterator<std::random_access_iterator_tag, int /*value*/, int /*distance*/,
                                              const int * /*pointer*/, const int & /*ref*/> {
      int fCursor = 0; ///< Current iteration position

   public:
      const_iterator() = default;

      /// Initialize a const_iterator with its position
      explicit const_iterator(int cursor) noexcept: fCursor(cursor) {}

      /// ++i
      const_iterator &operator++() noexcept
      {
         // Could check whether fCursor < fEnd - but what for?
         ++fCursor;
         return *this;
      }

      /// --i
      const_iterator &operator--() noexcept
      {
         // Could check whether fCursor > fBegin - but what for?
         --fCursor;
         return *this;
      }

      /// i++
      const_iterator operator++(int)noexcept
      {
         const_iterator old(*this);
         ++(*this);
         return old;
      }

      // i--
      const_iterator operator--(int)noexcept
      {
         const_iterator old(*this);
         --(*this);
         return old;
      }

      // i += 2
      const_iterator &operator+=(int d) noexcept
      {
         fCursor += d;
         return *this;
      }

      // i -= 2
      const_iterator &operator-=(int d) noexcept
      {
         fCursor -= d;
         return *this;
      }

      // i + 2
      const_iterator operator+(int d) noexcept
      {
         const_iterator ret(*this);
         ret += d;
         return ret;
      }
      friend const_iterator operator+(int d, const_iterator rhs) noexcept;

      // i - 2
      const_iterator operator-(int d) noexcept
      {
         const_iterator ret(*this);
         ret -= d;
         return ret;
      }

      // i - j
      int operator-(const const_iterator& j) noexcept
      {
         return fCursor - j.fCursor;
      }

      // i[2]
      int operator[](int d) noexcept
      {
         return fCursor + d;
      }

      // *i
      int operator*() const noexcept { return fCursor; }

      // i->
      const int *operator->() const noexcept { return &fCursor; }

      friend bool operator<(const_iterator lhs, const_iterator rhs) noexcept;
      friend bool operator>(const_iterator lhs, const_iterator rhs) noexcept;
      friend bool operator<=(const_iterator lhs, const_iterator rhs) noexcept;
      friend bool operator>=(const_iterator lhs, const_iterator rhs) noexcept;
      friend bool operator==(const_iterator lhs, const_iterator rhs) noexcept;
      friend bool operator!=(const_iterator lhs, const_iterator rhs) noexcept;
   };

   /// Special bin index returned to signify that no bin matches a request.
   constexpr static const int kInvalidBin = 0;

   /// Index of the underflow bin, if any.
   constexpr static const int kUnderflowBin = -1;

   /// Index of the overflow bin, if any.
   constexpr static const int kOverflowBin = -2;

   /// Get the axis's title
   const std::string &GetTitle() const { return fTitle; }

   /// Whether this axis can grow (and thus has no overflow bins).
   virtual bool CanGrow() const noexcept = 0;

   /// Get the number of bins, excluding under- and overflow.
   virtual int GetNBinsNoOver() const noexcept = 0;

   /// Get the number of bins, including under- and overflow.
   int GetNBins() const noexcept { return GetNBinsNoOver() + GetNOverflowBins(); }

   /// Get the number of over- and underflow bins: 0 for growable axes, 2 otherwise.
   int GetNOverflowBins() const noexcept
   {
      if (CanGrow())
         return 0;
      else
         return 2;
   };

   /// Get the bin index for the underflow bin (or `kInvalidBin`
   /// if CanGrow()).
   int GetUnderflowBin() const noexcept {
      if (CanGrow())
         return kInvalidBin;
      else
         return kUnderflowBin;
   }

   /// Get the bin index for the overflow bin (or `kInvalidBin`
   /// if CanGrow()).
   int GetOverflowBin() const noexcept {
      if (CanGrow())
         return kInvalidBin;
      else
         return kOverflowBin;
   }

   /// Get the bin index for the first bin of the axis
   int GetFirstBin() const noexcept { return 1; }

   /// Get the bin index for the last bin of the axis
   int GetLastBin() const noexcept { return GetNBinsNoOver(); }

   ///\name Iterator interfaces
   ///\{

   /// Get a const_iterator pointing to the first regular bin.
   const_iterator begin() const noexcept { return const_iterator{GetFirstBin()}; }

   /// Get a const_iterator pointing beyond the last regular bin
   const_iterator end() const noexcept { return const_iterator{GetLastBin() + 1}; }
   ///\}

   /// Find the adjusted bin index (returning `kUnderflowBin` for underflow and `kOverflowBin`
   /// for overflow) for the given coordinate.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   virtual int FindBin(double x) const noexcept = 0;

   /// Get the bin center for the given bin index.
   /// The result of this method on an overflow or underflow bin is unspecified.
   virtual double GetBinCenter(int bin) const = 0;

   /// Get the low bin border ("left edge") for the given bin index.
   /// The result of this method on an underflow bin is unspecified.
   virtual double GetBinFrom(int bin) const = 0;

   /// Get the high bin border ("right edge") for the given bin index.
   /// The result of this method on an overflow bin is unspecified.
   double GetBinTo(int bin) const {
      const double result = (bin == kUnderflowBin) ? GetMinimum() : GetBinFrom(bin + 1);
      return result;
   }

   /// Get the low end of the axis range.
   double GetMinimum() const { return GetBinFrom(GetFirstBin()); }

   /// Get the high end of the axis range.
   double GetMaximum() const { return GetBinTo(GetLastBin()); }

   /// Check if two axes use the same binning convention, i.e.
   ///
   /// - Either they are both growable or neither of them is growable.
   /// - Minimum, maximum, and all bin borders in the middle are the same.
   /// - Bin labels must match (exactly including order, for now).
   bool HasSameBinningAs(const RAxisBase& other) const;

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   virtual int GetBinIndexForLowEdge(double x) const noexcept = 0;

private:
   std::string fTitle;    ///< Title of this axis, used for graphics / text.
};

///\name RAxisBase::const_iterator external operators
///\{

/// 2 + i
inline RAxisBase::const_iterator operator+(int d, RAxisBase::const_iterator rhs) noexcept
{
   return rhs + d;
}

/// i < j
inline bool operator<(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor < rhs.fCursor;
}

/// i > j
inline bool operator>(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor > rhs.fCursor;
}

/// i <= j
inline bool operator<=(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor <= rhs.fCursor;
}

/// i >= j
inline bool operator>=(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor >= rhs.fCursor;
}

/// i == j
inline bool operator==(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor == rhs.fCursor;
}

/// i != j
inline bool operator!=(RAxisBase::const_iterator lhs, RAxisBase::const_iterator rhs) noexcept
{
   return lhs.fCursor != rhs.fCursor;
}
///\}

/**
 Axis with equidistant bin borders. Defined by lower l and upper u limit and
 the number of bins n. All bins have the same width (u-l)/n.

 This axis cannot grow; use `RAxisGrow` for that.
 */
class RAxisEquidistant: public RAxisBase {
protected:
   double fLow = 0.;          ///< The lower limit of the axis
   double fInvBinWidth = 0.;  ///< The inverse of the bin width
   unsigned int fNBinsNoOver; ///< Number of bins excluding under- and overflow.

   /// Determine the inverse bin width.
   /// \param nbinsNoOver - number of bins without unter-/overflow
   /// \param lowOrHigh - first axis boundary
   /// \param highOrLow - second axis boundary
   static double GetInvBinWidth(int nbinsNoOver, double lowOrHigh, double highOrLow)
   {
      return nbinsNoOver / std::fabs(highOrLow - lowOrHigh);
   }

   /// See RAxisBase::HasSameBinBordersAs
   bool HasSameBinBordersAs(const RAxisBase& other) const override;

   /// Find the raw bin index (not adjusted) for the given coordinate.
   /// The resulting raw bin is 0-based.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   double FindBinRaw(double x) const noexcept
   {
      return (x - fLow) * fInvBinWidth;
   }

public:
   RAxisEquidistant() = default;

   /// Initialize a RAxisEquidistant.
   /// \param[in] title - axis title used for graphics and text representation.
   /// \param nbinsNoOver - number of bins in the axis, excluding under- and overflow
   ///   bins.
   /// \param low - the low axis range. Any coordinate below that is considered
   ///   as underflow. The first bin's lower edge is at this value.
   /// \param high - the high axis range. Any coordinate above that is considered
   ///   as overflow. The last bin's higher edge is at this value.
   explicit RAxisEquidistant(std::string_view title, int nbinsNoOver, double low, double high) noexcept
      : RAxisBase(title)
      , fLow(low)
      , fInvBinWidth(GetInvBinWidth(nbinsNoOver, low, high))
      , fNBinsNoOver(nbinsNoOver)
   {}

   /// Initialize a RAxisEquidistant.
   /// \param nbinsNoOver - number of bins in the axis, excluding under- and overflow
   ///   bins.
   /// \param low - the low axis range. Any coordinate below that is considered
   ///   as underflow. The first bin's lower edge is at this value.
   /// \param high - the high axis range. Any coordinate above that is considered
   ///   as overflow. The last bin's higher edge is at this value.
   explicit RAxisEquidistant(int nbinsNoOver, double low, double high) noexcept
      : RAxisEquidistant("", nbinsNoOver, low, high)
   {}

   /// Convert to RAxisConfig.
   operator RAxisConfig() const { return RAxisConfig(GetTitle(), GetNBinsNoOver(), GetMinimum(), GetMaximum()); }

   /// Get the number of bins, excluding under- and overflow.
   int GetNBinsNoOver() const noexcept final override { return fNBinsNoOver; }

   /// Find the adjusted bin index (returning `kUnderflowBin` for underflow and
   /// `kOverflowBin` for overflow) for the given coordinate.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   int FindBin(double x) const noexcept final override
   {
      double rawbin = FindBinRaw(x);
      return AdjustOverflowBinNumber(rawbin);
   }

   /// This axis cannot grow.
   bool CanGrow() const noexcept override { return false; }

   /// Get the width of the bins.
   double GetBinWidth() const noexcept { return 1. / fInvBinWidth; }

   /// Get the inverse of the width of the bins.
   double GetInverseBinWidth() const noexcept { return fInvBinWidth; }

   /// Get the bin center for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.25.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(int bin) const final override { return fLow + (bin - GetFirstBin() + 0.5) / fInvBinWidth; }

   /// Get the low bin border for the given bin index.
   /// For the bin == 1 (the first bin) of 2 bins for an axis (0., 1.), this
   /// returns 0.
   /// The result of this method on an underflow bin is unspecified.
   double GetBinFrom(int bin) const final override {
      const double result = (bin == kOverflowBin) ? GetMaximum() : fLow + (bin - GetFirstBin()) / fInvBinWidth;
      return result;
   }

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   int GetBinIndexForLowEdge(double x) const noexcept final override;
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kEquidistant> {
   using Axis_t = RAxisEquidistant;

   Axis_t operator()(const RAxisConfig &cfg) noexcept
   {
      return RAxisEquidistant(cfg.GetTitle(), cfg.GetNBinsNoOver(), cfg.GetBinBorders()[0], cfg.GetBinBorders()[1]);
   }
};

} // namespace Internal

/** An axis that can extend its range, keeping the number of its bins unchanged.
 The axis is constructed with an initial range. Apart from its ability to
 grow, this axis behaves like a RAxisEquidistant.
 */
class RAxisGrow: public RAxisEquidistant {
public:
   /// Initialize a RAxisGrow.
   /// \param[in] title - axis title used for graphics and text representation.
   /// \param nbins - number of bins in the axis, excluding under- and overflow
   ///   bins. This value is fixed over the lifetime of the object.
   /// \param low - the initial value for the low axis range. Any coordinate
   ///   below that is considered as underflow. To trigger the growing of the
   ///   axis call `Grow()`.
   /// \param high - the initial value for the high axis range. Any coordinate
   ///   above that is considered as overflow. To trigger the growing of the
   ///   axis call `Grow()`.
   explicit RAxisGrow(std::string_view title, int nbins, double low, double high) noexcept
      : RAxisEquidistant(title, nbins, low, high)
   {}

   /// Initialize a RAxisGrow.
   /// \param nbins - number of bins in the axis, excluding under- and overflow
   ///   bins. This value is fixed over the lifetime of the object.
   /// \param low - the initial value for the low axis range. Any coordinate
   ///   below that is considered as underflow. To trigger the growing of the
   ///   axis call `Grow()`.
   /// \param high - the initial value for the high axis range. Any coordinate
   ///   above that is considered as overflow. To trigger the growing of the
   ///   axis call `Grow()`.
   explicit RAxisGrow(int nbins, double low, double high) noexcept: RAxisGrow("", nbins, low, high) {}

   /// Convert to RAxisConfig.
   operator RAxisConfig() const { return RAxisConfig(GetTitle(), RAxisConfig::Grow, GetNBinsNoOver(), GetMinimum(), GetMaximum()); }

   /// Grow this axis to make the "virtual bin" toBin in-range. This keeps the
   /// non-affected axis limit unchanged, and extends the other axis limit such
   /// that a number of consecutive bins are merged.
   ///
   /// Example, assuming an initial RAxisGrow with 10 bins from 0. to 1.:
   ///   - `Grow(0)`: that (virtual) bin spans from -0.1 to 0. To include it
   ///     in the axis range, the lower limit must be shifted. The minimal number
   ///     of bins that can be merged is 2, thus the new axis will span from
   ///     -1. to 1.
   ///   - `Grow(-1)`: that (virtual) bin spans from -0.2 to 0.1. To include it
   ///     in the axis range, the lower limit must be shifted. The minimal number
   ///     of bins that can be merged is 2, thus the new axis will span from
   ///     -1. to 1.
   ///   - `Grow(50)`: that (virtual) bin spans from 4.9 to 5.0. To include it
   ///     in the axis range, the higher limit must be shifted. Five bins need to
   ///     be merged, making the new axis range 0. to 5.0.
   ///
   /// \param toBin - the "virtual" bin number, as if the axis had an infinite
   ///   number of bins with the current bin width. For instance, for an axis
   ///   with ten bins in the range 0. to 1., the coordinate 2.05 has the virtual
   ///   bin index 20.
   /// \return Returns the number of bins that were merged to reach the value.
   ///   A value of 1 means that no bins were merged (toBin was in the original
   ///   axis range).
   int Grow(int toBin);

   /// This axis kind can increase its range.
   bool CanGrow() const noexcept final override { return true; }
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kGrow> {
   using Axis_t = RAxisGrow;

   Axis_t operator()(const RAxisConfig &cfg) noexcept
   {
      return RAxisGrow(cfg.GetTitle(), cfg.GetNBinsNoOver(), cfg.GetBinBorders()[0], cfg.GetBinBorders()[1]);
   }
};

} // namespace Internal

/**
 An axis with non-equidistant bins (also known as "variable binning"). It is
 defined by an array of bin borders - one more than the number of
 (non-overflow-) bins it has! As an example, an axis with two bin needs three
 bin borders:
   - lower edge of the first bin;
   - higher edge of the first bin, identical to the lower edge of the second
     bin;
   - higher edge of the second bin

 This axis cannot grow; the size of new bins would not be well defined.
 */
class RAxisIrregular: public RAxisBase {
private:
   /// Bin borders, one more than the number of regular bins.
   std::vector<double> fBinBorders;

protected:
   /// See RAxisBase::HasSameBinBordersAs
   bool HasSameBinBordersAs(const RAxisBase& other) const override;

   /// Find the raw bin index (not adjusted) for the given coordinate `x`.
   /// The resulting raw bin is 1-based.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   double FindBinRaw(double x) const noexcept
   {
      const auto bBegin = fBinBorders.begin();
      const auto bEnd = fBinBorders.end();
      // lower_bound finds the first bin border that is >= x.
      auto iNotLess = std::lower_bound(bBegin, bEnd, x);
      return iNotLess - bBegin;
   }

public:
   RAxisIrregular() = default;

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   explicit RAxisIrregular(const std::vector<double> &binborders)
      : RAxisBase(), fBinBorders(binborders)
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   /// Faster, noexcept version taking an rvalue of binborders. The compiler will
   /// know when it can take this one.
   explicit RAxisIrregular(std::vector<double> &&binborders) noexcept
      : RAxisBase(), fBinBorders(std::move(binborders))
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   explicit RAxisIrregular(std::string_view title, const std::vector<double> &binborders)
      : RAxisBase(title), fBinBorders(binborders)
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Construct a RAxisIrregular from a vector of bin borders.
   /// \note The bin borders must be sorted in increasing order!
   /// Faster, noexcept version taking an rvalue of binborders. The compiler will
   /// know when it can take this one.
   explicit RAxisIrregular(std::string_view title, std::vector<double> &&binborders) noexcept
      : RAxisBase(title), fBinBorders(std::move(binborders))
   {
#ifdef R__DO_RANGE_CHECKS
      if (!std::is_sorted(fBinBorders.begin(), fBinBorders.end()))
         R__LOG_ERROR("HIST") << "Bin borders must be sorted!";
#endif // R__DO_RANGE_CHECKS
   }

   /// Convert to RAxisConfig.
   operator RAxisConfig() const { return RAxisConfig(GetTitle(), GetBinBorders()); }

   /// Get the number of bins, excluding under- and overflow.
   int GetNBinsNoOver() const noexcept final override { return fBinBorders.size() - 1; }

   /// Find the bin index (adjusted with under- and overflow) for the given coordinate `x`.
   /// \note Passing a bin border coordinate can either return the bin above or
   /// below the bin border. I.e. don't do that for reliable results!
   int FindBin(double x) const noexcept final override
   {
      int rawbin = FindBinRaw(x);
      // No need for AdjustOverflowBinNumber(rawbin) here; lower_bound() is the
      // answer: e.g. for x < *bBegin, rawbin is -1.
      if (rawbin < GetFirstBin())
         return kUnderflowBin;
      if (rawbin >= GetLastBin() + 1)
         return kOverflowBin;
      return rawbin;
   }

   /// Get the bin center of the bin with the given index.
   /// The result of this method on an overflow or underflow bin is unspecified.
   double GetBinCenter(int bin) const final override { return 0.5 * (fBinBorders[bin - 1] + fBinBorders[bin]); }

   /// Get the lower bin border for a given bin index.
   /// The result of this method on an underflow bin is unspecified.
   double GetBinFrom(int bin) const final override
   {
      if (bin == kOverflowBin)
         return fBinBorders[GetLastBin()];
      return fBinBorders[bin - 1];
   }

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   int GetBinIndexForLowEdge(double x) const noexcept final override;

   /// This axis cannot be extended.
   bool CanGrow() const noexcept final override { return false; }

   /// Access to the bin borders used by this axis.
   const std::vector<double> &GetBinBorders() const noexcept { return fBinBorders; }
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kIrregular> {
   using Axis_t = RAxisIrregular;

   Axis_t operator()(const RAxisConfig &cfg) { return RAxisIrregular(cfg.GetTitle(), cfg.GetBinBorders()); }
};

} // namespace Internal

/**
 \class RAxisLabels
 A RAxisGrow that has a label assigned to each bin and a bin width of 1.

 While filling still works through coordinates (i.e. arrays of doubles),
 RAxisLabels allows to convert a string to a bin number or the bin's coordinate
 center. The number of labels and the number of bins reported by RAxisGrow might
 differ: the RAxisGrow will only grow when seeing a Fill(), while the RAxisLabels
 will add a new label whenever `GetBinCenter()` is called.

 Implementation details:
 Filling happens often; `GetBinCenter()` needs to be fast. Thus the unordered_map.
 The painter needs the reverse: it wants the label for bin 0, bin 1 etc. The axis
 should only store the bin labels once; referencing them is (due to re-allocation,
 hashing etc) non-trivial. So instead, build a `vector<string_view>` for the few
 times the axis needs to be painted.
 */
class RAxisLabels: public RAxisGrow {
private:
   /// Map of label (view on `fLabels`'s elements) to bin index
   std::unordered_map<std::string, int /*bin number*/> fLabelsIndex;

public:
   /// Construct a RAxisLables from a `vector` of `string_view`s, with title.
   explicit RAxisLabels(std::string_view title, const std::vector<std::string_view> &labels)
      : RAxisGrow(title, labels.size(), 0., static_cast<double>(labels.size()))
   {
      for (size_t i = 0, n = labels.size(); i < n; ++i)
         fLabelsIndex[std::string(labels[i])] = i;
   }

   /// Construct a RAxisLables from a `vector` of `string`s, with title.
   explicit RAxisLabels(std::string_view title, const std::vector<std::string> &labels)
      : RAxisGrow(title, labels.size(), 0., static_cast<double>(labels.size()))
   {
      for (size_t i = 0, n = labels.size(); i < n; ++i)
         fLabelsIndex[labels[i]] = i;
   }

   /// Construct a RAxisLables from a `vector` of `string_view`s
   explicit RAxisLabels(const std::vector<std::string_view> &labels): RAxisLabels("", labels) {}

   /// Construct a RAxisLables from a `vector` of `string`s
   explicit RAxisLabels(const std::vector<std::string> &labels): RAxisLabels("", labels) {}

   /// Convert to RAxisConfig.
   operator RAxisConfig() const { return RAxisConfig(GetTitle(), GetBinLabels()); }

   /// Get the bin index with label.
   int FindBinByName(const std::string &label)
   {
      auto insertResult = fLabelsIndex.insert({label, -1});
      if (insertResult.second) {
         // we have created a new label
         int idx = fLabelsIndex.size() - 1;
         insertResult.first->second = idx;
         return idx;
      }
      return insertResult.first->second;
   }

   /// Get the center of the bin with label.
   double GetBinCenterByName(const std::string &label)
   {
      return FindBinByName(label) + 0.5; // bin *center*
   }

   /// Build a vector of labels. The position in the vector defines the label's bin.
   std::vector<std::string_view> GetBinLabels() const
   {
      std::vector<std::string_view> vec(fLabelsIndex.size());
      for (const auto &kv: fLabelsIndex)
         vec.at(kv.second) = kv.first;
      return vec;
   }

   /// Result of an RAxisLabels label set comparison
   enum LabelsCmpFlags {
      /// Both axes have the same labels, mapping to the same bins
      kLabelsCmpSame = 0,

      /// The other axis doesn't have some labels from this axis
      kLabelsCmpSubset = 0b1,

      /// The other axis has some labels which this axis doesn't have
      kLabelsCmpSuperset = 0b10,

      /// The labels shared by both axes do not map into the same bins
      kLabelsCmpDisordered = 0b100,
   };

   /// Compare the labels of this axis with those of another axis
   LabelsCmpFlags CompareBinLabels(const RAxisLabels& other) const noexcept {
      // This will eventually contain the results of the labels comparison
      LabelsCmpFlags result = kLabelsCmpSame;
      size_t missing_in_other = 0;

      // First, check how this axis' labels map into the other axis
      for (const auto &kv: fLabelsIndex) {
         auto iter = other.fLabelsIndex.find(kv.first);
         if (iter == other.fLabelsIndex.cend()) {
            ++missing_in_other;
         } else if (iter->second != kv.second) {
            result = LabelsCmpFlags(result | kLabelsCmpDisordered);
         }
      }
      if (missing_in_other > 0)
         result = LabelsCmpFlags(result | kLabelsCmpSubset);

      // If this covered all labels in the other axis, we're done
      if (fLabelsIndex.size() == other.fLabelsIndex.size() + missing_in_other)
         return result;

      // Otherwise, we must check the labels of the other axis too
      for (const auto &kv: other.fLabelsIndex)
         if (fLabelsIndex.find(kv.first) == fLabelsIndex.cend())
            return LabelsCmpFlags(result | kLabelsCmpSuperset);
      return result;
   }
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kLabels> {
   using Axis_t = RAxisLabels;

   Axis_t operator()(const RAxisConfig &cfg) { return RAxisLabels(cfg.GetTitle(), cfg.GetBinLabels()); }
};

} // namespace Internal

///\name Axis Compatibility
///\{
enum class EAxisCompatibility {
   kIdentical, ///< Source and target axes are identical

   kContains, ///< The source is a subset of bins of the target axis

   /// The bins of the source axis have finer granularity, but the bin borders
   /// are compatible. Example:
   /// source: 0., 1., 2., 3., 4., 5., 6.; target: 0., 2., 5., 6.
   /// Note that this is *not* a symmetrical property: only one of
   /// CanMerge(source, target), CanMap(target, source) can return kContains.
   kSampling,

   /// The source axis and target axis have different binning. Example:
   /// source: 0., 1., 2., 3., 4., target: 0., 0.1, 0.2, 0.3, 0.4
   kIncompatible
};

/// Whether (and how) the source axis can be merged into the target axis.
EAxisCompatibility CanMap(const RAxisEquidistant &target, const RAxisEquidistant &source) noexcept;
///\}

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAxis header guard
