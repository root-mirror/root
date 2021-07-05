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
#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>
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

   /// Placeholder bin width to be used in ComparePosToBinBorder when there is
   /// no bin or the width of that bin is irrelevant.
   constexpr static const double kNoBinWidth = -1.;

   /// Compare an axis position with an axis bin border
   ///
   /// Given a position of interest on the axis, a bin border position, and the
   /// bin width on both sides of said bin border, tell if the position of
   /// interest should be considered to be located before (-1), at the same
   /// position (0), or after (+1) the bin border of interest.
   ///
   /// If there is no bin on one side of the target axis bin border, if it is an
   /// under/overflow bin, or if you do not care if the result is zero or
   /// nonzero on that side of the target bin border, please set the
   /// corresponding bin width to kNoBinWidth.
   ///
   static int ComparePosToBinBorder(double x,
                                    double binBorder,
                                    double leftBinWidth,
                                    double rightBinWidth) {
      // Current tolerance policy when there is no bin width on one side is to
      // treat the unknown bin width as a bin width of 1.
      if (leftBinWidth == kNoBinWidth) leftBinWidth = 1.;
      if (rightBinWidth == kNoBinWidth) rightBinWidth = 1.;

      // Perform an approximate bin border comparison
      const double offset = x - binBorder;
      const double tolerance = 1e-6;
      if (offset < 0.) {
         return -static_cast<int>(offset < -leftBinWidth*tolerance);
      } else {
         return static_cast<int>(offset > rightBinWidth*tolerance);
      }
   }

   /// Enum for specifying a side of an axis bin
   enum class BinSide {
      /// Left side of a bin
      kFrom,

      /// Right side of a bin
      kTo,
   };

   /// Compare an axis position with a bin border of this axis
   ///
   /// Given a position of interest, a bin number, and a BinSide indicating
   /// which side of the bin we're talking about, tell if the position of
   /// interest is located strictly before (-1), at about the same position (0)
   /// or strictly after (+1) the bin border of interest.
   ///
   /// This is a higher-level alternative to the static method above.
   ///
   int ComparePosToBinBorder(double x, int bin, BinSide side) const noexcept {
      // Handle underflow bin edge case
      if (bin == kUnderflowBin) {
         assert(!CanGrow());
         switch (side) {
            case BinSide::kFrom:
               // Everything is after the underflow bin's "left edge"
               return 1;
            case BinSide::kTo:
               return ComparePosToBinBorder(x, GetFirstBin(), BinSide::kFrom);
         }
      }

      // Handle overflow bin edge case
      if (bin == kOverflowBin) {
         assert(!CanGrow());
         switch (side) {
            case BinSide::kFrom:
               return ComparePosToBinBorder(x, GetLastBin(), BinSide::kTo);
            case BinSide::kTo:
               // Everything is before the overflow bin's "right edge"
               return -1;
         }
      }

      // Get regular bin border comparison parameters
      double borderPos = 0., leftBinWidth = 0., rightBinWidth = 0.;
      switch (side) {
         case BinSide::kFrom:
            borderPos = GetBinFrom(bin);
            leftBinWidth =
               (bin > GetFirstBin()) ? (borderPos - GetBinFrom(bin-1))
                                     : kNoBinWidth;
            rightBinWidth = GetBinTo(bin) - borderPos;
            break;
         case BinSide::kTo:
            borderPos = GetBinTo(bin);
            leftBinWidth = borderPos - GetBinFrom(bin);
            rightBinWidth =
               (bin < GetLastBin()) ? (GetBinTo(bin+1) - borderPos)
                                    : kNoBinWidth;
      }

      // Perform an approximate bin border comparison
      return ComparePosToBinBorder(x, borderPos, leftBinWidth, rightBinWidth);
   }

public:
   /// Result of comparing two axes with numerical bin borders for histogram
   /// merging
   class NumericBinningCompatibility {
   public:
      /// Flags representing various properties that can emerge from the
      /// comparison of two numerical axes (see methods for a full description)
      enum Flags {
         // The mapping from source to target regular bin indices is trivial
         kTrivialRegularBinMapping = 1 << 0,

         // The mapping between source and target regular bins is bijective
         kRegularBinBijection = 1 << 1,

         // The mapping between all source and target bins is bijective
         kFullBinBijection = 1 << 2,

         // Some bins from the source map to target bins that span extra range
         kMergingIsLossy = 1 << 3,

         // Some regular bins from the source axis map to multiple target bins
         kRegularBinAliasing = 1 << 4,

         // The source underflow bin must be empty to allow histogram merging
         kNeedEmptyUnderflow = 1 << 5,

         // The source overflow bin must be empty to allow histogram merging
         kNeedEmptyOverflow = 1 << 6,

         // The target axis must grow in order to span the full source range
         kTargetMustGrow = 1 << 7,
      };

      /// Build a numerical binning comparison result
      ///
      /// See the methods of this class for a more detailed description of what
      /// each of these flags mean.
      ///
      NumericBinningCompatibility(Flags flags)
         : fFlags(flags)
      {}

      // Check against another comparison result
      bool operator==(const NumericBinningCompatibility& other) const {
         return fFlags == other.fFlags;
      }

      /// Truth that there is a trivial mapping from the indices of the source
      /// axis's regular bins to those of the target axis
      ///
      /// If this property is true, then every regular source axis bin maps into
      /// a target axis bin which has the same bin index.
      ///
      /// For source axis bins which map into multiple target axis bins (see
      /// `HasRegularBinAliasing()`), the target bin which has the same index as
      /// the source bin must be the _first_ matching bin in target axis order.
      ///
      // NOTE: This property can be leveraged to avoid local bin index
      //       conversions in an histogram merging implementation.
      //
      bool HasTrivialRegularBinMapping() const {
         return fFlags & kTrivialRegularBinMapping;
      }

      /// Truth that there is a bijective mapping between the regular bins
      /// of the source and target axis
      ///
      /// This property, which implies `HasTrivialRegularBinMapping()`,
      /// indicates that every regular source axis bin maps into a regular
      /// target axis bin with the same index and vice versa.
      ///
      /// If this property (or its `!TargetWillHaveExtraBins()` cousin for
      /// labeled axes) is true for every dimension of a source and target
      /// histogram, then every regular source histogram bin maps into a target
      /// histogram bin with the same global bin index.
      ///
      // NOTE: This property can be leveraged to avoid local<->global regular
      //       bin index conversions in the histogram merging implementation.
      //
      bool HasRegularBinBijection() const {
         if (fFlags & kRegularBinBijection)
            assert(HasTrivialRegularBinMapping());
         return fFlags & kRegularBinBijection;
      }

      /// Truth that there is a bijective mapping between all bins of the source
      /// and target axis, both regular and under/overflow
      ///
      /// This property, which implies `HasBijectiveRegularBinMapping()`,
      /// indicates that every bin of the source axis, regular or
      /// under/overflow, maps into a target axis bin with the same index and
      /// vice versa.
      ///
      /// The mapping of the source overflow bin to the target overflow bin
      /// differs from the definition of `HasTrivialRegularBinMapping()` in that
      /// the target overflow bin is guaranteed to be the _last_ target axis bin
      /// which the source overflow bin maps into, not the first.
      ///
      /// If this property (or its `!TargetWillHaveExtraBins()` cousin for
      /// labeled axes) is true for every dimension of a source and target
      /// histogram, then every source histogram bin, whether regular or
      /// under/overflow, maps into a target histogram bin with the same global
      /// bin index.
      ///
      // NOTE: This property allows applying the same sort of histogram merging
      //       optimizations as RegularBinBijection(), but for overflow bins
      //       too in addition to regular bins.
      //
      //       Regular bin bijection is treated as a prerequisite of overflow
      //       bin bijection because in a multi-dimensional histogram, regular
      //       bins of one axis will be part of the under/overflow hyperplane of
      //       other axes. That condition is slightly more pessimistic than what
      //       we actually need, but easier to check and good enough for the
      //       common case of merging two histograms with identical axis config.
      //
      bool HasFullBinBijection() const {
         if (fFlags & kFullBinBijection)
            assert(HasRegularBinBijection());
         return fFlags & kFullBinBijection;
      }

      /// Truth that some bins from the source axis map into target axis bins
      /// that span some extra axis range
      ///
      /// From the perspective of histogram merging, this property means that
      /// some information about the location of previous histogram fills may be
      /// lost in the histogram merging process.
      ///
      /// If that happens, the histogram merge is irreversible: trying to merge
      /// data back from the target histogram to another histogram which has the
      /// same binning as the source histogram will lead to bin aliasing issues
      /// (see below).
      ///
      // NOTE: This property does not affect the histogram merging
      //       implementation, but should be reported as a warning to its user.
      //
      // NOTE: This does _not_ imply that the target bins are bigger, for a
      //       counter-example this flag would be set in the example below:
      //
      //           Source axis bins:   |---|
      //           Target axis bins: |---|---|
      //
      bool MergingIsLossy() const {
         return fFlags & kMergingIsLossy;
      }

      /// Truth that some regular bins from the source axis map into multiple
      /// bins (regular or under/overflow) on the target axis
      ///
      /// These bins must be empty for histogram merging to be possible, because
      /// if they had some content, no automated histogram merging routine would
      /// be able to correctly split this content across matching target
      /// histogram bins. The required information about how fills were
      /// distributed inside of the source bins simply isn't there.
      ///
      // NOTE: If this property is true, then for each source axis bin, the
      //       histogram merging implementation must locate a matching target
      //       axis bin (possibly using the optimizations permitted by
      //       `HasTrivialRegularBinMapping()` and `HasRegularBinBijection()`),
      //       and check if its borders match those of the source. If not, the
      //       source bin in question must be empty or else the histogram
      //       merging implementation will have to error out.
      //
      bool HasRegularBinAliasing() const {
         return fFlags & kRegularBinAliasing;
      }

      /// Truth that the source axis' underflow bin(s) must be empty for
      /// histogram merging to be possible
      ///
      /// This property may only be true if the source axis has underflow bins.
      ///
      bool MergingNeedsEmptyUnderflow() const {
         return fFlags & kNeedEmptyUnderflow;
      }

      /// Truth that the source axis' overflow bin(s) must be empty for
      /// histogram merging to be possible
      ///
      /// This property may only be true if the source axis has overflow bins.
      ///
      bool MergingNeedsEmptyOverflow() const {
         return fFlags & kNeedEmptyOverflow;
      }

      /// Truth that the target axis must grow to span the source axis range
      ///
      /// This property may only be true if the target axis is growable.
      ///
      /// All other properties are computed under the assumption that such
      /// growth has indeed occurred.
      ///
      // NOTE: Such growth may not actually be necessary if the source axis bins
      //       which are not covered by the target axis are empty. However,
      //       figuring this out requires a scan of the source histogram bins,
      //       so it's dubious whether it's a good idea to check for it from a
      //       performance point of view.
      //
      bool MergingNeedsTargetGrowth() const {
         return fFlags & kTargetMustGrow;
      }

   private:
      Flags fFlags;
   };

protected:
   /// Compare the numerical bin borders of one axis with that of another for
   /// the purpose of evaluating an histogram merging scenario
   ///
   /// Since histogram merging has asymmetric properties, this axis is the
   /// target axis, and the other axis is the source axis.
   ///
   /// Growable RAxis subclasses **must** override this method in a manner which
   /// evaluates possible axis growth scenarios, before calling back
   /// `CheckFixedNumericalBinningCompat()`, which implements the actual
   /// binning comparison, on a simulated grown axis.
   ///
   // FIXME: As of 2020-05-27, every RAxis type is considered to have numerical
   //        bin borders, but that doesn't make sense for RAxisLabels. This will
   //        be fixed by adding a lower-level base class which doesn't assume
   //        the existence of floating-point bin borders, which will also lay
   //        out the groundwork for supporting boolean/integer bin borders.
   //
   virtual NumericBinningCompatibility
   CheckNumericBinningCompat(const RAxisBase& source) const {
      assert(!CanGrow());
      return CheckFixedNumericBinningCompat(source, false);
   }

   /// Callback of `CheckNumericalBinningCompat()` containing the actual bin
   /// border comparison logic, to be invoked after simulating any required axis
   /// growth scenarios, on the grown target axis.
   ///
   /// This method can be overriden for performance optimization purposes.
   ///
   virtual NumericBinningCompatibility
   CheckFixedNumericBinningCompat(const RAxisBase& source,
                                  bool growthOccured) const;

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

   /// If the coordinate `x` is within 10 ULPs of a bin low edge coordinate,
   /// return the bin for which this is a low edge. If it's not a bin edge,
   /// return `kInvalidBin`.
   virtual int GetBinIndexForLowEdge(double x) const noexcept = 0;

   /// Result of comparing two labeled axis for histogram merging
   class LabeledBinningCompatibility {
   public:
      /// Flags representing various properties that can emerge from the
      /// comparison of two labeled axes (see methods for a full description)
      enum Flags {
         // The target axis must grow to cover all source axis labels w/ bins
         kTargetMustGrow = 1 << 0,

         // Need target bin reordering or an order-insensitive merge algorithm
         kLabelOrderDiffers = 1 << 1,

         // After growing, the target axis will have more bins than the source
         kExtraTargetBins = 1 << 2,
      };

      /// Build a labeled axis comparison result
      ///
      /// See the methods of this class for a more detailed description of what
      /// each of these flags mean.
      ///
      LabeledBinningCompatibility(Flags flags,
                                  std::vector<std::string_view>&& sourceOnlyLabels)
         : fFlags(flags)
         , fSourceOnlyLabels(std::move(sourceOnlyLabels))
      {}

      // Check against another comparison result
      bool operator==(const LabeledBinningCompatibility& other) const {
         return (fFlags == other.fFlags)
                && (fSourceOnlyLabels == other.fSourceOnlyLabels);
      }

      /// List of source axis labels which the target axis does not have
      ///
      /// These labels will need to be added to the target axis before histogram
      /// data merging becomes possible.
      ///
      /// Labels are listed in the order in which they appear in the source
      /// axis, and should be added to the target axis in this order, both for
      /// good user experience and for other `LabeledBinningCompatibility`
      /// properties to be usable.
      ///
      /// Labels with have no associated bins are not accounted for, since they
      /// contain no useful data to be transferred to the target axis.
      ///
      /// This property implies `TargetMustGrow()`.
      ///
      const std::vector<std::string_view>& SourceOnlyLabels() const {
         if (!fSourceOnlyLabels.empty()) {
            assert(TargetMustGrow());
         }
         return fSourceOnlyLabels;
      }

      /// The target axis must grow to encompass all source bins
      ///
      /// This can happen either because more source labels need to be added
      /// to the target axis (see above) or because the target axis already has
      /// the proper labels, but they don't have associated bins yet.
      ///
      /// Target axis growth should be carried out after adding missing labels.
      ///
      bool TargetMustGrow() const {
         return fFlags & kTargetMustGrow;
      }

      /// Source and target labels are not (or won't be) in the same order
      ///
      /// After adding any missing bins and labels to the target axis, the
      /// histogram merging algorithm will need to either reorder the target
      /// axis' labels (and associated histogram bins) or support merge
      /// scenarios in which source->target bin correspondence is nontrivial.
      ///
      /// If needed, bin reordering must be carried out after adding any missing
      /// labels and growing the target axis.
      ///
      bool LabelOrderDiffers() const {
         return fFlags & kLabelOrderDiffers;
      }

      /// After growing, the target axis will have more bins than the source
      ///
      /// This means that even if `LabelOrderDiffers()` is false or bin
      /// reordering is applied, global bin indices will not match between the
      /// source and target histogram, and local<->global bin index conversions
      /// will be necessary during histogram merging.
      ///
      /// This is the moral equivalent of
      /// `NumericalBinCompatibility::FullBinBijection()` for labeled axes.
      ///
      // NOTE: Although this technically applies to regular bins only, there is
      //       no need for a regular/full bin bijection distinction here because
      //       labeled axes may only be successfully compared with other labeled
      //       axes, and therefore the number of under/overflow bins is
      //       guaranteed to be the same on both sides.
      //
      bool TargetWillHaveExtraBins() const {
         return fFlags & kExtraTargetBins;
      }

   private:
      Flags fFlags;
      std::vector<std::string_view> fSourceOnlyLabels;
   };

   /// Result of comparing two axes for histogram merging
   //
   // TODO: Replace with std::variant once RHist goes C++17
   //
   class BinningCompatibility {
   public:
      /// Case where two axes of incompatible types were compared
      BinningCompatibility()
         : fKind(CompatKind::kIncompatible)
      {}

      /// Case where two axes using numerical bin borders were compared
      ///
      /// See the methods of this class for a more detailed description of what
      /// each of these flags mean.
      ///
      explicit BinningCompatibility(NumericBinningCompatibility numeric)
         : fKind(CompatKind::kNumeric)
         , fNumeric(numeric)
      {}

      /// Case where two RAxisLabels were compared
      ///
      /// See the methods of this class for a more detailed description of what
      /// each of these flags mean.
      ///
      explicit BinningCompatibility(LabeledBinningCompatibility labeled)
         : fKind(CompatKind::kLabeled)
         , fLabeled(labeled)
      {}

      // Handle copies and moves
      BinningCompatibility(const BinningCompatibility& other) {
         *this = other;
      }
      BinningCompatibility(BinningCompatibility&& other) {
         *this = std::move(other);
      }
      BinningCompatibility& operator=(const BinningCompatibility& other) {
         switch (other.fKind) {
            case CompatKind::kIncompatible:
               rebuild();
               break;

            case CompatKind::kNumeric:
               rebuild(other.fNumeric);
               break;

            case CompatKind::kLabeled:
               rebuild(other.fLabeled);
               break;
         }
         return *this;
      }
      BinningCompatibility& operator=(BinningCompatibility&& other) {
         switch (other.fKind) {
            case CompatKind::kIncompatible:
               rebuild();
               break;

            case CompatKind::kNumeric:
               rebuild(std::move(other.fNumeric));
               break;

            case CompatKind::kLabeled:
               rebuild(std::move(other.fLabeled));
               break;
         }
         return *this;
      }

      /// Destroy any inner data on destruction
      ~BinningCompatibility() {
         switch (fKind) {
            case CompatKind::kIncompatible:
               // No inner data to be destroyed
               break;

            case CompatKind::kNumeric:
               fNumeric.~NumericBinningCompatibility();
               break;

            case CompatKind::kLabeled:
               fLabeled.~LabeledBinningCompatibility();
               break;
         }
      }

      /// Check if two comparisons led to the same result
      bool operator==(const BinningCompatibility& other) const {
         if (other.fKind != fKind) return false;
         switch (fKind) {
            case CompatKind::kIncompatible:
               return true;

            case CompatKind::kNumeric:
               return fNumeric == other.fNumeric;

            case CompatKind::kLabeled:
               return fLabeled == other.fLabeled;
         }
         return false;
      }

      /// Broad classification of possible axis comparisons
      enum class CompatKind {
         /// Two axes using a fundamentally incompatible binning scheme (e.g.
         /// `RAxisIrregular` vs `RAxisLabels`) were compared
         ///
         /// It is impossible to automatically merge two histograms when the
         /// axis types for some dimension differ so much.
         ///
         kIncompatible,

         /// Two axes using numerical bin borders (e.g. `RAxisIrregular` vs
         /// `RAxisEquidistant`) were compared
         ///
         /// In this case, bin borders are compared up to a certain tolerance,
         /// and differences are only reported when the source axis bin borders
         /// are not within that tolerance of the closest target bin border.
         ///
         /// The ability to automatically merge two histograms with numerical
         /// bin borders depends on many factors, and may not decidable without
         /// looking at the details of the source histogram's bin contents.
         ///
         kNumeric,

         /// Two `RAxisLabels` were compared
         ///
         /// It is always possible to automatically merge two histograms if all
         /// of their axes use labeled bins.
         ///
         kLabeled,
      };

      /// Kind of axis comparison that was carried out
      CompatKind Kind() const noexcept { return fKind; }

      /// Get the detailed result of a numerical axis comparison
      const NumericBinningCompatibility& GetNumeric() const {
         CheckKind(CompatKind::kNumeric);
         return fNumeric;
      }

      /// Get the detailed result of a labeled axis comparison
      const LabeledBinningCompatibility& GetLabeled() const {
         CheckKind(CompatKind::kLabeled);
         return fLabeled;
      }

   private:
      /// Kind of axis comparison that was carried out
      CompatKind fKind;

      /// Check that the axis comparison kind is correct in preparation to
      /// querying a binning property which is specific to that axis kind
      void CheckKind(CompatKind expectedKind) const;

      /// Details of the axis comparison results
      union {
         // Valid if fKind is CompatKind::kNumeric
         NumericBinningCompatibility fNumeric;

         // Valid if fKind is CompatKind::kLabeled
         LabeledBinningCompatibility fLabeled;
      };

      /// Replace with another instance of BinningCompatibility, constructed
      /// using the provided arguments.
      template <typename... Args>
      void rebuild(Args&&... args) {
         this->~BinningCompatibility();
         new(this) BinningCompatibility(std::forward<Args>(args)...);
      }
   };

   /// Compare the binning of this axis with that of another axis for the
   /// purpose of evaluating an histogram merging scenario
   ///
   /// Since histogram merging has asymmetric properties, this axis is the
   /// target axis, and the other axis is the source axis.
   ///
   BinningCompatibility CheckBinningCompat(const RAxisBase& source) const;

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
   RAxisGrow() = default;

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

   /// CheckNumericalBinningCompat must be overriden to handle axis growth
   NumericBinningCompatibility
   CheckNumericBinningCompat(const RAxisBase& source) const final override;
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

   /// Compare the labels of this axis with those of another axis for the
   /// purpose of investigating a histogram merging scenario
   LabeledBinningCompatibility
   CheckLabeledBinningCompat(const RAxisLabels& source) const noexcept {
      // For each axis, we must carefully distinguish the number of bins from
      // the number of labels. An RAxisLabels may have more bin labels than it
      // has bins if a label has been queried (which automatically allocates a
      // bin index) but the corresponding bin never registered any fill.
      const int numSourceBins = source.GetNBinsNoOver();
      const int numTargetLabels = fLabelsIndex.size();
      const int numTargetBins = GetNBinsNoOver();

      // Collect the set of source bin labels in the order where the user sees
      // them. Ignore uncommitted source labels: no bin = no data to be merged.
      std::vector<std::string_view> sourceLabels = source.GetBinLabels();
      sourceLabels.resize(numSourceBins);

      // Check how source _bins_ map into target bins and labels, simulating
      // label and bin creation as necessary.
      int newTargetLabelIdx = numTargetLabels;
      int numBinsAfterGrowth = numTargetBins;
      bool labelOrderDiffers = false;
      for (int sourceLabelIdx = 0; sourceLabelIdx < numSourceBins; ++sourceLabelIdx) {
         // Look for the the source bin's label in the target axis' label set
         //
         // FIXME: Remove string allocation once ROOT goes C++20 and
         //        unordered_map supports searching a hashmap without
         //        materializing an owned key...
         //
         auto iter = fLabelsIndex.find(std::string(sourceLabels[sourceLabelIdx]));

         // Use the bin index of any existing target label, or simulate creating
         // a new label in the target axis.
         //
         // Along the way, keep track of which labels only exist in the source
         // axis, and thus need to be added to the target axis, by clearing out
         // the labels which _do_ exist in the target axis within the list of
         // source labels. Obviously, the source label must not be used after
         // this operation has been performed.
         //
         int targetLabelIdx;
         if (iter != fLabelsIndex.cend()) {
            targetLabelIdx = iter->second;
            sourceLabels[sourceLabelIdx] = std::string_view();
         } else {
            targetLabelIdx = newTargetLabelIdx++;
         }

         // If the target label has no associated bin yet, simulate growing
         // the target axis to materialize that bin
         numBinsAfterGrowth = std::max(numBinsAfterGrowth, targetLabelIdx+1);

         // Check if label order is consistent in the source and target axes
         labelOrderDiffers |= (targetLabelIdx != sourceLabelIdx);
      }

      // At this point, sourceLabels contains non-empty string views for labels
      // which are unique to the source axis, and empty string views for labels
      // which exist in both the source and target axis. From this, we can
      // easily get a list of labels which only exist in the source axis.
      const auto sourceOnlyEnd =
         std::remove(sourceLabels.begin(), sourceLabels.end(), std::string_view());
      sourceLabels.resize(std::distance(sourceLabels.begin(), sourceOnlyEnd));
      // From this point on, the sourceLabels vectors contains the set of labels
      // which only exists in the source axis, and not on the target axis.

      // Figure if after performing any missing source label addition, the
      // histogram merging implementation will need to grow the target axis.
      const bool targetMustGrow = (numBinsAfterGrowth > numTargetBins);

      // Figure out if after growing the target axis as needed, it will feature
      // some labels which the source axis doesn't have.
      bool extraTargetBins = (numBinsAfterGrowth > numSourceBins);

      // Produce the results of the comparison
      using Flags = LabeledBinningCompatibility::Flags;
      return LabeledBinningCompatibility(
         Flags(targetMustGrow * Flags::kTargetMustGrow
               + labelOrderDiffers * Flags::kLabelOrderDiffers
               + extraTargetBins * Flags::kExtraTargetBins),
         std::move(sourceLabels)
      );
   }
};

namespace Internal {

template <>
struct AxisConfigToType<RAxisConfig::kLabels> {
   using Axis_t = RAxisLabels;

   Axis_t operator()(const RAxisConfig &cfg) { return RAxisLabels(cfg.GetTitle(), cfg.GetBinLabels()); }
};

} // namespace Internal

} // namespace Experimental
} // namespace ROOT

// Display operator for nicer test assertion errors
std::ostream& operator<<(
   std::ostream&,
   const ROOT::Experimental::RAxisBase::NumericBinningCompatibility&
);
std::ostream& operator<<(
   std::ostream&,
   const ROOT::Experimental::RAxisBase::LabeledBinningCompatibility&
);
std::ostream& operator<<(
   std::ostream&,
   const ROOT::Experimental::RAxisBase::BinningCompatibility&
);

#endif // ROOT7_RAxis header guard
