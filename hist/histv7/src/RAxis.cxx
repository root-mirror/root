/// \file RAxis.cxx
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-08-06
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RAxis.hxx"

#include "ROOT/RHistUtils.hxx"
#include <ROOT/RLogger.hxx>

#include <cmath>
#include <limits>
#include <stdexcept>

// From RHistUtils.hxx:
ROOT::Experimental::RLogChannel &ROOT::Experimental::HistLog() {
   static RLogChannel sLog("ROOT.Hist");
   return sLog;
}


ROOT::Experimental::RAxisBase::~RAxisBase() {}

std::ostream& operator<<(
   std::ostream& s,
   const ROOT::Experimental::RAxisBase::NumericBinningCompatibility& x
) {
   return s << "NumericBinningCompatibility { Flags("
            << ' ' << (x.HasTrivialRegularBinMapping() ? "TRIVIAL_REG" : "trivial_reg")
            << ' ' << (x.HasRegularBinBijection() ? "REG_BIJECTION" : "reg_bijection")
            << ' ' << (x.HasFullBinBijection() ? "FULL_BIJECTION" : "full_bijection")
            << ' ' << (x.MergingIsLossy() ? "LOSSY_MERGE" : "lossy_merge")
            << ' ' << (x.HasRegularBinAliasing() ? "REG_ALIASING" : "reg_aliasing")
            << ' ' << (x.MergingNeedsEmptyUnderflow() ? "NEED_EMPTY_UF" : "need_empty_uf")
            << ' ' << (x.MergingNeedsEmptyOverflow() ? "NEED_EMPTY_OF" : "need_empty_of")
            << ' ' << (x.MergingNeedsTargetGrowth() ? "NEED_GROWTH" : "need_growth")
            << " ) }";
}

ROOT::Experimental::RAxisBase::NumericBinningCompatibility
ROOT::Experimental::RAxisBase::CheckFixedNumericBinningCompat(
   const RAxisBase& source,
   bool growthOccured
) const {
   // Set some shorthands for frequently used quantities
   const double sourceMin = source.GetMinimum();
   const double sourceMax = source.GetMaximum();

   // Compare the positions of the minima/maxima of the source and target axes
   const int minComparison =
      ComparePosToBinBorder(sourceMin, GetFirstBin(), BinSide::kFrom);
   const int maxComparison =
      ComparePosToBinBorder(sourceMax, GetLastBin(), BinSide::kTo);

   // Check if the source underflow and overflow bins must be empty
   //
   // Presuming that the source does indeed have under/overflow bins, this can
   // happen in two different situations:
   //
   // - Target axis is growable, and therefore doesn't support spilling of
   //   under/overflow bins (it would require infinite growth, and also alias).
   //   In this case, both the source underflow and overflow bins must be empty.
   // - Either of these source bins map into multiple target bins, which in the
   //   presence of target under/overflow bins happens if they cover at least
   //   one target regular/under/overflow bin modulo bin comparison tolerance.
   //
   const bool sourceHasUnderOver = !source.CanGrow();
   const bool needEmptyUnderOver = sourceHasUnderOver && CanGrow();
   const bool sourceUnderflowAliasing = sourceHasUnderOver && (minComparison > 0);
   const bool sourceOverflowAliasing = sourceHasUnderOver && (maxComparison < 0);
   const bool needEmptyUnderflow = needEmptyUnderOver || sourceUnderflowAliasing;
   const bool needEmptyOverflow = needEmptyUnderOver || sourceOverflowAliasing;

   // Check for any lossy merge of source under/overflow bins
   //
   // We'll then need to update this variable to account for the lossiness of
   // merging source regular bins into target bins.
   //
   bool mergingIsLossy =
      sourceHasUnderOver && ((minComparison < 0) || (maxComparison > 0));

   // Now, time to look at regular bins
   //
   // The lambda is run immediately, and only once. Its purpose is to enable
   // exiting the regular bin comparison logic early as soon as all source
   // regular bins have been processed. You could call it the ultimate goto.
   //
   bool trivialRegularBinMapping = true;
   bool regularBinAliasing = false;
   [&]{
      // Handle the edge case where the source axis has no regular bin
      if (source.GetNBinsNoOver() == 0) {
         return;
      }

      // Handle the edge case where all regular source axis bins are located
      // either before or after the end of the target axis.
      const bool sourceBeforeTarget =
         (ComparePosToBinBorder(sourceMax, GetFirstBin(), BinSide::kFrom) <= 0);
      const bool sourceAfterTarget =
         (ComparePosToBinBorder(sourceMin, GetLastBin(), BinSide::kTo) >= 0);
      if (sourceBeforeTarget || sourceAfterTarget) {
         // The source axis has at least one regular bin, and it will be merged
         // into a conceptually infinite target under/overflow bin, so an
         // histogram merge from the source to the target loses information.
         mergingIsLossy = true;

         // Unless the source axis has no bin, it is pretty clear that its first
         // regular bin will not map into the first regular target bin, so the
         // mapping from source to target regular bin indices cannot be trivial.
         trivialRegularBinMapping = false;

         // On the other hand, since all source regular bins map into the target
         // overflow bin, it is clear that no source regular bin maps into
         // multiple target bins, so we can leave regularBinAliasing at `false`.
         return;
      }

      // Find the first source bin which doesn't completely map into the target
      // underflow bin.
      //
      // We know that there is one such bin, as we have checked that there is
      // at least one source bin and that not all source bins are fully in the
      // underflow range of the target axis.
      //
      int sourceBin = source.GetFirstBin();
      while (ComparePosToBinBorder(source.GetBinTo(sourceBin),
                                   GetFirstBin(),
                                   BinSide::kFrom) <= 0) {
         ++sourceBin;
      }

      // If any source bin mapped into the target underflow bin like this...
      if (sourceBin != source.GetFirstBin()) {
         // Then those bins don't map into identically numbered target bin
         trivialRegularBinMapping = false;

         // And they map into a larger bin, the infinitely large underflow bin
         mergingIsLossy = true;
      }

      // If the selected bin partially maps into the target underflow bin...
      if (ComparePosToBinBorder(source.GetBinFrom(sourceBin),
                                GetFirstBin(),
                                BinSide::kFrom) < 0) {
         // Then the first bin it maps to is the target underflow bin
         trivialRegularBinMapping = false;

         // And it maps into two bins, underflow and first regular or overflow
         regularBinAliasing = true;
      }
      // At this point, we have taken care of mappings from the first source
      // regular bins to the target underflow bin. What we have left to do is to
      // handle the mapping of remaining source bins, starting from the current
      // one, into target regular and overflow bins.

      // Handle the edge case where the target axis has no regular bin
      if (GetNBinsNoOver() == 0) {
         // There is at least one regular source bin, and the only target bins
         // that it can map into are the infinite underflow and overflow bins.
         // Therefore, this histogram merge is lossy.
         mergingIsLossy = true;

         // The mapping from source to target bins is obviously nontrivial,
         // since the first source bin (which is known to exist) cannot map into
         // the nonexistent first target bin.
         trivialRegularBinMapping = false;

         // Whether there is regular bin aliasing is fully determined by the
         // computation that was carried out above, since the boundary between
         // target underflow and overflow is the only place on the target axis
         // where a source bin can map into multiple target bins.
         return;
      }

      // Find the first target bin which a source bin maps into
      //
      // We know that there will be one such bin, as we have checked that there
      // is at least one target bin and that the source bins are not all located
      // in the target overflow range.
      //
      int targetBin = GetFirstBin();
      while (ComparePosToBinBorder(sourceMin, targetBin, BinSide::kTo) >= 0) {
         ++targetBin;
      }
      // At this point, we know that sourceBin maps into targetBin, and that
      // targetBin is the first bin on the target axis which sourceBin maps to.

      // Iterate over source bins, advancing the target bin index as needed,
      // until either axis has been fully covered
      //
      // The key loop invariant here is that anytime a loop iteration begins,
      // sourceBin designates a source bin which we haven't studied (underflow
      // bin mapping aside), and targetBin designates the first target axis bin
      // which sourceBin maps into.
      //
      for (; sourceBin <= source.GetNBinsNoOver(); ++sourceBin) {
         // Get the source bin's limits
         const double sourceFrom = source.GetBinFrom(sourceBin);
         const double sourceTo = source.GetBinTo(sourceBin);

         // Does the source->target bin mapping remain trivial so far?
         if (targetBin != sourceBin) {
            trivialRegularBinMapping = false;
         }

         // Does the first target bin cover nontrivial extra range on the left
         // of the source bin?
         if (ComparePosToBinBorder(sourceFrom, targetBin, BinSide::kFrom) > 0) {
            // If so, some information about the position of past source
            // histogram fills will be lost upon merging
            mergingIsLossy = true;
         }

         // Next, iterate over target bins until we find a target bin which
         // extends beyond the end of the current source bin (and therefore
         // into the next source bin, if any) or we reach the end of the target
         // axis in attempting to do so.
         bool endOfTargetAxis = false;
         const int firstTargetBin = targetBin;
         while (ComparePosToBinBorder(sourceTo, targetBin, BinSide::kTo) >= 0) {
            if (targetBin < GetNBinsNoOver()) {
               ++targetBin;
            } else {
               endOfTargetAxis = true;
               break;
            }
         }

         // Whether iteration succeeded or failed, we know that every targetBin
         // that was covered by iteration, with the possible exception of the
         // current target, is a bin that the source bin maps into.
         int numCoveredBins = targetBin - firstTargetBin;

         // Next, we need to tell which other bins the source bin maps into
         int lastBinCmpResult;
         if (endOfTargetAxis) {
            // If the end of the target axis was reached, then we know that
            // sourceBin maps into the current targetBin, because we didn't
            // manage to find a targetBin which even extends beyond the end of
            // the current sourceBin.
            ++numCoveredBins;

            // In that case, however, we need to check if the current source bin
            // maps into the target overflow bin.
            lastBinCmpResult =
               ComparePosToBinBorder(sourceTo, GetLastBin(), BinSide::kTo);
         } else {
            // If we managed to find a targetBin which extends beyond the end of
            // the current sourceBin, then we must check if this bin still
            // covers some of the current sourceBin range on the left.
            lastBinCmpResult =
               ComparePosToBinBorder(sourceTo, targetBin, BinSide::kFrom);
         }

         // Does the current source bin map into the current targetBin or into
         // the target overflow bin?
         if (lastBinCmpResult > 0) {
            // If so, then that's one more covered bin...
            ++numCoveredBins;

            // ...and we map into another bin that, by definition, spans some
            // extra range, so the merge loses fill location information.
            mergingIsLossy = true;
         }

         // If the current source bin maps into multiple target bins, then it
         // must be empty for histogram merging to succeed.
         if (numCoveredBins > 1) {
            regularBinAliasing = true;
         }

         // If the end of the target axis was reached, then we must abort the
         // loop, as we cannot maintain the loop invariant that at the beginning
         // of a loop iteration, `targetBin` must be the first bin which the
         // active `sourceBin` maps into.
         if (endOfTargetAxis) break;
      }

      // Was the end of the target axis reached w/o covering all source bins?
      if (sourceBin < source.GetNBinsNoOver()) {
         // In that case, the extra source bins map into the infinite target
         // overflow bin, so the merge loses information...
         mergingIsLossy = true;

         // ...and these source bins do not map into target bins with the same
         // indices, so the bin index mapping is nontrivial
         trivialRegularBinMapping = false;
      }
   }();

   // Compute the remaining properties that we need
   const bool regularBinBijection =
      trivialRegularBinMapping && (GetNBinsNoOver() == source.GetNBinsNoOver());
   const bool fullBinBijection =
      regularBinBijection && (source.CanGrow() == CanGrow());

   // Produce the final result of the numerical axis binning comparison
   using Flags = NumericBinningCompatibility::Flags;
   return NumericBinningCompatibility(
      Flags(trivialRegularBinMapping * Flags::kTrivialRegularBinMapping
            + regularBinBijection * Flags::kRegularBinBijection
            + fullBinBijection * Flags::kFullBinBijection
            + mergingIsLossy * Flags::kMergingIsLossy
            + regularBinAliasing * Flags::kRegularBinAliasing
            + needEmptyUnderflow * Flags::kNeedEmptyUnderflow
            + needEmptyOverflow * Flags::kNeedEmptyOverflow
            + growthOccured * Flags::kTargetMustGrow)
   );
}

std::ostream& operator<<(
   std::ostream& s,
   const ROOT::Experimental::RAxisBase::LabeledBinningCompatibility& x
) {
   s << "LabeledBinningCompatibility { Flags("
     << ' ' << (x.TargetMustGrow() ? "MUST_GROW" : "must_grow")
     << ' ' << (x.LabelOrderDiffers() ? "BAD_ORDER" : "bad_order")
     << ' ' << (x.TargetWillHaveExtraBins() ? "EXTRA_BINS" : "extra_bins")
     << " ), SourceOnlyLabels { ";
   for (const auto& label: x.SourceOnlyLabels()) {
      s << '"' << label << "\" ";
   }
   return s << "} }";
}

std::ostream& operator<<(
   std::ostream& s,
   const ROOT::Experimental::RAxisBase::BinningCompatibility& x
) {
   s << "BinningCompatibility { ";
   using CompatKind = ROOT::Experimental::RAxisBase::BinningCompatibility::CompatKind;
   switch (x.Kind()) {
      case CompatKind::kIncompatible:
         s << "Incompatible";
         break;

      case CompatKind::kNumeric:
         s << x.GetNumeric();
         break;

      case CompatKind::kLabeled:
         s << x.GetLabeled();
         break;
   };
   return s << " }";
}

void ROOT::Experimental::RAxisBase::BinningCompatibility::CheckKind(CompatKind expectedKind) const {
   if (fKind != expectedKind) {
      throw std::runtime_error("The queried property is invalid for this "
         "kind of axis binning comparison");
   }
}

ROOT::Experimental::RAxisBase::BinningCompatibility
ROOT::Experimental::RAxisBase::CheckBinningCompat(const RAxisBase& source) const {
   // Handle labeled axis edge case
   //
   // NOTE: This must be handled at the axis base class level, because C++ does
   //       not provide a way to dispatch at runtime based on _both_ types of
   //       the `this` and `source` axes.
   //
   const auto target_lbl_ptr = dynamic_cast<const RAxisLabels*>(this);
   const auto source_lbl_ptr = dynamic_cast<const RAxisLabels*>(&source);
   if (bool(target_lbl_ptr) != bool(source_lbl_ptr)) {
      return BinningCompatibility();
   } else if (target_lbl_ptr) {
      return BinningCompatibility(
         target_lbl_ptr->CheckLabeledBinningCompat(*source_lbl_ptr)
      );
   }

   // If control reached this point, then we know that both the source and the
   // target axis use numerical bin borders
   return BinningCompatibility(CheckNumericBinningCompat(source));
}

int ROOT::Experimental::RAxisEquidistant::GetBinIndexForLowEdge(double x) const noexcept
{
   // fracBinIdx is the fractional bin index of x in this axis. It's (close to)
   // an integer if it's an axis border.
   double fracBinIdx = GetFirstBin() + FindBinRaw(x);

   // fracBinIdx might be 12.99999999. It's a bin border if the deviation from
   // an regular bin border is "fairly small".
   int binIdx = std::round(fracBinIdx);
   double binOffset = fracBinIdx - binIdx;
   if (std::fabs(binOffset) > 10 * std::numeric_limits<double>::epsilon())
      return RAxisBase::kInvalidBin;

   // If the bin index is below the first bin (i.e. x is the lower edge of the
   // underflow bin) then it's out of range.
   if (binIdx < GetFirstBin())
      return RAxisBase::kInvalidBin;
   // If x is the lower edge of the overflow bin then that's still okay - but if
   // even the bin before binIdx is an overflow it's out of range.
   if (binIdx > GetLastBin() + 1)
      return RAxisBase::kInvalidBin;

   return binIdx;
}

ROOT::Experimental::RAxisBase::NumericBinningCompatibility
ROOT::Experimental::RAxisGrow::CheckNumericBinningCompat(
   const RAxisBase& source
) const {
   // Convenience shorthands
   const double sourceMin = source.GetMinimum();
   const double sourceMax = source.GetMaximum();

   // Prepare to simulate axis growth if need be
   RAxisGrow targetAfterGrowth;
   const RAxisGrow* targetPtr = this;

   // Check if axis growth is needed
   const bool growLeft =
      (ComparePosToBinBorder(sourceMin, GetFirstBin(), BinSide::kFrom) < 0);
   const bool growRight =
      (ComparePosToBinBorder(sourceMax, GetLastBin(), BinSide::kTo) > 0);
   const bool mustGrow = growLeft || growRight;
   if (mustGrow) {
      // Simulate axis growth on the left-hand side
      const double leftGrowth =
         static_cast<double>(growLeft) * (GetMinimum() - sourceMin);
      int leftBins = std::floor(leftGrowth / GetBinWidth());
      double leftBorder = GetMinimum() - leftBins*GetBinWidth();
      if (ComparePosToBinBorder(sourceMin,
                                leftBorder,
                                kNoBinWidth,
                                GetBinWidth()) < 0) {
         ++leftBins;
         leftBorder -= GetBinWidth();
      }

      // Simulate axis growth on the right-hand side
      const double rightGrowth =
         static_cast<double>(growRight) * (sourceMax - GetMaximum());
      int rightBins = std::floor(rightGrowth / GetBinWidth());
      double rightBorder = GetMaximum() + rightBins*GetBinWidth();
      if (ComparePosToBinBorder(sourceMax,
                                rightBorder,
                                GetBinWidth(),
                                kNoBinWidth) > 0) {
         ++rightBins;
         rightBorder += GetBinWidth();
      }

      // Build simulated grown axis and inject it into targetPtr
      targetAfterGrowth =
         RAxisGrow(GetNBinsNoOver() + leftBins + rightBins,
                   leftBorder,
                   rightBorder);
      targetPtr = &targetAfterGrowth;
   }

   // Call back binning comparison hook on the possibly grown axis
   return targetPtr->CheckFixedNumericBinningCompat(source, mustGrow);
}

int ROOT::Experimental::RAxisIrregular::GetBinIndexForLowEdge(double x) const noexcept
{
   // Check in which bin `x` resides
   double fracBinIdx = FindBinRaw(x);
   const int binIdx = fracBinIdx;

   // Are we close to the lower and upper bin boundaries, if any?
   constexpr double tol = 10 * std::numeric_limits<double>::epsilon();
   if (binIdx >= GetFirstBin()) {
      const double lowBound = GetBinFrom(binIdx);
      if (std::fabs(x - lowBound) < tol * std::fabs(lowBound))
         return binIdx;
   }
   if (binIdx <= GetLastBin()) {
      const double upBound = GetBinTo(binIdx);
      if (std::fabs(x - upBound) < tol * std::fabs(upBound))
         return binIdx + 1;
   }

   // If not, report failure
   return RAxisBase::kInvalidBin;
}
