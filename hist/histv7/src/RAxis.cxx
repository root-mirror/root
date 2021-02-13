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

// From RHistUtils.hxx:
ROOT::Experimental::RLogChannel &ROOT::Experimental::HistLog() {
   static RLogChannel sLog("ROOT.Hist");
   return sLog;
}


ROOT::Experimental::RAxisBase::~RAxisBase() {}

bool ROOT::Experimental::RAxisBase::HasSameBinningAs(const RAxisBase& other) const {
   // Bin borders must match
   if (!HasSameBinBordersAs(other))
      return false;

   // Bin labels must match
   auto lbl_ptr = dynamic_cast<const RAxisLabels*>(this);
   auto other_lbl_ptr = dynamic_cast<const RAxisLabels*>(&other);
   if (bool(lbl_ptr) != bool(other_lbl_ptr)) {
      return false;
   } else if (lbl_ptr) {
      auto lbl_cmp = lbl_ptr->CompareBinLabels(*other_lbl_ptr);
      return (lbl_cmp == RAxisLabels::kLabelsCmpSame);
   } else {
      return true;
   }
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

bool ROOT::Experimental::RAxisEquidistant::HasSameBinBordersAs(const RAxisBase& other) const {
   // This is an optimized override for the equidistant-equidistant case,
   // fall back to the default implementation if we're not in that case.
   auto other_eq_ptr = dynamic_cast<const RAxisEquidistant*>(&other);
   if (!other_eq_ptr)
      return RAxisBase::HasSameBinBordersAs(other);
   const RAxisEquidistant& other_eq = *other_eq_ptr;

   // Can directly compare equidistant/growable axis properties in this case
   return fInvBinWidth == other_eq.fInvBinWidth &&
          fLow == other_eq.fLow &&
          fNBinsNoOver == other_eq.fNBinsNoOver &&
          CanGrow() == other_eq.CanGrow();
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

bool ROOT::Experimental::RAxisIrregular::HasSameBinBordersAs(const RAxisBase& other) const {
   // This is an optimized override for the irregular-irregular case,
   // fall back to the default implementation if we're not in that case.
   auto other_irr_ptr = dynamic_cast<const RAxisIrregular*>(&other);
   if (!other_irr_ptr)
      return RAxisBase::HasSameBinBordersAs(other);
   const RAxisIrregular& other_irr = *other_irr_ptr;

   // Only need to compare bin borders in this specialized case
   return fBinBorders == other_irr.fBinBorders;
}

ROOT::Experimental::EAxisCompatibility ROOT::Experimental::CanMap(const RAxisEquidistant &target,
                                                                  const RAxisEquidistant &source) noexcept
{
   // First, let's get the common "all parameters are equal" case out of the way
   if (source.HasSameBinningAs(target))
      return EAxisCompatibility::kIdentical;

   // Do the source min/max boundaries correspond to target bin boundaries?
   int idxTargetLow = target.GetBinIndexForLowEdge(source.GetMinimum());
   int idxTargetHigh = target.GetBinIndexForLowEdge(source.GetMaximum());
   if (idxTargetLow < 0 || idxTargetHigh < 0)
      // If not, the source is incompatible with the target since the first or
      // last source bin does not map into a target axis bin.
      return EAxisCompatibility::kIncompatible;

   // If so, and if the bin width is the same, then since we've eliminated the
   // care where min/max/width are equal, source must be a subset of target.
   if (source.GetInverseBinWidth() == target.GetInverseBinWidth())
      return EAxisCompatibility::kContains;

   // Now we are left with the case
   //   source: 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6
   //   target: ...0.0, 0.3, 0.6...
   // The question is: is the ratio of the bin width identical to the ratio of
   // the number of bin?
   if (std::fabs(target.GetInverseBinWidth() * source.GetNBinsNoOver() -
                 source.GetInverseBinWidth() * (idxTargetHigh - idxTargetLow)) > 1E-6 * target.GetInverseBinWidth())
      return EAxisCompatibility::kIncompatible;

   // source is a fine-grained version of target.
   return EAxisCompatibility::kSampling;
}

// TODO: the other CanMap() overloads
