// @(#)root/hist:$Id$
// Author: Filip Ilic

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TProfile2Poly.h"
#include "TProfileHelper.h"
#include "TMultiGraph.h"
#include "TList.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <set>

/** \class TProfile2Poly
    \ingroup Hist
2D Profile Histogram with Polygonal Bins.

tprofile2polyRealisticModuleError.C and tprofile2polyRealistic.C illustrate how
to use this class.
*/

/** \class TProfile2PolyBin
    \ingroup Hist
Helper class to represent a bin in the TProfile2Poly histogram
*/

ClassImp(TProfile2Poly);

////////////////////////////////////////////////////////////////////////////////
/// TProfile2PolyBin constructor.

TProfile2PolyBin::TProfile2PolyBin()
{
   fSumw = 0;
   fSumvw = 0;
   fSumw2 = 0;
   fSumwv2 = 0;
   fError = 0;
   fAverage = 0;
   fErrorMode = kERRORMEAN;
}

////////////////////////////////////////////////////////////////////////////////
/// TProfile2PolyBin constructor.

TProfile2PolyBin::TProfile2PolyBin(TObject *poly, Int_t bin_number) : TH2PolyBin(poly, bin_number)
{
   fSumw = 0;
   fSumvw = 0;
   fSumw2 = 0;
   fSumwv2 = 0;
   fError = 0;
   fAverage = 0;
   fErrorMode = kERRORMEAN;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge.

void TProfile2PolyBin::Merge(const TProfile2PolyBin *toMerge)
{
   this->fSumw += toMerge->fSumw;
   this->fSumvw += toMerge->fSumvw;
   this->fSumw2 += toMerge->fSumw2;
   this->fSumwv2 += toMerge->fSumwv2;
}

////////////////////////////////////////////////////////////////////////////////
/// Update.

void TProfile2PolyBin::Update()
{
   UpdateAverage();
   UpdateError();
   SetChanged(true);
}

////////////////////////////////////////////////////////////////////////////////
/// Update average.

void TProfile2PolyBin::UpdateAverage()
{
   if (fSumw != 0) fAverage = fSumvw / fSumw;
}

////////////////////////////////////////////////////////////////////////////////
/// Update error.

void TProfile2PolyBin::UpdateError()
{
   Double_t tmp = 0;
   if (fSumw != 0) tmp = std::sqrt((fSumwv2 / fSumw) - (fAverage * fAverage));

   fError = tmp;

   return;

}

////////////////////////////////////////////////////////////////////////////////
/// Clear statistics.

void TProfile2PolyBin::ClearStats()
{
   fSumw = 0;
   fSumvw = 0;
   fSumw2 = 0;
   fSumwv2 = 0;
   fError = 0;
   fAverage = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill.

void TProfile2PolyBin::Fill(Double_t value, Double_t weight)
{
   fSumw += weight;
   fSumvw += value * weight;
   fSumw2 += weight * weight;
   fSumwv2 += weight * value * value;
   this->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// TProfile2Poly constructor.

TProfile2Poly::TProfile2Poly(const char *name, const char *title, Double_t xlow, Double_t xup, Double_t ylow,
                             Double_t yup)
   : TH2Poly(name, title, xlow, xup, ylow, yup)
{
}

////////////////////////////////////////////////////////////////////////////////
/// TProfile2Poly constructor.

TProfile2Poly::TProfile2Poly(const char *name, const char *title, Int_t nX, Double_t xlow, Double_t xup, Int_t nY,
                             Double_t ylow, Double_t yup)
   : TH2Poly(name, title, nX, xlow, xup, nY, ylow, yup)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create bin.

TProfile2PolyBin *TProfile2Poly::CreateBin(TObject *poly)
{
   if (!poly) return 0;

   if (fBins == 0) {
      fBins = new TList();
      fBins->SetOwner();
   }

   fNcells++;
   Int_t ibin = fNcells - kNOverflow;
   return new TProfile2PolyBin(poly, ibin);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill

Int_t TProfile2Poly::Fill(Double_t xcoord, Double_t ycoord, Double_t value)
{
   return Fill(xcoord, ycoord, value, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill

Int_t TProfile2Poly::Fill(Double_t xcoord, Double_t ycoord, Double_t value, Double_t weight)
{
   // Find region in which the hit occurred
   Int_t tmp = GetOverflowRegionFromCoordinates(xcoord, ycoord);
   if (tmp < 0) {
      Int_t overflow_idx = OverflowIdxToArrayIdx(tmp);
      fOverflowBins[overflow_idx].Fill(value, weight);
      fOverflowBins[overflow_idx].SetContent(fOverflowBins[overflow_idx].fAverage );
   }

   // Find the cell to which (x,y) coordinates belong to
   Int_t n = (Int_t)(floor((xcoord - fXaxis.GetXmin()) / fStepX));
   Int_t m = (Int_t)(floor((ycoord - fYaxis.GetXmin()) / fStepY));

   // Make sure the array indices are correct.
   if (n >= fCellX) n = fCellX - 1;
   if (m >= fCellY) m = fCellY - 1;
   if (n < 0) n = 0;
   if (m < 0) m = 0;

   // ------------ Update global (per histo) statistics
   fTsumw += weight;
   fTsumw2 += weight * weight;
   fTsumwx += weight * xcoord;
   fTsumwx2 += weight * xcoord * xcoord;
   fTsumwy += weight * ycoord;
   fTsumwy2 += weight * ycoord * ycoord;
   fTsumwxy += weight * xcoord * ycoord;
   fTsumwz += weight * value;
   fTsumwz2 += weight * value * value;

   // ------------ Update local (per bin) statistics
   TProfile2PolyBin *bin;
   TIter next(&fCells[n + fCellX * m]);
   TObject *obj;
   while ((obj = next())) {
      bin = (TProfile2PolyBin *)obj;
      if (bin->IsInside(xcoord, ycoord)) {
         fEntries++;
         bin->Fill(value, weight);
         bin->Update();
         bin->SetContent(bin->fAverage);
      }
   }

   return tmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge

Long64_t TProfile2Poly::Merge(TCollection *in)
{
   Int_t size = in->GetSize();

   std::vector<TProfile2Poly *> list;
   list.reserve(size);

   for (int i = 0; i < size; i++) {
      list.push_back((TProfile2Poly *)((TList *)in)->At(i));
   }
   return this->Merge(list);
}

////////////////////////////////////////////////////////////////////////////////
/// Merge

Long64_t TProfile2Poly::Merge(const std::vector<TProfile2Poly *> &list)
{
   if (list.size() == 0) {
      std::cout << "[FAIL] TProfile2Poly::Merge: No objects to be merged " << std::endl;
      return -1;
   }

   // ------------ Check that bin numbers of TP2P's to be merged are equal
   std::set<Int_t> numBinUnique;
   for (const auto &histo : list) {
      if (histo->fBins) numBinUnique.insert(histo->fBins->GetSize());
   }
   if (numBinUnique.size() != 1) {
      std::cout << "[FAIL] TProfile2Poly::Merge: Bin numbers of TProfile2Polys to be merged differ!" << std::endl;
      return -1;
   }
   Int_t nbins = *numBinUnique.begin();

   // ------------ Update global (per histo) statistics
   for (const auto &histo : list) {
      this->fEntries += histo->fEntries;
      this->fTsumw += histo->fTsumw;
      this->fTsumw2 += histo->fTsumw2;
      this->fTsumwx += histo->fTsumwx;
      this->fTsumwx2 += histo->fTsumwx2;
      this->fTsumwy += histo->fTsumwy;
      this->fTsumwy2 += histo->fTsumwy2;
      this->fTsumwxy += histo->fTsumwxy;
      this->fTsumwz += histo->fTsumwz;
      this->fTsumwz2 += histo->fTsumwz2;

      // Merge overflow bins
      for (Int_t i = 0; i < kNOverflow; ++i) {
         this->fOverflowBins[i].Merge(&histo->fOverflowBins[i]);
      }
   }

   // ------------ Update local (per bin) statistics
   TProfile2PolyBin *dst = nullptr;
   TProfile2PolyBin *src = nullptr;
   for (Int_t i = 0; i < nbins; i++) {
      dst = (TProfile2PolyBin *)fBins->At(i);

      for (const auto &e : list) {
         src = (TProfile2PolyBin *)e->fBins->At(i);
         dst->Merge(src);
      }

      dst->Update();
   }

   this->SetContentToAverage();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set content to average.

void TProfile2Poly::SetContentToAverage()
{
   Int_t nbins = fBins ? fBins->GetSize() : 0;
   for (Int_t i = 0; i < nbins; i++) {
      TProfile2PolyBin *bin = (TProfile2PolyBin *)fBins->At(i);
      bin->Update();
      bin->SetContent(bin->fAverage);
   }
   for (Int_t i = 0; i < kNOverflow; ++i) {
      TProfile2PolyBin & bin = fOverflowBins[i];
      bin.Update();
      bin.SetContent(bin.fAverage);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set content to error.

void TProfile2Poly::SetContentToError()
{
   Int_t nbins = fBins ? fBins->GetSize() : 0;
   for (Int_t i = 0; i < nbins; i++) {
      TProfile2PolyBin *bin = (TProfile2PolyBin *)fBins->At(i);
      bin->Update();
      bin->SetContent(bin->fError);
   }
   for (Int_t i = 0; i < kNOverflow; ++i) {
      TProfile2PolyBin & bin = fOverflowBins[i];
      bin.Update();
      bin.SetContent(bin.fError);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin content.

Double_t TProfile2Poly::GetBinContent(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin<0) return fOverflowBins[-bin - 1].GetContent();
   return ((TProfile2PolyBin*) fBins->At(bin-1))->GetContent();
}


////////////////////////////////////////////////////////////////////////////////
/// Get bin effective entries.

Double_t TProfile2Poly::GetBinEffectiveEntries(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0) return fOverflowBins[-bin - 1].GetEffectiveEntries();
   return ((TProfile2PolyBin *)fBins->At(bin - 1))->GetEffectiveEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin entries.

Double_t TProfile2Poly::GetBinEntries(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0) return fOverflowBins[-bin - 1].GetEntries();
   return ((TProfile2PolyBin *)fBins->At(bin - 1))->GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin entries W2.

Double_t TProfile2Poly::GetBinEntriesW2(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0) return fOverflowBins[-bin - 1].GetEntriesW2();
   return ((TProfile2PolyBin *)fBins->At(bin - 1))->GetEntriesW2();
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin entries VW.

Double_t TProfile2Poly::GetBinEntriesVW(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0) return fOverflowBins[-bin - 1].GetEntriesVW();
   return ((TProfile2PolyBin *)fBins->At(bin - 1))->GetEntriesVW();
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin entries WV2.

Double_t TProfile2Poly::GetBinEntriesWV2(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0) return fOverflowBins[-bin - 1].GetEntriesWV2();
   return ((TProfile2PolyBin *)fBins->At(bin - 1))->GetEntriesWV2();
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin error.

Double_t TProfile2Poly::GetBinError(Int_t bin) const
{
   Double_t tmp = 0;
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin < 0)
      tmp =  fOverflowBins[-bin - 1].GetError();
   else
      tmp = ((TProfile2PolyBin *)fBins->At(bin - 1))->GetError();

   return (fErrorMode == kERRORSPREAD) ?  tmp : tmp / std::sqrt(GetBinEffectiveEntries(bin));

}

////////////////////////////////////////////////////////////////////////////////
/// Fill the array stats from the contents of this profile.
/// The array stats must be correctly dimensioned in the calling program.
///
/// - stats[0] = sumw
/// - stats[1] = sumw2
/// - stats[2] = sumwx
/// - stats[3] = sumwx2
/// - stats[4] = sumwy
/// - stats[5] = sumwy2
/// - stats[6] = sumwxy
/// - stats[7] = sumwz
/// - stats[8] = sumwz2
///
/// If no axis-subrange is specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If a sub-range is specified, the function recomputes these quantities
/// from the bin contents in the current axis range.

void TProfile2Poly::GetStats(Double_t *stats) const
{
   stats[0] = fTsumw;
   stats[1] = fTsumw2;
   stats[2] = fTsumwx;
   stats[3] = fTsumwx2;
   stats[4] = fTsumwy;
   stats[5] = fTsumwy2;
   stats[6] = fTsumwxy;
   stats[7] = fTsumwz;
   stats[8] = fTsumwz2;
}

////////////////////////////////////////////////////////////////////////////////
/// Print overflow regions.

void TProfile2Poly::PrintOverflowRegions()
{
   Double_t total = 0;
   Double_t cont = 0;
   for (Int_t i = 0; i < kNOverflow; ++i) {
      cont = GetOverflowContent(i);
      total += cont;
      std::cout << "\t" << cont << "\t";
      if ((i + 1) % 3 == 0) std::cout << std::endl;
   }

   std::cout << "Total: " << total << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset

void TProfile2Poly::Reset(Option_t *opt)
{
   TIter next(fBins);
   TObject *obj;
   TProfile2PolyBin *bin;

   // Clears bin contents
   while ((obj = next())) {
      bin = (TProfile2PolyBin *)obj;
      bin->ClearContent();
      bin->ClearStats();
   }
   TH2::Reset(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// The overflow regions are calculated by considering x, y coordinates.
/// The Middle bin at -5 contains all the TProfile2Poly bins.
///
/// ~~~ {.cpp}
///           -0 -1 -2
///           ________
///    -1:   |__|__|__|
///    -4:   |__|__|__|
///    -7:   |__|__|__|
/// ~~~

Int_t TProfile2Poly::GetOverflowRegionFromCoordinates(Double_t x, Double_t y)
{


   Int_t region = 0;

   if (fNcells <= kNOverflow) return 0;

   // --- y offset
   if (y > fYaxis.GetXmax())
      region += -1;
   else if (y > fYaxis.GetXmin())
      region += -4;
   else
      region += -7;

   // --- x offset
   if (x > fXaxis.GetXmax())
      region += -2;
   else if (x > fXaxis.GetXmin())
      region += -1;
   else
      region += 0;

   return region;
}

////////////////////////////////////////////////////////////////////////////////
/// Set error option.

void TProfile2Poly::SetErrorOption(EErrorType type)
{
   fErrorMode = type;
}
