// @(#)root/hist:$Id$
// Author: Lorenzo MOneta 11/2020

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THistRange.h"

#include "TH1.h"
#include "TH2Poly.h"
#include "TProfile2Poly.h"

THistRange::THistRange(const TH1 *h1, RangeType type)
{
   fBegin = TBinIterator(h1, type);
   fEnd = TBinIterator::End();
}

// implementation of TBinIterator class

/// constructor of  TBInIterator taking as input an histogram pointer
/// This constructors set iterator, the global bin,  to its first value (begin)
TBinIterator::TBinIterator(const TH1 *h, ERangeType type)
   : fNx(0), fNy(0), fNz(0),
     fXmin(0), fXmax(0), fYmin(0), fYmax(0), fZmin(0), fZmax(0)
{
   // deal with special cases (e.g. TH2Poly)
   if (h->IsA() == TH2Poly::Class() || h->IsA() == TProfile2Poly::Class()) {
      const TH2Poly *hpoly = static_cast<const TH2Poly *>(h);
      R__ASSERT(hpoly);
      if (type == TBinIterator::kAllBins) {
         // this will loop on all bins (one should exclude 0)
         fXmin = -9;
         fXmax = hpoly->GetNumberOfBins();
      } else if (type == TBinIterator::kUnOfBins) {
         // overflow bins in TH2Poly are from -9 to -1
         fXmin = -9;
         fXmax = -1;
      } else {
         // standard bin loop
         fXmin = 1;
         fXmax = hpoly->GetNumberOfBins();
      }

      fYmin = 0;
      fYmax = 0;
      fZmin = 0;
      fZmax = 0;
      fDim = 1; // this case is equivalent to have one dimension
   }
   // general case for TH1,TH2, TH3 and prfile classes
   else {

      fNx = h->GetNbinsX() + 2;
      fNy = h->GetNbinsY() + 2;
      fNz = h->GetNbinsZ() + 2;
      fDim = h->GetDimension();

      if (type == TBinIterator::kHistRange) {
         fXmin = h->GetXaxis()->GetFirst();
         fXmax = h->GetXaxis()->GetLast();
         fYmin = h->GetYaxis()->GetFirst();
         fYmax = h->GetYaxis()->GetLast();
         fZmin = h->GetZaxis()->GetFirst();
         fZmax = h->GetZaxis()->GetLast();
      } else if (type == TBinIterator::kAxisBins) {
         fXmin = 1;
         fXmax = h->GetNbinsX();
         fYmin = 1;
         fYmax = h->GetNbinsY();
         fZmin = 1;
         fZmax = h->GetNbinsZ();
      } else if (type == TBinIterator::kAllBins || type == TBinIterator::kUnOfBins) {
         fXmin = 0;
         fXmax = h->GetNbinsX() + 1;
         fYmin = 0;
         fYmax = h->GetNbinsY() + 1;
         fZmin = 0;
         fZmax = h->GetNbinsZ() + 1;
      }
   }

   // set bins to initial value
   fXbin = fXmin;
   fYbin = fYmin;
   fZbin = fZmin;

   SetGlobalBin();
}
