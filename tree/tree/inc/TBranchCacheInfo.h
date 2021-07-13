// @(#)root/tree:$Id$
// Author: Philippe Canal, 2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchCacheInfo
#define ROOT_TBranchCacheInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBranchCacheInfo                                                     //
//                                                                      //
// Hold info about which basket are in the cache and if they            //
// have been retrieved from the cache.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include "TBits.h"
#include "TString.h" // for Printf

#include <vector>

class TBranch;

namespace ROOT {
namespace Internal {

class TBranchCacheInfo {

   enum EStates {
      kLoaded = 0,
      kUsed = 1,
      kVetoed = 2,
      kSize = 3
   };

   Int_t fBasketPedestal{-1}; // Lowest basket we have information for.  Its data is in bits [0-3[.
   TBits fInfo;               // kSize bits per baskets (loaded, used, vetoed)

   /// Update the pedestal to be less or equal to basketNumber, shift the bits if needed.
   void UpdatePedestal(Int_t basketNumber)
   {
      if (fBasketPedestal == -1) {
         fBasketPedestal = basketNumber;
      } else if (basketNumber < fBasketPedestal) {
         auto delta = fBasketPedestal - basketNumber;
         fBasketPedestal = basketNumber;
         fInfo <<= (delta * kSize);
      }
   }

   /// Return true if the basket has been marked as having the 'what' state.
   Bool_t TestState(Int_t basketNumber, EStates what) const
   {
      if (basketNumber < fBasketPedestal)
         return kFALSE;
      return fInfo.TestBitNumber(kSize * (basketNumber - fBasketPedestal) + what);
   }

   /// Mark if the basket has been marked has the 'what' state.
   void SetState(Int_t basketNumber, EStates what)
   {
      if (fBasketPedestal <= basketNumber)
         fInfo.SetBitNumber(kSize * (basketNumber - fBasketPedestal) + what, true);
   }

public:
   /// Return true if the basket has been marked as 'used'
   Bool_t HasBeenUsed(Int_t basketNumber) const { return TestState(basketNumber, kUsed); }

   /// Mark if the basket has been marked as 'used'
   void SetUsed(Int_t basketNumber)
   {
      UpdatePedestal(basketNumber);
      SetState(basketNumber, kUsed);
   }

   /// Return true if the basket is currently in the cache.
   Bool_t IsInCache(Int_t basketNumber) const { return TestState(basketNumber, kLoaded); }

   /// Mark if the basket is currently in the cache.
   void SetIsInCache(Int_t basketNumber)
   {
      UpdatePedestal(basketNumber);
      SetState(basketNumber, kLoaded);
   }

   /// Mark if the basket should be vetoed in the next round.
   /// This happens when the basket was loaded in the previous round
   /// and was not used and is overlapping to the next round/cluster
   void Veto(Int_t basketNumber)
   {
      UpdatePedestal(basketNumber);
      SetState(basketNumber, kVetoed);
   }

   /// Return true if the basket is currently vetoed.
   Bool_t IsVetoed(Int_t basketNumber) const { return TestState(basketNumber, kVetoed); }

   /// Return true if all the baskets that are marked loaded are also
   /// mark as used.
   Bool_t AllUsed() const
   {
      auto len = fInfo.GetNbits() / kSize + 1;
      for (UInt_t b = 0; b < len; ++b) {
         if (fInfo[kSize * b + kLoaded] && !fInfo[kSize * b + kUsed]) {
            // Not loaded or (loaded and used)
            return kFALSE;
         }
      }
      return kTRUE;
   }

   /// Return a set of unused basket, let's not re-read them.
   void GetUnused(std::vector<Int_t> &unused)
   {
      unused.clear();
      auto len = fInfo.GetNbits() / kSize + 1;
      for (UInt_t b = 0; b < len; ++b) {
         if (fInfo[kSize * b + kLoaded] && !fInfo[kSize * b + kUsed]) {
            unused.push_back(fBasketPedestal + b);
         }
      }
   }

   /// Reset all info.
   void Reset()
   {
      fBasketPedestal = -1;
      fInfo.ResetAllBits();
   }

   /// Print the info we have for the baskets.
   void Print(const char *owner, Long64_t *entries, Int_t nBaskets, Long64_t nEntries) const
   {
      if (!owner || !entries)
         return;
      auto len = fInfo.GetNbits() / kSize + 1;
      if ( (fBasketPedestal + len) > (UInt_t)nBaskets)
         len = nBaskets - fBasketPedestal;
      if (fBasketPedestal >= 0 && len) {
         for (UInt_t b = 0; b < len; ++b) {
            Printf("Branch %s : basket %d loaded=%d used=%d start entry=%lld", owner, b + fBasketPedestal,
                   (bool)fInfo[kSize * b + kLoaded], (bool)fInfo[kSize * b + kUsed], entries[fBasketPedestal + b]);
         }
         Printf("Branch %s : entry range: [ %lld, %lld [", owner,
                entries[fBasketPedestal],
                (fBasketPedestal + len) >= (UInt_t)nBaskets ? nEntries :  entries[fBasketPedestal + len]);
      }
   }
};

} // Internal
} // ROOT

#endif // ROOT_TBranchCacheInfo
