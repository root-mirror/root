// @(#)root/hist:$Id$
// Author: Lorenzo MOneta 11/2020

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THistRange
#define ROOT_THistRange

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THistRange                                                           //
//                                                                      //
// Class defining a generic range of the histogram                      //
// Used to iterated between bins                                        //
//////////////////////////////////////////////////////////////////////////

#include "TError.h"  // for R__ASSERT

class TH1;



class TBinIterator {
private:
   int fBin;  // global bin number used to advanced
   int fXbin; //  bin X number
   int fYbin; // bin y number
   int fZbin; // bin Z number

   int fNx;   // total x size (nbins+2)
   int fNy;   //  y size
   int fNz;   // z size
   int fXmin; // min x value
   int fXmax; // max x value
   int fYmin; // min y value
   int fYmax; // max y value
   int fZmin; // min z value
   int fZmax; // max z value

   int fDim; // histogram dimension

   // compute global bin number given x,y,x bin numbers
   void SetGlobalBin()
   {
      if (fDim == 1)
         fBin = fXbin;
      else if (fDim == 2)
         fBin = fXbin + fNx * fYbin;
      else
         fBin = fXbin + fNx * (fYbin + fNy * fZbin);
   }

   // private default ctor (used by THistRange)
   TBinIterator()
      : fBin(0), fXbin(0), fYbin(0), fZbin(0), fNx(0), fNy(0), fNz(0), fXmin(0), fXmax(0), fYmin(0), fYmax(0), fZmin(0),
        fZmax(0), fDim(0)
   {
   }

public:
   friend class THistRange;

   // enum defining option range type:
   enum ERangeType {
      kHistRange, // use range provided by histogram
      kAxisBins,  // use allbins within axis limits (no underflow/overflows)
      kAllBins,   // use all bins including underflows/overflows
      kUnOfBins   // collection of all underflow/overflow bins
   };

   TBinIterator(const TH1 *h, ERangeType type);
   // TBinIterator(TBinIterator &rhs) = default;
   // TBinIterator &operator=(const TBinIterator &) = default;

   // implement end
   static TBinIterator End()
   {
      TBinIterator end;
      end.fBin = -1;
      return end;
   }

   // keep inline to be faster
   TBinIterator &operator++()
   {
      if (fXbin < fXmax)
         fXbin++;
      else if (fDim > 1) {
         fXbin = fXmin;
         if (fYbin < fYmax)
            fYbin++;
         else if (fDim > 2) {
            fYbin = fYmin;
            if (fZbin < fZmax)
               fZbin++;
            else {
               fBin = -1;
               return *this;
            }
         } else {
            R__ASSERT(fDim == 2);
            fBin = -1;
            return *this;
         }
      } else {
         R__ASSERT(fDim == 1);
         fBin = -1;
         return *this;
      }
      // fXbin can be incremented to zero only in case of TH2Poly
      // where you can start from negative bin numbers
      // In that case fXbin = 0 will be excluded
      if (fXbin == 0) {
         R__ASSERT(fXmin < 0 && fDim == 1);
         fXbin = 1; // this happens in case of TH2Poly
      }
      SetGlobalBin();
      return *this;
   }

   TBinIterator operator++(int)
   {
      TBinIterator tmp(*this);
      operator++();
      return tmp;
   }

   bool operator==(const TBinIterator &rhs) const { return fBin == rhs.fBin; }
   bool operator!=(const TBinIterator &rhs) const { return fBin != rhs.fBin; }
   int &operator*() { return fBin; }
};

class THistRange {

public:
   typedef TBinIterator iterator;
   typedef TBinIterator::ERangeType RangeType;

   iterator begin() { return fBegin; }

   iterator end() { return fEnd; }

   THistRange(const TH1 *h1, TBinIterator::ERangeType type = TBinIterator::kHistRange);

protected:

   TBinIterator fBegin;
   TBinIterator fEnd;
};


#endif
