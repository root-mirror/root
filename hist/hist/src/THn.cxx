// @(#)root/hist:$Id$
// Author: Axel Naumann (2011-12-13)

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THn.h"

namespace {
   //______________________________________________________________________________
   //
   // Helper struct to hold one dimension's bin range for THnBinIter.
   /////////////////////////////////////////////////////////////////////////////

   struct CounterRange_t {
      Int_t i;
      Int_t first;
      Int_t last;
      Int_t len;
      Long64_t cellSize;
   };

   //______________________________________________________________________________
   //
   // THnBinIter iterates over all bins of a THn, recursing over all dimensions.
   /////////////////////////////////////////////////////////////////////////////

   class THnBinIter: public ROOT::Internal::THnBaseBinIter {
   public:
      THnBinIter(Int_t dim, const TObjArray* axes, const TNDArray* arr,
                 Bool_t respectAxisRange);
      ~THnBinIter() { delete [] fCounter; }

      Long64_t Next(Int_t* coord = 0);
      Int_t GetCoord(Int_t dim) const { return fCounter[dim].i; }
   private:
      THnBinIter(const THnBinIter&); // intentionally unimplemented
      THnBinIter& operator=(const THnBinIter&); // intentionally unimplemented

   public:
      Int_t fNdimensions;
      Long64_t fIndex;
      const TNDArray* fArray;
      CounterRange_t* fCounter;
   };


   /////////////////////////////////////////////////////////////////////////////
   /// Construct a THnBinIter.

   THnBinIter::THnBinIter(Int_t dim, const TObjArray* axes,
                              const TNDArray* arr, Bool_t respectAxisRange):
      ROOT::Internal::THnBaseBinIter(respectAxisRange),
      fNdimensions(dim), fIndex(-1), fArray(arr) {
      fCounter = new CounterRange_t[dim]();
      for (Int_t i = 0; i < dim; ++i) {
         TAxis *axis = (TAxis*) axes->At(i);
         fCounter[i].len  = axis->GetNbins() + 2;
         fCounter[i].cellSize  = arr->GetCellSize(i);
         if (!respectAxisRange || !axis->TestBit(TAxis::kAxisRange)) {
            fCounter[i].first = 0;
            fCounter[i].last  = fCounter[i].len - 1;
            fCounter[i].i     = 0;
            continue;
         }
         fHaveSkippedBin = kTRUE;
         Int_t min = axis->GetFirst();
         Int_t max = axis->GetLast();
         if (min == 0 && max == 0) {
            // special case where TAxis::SetBit(kAxisRange) and
            // over- and underflow bins are de-selected.
            // first and last are == 0 due to axis12->SetRange(1, axis12->GetNbins());
            min = 1;
            max = axis->GetNbins();
         }
         fCounter[i].first = min;
         fCounter[i].last  = max;
         fCounter[i].i     = min;
         fIndex += fCounter[i].first * fCounter[i].cellSize;
      }
      // First Next() will increment it:
      --fCounter[dim - 1].i;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Return the current linear bin index (in range), then go to the next bin.
   /// If all bins have been visited, return -1.

   Long64_t THnBinIter::Next(Int_t* coord /*= 0*/) {
      if (fNdimensions < 0) return -1; // end
      ++fCounter[fNdimensions - 1].i;
      ++fIndex;
      // Wrap around if needed
      for (Int_t d = fNdimensions - 1; d > 0 && fCounter[d].i > fCounter[d].last; --d) {
         // We skip last + 1..size and 0..first - 1, adjust fIndex
         Int_t skippedCells = fCounter[d].len - (fCounter[d].last + 1);
         skippedCells += fCounter[d].first;
         fIndex += skippedCells * fCounter[d].cellSize;
         fCounter[d].i = fCounter[d].first;
         ++fCounter[d - 1].i;
      }
      if (fCounter[0].i > fCounter[0].last) {
         fNdimensions = -1;
         return -1;
      }
      if (coord) {
         for (Int_t d = 0; d < fNdimensions; ++d) {
            coord[d] = fCounter[d].i;
         }
      }
      return fIndex;
   }
} // unnamed namespce



/** \class THn
    \ingroup Hist
Multidimensional histogram.

Use a THn if you really, really have to store more than three dimensions,
and if a large fraction of all bins are filled.
Better alternatives are
  -  THnSparse if a fraction of all bins are filled
  - TTree

The major problem of THn is the memory use caused by n-dimensional
histogramming: a THnD with 8 dimensions and 100 bins per dimension needs
more than 2.5GB of RAM!

To construct a THn object you must use one of its templated, derived
classes:

    THnD (typedef for THnT<Double_t>): bin content held by a Double_t,
    THnF (typedef for THnT<Float_t>): bin content held by a Float_t,
    THnL (typedef for THnT<Long_t>): bin content held by a Long_t,
    THnI (typedef for THnT<Int_t>): bin content held by an Int_t,
    THnS (typedef for THnT<Short_t>): bin content held by a Short_t,
    THnC (typedef for THnT<Char_t>): bin content held by a Char_t,

They take name and title, the number of dimensions, and for each dimension
the number of bins, the minimal, and the maximal value on the dimension's
axis. A TH2F h("h","h",10, 0., 10., 20, -5., 5.) would correspond to

    Int_t bins[2] = {10, 20};
    Double_t xmin[2] = {0., -5.};
    Double_t xmax[2] = {10., 5.};
    THnF hn("hn", "hn", 2, bins, xmin, xmax);

## Filling
A THn is filled just like a regular histogram, using
THn::Fill(x, weight), where x is a n-dimensional Double_t value.
To take errors into account, Sumw2() must be called before filling the
histogram.
Storage is allocated when the first bin content is stored.

## Projections
The dimensionality of a THn can be reduced by projecting it to
1, 2, 3, or n dimensions, which can be represented by a TH1, TH2, TH3, or
a THn. See the Projection() members. To only project parts of the
histogram, call

    hn->GetAxis(12)->SetRange(from_bin, to_bin);

## Conversion from other histogram classes
The static factory function THn::CreateHn() can be used to create a THn
from a TH1, TH2, TH3, THnSparse and (for copying) even from a THn. The
created THn will have compatble storage type, i.e. calling CreateHn() on
a TH2F will create a THnF.
*/

ClassImp(THn);

////////////////////////////////////////////////////////////////////////////////
/// Construct a THn.

THn::THn(const char* name, const char* title,
         Int_t dim, const Int_t* nbins,
         const Double_t* xmin, const Double_t* xmax):
   THnBase(name, title, dim, nbins, xmin, xmax),
   fSumw2(dim, nbins, kTRUE /*overflow*/),
   fCoordBuf() {
}

////////////////////////////////////////////////////////////////////////////////
/// Destruct a THn

THn::~THn()
{
   delete [] fCoordBuf;
}


////////////////////////////////////////////////////////////////////////////////
/// Create an iterator over all bins. Public interface is THnIter.

ROOT::Internal::THnBaseBinIter* THn::CreateIter(Bool_t respectAxisRange) const
{
   return new THnBinIter(GetNdimensions(), GetListOfAxes(), &GetArray(),
                         respectAxisRange);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable calculation of errors

void THn::Sumw2() {
   if (!GetCalculateErrors()) {
      fTsumw2 = 0.;
   }
   // fill sumw2 array with current content
   TNDArray & content = GetArray();
   Long64_t nbins = GetNbins();
   for (Long64_t ibin = 0; ibin < nbins; ++ibin)
      fSumw2.At(ibin) = content.AtAsDouble(ibin);
}


////////////////////////////////////////////////////////////////////////////////
/// Create the coordinate buffer. Outlined to hide allocation
/// from inlined functions.

void THn::AllocCoordBuf() const
{
   fCoordBuf = new Int_t[fNdimensions]();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the storage of a histogram created via Init()

void THn::InitStorage(Int_t* nbins, Int_t /*chunkSize*/)
{
   fCoordBuf = new Int_t[fNdimensions]();
   GetArray().Init(fNdimensions, nbins, true /*addOverflow*/);
   fSumw2.Init(fNdimensions, nbins, true /*addOverflow*/);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the contents of a THn.

void THn::Reset(Option_t* option /*= ""*/)
{
   GetArray().Reset(option);
   fSumw2.Reset(option);
}
