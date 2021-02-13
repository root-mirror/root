// @(#)root/hist:$Id$
// TH2Poly v2.1
// Author: Olivier Couet, Deniz Gunceler, Danilo Piparo

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TH2Poly.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "Riostream.h"
#include "TList.h"
#include "TMath.h"
#include <cassert>

ClassImp(TH2Poly);

/** \class TH2Poly
    \ingroup Hist
2D Histogram with Polygonal Bins

## Overview
`TH2Poly` is a 2D Histogram class (TH2) allowing to define polygonal
bins of arbitrary shape.

Each bin in the `TH2Poly` histogram is a `TH2PolyBin` object.
`TH2PolyBin` is a very simple class containing the vertices (stored
as `TGraph`s or `TMultiGraph`s ) and contents of the polygonal
bin as well as several related functions.

Essentially, a `TH2Poly` is a TList of `TH2PolyBin` objects
with methods to manipulate them.

Bins are defined using one of the `AddBin()` methods. The bin definition
should be done before filling.

The histogram can be filled with `Fill(Double_t x, Double_t y, Double_t w)
`. `w` is the weight.
If no weight is specified, it is assumed to be 1.

Not all histogram's area need to be binned. Filling an area without bins,
will falls into the overflows. Adding a bin is not retroactive; it doesn't
affect previous fillings. A `Fill()` call, that
was previously ignored due to the lack of a bin at the specified location, is
not reconsidered when that location is binned later.

If there are two overlapping bins, the first one in the list will be incremented
by `Fill()`.

The histogram may automatically extends its limits if a bin outside the
histogram limits is added. This is done when the default constructor (with no
arguments) is used. It generates a histogram with no limits along the X and Y
axis. Adding bins to it will extend it up to a proper size.

`TH2Poly` implements a partitioning algorithm to speed up bins' filling
(see the "Partitioning Algorithm" section for details).
The partitioning algorithm divides the histogram into regions called cells.
The bins that each cell intersects are recorded in an array of `TList`s.
When a coordinate in the histogram is to be filled; the method (quickly) finds
which cell the coordinate belongs.  It then only loops over the bins
intersecting that cell to find the bin the input coordinate corresponds to.
The partitioning of the histogram is updated continuously as each bin is added.
The default number of cells on each axis is 25. This number could be set to
another value in the constructor or adjusted later by calling the
`ChangePartition(Int_t, Int_t)` method. The partitioning algorithm is
considerably faster than the brute force algorithm (i.e. checking if each bin
contains the input coordinates), especially if the histogram is to be filled
many times.

The following very simple macro shows how to build and fill a `TH2Poly`:
~~~ {.cpp}
{
    TH2Poly *h2p = new TH2Poly();

    Double_t x1[] = {0, 5, 6};
    Double_t y1[] = {0, 0, 5};
    Double_t x2[] = {0, -1, -1, 0};
    Double_t y2[] = {0, 0, -1, 3};
    Double_t x3[] = {4, 3, 0, 1, 2.4};
    Double_t y3[] = {4, 3.7, 1, 3.7, 2.5};

    h2p->AddBin(3, x1, y1);
    h2p->AddBin(4, x2, y2);
    h2p->AddBin(5, x3, y3);

    h2p->Fill(0.1, 0.01, 3);
    h2p->Fill(-0.5, -0.5, 7);
    h2p->Fill(-0.7, -0.5, 1);
    h2p->Fill(1, 3, 1.5);
}
~~~

More examples can be found in th2polyBoxes.C, th2polyEurope.C, th2polyHoneycomb.C
and th2polyUSA.C.

## Partitioning Algorithm
The partitioning algorithm forms an essential part of the `TH2Poly`
class. It is implemented to speed up the filling of bins.

With the brute force approach, the filling is done in the following way:  An
iterator loops over all bins in the `TH2Poly` and invokes the
method `IsInside()` for each of them.
This method checks if the input location is in that bin. If the filling
coordinate is inside, the bin is filled. Looping over all the bin is
very slow.

The alternative is to divide the histogram into virtual rectangular regions
called "cells". Each cell stores the pointers of the bins intersecting it.
When a coordinate is to be filled, the method finds which cell the coordinate
falls into. Since the cells are rectangular, this can be done very quickly.
It then only loops over the bins associated with that cell and calls `IsInside()`
only on that bins. This reduces considerably the number of bins on which `IsInside()`
is called and therefore speed up by a huge factor the filling compare to the brute force
approach where `IsInside()` is called for all bins.

The addition of bins to the appropriate cells is done when the bin is added
to the histogram. To do this, `AddBin()` calls the
`AddBinToPartition()` method.
This method adds the input bin to the partitioning matrix.

The number of partition cells per axis can be specified in the constructor.
If it is not specified, the default value of 25 along each axis will be
assigned. This value was chosen because it is small enough to avoid slowing
down AddBin(), while being large enough to enhance Fill() by a considerable
amount. Regardless of how it is initialized at construction time, it can be
changed later with the `ChangePartition()` method.
`ChangePartition()` deletes the
old partition matrix and generates a new one with the specified number of cells
on each axis.

The optimum number of partition cells per axis changes with the number of
times `Fill()` will be called.  Although partitioning greatly speeds up
filling, it also adds a constant time delay into the code. When `Fill()`
is to be called many times, it is more efficient to divide the histogram into
a large number cells. However, if the histogram is to be filled only a few
times, it is better to divide into a small number of cells.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default Constructor. No boundaries specified.

TH2Poly::TH2Poly()
{
   Initialize(0., 0., 0., 0., 25, 25);
   SetName("NoName");
   SetTitle("NoTitle");
   SetFloat();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with specified name and boundaries,
/// but no partition cell number.

TH2Poly::TH2Poly(const char *name,const char *title, Double_t xlow,Double_t xup
                                             , Double_t ylow,Double_t yup)
{
   Initialize(xlow, xup, ylow, yup, 25, 25);
   SetName(name);
   SetTitle(title);
   SetFloat(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with specified name and boundaries and partition cell number.

TH2Poly::TH2Poly(const char *name,const char *title,
           Int_t nX, Double_t xlow, Double_t xup,
           Int_t nY, Double_t ylow, Double_t yup)
{
   Initialize(xlow, xup, ylow, yup, nX, nY);
   SetName(name);
   SetTitle(title);
   SetFloat(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2Poly::~TH2Poly()
{
   delete[] fCells;
   delete[] fIsEmpty;
   delete[] fCompletelyInside;
   // delete at the end the bin List since it owns the objects
   delete fBins;
}

////////////////////////////////////////////////////////////////////////////////
/// Create appropriate histogram bin.
///  e.g. TH2Poly        creates TH2PolyBin,
///       TProfile2Poly  creates TProfile2PolyBin
/// This is done so that TH2Poly::AddBin does not have to be duplicated,
/// but only create needs to be reimplemented for additional histogram types

TH2PolyBin *TH2Poly::CreateBin(TObject *poly)
{
   if (!poly) return 0;

   if (fBins == 0) {
      fBins = new TList();
      fBins->SetOwner();
   }

   fNcells++;
   Int_t ibin = fNcells - kNOverflow;
   // if structure fsumw2 is created extend it
   if (fSumw2.fN) fSumw2.Set(fNcells);
   return new TH2PolyBin(poly, ibin);
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a new bin to the histogram. It can be any object having the method
/// IsInside(). It returns the bin number in the histogram. It returns 0 if
/// it failed to add. To allow the histogram limits to expand when a bin
/// outside the limits is added, call SetFloat() before adding the bin.

Int_t TH2Poly::AddBin(TObject *poly)
{
   auto *bin = CreateBin(poly);
   Int_t ibin = fNcells-kNOverflow;
   if(!bin) return 0;

   // If the bin lies outside histogram boundaries, then extends the boundaries.
   // Also changes the partition information accordingly
   Bool_t flag = kFALSE;
   if (fFloat) {
      if (fXaxis.GetXmin() > bin->GetXMin()) {
         fXaxis.Set(100, bin->GetXMin(), fXaxis.GetXmax());
         flag = kTRUE;
      }
      if (fXaxis.GetXmax() < bin->GetXMax()) {
         fXaxis.Set(100, fXaxis.GetXmin(), bin->GetXMax());
         flag = kTRUE;
      }
      if (fYaxis.GetXmin() > bin->GetYMin()) {
         fYaxis.Set(100, bin->GetYMin(), fYaxis.GetXmax());
         flag = kTRUE;
      }
      if (fYaxis.GetXmax() < bin->GetYMax()) {
         fYaxis.Set(100, fYaxis.GetXmin(), bin->GetYMax());
         flag = kTRUE;
      }
      if (flag) ChangePartition(fCellX, fCellY);
   } else {
      /*Implement polygon clipping code here*/
   }

   fBins->Add((TObject*) bin);
   SetNewBinAdded(kTRUE);

   // Adds the bin to the partition matrix
   AddBinToPartition(bin);

   return ibin;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a new bin to the histogram. The number of vertices and their (x,y)
/// coordinates are required as input. It returns the bin number in the
/// histogram.

Int_t TH2Poly::AddBin(Int_t n, const Double_t *x, const Double_t *y)
{
   TGraph *g = new TGraph(n, x, y);
   Int_t bin = AddBin(g);
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new bin to the histogram. The bin shape is a rectangle.
/// It returns the bin number of the bin in the histogram.

Int_t TH2Poly::AddBin(Double_t x1, Double_t y1, Double_t x2, Double_t  y2)
{
   Double_t x[] = {x1, x1, x2, x2, x1};
   Double_t y[] = {y1, y2, y2, y1, y1};
   TGraph *g = new TGraph(5, x, y);
   Int_t bin = AddBin(g);
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this + c1*h1.

Bool_t TH2Poly::Add(const TH1 *h1, Double_t c1)
{
   Int_t bin;

   TH2Poly *h1p = (TH2Poly *)h1;

   // Check if number of bins is the same.
   if (h1p->GetNumberOfBins() != GetNumberOfBins()) {
      Error("Add", "Attempt to add histograms with different number of bins");
      return kFALSE;
   }

   // Check if the bins are the same.
   TList *h1pBins = h1p->GetBins();
   TH2PolyBin *thisBin, *h1pBin;
   for (bin = 1; bin <= GetNumberOfBins(); bin++) {
      thisBin = (TH2PolyBin *)fBins->At(bin - 1);
      h1pBin  = (TH2PolyBin *)h1pBins->At(bin - 1);
      if (thisBin->GetXMin() != h1pBin->GetXMin() ||
            thisBin->GetXMax() != h1pBin->GetXMax() ||
            thisBin->GetYMin() != h1pBin->GetYMin() ||
            thisBin->GetYMax() != h1pBin->GetYMax()) {
         Error("Add", "Attempt to add histograms with different bin limits");
         return kFALSE;
      }
   }


   // Create Sumw2 if h1p has Sumw2 set
   if (fSumw2.fN == 0 && h1p->GetSumw2N() != 0) Sumw2();

   // statistics can be preserved only in case of positive coefficients
   // otherwise with negative c1 (histogram subtraction) one risks to get negative variances
   Bool_t resetStats = (c1 < 0);
   Double_t s1[kNstat] = {0};
   Double_t s2[kNstat] = {0};
   if (!resetStats) {
      // need to initialize to zero s1 and s2 since
      // GetStats fills only used elements depending on dimension and type
      GetStats(s1);
      h1->GetStats(s2);
   }
   //   get number of entries now because afterwards UpdateBinContent will change it
   Double_t entries = TMath::Abs( GetEntries() + c1 * h1->GetEntries() );


   // Perform the Add.
   Double_t factor = 1;
   if (h1p->GetNormFactor() != 0)
      factor = h1p->GetNormFactor() / h1p->GetSumOfWeights();
   for (bin = 0; bin < fNcells; bin++) {
      Double_t y = this->RetrieveBinContent(bin) + c1 * h1p->RetrieveBinContent(bin);
      UpdateBinContent(bin, y);
      if (fSumw2.fN) {
         Double_t esq = factor * factor * h1p->GetBinErrorSqUnchecked(bin);
         fSumw2.fArray[bin] += c1 * c1 * factor * factor * esq;
      }
   }

   // update statistics (do here to avoid changes by SetBinContent)
   if (resetStats)  {
      // statistics need to be reset in case coefficient are negative
      ResetStats();
   } else {
      for (Int_t i = 0; i < kNstat; i++) {
         if (i == 1) s1[i] += c1 * c1 * s2[i];
         else        s1[i] += c1 * s2[i];
      }
      PutStats(s1);
      SetEntries(entries);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds the input bin into the partition cell matrix. This method is called
/// in AddBin() and ChangePartition().

void TH2Poly::AddBinToPartition(TH2PolyBin *bin)
{
   // Cell Info
   Int_t nl, nr, mb, mt; // Max/min indices of the cells that contain the bin
   Double_t xclipl, xclipr, yclipb, yclipt; // x and y coordinates of a cell
   Double_t binXmax, binXmin, binYmax, binYmin; // The max/min bin coordinates

   binXmax = bin->GetXMax();
   binXmin = bin->GetXMin();
   binYmax = bin->GetYMax();
   binYmin = bin->GetYMin();
   nl = (Int_t)(floor((binXmin - fXaxis.GetXmin())/fStepX));
   nr = (Int_t)(floor((binXmax - fXaxis.GetXmin())/fStepX));
   mb = (Int_t)(floor((binYmin - fYaxis.GetXmin())/fStepY));
   mt = (Int_t)(floor((binYmax - fYaxis.GetXmin())/fStepY));

   // Make sure the array indices are correct.
   if (nr>=fCellX) nr = fCellX-1;
   if (mt>=fCellY) mt = fCellY-1;
   if (nl<0)       nl = 0;
   if (mb<0)       mb = 0;

   // number of cells in the grid
   //N.B. not to be confused with fNcells (the number of bins) !
   fNCells = fCellX*fCellY;

   // Loop over all cells
   for (int i = nl; i <= nr; i++) {
      xclipl = fXaxis.GetXmin() + i*fStepX;
      xclipr = xclipl + fStepX;
      for (int j = mb; j <= mt; j++) {
         yclipb = fYaxis.GetXmin() + j*fStepY;
         yclipt = yclipb + fStepY;

         // If the bin is completely inside the cell,
         // add that bin to the cell then return
         if ((binXmin >= xclipl) && (binXmax <= xclipr) &&
             (binYmax <= yclipt) && (binYmin >= yclipb)){
            fCells[i + j*fCellX].Add((TObject*) bin);
            fIsEmpty[i + j*fCellX] = kFALSE;  // Makes the cell non-empty
            return;
         }

         // If any of the sides of the cell intersect with any side of the bin,
         // add that bin then continue
         if (IsIntersecting(bin, xclipl, xclipr, yclipb, yclipt)) {
            fCells[i + j*fCellX].Add((TObject*) bin);
            fIsEmpty[i + j*fCellX] = kFALSE;  // Makes the cell non-empty
            continue;
         }
         // If a corner of the cell is inside the bin and since there is no
         // intersection, then that cell completely inside the bin.
         if((bin->IsInside(xclipl,yclipb)) || (bin->IsInside(xclipl,yclipt))){
            fCells[i + j*fCellX].Add((TObject*) bin);
            fIsEmpty[i + j*fCellX] = kFALSE;  // Makes the cell non-empty
            fCompletelyInside[i + fCellX*j] = kTRUE;
            continue;
         }
         if((bin->IsInside(xclipr,yclipb)) || (bin->IsInside(xclipr,yclipt))){
            fCells[i + j*fCellX].Add((TObject*) bin);
            fIsEmpty[i + j*fCellX] = kFALSE;  // Makes the cell non-empty
            fCompletelyInside[i + fCellX*j] = kTRUE;
            continue;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the number of partition cells in the histogram.
/// Deletes the old partition and constructs a new one.

void TH2Poly::ChangePartition(Int_t n, Int_t m)
{
   fCellX = n;                          // Set the number of cells
   fCellY = m;                          // Set the number of cells

   delete [] fCells;                    // Deletes the old partition

   // number of cells in the grid
   //N.B. not to be confused with fNcells (the number of bins) !
   fNCells = fCellX*fCellY;
   fCells  = new TList [fNCells];  // Sets an empty partition

   fStepX = (fXaxis.GetXmax() - fXaxis.GetXmin())/fCellX;
   fStepY = (fYaxis.GetXmax() - fYaxis.GetXmin())/fCellY;

   delete [] fIsEmpty;
   delete [] fCompletelyInside;
   fIsEmpty = new Bool_t [fNCells];
   fCompletelyInside = new Bool_t [fNCells];

   // Initializes the flags
   for (int i = 0; i<fNCells; i++) {
      fIsEmpty[i]          = kTRUE;
      fCompletelyInside[i] = kFALSE;
   }

   // TList iterator
   TIter    next(fBins);
   TObject  *obj;

   while((obj = next())){   // Loop over bins and add them to the partition
      AddBinToPartition((TH2PolyBin*) obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make a complete copy of the underlying object.  If 'newname' is set,
/// the copy's name will be set to that name.

TObject* TH2Poly::Clone(const char* newname) const
{
   // TH1::Clone relies on ::Copy to implemented by the derived class.
   // Until this is implemented, revert to the much slower default version
   // (and possibly non-thread safe).

   return TNamed::Clone(newname);
}

////////////////////////////////////////////////////////////////////////////////
/// Clears the contents of all bins in the histogram.

void TH2Poly::ClearBinContents()
{
   TIter next(fBins);
   TObject *obj;
   TH2PolyBin *bin;

   // Clears the bin contents
   while ((obj = next())) {
      bin = (TH2PolyBin*) obj;
      bin->ClearContent();
   }

   // Clears the statistics
   fTsumw   = 0;
   fTsumw2  = 0;
   fTsumwx  = 0;
   fTsumwx2 = 0;
   fTsumwy  = 0;
   fTsumwy2 = 0;
   fEntries = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2Poly::Reset(Option_t *opt)
{
   TIter next(fBins);
   TObject *obj;
   TH2PolyBin *bin;

   // Clears the bin contents
   while ((obj = next())) {
      bin = (TH2PolyBin*) obj;
      bin->ClearContent();
   }

   TH2::Reset(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the bin number of the bin at the given coordinate. -1 to -9 are
/// the overflow and underflow bins.  overflow bin -5 is the unbinned areas in
/// the histogram (also called "the sea"). The third parameter can be left
/// blank.
/// The overflow/underflow bins are:
///~~~ {.cpp}
/// -1 | -2 | -3
/// -------------
/// -4 | -5 | -6
/// -------------
/// -7 | -8 | -9
///~~~
/// where -5 means is the "sea" bin (i.e. unbinned areas)

Int_t TH2Poly::FindBin(Double_t x, Double_t y, Double_t)
{

   // Checks for overflow/underflow
   Int_t overflow = 0;
   if      (y > fYaxis.GetXmax()) overflow += -1;
   else if (y > fYaxis.GetXmin()) overflow += -4;
   else                           overflow += -7;
   if      (x > fXaxis.GetXmax()) overflow += -2;
   else if (x > fXaxis.GetXmin()) overflow += -1;
   if (overflow != -5) return overflow;

   // Finds the cell (x,y) coordinates belong to
   Int_t n = (Int_t)(floor((x-fXaxis.GetXmin())/fStepX));
   Int_t m = (Int_t)(floor((y-fYaxis.GetXmin())/fStepY));

   // Make sure the array indices are correct.
   if (n>=fCellX) n = fCellX-1;
   if (m>=fCellY) m = fCellY-1;
   if (n<0)       n = 0;
   if (m<0)       m = 0;

   if (fIsEmpty[n+fCellX*m]) return -5;

   TH2PolyBin *bin;

   TIter next(&fCells[n+fCellX*m]);
   TObject *obj;

   // Search for the bin in the cell
   while ((obj=next())) {
      bin  = (TH2PolyBin*)obj;
      if (bin->IsInside(x,y)) return bin->GetBinNumber();
   }

   // If the search has not returned a bin, the point must be on "the sea"
   return -5;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment the bin containing (x,y) by 1.
/// Uses the partitioning algorithm.

Int_t TH2Poly::Fill(Double_t x, Double_t y)
{
   return Fill(x, y, 1.0);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment the bin containing (x,y) by w.
/// Uses the partitioning algorithm.

Int_t TH2Poly::Fill(Double_t x, Double_t y, Double_t w)
{
   // see GetBinCOntent for definition of overflow bins
   // in case of weighted events store weight square in fSumw2.fArray
   // but with this indexing:
   // fSumw2.fArray[0:kNOverflow-1] : sum of weight squares for the overflow bins
   // fSumw2.fArray[kNOverflow:fNcells] : sum of weight squares for the standard bins
   // where fNcells = kNOverflow + Number of bins. kNOverflow=9

   if (fNcells <= kNOverflow) return 0;

   // create sum of weight square array if weights are different than 1
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW) )  Sumw2();

   Int_t overflow = 0;
   if      (y > fYaxis.GetXmax()) overflow += -1;
   else if (y > fYaxis.GetXmin()) overflow += -4;
   else                           overflow += -7;
   if      (x > fXaxis.GetXmax()) overflow += -2;
   else if(x > fXaxis.GetXmin())  overflow += -1;
   if (overflow != -5) {
      fOverflow[-overflow - 1]+= w;
      if (fSumw2.fN) fSumw2.fArray[-overflow - 1] += w*w;
      return overflow;
   }

   // Finds the cell (x,y) coordinates belong to
   Int_t n = (Int_t)(floor((x-fXaxis.GetXmin())/fStepX));
   Int_t m = (Int_t)(floor((y-fYaxis.GetXmin())/fStepY));

   // Make sure the array indices are correct.
   if (n>=fCellX) n = fCellX-1;
   if (m>=fCellY) m = fCellY-1;
   if (n<0)       n = 0;
   if (m<0)       m = 0;

   if (fIsEmpty[n+fCellX*m]) {
      fOverflow[4]+= w;
      if (fSumw2.fN) fSumw2.fArray[4] += w*w;
      return -5;
   }

   TH2PolyBin *bin;
   Int_t bi;

   TIter next(&fCells[n+fCellX*m]);
   TObject *obj;

   while ((obj=next())) {
      bin  = (TH2PolyBin*)obj;
      // needs to account offset in array for overflow bins
      bi = bin->GetBinNumber()-1+kNOverflow;
      if (bin->IsInside(x,y)) {
         bin->Fill(w);

         // Statistics
         fTsumw   = fTsumw + w;
         fTsumw2  = fTsumw2 + w*w;
         fTsumwx  = fTsumwx + w*x;
         fTsumwx2 = fTsumwx2 + w*x*x;
         fTsumwy  = fTsumwy + w*y;
         fTsumwy2 = fTsumwy2 + w*y*y;
         if (fSumw2.fN) {
            assert(bi < fSumw2.fN);
            fSumw2.fArray[bi] += w*w;
         }
         fEntries++;

         SetBinContentChanged(kTRUE);

         return bin->GetBinNumber();
      }
   }

   fOverflow[4]+= w;
   if (fSumw2.fN) fSumw2.fArray[4] += w*w;
   return -5;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment the bin named "name" by w.

Int_t TH2Poly::Fill(const char* name, Double_t w)
{
   TString sname(name);

   TIter    next(fBins);
   TObject  *obj;
   TH2PolyBin *bin;

   while ((obj = next())) {
      bin = (TH2PolyBin*) obj;
      if (!sname.CompareTo(bin->GetPolygon()->GetName())) {
         bin->Fill(w);
         fEntries++;
         SetBinContentChanged(kTRUE);
         return bin->GetBinNumber();
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a 2-D histogram with an array of values and weights.
///
/// \param [in] ntimes:  number of entries in arrays x and w
///                      (array size must be ntimes*stride)
/// \param [in] x:       array of x values to be histogrammed
/// \param [in] y:       array of y values to be histogrammed
/// \param [in] w:       array of weights
/// \param [in] stride:  step size through arrays x, y and w

void TH2Poly::FillN(Int_t ntimes, const Double_t* x, const Double_t* y,
                               const Double_t* w, Int_t stride)
{
   for (int i = 0; i < ntimes; i += stride) {
      Fill(x[i], y[i], w[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the integral of bin contents.
/// By default the integral is computed as the sum of bin contents.
/// If option "width" or "area" is specified, the integral is the sum of
/// the bin contents multiplied by the area of the bin.

Double_t TH2Poly::Integral(Option_t* option) const
{
   TString opt = option;
   opt.ToLower();

   Double_t w;
   Double_t integral = 0.;

   TIter next(fBins);
   TObject *obj;
   TH2PolyBin *bin;
   if ((opt.Contains("width")) || (opt.Contains("area"))) {
      while ((obj = next())) {
         bin = (TH2PolyBin *)obj;
         w = bin->GetArea();
         integral += w * (bin->GetContent());
      }
   } else {
      // need to recompute integral in case somebith called SetBinContent and
      // cannot use fTsumw since it is not updated in that case
      while ((obj = next())) {
         bin = (TH2PolyBin *)obj;
         integral += (bin->GetContent());
      }
   }
   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the content of the input bin
/// For the overflow/underflow/sea bins:
///~~~ {.cpp}
/// -1 | -2 | -3
/// ---+----+----
/// -4 | -5 | -6
/// ---+----+----
/// -7 | -8 | -9
///~~~
/// where -5 is the "sea" bin (i.e. unbinned areas)

Double_t TH2Poly::GetBinContent(Int_t bin) const
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -kNOverflow) return 0;
   if (bin<0) return fOverflow[-bin - 1];
   return ((TH2PolyBin*) fBins->At(bin-1))->GetContent();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the value of error associated to bin number bin.
/// If the sum of squares of weights has been defined (via Sumw2),
/// this function returns the sqrt(sum of w2).
/// otherwise it returns the sqrt(contents) for this bin.
/// Bins are in range [1:nbins] and for bin < 0 in range [-9:-1] it returns errors for overflow bins.
/// See also TH2Poly::GetBinContent

Double_t TH2Poly::GetBinError(Int_t bin) const
{
   if (bin == 0 || bin > GetNumberOfBins() || bin < - kNOverflow) return 0;
   if (fBuffer) ((TH1*)this)->BufferEmpty();
   // in case of weighted events the sum of the weights are stored in a different way than
   // a normal histogram
   // fSumw2.fArray[0:kNOverflow-1] : sum of weight squares for the overflow bins (
   // fSumw2.fArray[kNOverflow:fNcells] : sum of weight squares for the standard bins
   //  fNcells = kNOverflow (9) + Number of bins
   if (fSumw2.fN) {
      Int_t binIndex = (bin > 0) ? bin+kNOverflow-1 : -(bin+1);
      Double_t err2 = fSumw2.fArray[binIndex];
      return TMath::Sqrt(err2);
   }
   Double_t error2 = TMath::Abs(GetBinContent(bin));
   return TMath::Sqrt(error2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bin Error.
/// Re-implementation for TH2Poly given the different bin indexing in the
/// stored squared error array.
/// See also notes in TH1::SetBinError
///
/// Bins are in range [1:nbins] and for bin < 0 in the range [-9:-1] the  errors is set for the overflow bins


void TH2Poly::SetBinError(Int_t bin, Double_t error)
{
   if (bin == 0 || bin > GetNumberOfBins() || bin < - kNOverflow) return;
   if (!fSumw2.fN) Sumw2();
   SetBinErrorOption(kNormal);
   // see comment in GetBinError for special convention of bin index in fSumw2 array
   Int_t binIndex = (bin > 0) ? bin+kNOverflow-1 : -(bin+1);
   fSumw2.fArray[binIndex] = error * error;
}



////////////////////////////////////////////////////////////////////////////////
/// Returns the bin name.

const char *TH2Poly::GetBinName(Int_t bin) const
{
   if (bin > GetNumberOfBins())  return "";
   if (bin < 0)          return "";
   return ((TH2PolyBin*) fBins->At(bin-1))->GetPolygon()->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the bin title.

const char *TH2Poly::GetBinTitle(Int_t bin) const
{
   if (bin > GetNumberOfBins())  return "";
   if (bin < 0)          return "";
   return ((TH2PolyBin*) fBins->At(bin-1))->GetPolygon()->GetTitle();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value of the histogram.

Double_t TH2Poly::GetMaximum() const
{
   if (fNcells <= kNOverflow) return 0;
   if (fMaximum != -1111) return fMaximum;

   TH2PolyBin  *b;

   TIter next(fBins);
   TObject *obj;
   Double_t max,c;

   max = ((TH2PolyBin*) next())->GetContent();

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      c = b->GetContent();
      if (c>max) max = c;
   }
   return max;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value of the histogram that is less than maxval.

Double_t TH2Poly::GetMaximum(Double_t maxval) const
{
   if (fNcells <= kNOverflow) return 0;
   if (fMaximum != -1111) return fMaximum;

   TH2PolyBin  *b;

   TIter next(fBins);
   TObject *obj;
   Double_t max,c;

   max = ((TH2PolyBin*) next())->GetContent();

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      c = b->GetContent();
      if (c>max && c<maxval) max=c;
   }
   return max;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value of the histogram.

Double_t TH2Poly::GetMinimum() const
{
   if (fNcells <= kNOverflow) return 0;
   if (fMinimum != -1111) return fMinimum;

   TH2PolyBin  *b;

   TIter next(fBins);
   TObject *obj;
   Double_t min,c;

   min = ((TH2PolyBin*) next())->GetContent();

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      c = b->GetContent();
      if (c<min) min=c;
   }
   return min;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value of the histogram that is greater than minval.

Double_t TH2Poly::GetMinimum(Double_t minval) const
{
   if (fNcells <= kNOverflow) return 0;
   if (fMinimum != -1111) return fMinimum;

   TH2PolyBin  *b;

   TIter next(fBins);
   TObject *obj;
   Double_t min,c;

   min = ((TH2PolyBin*) next())->GetContent();

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      c = b->GetContent();
      if (c<min && c>minval) min=c;
   }
   return min;
}

////////////////////////////////////////////////////////////////////////////////
/// Bins the histogram using a honeycomb structure

void TH2Poly::Honeycomb(Double_t xstart, Double_t ystart, Double_t a,
                     Int_t k, Int_t s)
{
   // Add the bins
   Double_t numberOfHexagonsInTheRow;
   Double_t x[6], y[6];
   Double_t xloop, yloop, xtemp;
   xloop = xstart; yloop = ystart + a/2.0;
   for (int sCounter = 0; sCounter < s; sCounter++) {

      xtemp = xloop; // Resets the temp variable

      // Determine the number of hexagons in that row
      if(sCounter%2 == 0){numberOfHexagonsInTheRow = k;}
      else{numberOfHexagonsInTheRow = k - 1;}

      for (int kCounter = 0; kCounter <  numberOfHexagonsInTheRow; kCounter++) {

         // Go around the hexagon
         x[0] = xtemp;
         y[0] = yloop;
         x[1] = x[0];
         y[1] = y[0] + a;
         x[2] = x[1] + a*TMath::Sqrt(3)/2.0;
         y[2] = y[1] + a/2.0;
         x[3] = x[2] + a*TMath::Sqrt(3)/2.0;
         y[3] = y[1];
         x[4] = x[3];
         y[4] = y[0];
         x[5] = x[2];
         y[5] = y[4] - a/2.0;

         this->AddBin(6, x, y);

         // Go right
         xtemp += a*TMath::Sqrt(3);
      }

      // Increment the starting position
      if (sCounter%2 == 0) xloop += a*TMath::Sqrt(3)/2.0;
      else                 xloop -= a*TMath::Sqrt(3)/2.0;
      yloop += 1.5*a;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initializes the TH2Poly object.  This method is called by the constructor.

void TH2Poly::Initialize(Double_t xlow, Double_t xup,
                      Double_t ylow, Double_t yup, Int_t n, Int_t m)
{
   Int_t i;
   fDimension = 2;  //The dimension of the histogram

   fBins   = 0;
   fNcells = kNOverflow;

   // Sets the boundaries of the histogram
   fXaxis.Set(100, xlow, xup);
   fYaxis.Set(100, ylow, yup);

   for (i=0; i<9; i++) fOverflow[i] = 0.;

   // Statistics
   fEntries = 0;   // The total number of entries
   fTsumw   = 0.;  // Total amount of content in the histogram
   fTsumw2  = 0.;  // Sum square of the weights
   fTsumwx  = 0.;  // Weighted sum of x coordinates
   fTsumwx2 = 0.;  // Weighted sum of the squares of x coordinates
   fTsumwy2 = 0.;  // Weighted sum of the squares of y coordinates
   fTsumwy  = 0.;  // Weighted sum of y coordinates

   fCellX = n; // Set the number of cells to default
   fCellY = m; // Set the number of cells to default

   // number of cells in the grid
   //N.B. not to be confused with fNcells (the number of bins) !
   fNCells = fCellX*fCellY;
   fCells  = new TList [fNCells];  // Sets an empty partition
   fStepX  = (fXaxis.GetXmax() - fXaxis.GetXmin())/fCellX; // Cell width
   fStepY  = (fYaxis.GetXmax() - fYaxis.GetXmin())/fCellY; // Cell height

   fIsEmpty = new Bool_t [fNCells]; // Empty partition
   fCompletelyInside = new Bool_t [fNCells]; // Cell is completely inside bin

   for (i = 0; i<fNCells; i++) {   // Initializes the flags
      fIsEmpty[i] = kTRUE;
      fCompletelyInside[i] = kFALSE;
   }

   // 3D Painter flags
   SetNewBinAdded(kFALSE);
   SetBinContentChanged(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the input bin is intersecting with the
/// input rectangle (xclipl, xclipr, yclipb, yclipt)

Bool_t TH2Poly::IsIntersecting(TH2PolyBin *bin,
                               Double_t xclipl, Double_t xclipr,
                               Double_t yclipb, Double_t yclipt)
{
   Int_t     gn;
   Double_t *gx;
   Double_t *gy;
   Bool_t inter = kFALSE;
   TObject *poly = bin->GetPolygon();

   if (poly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)poly;
      gx = g->GetX();
      gy = g->GetY();
      gn = g->GetN();
      inter = IsIntersectingPolygon(gn, gx, gy, xclipl, xclipr, yclipb, yclipt);
   }

   if (poly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)poly;
      TList *gl = mg->GetListOfGraphs();
      if (!gl) return inter;
      TGraph *g;
      TIter next(gl);
      while ((g = (TGraph*) next())) {
         gx = g->GetX();
         gy = g->GetY();
         gn = g->GetN();
         inter = IsIntersectingPolygon(gn, gx, gy, xclipl, xclipr,
                                                   yclipb, yclipt);
         if (inter) return inter;
      }
   }

   return inter;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the input polygon (bn, x, y) is intersecting with the
/// input rectangle (xclipl, xclipr, yclipb, yclipt)

Bool_t TH2Poly::IsIntersectingPolygon(Int_t bn, Double_t *x, Double_t *y,
                                      Double_t xclipl, Double_t xclipr,
                                      Double_t yclipb, Double_t yclipt)
{
   Bool_t p0R, p0L, p0T, p0B, p0xM, p0yM, p1R, p1L, p1T;
   Bool_t p1B, p1xM, p1yM, p0In, p1In;

   for (int counter = 0; counter < (bn-1); counter++) {
      // If both are on the same side, return kFALSE
      p0L = x[counter]     <= xclipl; // Point 0 is on the left
      p1L = x[counter + 1] <= xclipl; // Point 1 is on the left
      if (p0L && p1L) continue;
      p0R = x[counter]     >= xclipr; // Point 0 is on the right
      p1R = x[counter + 1] >= xclipr; // Point 1 is on the right
      if (p0R && p1R) continue;
      p0T = y[counter]     >= yclipt; // Point 0 is at the top
      p1T = y[counter + 1] >= yclipt; // Point 1 is at the top
      if (p0T && p1T) continue;
      p0B = y[counter]     <= yclipb; // Point 0 is at the bottom
      p1B = y[counter + 1] <= yclipb; // Point 1 is at the bottom
      if (p0B && p1B) continue;

      // Checks to see if any are inside
      p0xM = !p0R && !p0L; // Point 0 is inside along x
      p0yM = !p0T && !p0B; // Point 1 is inside along x
      p1xM = !p1R && !p1L; // Point 0 is inside along y
      p1yM = !p1T && !p1B; // Point 1 is inside along y
      p0In = p0xM && p0yM; // Point 0 is inside
      p1In = p1xM && p1yM; // Point 1 is inside
      if (p0In) {
         if (p1In) continue;
         return kTRUE;
      } else {
         if (p1In) return kTRUE;
      }

      // We know by now that the points are not in the same side and not inside.

      // Checks to see if they are opposite

      if (p0xM && p1xM) return kTRUE;
      if (p0yM && p1yM) return kTRUE;

      // We now know that the points are in different x and y indices

      Double_t xcoord[3], ycoord[3];
      xcoord[0] = x[counter];
      xcoord[1] = x[counter + 1];
      ycoord[0] = y[counter];
      ycoord[1] = y[counter + 1];

      if (p0L) {
         if(p1T){
            xcoord[2] = xclipl;
            ycoord[2] = yclipb;
            if((TMath::IsInside(xclipl, yclipt, 3, xcoord, ycoord)) ||
               (TMath::IsInside(xclipr, yclipb, 3, xcoord, ycoord))) continue;
            else return kTRUE;
         } else if (p1B) {
            xcoord[2] = xclipl;
            ycoord[2] = yclipt;
            if((TMath::IsInside(xclipl, yclipb, 3, xcoord, ycoord)) ||
               (TMath::IsInside(xclipr, yclipt, 3, xcoord, ycoord))) continue;
            else return kTRUE;
         } else { // p1yM
            if (p0T) {
               xcoord[2] = xclipl;
               ycoord[2] = yclipb;
               if (TMath::IsInside(xclipr, yclipt, 3, xcoord, ycoord)) continue;
               else return kTRUE;
            }
            if (p0B) {
               xcoord[2] = xclipl;
               ycoord[2] = yclipt;
               if (TMath::IsInside(xclipr, yclipb, 3, xcoord, ycoord)) continue;
               else return kTRUE;
            }
         }
      } else if (p0R) {
         if (p1T) {
            xcoord[2] = xclipl;
            ycoord[2] = yclipb;
            if ((TMath::IsInside(xclipr, yclipb, 3, xcoord, ycoord)) ||
                (TMath::IsInside(xclipl, yclipt, 3, xcoord, ycoord))) continue;
            else return kTRUE;
         } else if (p1B) {
            xcoord[2] = xclipl;
            ycoord[2] = yclipt;
            if ((TMath::IsInside(xclipl, yclipb, 3, xcoord, ycoord)) ||
                (TMath::IsInside(xclipr, yclipt, 3, xcoord, ycoord))) continue;
            else return kTRUE;
         } else{ // p1yM
            if (p0T) {
               xcoord[2] = xclipr;
               ycoord[2] = yclipb;
               if (TMath::IsInside(xclipl, yclipt, 3, xcoord, ycoord)) continue;
               else return kTRUE;
            }
            if (p0B) {
               xcoord[2] = xclipr;
               ycoord[2] = yclipt;
               if (TMath::IsInside(xclipl, yclipb, 3, xcoord, ycoord)) continue;
               else return kTRUE;
            }
         }
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge TH2Polys
/// Given the special nature of the TH2Poly, the merge is implemented in
/// terms of subsequent TH2Poly::Add calls.
Long64_t TH2Poly::Merge(TCollection *coll)
{
   for (auto h2pAsObj : *coll) {
      if (!Add((TH1*)h2pAsObj, 1.)) {
         Warning("Merge", "An issue was encountered during the merge operation.");
         return 0L;
      }
   }
   return GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TH2Poly::SavePrimitive(std::ostream &out, Option_t *option)
{
   out <<"   "<<std::endl;
   out <<"   "<< ClassName() <<" *";

   //histogram pointer has by default the histogram name.
   //however, in case histogram has no directory, it is safer to add a
   //incremental suffix
   static Int_t hcounter = 0;
   TString histName = GetName();
   if (!fDirectory && !histName.Contains("Graph")) {
      hcounter++;
      histName += "__";
      histName += hcounter;
   }
   const char *hname = histName.Data();

   //Construct the class initialization
   out << hname << " = new " << ClassName() << "(\"" << hname << "\", \""
       << GetTitle() << "\", " << fCellX << ", " << fXaxis.GetXmin()
       << ", " << fXaxis.GetXmax()
       << ", " << fCellY << ", " << fYaxis.GetXmin() << ", "
       << fYaxis.GetXmax() << ");" << std::endl;

   // Save Bins
   TIter       next(fBins);
   TObject    *obj;
   TH2PolyBin *th2pBin;

   while((obj = next())){
      th2pBin = (TH2PolyBin*) obj;
      th2pBin->GetPolygon()->SavePrimitive(out,
                                           TString::Format("th2poly%s",histName.Data()));
   }

   // save bin contents
   out<<"   "<<std::endl;
   Int_t bin;
   for (bin=1;bin<=GetNumberOfBins();bin++) {
      Double_t bc = GetBinContent(bin);
      if (bc) {
         out<<"   "<<hname<<"->SetBinContent("<<bin<<","<<bc<<");"<<std::endl;
      }
   }

   // save bin errors
   if (fSumw2.fN) {
      for (bin=1;bin<=GetNumberOfBins();bin++) {
         Double_t be = GetBinError(bin);
         if (be) {
            out<<"   "<<hname<<"->SetBinError("<<bin<<","<<be<<");"<<std::endl;
         }
      }
   }
   TH1::SavePrimitiveHelp(out, hname, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this histogram by a constant c1.

void TH2Poly::Scale(Double_t c1, Option_t*)
{
   for( int i = 0; i < this->GetNumberOfBins(); i++ ) {
      this->SetBinContent(i+1, c1*this->GetBinContent(i+1));
   }
   for( int i = 0; i < kNOverflow; i++ ) {
      this->SetBinContent(-i-1, c1*this->GetBinContent(-i-1) );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the contents of the input bin to the input content
/// Negative values between -1 and -9 are for the overflows and the sea

void TH2Poly::SetBinContent(Int_t bin, Double_t content)
{
   if (bin > GetNumberOfBins() || bin == 0 || bin < -9 ) return;
   if (bin > 0) {
      ((TH2PolyBin*) fBins->At(bin-1))->SetContent(content);
   }
   else
      fOverflow[-bin - 1] = content;

   SetBinContentChanged(kTRUE);
   fEntries++;
}

////////////////////////////////////////////////////////////////////////////////
/// When set to kTRUE, allows the histogram to expand if a bin outside the
/// limits is added.

void TH2Poly::SetFloat(Bool_t flag)
{
   fFloat = flag;
}

////////////////////////////////////////////////////////////////////////////////
/// Return "true" if the point (x,y) is inside the bin of binnr.

Bool_t TH2Poly::IsInsideBin(Int_t binnr, Double_t x, Double_t y)
{
   if (!fBins) return false;
   TH2PolyBin* bin = (TH2PolyBin*)fBins->At(binnr);
   if (!bin) return false;
   return bin->IsInside(x,y);
}

void TH2Poly::GetStats(Double_t *stats) const
{
   stats[0] = fTsumw;
   stats[1] = fTsumw2;
   stats[2] = fTsumwx;
   stats[3] = fTsumwx2;
   stats[4] = fTsumwy;
   stats[5] = fTsumwy2;
   stats[6] = fTsumwxy;
}

/** \class TH2PolyBin
    \ingroup Hist
Helper class to represent a bin in the TH2Poly histogram
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TH2PolyBin::TH2PolyBin()
{
   fPoly    = 0;
   fContent = 0.;
   fNumber  = 0;
   fXmax    = -1111;
   fXmin    = -1111;
   fYmax    = -1111;
   fYmin    = -1111;
   fArea    = 0;
   SetChanged(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TH2PolyBin::TH2PolyBin(TObject *poly, Int_t bin_number)
{
   fContent = 0.;
   fNumber  = bin_number;
   fArea    = 0.;
   fPoly    = poly;
   fXmax    = -1111;
   fXmin    = -1111;
   fYmax    = -1111;
   fYmin    = -1111;
   SetChanged(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2PolyBin::~TH2PolyBin()
{
   if (fPoly) delete fPoly;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the area of the bin.

Double_t TH2PolyBin::GetArea()
{
   Int_t     bn;

   if (fArea == 0) {
      if (fPoly->IsA() == TGraph::Class()) {
         TGraph *g = (TGraph*)fPoly;
         bn    = g->GetN();
         fArea = g->Integral(0,bn-1);
      }

      if (fPoly->IsA() == TMultiGraph::Class()) {
         TMultiGraph *mg = (TMultiGraph*)fPoly;
         TList *gl = mg->GetListOfGraphs();
         if (!gl) return fArea;
         TGraph *g;
         TIter next(gl);
         while ((g = (TGraph*) next())) {
            bn    = g->GetN();
            fArea = fArea + g->Integral(0,bn-1);
         }
      }
   }

   return fArea;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value for the x coordinates of the bin.

Double_t TH2PolyBin::GetXMax()
{
   if (fXmax != -1111) return fXmax;

   Int_t     bn,i;
   Double_t *bx;

   if (fPoly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)fPoly;
      bx    = g->GetX();
      bn    = g->GetN();
      fXmax = bx[0];
      for (i=1; i<bn; i++) {if (fXmax < bx[i]) fXmax = bx[i];}
   }

   if (fPoly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)fPoly;
      TList *gl = mg->GetListOfGraphs();
      if (!gl) return fXmax;
      TGraph *g;
      TIter next(gl);
      Bool_t first = kTRUE;
      while ((g = (TGraph*) next())) {
         bx = g->GetX();
         bn = g->GetN();
         if (first) {fXmax = bx[0]; first = kFALSE;}
         for (i=0; i<bn; i++) {if (fXmax < bx[i]) fXmax = bx[i];}
      }
   }

   return fXmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value for the x coordinates of the bin.

Double_t TH2PolyBin::GetXMin()
{
   if (fXmin != -1111) return fXmin;

   Int_t     bn,i;
   Double_t *bx;

   if (fPoly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)fPoly;
      bx    = g->GetX();
      bn    = g->GetN();
      fXmin = bx[0];
      for (i=1; i<bn; i++) {if (fXmin > bx[i]) fXmin = bx[i];}
   }

   if (fPoly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)fPoly;
      TList *gl = mg->GetListOfGraphs();
      if (!gl) return fXmin;
      TGraph *g;
      TIter next(gl);
      Bool_t first = kTRUE;
      while ((g = (TGraph*) next())) {
         bx = g->GetX();
         bn = g->GetN();
         if (first) {fXmin = bx[0]; first = kFALSE;}
         for (i=0; i<bn; i++) {if (fXmin > bx[i]) fXmin = bx[i];}
      }
   }

   return fXmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value for the y coordinates of the bin.

Double_t TH2PolyBin::GetYMax()
{
   if (fYmax != -1111) return fYmax;

   Int_t     bn,i;
   Double_t *by;

   if (fPoly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)fPoly;
      by    = g->GetY();
      bn    = g->GetN();
      fYmax = by[0];
      for (i=1; i<bn; i++) {if (fYmax < by[i]) fYmax = by[i];}
   }

   if (fPoly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)fPoly;
      TList *gl = mg->GetListOfGraphs();
      if (!gl) return fYmax;
      TGraph *g;
      TIter next(gl);
      Bool_t first = kTRUE;
      while ((g = (TGraph*) next())) {
         by = g->GetY();
         bn = g->GetN();
         if (first) {fYmax = by[0]; first = kFALSE;}
         for (i=0; i<bn; i++) {if (fYmax < by[i]) fYmax = by[i];}
      }
   }

   return fYmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value for the y coordinates of the bin.

Double_t TH2PolyBin::GetYMin()
{
   if (fYmin != -1111) return fYmin;

   Int_t     bn,i;
   Double_t *by;

   if (fPoly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)fPoly;
      by    = g->GetY();
      bn    = g->GetN();
      fYmin = by[0];
      for (i=1; i<bn; i++) {if (fYmin > by[i]) fYmin = by[i];}
   }

   if (fPoly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)fPoly;
      TList *gl = mg->GetListOfGraphs();
      if (!gl) return fYmin;
      TGraph *g;
      TIter next(gl);
      Bool_t first = kTRUE;
      while ((g = (TGraph*) next())) {
         by = g->GetY();
         bn = g->GetN();
         if (first) {fYmin = by[0]; first = kFALSE;}
         for (i=0; i<bn; i++) {if (fYmin > by[i]) fYmin = by[i];}
      }
   }

   return fYmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Return "true" if the point (x,y) is inside the bin.

Bool_t TH2PolyBin::IsInside(Double_t x, Double_t y) const
{
   Int_t in=0;

   if (fPoly->IsA() == TGraph::Class()) {
      TGraph *g = (TGraph*)fPoly;
      in = g->IsInside(x, y);
   }

   if (fPoly->IsA() == TMultiGraph::Class()) {
      TMultiGraph *mg = (TMultiGraph*)fPoly;
      in = mg->IsInside(x, y);
   }

   return in;
}

////////////////////////////////////////////////////////////////////////
/// RE-implement dummy functions to avoid users calling the
/// corresponding implemntations in TH1 or TH2
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this + c1*f1. NOT IMPLEMENTED for TH2Poly
Bool_t TH2Poly::Add(TF1 *, Double_t, Option_t *)
{
   Error("Add","Not implement for TH2Poly");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this histogram by the addition of h1 and h2. NOT IMPLEMENTED for TH2Poly
Bool_t TH2Poly::Add(const TH1 *, const TH1 *, Double_t, Double_t)
{
   Error("Add","Not implement for TH2Poly");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this / c1*f1. NOT IMPLEMENTED for TH2Poly
Bool_t TH2Poly::Divide(TF1 *, Double_t)
{
   Error("Divide","Not implement for TH2Poly");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED for TH2Poly
Bool_t TH2Poly::Multiply(TF1 *, Double_t)
{
   Error("Multiply","Not implement for TH2Poly");
   return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED for TH2Poly
Double_t TH2Poly::ComputeIntegral(Bool_t )
{
   Error("ComputeIntegral","Not implement for TH2Poly");
   return TMath::QuietNaN();
}
////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED for TH2Poly
TH1 * TH2Poly::FFT(TH1*, Option_t * )
{
   Error("FFT","Not implement for TH2Poly");
   return nullptr;
}
////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED for TH2Poly
TH1 * TH2Poly::GetAsymmetry(TH1* , Double_t,  Double_t)
{
   Error("GetAsymmetry","Not implement for TH2Poly");
   return nullptr;
}
////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED for TH2Poly
Double_t TH2Poly::Interpolate(Double_t, Double_t)
{
   Error("Interpolate","Not implement for TH2Poly");
   return TMath::QuietNaN();
}
