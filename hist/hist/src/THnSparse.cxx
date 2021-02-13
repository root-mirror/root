// @(#)root/hist:$Id$
// Author: Axel Naumann (2007-09-11)

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THnSparse.h"

#include "TAxis.h"
#include "TClass.h"
#include "TDataMember.h"
#include "TDataType.h"

namespace {
//______________________________________________________________________________
//
// THnSparseBinIter iterates over all filled bins of a THnSparse.
//______________________________________________________________________________

   class THnSparseBinIter: public ROOT::Internal::THnBaseBinIter {
   public:
      THnSparseBinIter(Bool_t respectAxisRange, const THnSparse* hist):
         ROOT::Internal::THnBaseBinIter(respectAxisRange), fHist(hist),
         fNbins(hist->GetNbins()), fIndex(-1) {
         // Construct a THnSparseBinIter
         fCoord = new Int_t[hist->GetNdimensions()];
         fCoord[0] = -1;
      }
      virtual ~THnSparseBinIter() { delete [] fCoord; }

      virtual Int_t GetCoord(Int_t dim) const;
      virtual Long64_t Next(Int_t* coord = 0);

   private:
      THnSparseBinIter(const THnSparseBinIter&); // intentionally unimplemented
      THnSparseBinIter& operator=(const THnSparseBinIter&); // intentionally unimplemented

      const THnSparse* fHist;
      Int_t* fCoord; // coord buffer for fIndex; fCoord[0] == -1 if not yet calculated
      Long64_t fNbins; // number of bins to iterate over
      Long64_t fIndex; // current bin index
   };
}

Int_t THnSparseBinIter::GetCoord(Int_t dim) const
{
   if (fCoord[0] == -1) {
      fHist->GetBinContent(fIndex, fCoord);
   }
   return fCoord[dim];
}

Long64_t THnSparseBinIter::Next(Int_t* coord /*= 0*/)
{
   // Get next bin index (in range if RespectsAxisRange()).
   // If coords != 0, set it to the index's axis coordinates
   // (i.e. coord must point to an array of Int_t[fNdimension]
   if (!fHist) return -1;

   fCoord[0] = -1;
   Int_t* useCoordBuf = fCoord;
   if (coord) {
      useCoordBuf = coord;
      coord[0] = -1;
   }

   do {
      ++fIndex;
      if (fIndex >= fHist->GetNbins()) {
         fHist = 0;
         return -1;
      }
      if (RespectsAxisRange()) {
         fHist->GetBinContent(fIndex, useCoordBuf);
      }
   } while (RespectsAxisRange()
            && !fHist->IsInRange(useCoordBuf)
            && (fHaveSkippedBin = kTRUE /* assignment! */));

   if (coord && coord[0] == -1) {
      if (fCoord[0] == -1) {
         fHist->GetBinContent(fIndex, coord);
      } else {
         memcpy(coord, fCoord, fHist->GetNdimensions() * sizeof(Int_t));
      }
   }

   return fIndex;
}



/** \class THnSparseCoordCompression
THnSparseCoordCompression is a class used by THnSparse internally. It
represents a compacted n-dimensional array of bin coordinates (indices).
As the total number of bins in each dimension is known by THnSparse, bin
indices can be compacted to only use the amount of bins needed by the total
number of bins in each dimension. E.g. for a THnSparse with
{15, 100, 2, 20, 10, 100} bins per dimension, a bin index will only occupy
28 bits (4+7+1+5+4+7), i.e. less than a 32bit integer. The tricky part is
the fast compression and decompression, the platform-independent storage
(think of endianness: the bits of the number 0x123456 depend on the
platform), and the hashing needed by THnSparseArrayChunk.
*/


class THnSparseCoordCompression {
public:
   THnSparseCoordCompression(Int_t dim, const Int_t* nbins);
   THnSparseCoordCompression(const THnSparseCoordCompression& other);
   ~THnSparseCoordCompression();

   THnSparseCoordCompression& operator=(const THnSparseCoordCompression& other);

   ULong64_t GetHashFromBuffer(const Char_t* buf) const;
   Int_t     GetBufferSize() const { return fCoordBufferSize; }
   Int_t     GetNdimensions() const { return fNdimensions; }
   void      SetCoordFromBuffer(const Char_t* buf_in, Int_t* coord_out) const;
   ULong64_t SetBufferFromCoord(const Int_t* coord_in, Char_t* buf_out) const;

protected:
   Int_t     GetNumBits(Int_t n) const {
      // return the number of bits allocated by the number "n"
      Int_t r = (n > 0);
      while (n/=2) ++r;
      return r;
   }
private:
   Int_t  fNdimensions;     // number of dimensions
   Int_t  fCoordBufferSize; // size of coordbuf
   Int_t *fBitOffsets;      //[fNdimensions + 1] bit offset of each axis index
};


//______________________________________________________________________________
//______________________________________________________________________________


////////////////////////////////////////////////////////////////////////////////
/// Initialize a THnSparseCoordCompression object with "dim" dimensions
/// and "bins" holding the number of bins for each dimension; it
/// stores the

THnSparseCoordCompression::THnSparseCoordCompression(Int_t dim, const Int_t* nbins):
   fNdimensions(dim), fCoordBufferSize(0), fBitOffsets(0)
{
   fBitOffsets = new Int_t[dim + 1];

   int shift = 0;
   for (Int_t i = 0; i < dim; ++i) {
      fBitOffsets[i] = shift;
      shift += GetNumBits(nbins[i] + 2);
   }
   fBitOffsets[dim] = shift;
   fCoordBufferSize = (shift + 7) / 8;
}


////////////////////////////////////////////////////////////////////////////////
/// Construct a THnSparseCoordCompression from another one

THnSparseCoordCompression::THnSparseCoordCompression(const THnSparseCoordCompression& other)
{
   fNdimensions = other.fNdimensions;
   fCoordBufferSize = other.fCoordBufferSize;
   fBitOffsets = new Int_t[fNdimensions + 1];
   memcpy(fBitOffsets, other.fBitOffsets, sizeof(Int_t) * fNdimensions);
}


////////////////////////////////////////////////////////////////////////////////
/// Set this to other if different.

THnSparseCoordCompression& THnSparseCoordCompression::operator=(const THnSparseCoordCompression& other)
{
   if (&other == this) return *this;

   fNdimensions = other.fNdimensions;
   fCoordBufferSize = other.fCoordBufferSize;
   delete [] fBitOffsets;
   fBitOffsets = new Int_t[fNdimensions + 1];
   memcpy(fBitOffsets, other.fBitOffsets, sizeof(Int_t) * fNdimensions);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// destruct a THnSparseCoordCompression

THnSparseCoordCompression::~THnSparseCoordCompression()
{
   delete [] fBitOffsets;
}


////////////////////////////////////////////////////////////////////////////////
/// Given the compressed coordinate buffer buf_in, calculate ("decompact")
/// the bin coordinates and return them in coord_out.

void THnSparseCoordCompression::SetCoordFromBuffer(const Char_t* buf_in,
                                                  Int_t* coord_out) const
{
   for (Int_t i = 0; i < fNdimensions; ++i) {
      const Int_t offset = fBitOffsets[i] / 8;
      Int_t shift = fBitOffsets[i] % 8;
      Int_t nbits = fBitOffsets[i + 1] - fBitOffsets[i];
      const UChar_t* pbuf = (const UChar_t*) buf_in + offset;
      coord_out[i] = *pbuf >> shift;
      Int_t subst = (Int_t) -1;
      subst = subst << nbits;
      nbits -= (8 - shift);
      shift = 8 - shift;
      for (Int_t n = 0; n * 8 < nbits; ++n) {
         ++pbuf;
         coord_out[i] += *pbuf << shift;
         shift += 8;
      }
      coord_out[i] &= ~subst;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Given the cbin coordinates coord_in, calculate ("compact")
/// the bin coordinates and return them in buf_in.
/// Return the hash value.

ULong64_t THnSparseCoordCompression::SetBufferFromCoord(const Int_t* coord_in,
                                                       Char_t* buf_out) const
{
   if (fCoordBufferSize <= 8) {
      ULong64_t l64buf = 0;
      for (Int_t i = 0; i < fNdimensions; ++i) {
         l64buf += ((ULong64_t)((UInt_t)coord_in[i])) << fBitOffsets[i];
      }
      memcpy(buf_out, &l64buf, sizeof(Long64_t));
      return l64buf;
   }

   // else: doesn't fit into a Long64_t:
   memset(buf_out, 0, fCoordBufferSize);
   for (Int_t i = 0; i < fNdimensions; ++i) {
      const Int_t offset = fBitOffsets[i] / 8;
      const Int_t shift = fBitOffsets[i] % 8;
      ULong64_t val = coord_in[i];

      Char_t* pbuf = buf_out + offset;
      *pbuf += 0xff & (val << shift);
      val = val >> (8 - shift);
      while (val) {
         ++pbuf;
         *pbuf += 0xff & val;
         val = val >> 8;
      }
   }

   return GetHashFromBuffer(buf_out);
}

/*
////////////////////////////////////////////////////////////////////////////////
/// Calculate hash from bin indexes.

ULong64_t THnSparseCoordCompression::GetHashFromCoords(const Int_t* coord) const
{
   // Bins are addressed in two different modes, depending
   // on whether the compact bin index fits into a Long64_t or not.
   // If it does, we can use it as a "perfect hash" for the TExMap.
   // If not we build a hash from the compact bin index, and use that
   // as the TExMap's hash.

   if (fCoordBufferSize <= 8) {
      // fits into a Long64_t
      ULong64_t hash1 = 0;
      for (Int_t i = 0; i < fNdimensions; ++i) {
         hash1 += coord[i] << fBitOffsets[i];
      }
      return hash1;
   }

   // else: doesn't fit into a Long64_t:
   memset(coord, 0, fCoordBufferSize);
   for (Int_t i = 0; i < fNdimensions; ++i) {
      const Int_t offset = fBitOffsets[i] / 8;
      const Int_t shift = fBitOffsets[i] % 8;
      ULong64_t val = coord[i];

      Char_t* pbuf = fCoordBuffer + offset;
      *pbuf += 0xff & (val << shift);
      val = val >> (8 - shift);
      while (val) {
         ++pbuf;
         *pbuf += 0xff & val;
         val = val >> 8;
      }
   }

   ULong64_t hash = 5381;
   Char_t* str = fCoordBuffer;
   while (str - fCoordBuffer < fCoordBufferSize) {
      hash *= 5;
      hash += *(str++);
   }
   return hash;
}
*/


////////////////////////////////////////////////////////////////////////////////
/// Calculate hash from compact bin index.

ULong64_t THnSparseCoordCompression::GetHashFromBuffer(const Char_t* buf) const
{
   // Bins are addressed in two different modes, depending
   // on whether the compact bin index fits into a Long64_t or not.
   // If it does, we can use it as a "perfect hash" for the TExMap.
   // If not we build a hash from the compact bin index, and use that
   // as the TExMap's hash.

   if (fCoordBufferSize <= 8) {
      // fits into a Long64_t
      ULong64_t hash1 = 0;
      memcpy(&hash1, buf, fCoordBufferSize);
      return hash1;
   }

   // else: doesn't fit into a Long64_t:
   ULong64_t hash = 5381;
   const Char_t* str = buf;
   while (str - buf < fCoordBufferSize) {
      hash *= 5;
      hash += *(str++);
   }
   return hash;
}




/** \class THnSparseCompactBinCoord
THnSparseCompactBinCoord is a class used by THnSparse internally. It
maps between an n-dimensional array of bin coordinates (indices) and
its compact version, the THnSparseCoordCompression.
*/

class THnSparseCompactBinCoord: public THnSparseCoordCompression {
public:
   THnSparseCompactBinCoord(Int_t dim, const Int_t* nbins);
   ~THnSparseCompactBinCoord();
   Int_t*    GetCoord() { return fCurrentBin; }
   const Char_t*   GetBuffer() const { return fCoordBuffer; }
   ULong64_t GetHash() const { return fHash; }
   void UpdateCoord() {
      fHash = SetBufferFromCoord(fCurrentBin, fCoordBuffer);
   }
   void      SetCoord(const Int_t* coord) {
      memcpy(fCurrentBin, coord, sizeof(Int_t) * GetNdimensions());
      fHash = SetBufferFromCoord(coord, fCoordBuffer);
   }
   void      SetBuffer(const Char_t* buf) {
      memcpy(fCoordBuffer, buf, GetBufferSize());
      fHash = GetHashFromBuffer(fCoordBuffer);
   }

private:
   // intentionally not implemented
   THnSparseCompactBinCoord(const THnSparseCompactBinCoord&);
   // intentionally not implemented
   THnSparseCompactBinCoord& operator=(const THnSparseCompactBinCoord&);

private:
   ULong64_t fHash;      // hash for current coordinates; 0 if not calculated
   Char_t *fCoordBuffer; // compact buffer of coordinates
   Int_t *fCurrentBin;   // current coordinates
};


//______________________________________________________________________________
//______________________________________________________________________________


////////////////////////////////////////////////////////////////////////////////
/// Initialize a THnSparseCompactBinCoord object with "dim" dimensions
/// and "bins" holding the number of bins for each dimension.

THnSparseCompactBinCoord::THnSparseCompactBinCoord(Int_t dim, const Int_t* nbins):
   THnSparseCoordCompression(dim, nbins),
   fHash(0), fCoordBuffer(0), fCurrentBin(0)
{
   fCurrentBin = new Int_t[dim];
   size_t bufAllocSize = GetBufferSize();
   if (bufAllocSize < sizeof(Long64_t))
      bufAllocSize = sizeof(Long64_t);
   fCoordBuffer = new Char_t[bufAllocSize];
}


////////////////////////////////////////////////////////////////////////////////
/// destruct a THnSparseCompactBinCoord

THnSparseCompactBinCoord::~THnSparseCompactBinCoord()
{
   delete [] fCoordBuffer;
   delete [] fCurrentBin;
}

/** \class THnSparseArrayChunk
THnSparseArrayChunk is used internally by THnSparse.
THnSparse stores its (dynamic size) array of bin coordinates and their
contents (and possibly errors) in a TObjArray of THnSparseArrayChunk. Each
of the chunks holds an array of THnSparseCompactBinCoord and the content
(a TArray*), which is created outside (by the templated derived classes of
THnSparse) and passed in at construction time.
*/

ClassImp(THnSparseArrayChunk);

////////////////////////////////////////////////////////////////////////////////
/// (Default) initialize a chunk. Takes ownership of cont (~THnSparseArrayChunk deletes it),
/// and create an ArrayF for errors if "errors" is true.

THnSparseArrayChunk::THnSparseArrayChunk(Int_t coordsize, bool errors, TArray* cont):
      fCoordinateAllocationSize(-1), fSingleCoordinateSize(coordsize), fCoordinatesSize(0),
      fCoordinates(0), fContent(cont),
      fSumw2(0)
{
   fCoordinateAllocationSize = fSingleCoordinateSize * cont->GetSize();
   fCoordinates = new Char_t[fCoordinateAllocationSize];
   if (errors) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

THnSparseArrayChunk::~THnSparseArrayChunk()
{
   delete fContent;
   delete [] fCoordinates;
   delete fSumw2;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new bin in this chunk

void THnSparseArrayChunk::AddBin(Int_t idx, const Char_t* coordbuf)
{
   // When streaming out only the filled chunk is saved.
   // When reading back only the memory needed for that filled part gets
   // allocated. We need to check whether the allowed chunk size is
   // bigger than the allocated size. If fCoordinateAllocationSize is
   // set to -1 this chunk has been allocated by the  streamer and the
   // buffer allocation size is defined by [fCoordinatesSize]. In that
   // case we need to compare fCoordinatesSize to
   // fSingleCoordinateSize * fContent->GetSize()
   // to determine whether we need to expand the buffer.
   if (fCoordinateAllocationSize == -1 && fContent) {
      Int_t chunksize = fSingleCoordinateSize * fContent->GetSize();
      if (fCoordinatesSize < chunksize) {
         // need to re-allocate:
         Char_t *newcoord = new Char_t[chunksize];
         memcpy(newcoord, fCoordinates, fCoordinatesSize);
         delete [] fCoordinates;
         fCoordinates = newcoord;
      }
      fCoordinateAllocationSize = chunksize;
   }

   memcpy(fCoordinates + idx * fSingleCoordinateSize, coordbuf, fSingleCoordinateSize);
   fCoordinatesSize += fSingleCoordinateSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on support of errors

void THnSparseArrayChunk::Sumw2()
{
   if (!fSumw2)
      fSumw2 = new TArrayD(fContent->GetSize());
   // fill the structure with the current content
   for (Int_t bin=0; bin < fContent->GetSize(); bin++) {
      fSumw2->fArray[bin] = fContent->GetAt(bin);
   }

}


/** \class THnSparse
    \ingroup Hist

Efficient multidimensional histogram.

Use a THnSparse instead of TH1 / TH2 / TH3 / array for histogramming when
only a small fraction of bins is filled. A 10-dimensional histogram with 10
bins per dimension has 10^10 bins; in a naive implementation this will not
fit in memory. THnSparse only allocates memory for the bins that have
non-zero bin content instead, drastically reducing both the memory usage
and the access time.

To construct a THnSparse object you must use one of its templated, derived
classes:
- THnSparseD (typedef for THnSparseT<ArrayD>): bin content held by a Double_t,
- THnSparseF (typedef for THnSparseT<ArrayF>): bin content held by a Float_t,
- THnSparseL (typedef for THnSparseT<ArrayL>): bin content held by a Long_t,
- THnSparseI (typedef for THnSparseT<ArrayI>): bin content held by an Int_t,
- THnSparseS (typedef for THnSparseT<ArrayS>): bin content held by a Short_t,
- THnSparseC (typedef for THnSparseT<ArrayC>): bin content held by a Char_t,

They take name and title, the number of dimensions, and for each dimension
the number of bins, the minimal, and the maximal value on the dimension's
axis. A TH2 h("h","h",10, 0., 10., 20, -5., 5.) would correspond to

    Int_t bins[2] = {10, 20};
    Double_t xmin[2] = {0., -5.};
    Double_t xmax[2] = {10., 5.};
    THnSparseD hs("hs", "hs", 2, bins, xmin, xmax);

## Filling
A THnSparse is filled just like a regular histogram, using
THnSparse::Fill(x, weight), where x is a n-dimensional Double_t value.
To take errors into account, Sumw2() must be called before filling the
histogram.

Bins are allocated as needed; the status of the allocation can be observed
by GetSparseFractionBins(), GetSparseFractionMem().

## Fast Bin Content Access
When iterating over a THnSparse one should only look at filled bins to save
processing time. The number of filled bins is returned by
THnSparse::GetNbins(); the bin content for each (linear) bin number can
be retrieved by THnSparse::GetBinContent(linidx, (Int_t*)coord).
After the call, coord will contain the bin coordinate of each axis for the bin
with linear index linidx. A possible call would be

   std::cout << hs.GetBinContent(0, coord);
   std::cout <<" is the content of bin [x = " << coord[0] "
        << " | y = " << coord[1] << "]" << std::endl;

## Efficiency
TH1 and TH2 are generally faster than THnSparse for one and two dimensional
distributions. THnSparse becomes competitive for a sparsely filled TH3
with large numbers of bins per dimension. The tutorial sparsehist.C
shows the turning point. On a AMD64 with 8GB memory, THnSparse "wins"
starting with a TH3 with 30 bins per dimension. Using a THnSparse for a
one-dimensional histogram is only reasonable if it has a huge number of bins.

## Projections
The dimensionality of a THnSparse can be reduced by projecting it to
1, 2, 3, or n dimensions, which can be represented by a TH1, TH2, TH3, or
a THnSparse. See the Projection() members. To only project parts of the
histogram, call

    THnSparse::GetAxis(12)->SetRange(from_bin, to_bin);

## Internal Representation
An entry for a filled bin consists of its n-dimensional coordinates and
its bin content. The coordinates are compacted to use as few bits as
possible; e.g. a histogram with 10 bins in x and 20 bins in y will only
use 4 bits for the x representation and 5 bits for the y representation.
This is handled by the internal class THnSparseCompactBinCoord.
Bin data (content and coordinates) are allocated in chunks of size
fChunkSize; this parameter can be set when constructing a THnSparse. Each
chunk is represented by an object of class THnSparseArrayChunk.

Translation from an n-dimensional bin coordinate to the linear index within
the chunks is done by GetBin(). It creates a hash from the compacted bin
coordinates (the hash of a bin coordinate is the compacted coordinate itself
if it takes less than 8 bytes, the size of a Long64_t.
This hash is used to lookup the linear index in the TExMap member fBins;
the coordinates of the entry fBins points to is compared to the coordinates
passed to GetBin(). If they do not match, these two coordinates have the same
hash - which is extremely unlikely but (for the case where the compact bin
coordinates are larger than 4 bytes) possible. In this case, fBinsContinued
contains a chain of linear indexes with the same hash. Iterating through this
chain and comparing each bin coordinates with the one passed to GetBin() will
retrieve the matching bin.
*/


ClassImp(THnSparse);

////////////////////////////////////////////////////////////////////////////////
/// Construct an empty THnSparse.

THnSparse::THnSparse():
   fChunkSize(1024), fFilledBins(0), fCompactCoord(0)
{
   fBinContent.SetOwner();
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a THnSparse with "dim" dimensions,
/// with chunksize as the size of the chunks.
/// "nbins" holds the number of bins for each dimension;
/// "xmin" and "xmax" the minimal and maximal value for each dimension.
/// The arrays "xmin" and "xmax" can be NULL; in that case SetBinEdges()
/// must be called for each dimension.

THnSparse::THnSparse(const char* name, const char* title, Int_t dim,
                     const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
                     Int_t chunksize):
   THnBase(name, title, dim, nbins, xmin, xmax),
   fChunkSize(chunksize), fFilledBins(0), fCompactCoord(0)
{
   fCompactCoord = new THnSparseCompactBinCoord(dim, nbins);
   fBinContent.SetOwner();
}

////////////////////////////////////////////////////////////////////////////////
/// Destruct a THnSparse

THnSparse::~THnSparse() {
   delete fCompactCoord;
}

////////////////////////////////////////////////////////////////////////////////
/// Add "v" to the content of bin with index "bin"

void THnSparse::AddBinContent(Long64_t bin, Double_t v)
{
   THnSparseArrayChunk* chunk = GetChunk(bin / fChunkSize);
   bin %= fChunkSize;
   v += chunk->fContent->GetAt(bin);
   return chunk->fContent->SetAt(v, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new chunk of bin content

THnSparseArrayChunk* THnSparse::AddChunk()
{
   THnSparseArrayChunk* chunk =
      new THnSparseArrayChunk(GetCompactCoord()->GetBufferSize(),
                              GetCalculateErrors(), GenerateArray());
   fBinContent.AddLast(chunk);
   return chunk;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the storage of a histogram created via Init()

void THnSparse::InitStorage(Int_t* nbins, Int_t chunkSize)
{
   fChunkSize = chunkSize;
   fCompactCoord = new THnSparseCompactBinCoord(fNdimensions, nbins);
}

////////////////////////////////////////////////////////////////////////////////
///We have been streamed; set up fBins

void THnSparse::FillExMap()
{
   TIter iChunk(&fBinContent);
   THnSparseArrayChunk* chunk = 0;
   THnSparseCoordCompression compactCoord(*GetCompactCoord());
   Long64_t idx = 0;
   if (2 * GetNbins() > fBins.Capacity())
      fBins.Expand(3 * GetNbins());
   while ((chunk = (THnSparseArrayChunk*) iChunk())) {
      const Int_t chunkSize = chunk->GetEntries();
      Char_t* buf = chunk->fCoordinates;
      const Int_t singleCoordSize = chunk->fSingleCoordinateSize;
      const Char_t* endbuf = buf + singleCoordSize * chunkSize;
      for (; buf < endbuf; buf += singleCoordSize, ++idx) {
         Long64_t hash = compactCoord.GetHashFromBuffer(buf);
         Long64_t linidx = fBins.GetValue(hash);
         if (linidx) {
            Long64_t nextidx = linidx;
            while (nextidx) {
               // must be a collision, so go to fBinsContinued.
               linidx = nextidx;
               nextidx = fBinsContinued.GetValue(linidx);
            }
            fBinsContinued.Add(linidx, idx + 1);
         } else {
            fBins.Add(hash, idx + 1);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize storage for nbins

void THnSparse::Reserve(Long64_t nbins) {
   if (!fBins.GetSize() && fBinContent.GetSize()) {
      FillExMap();
   }
   if (2 * nbins > fBins.Capacity()) {
      fBins.Expand(3 * nbins);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the bin index for the n dimensional tuple x,
/// allocate one if it doesn't exist yet and "allocate" is true.

Long64_t THnSparse::GetBin(const Double_t* x, Bool_t allocate /* = kTRUE */)
{
   THnSparseCompactBinCoord* cc = GetCompactCoord();
   Int_t *coord = cc->GetCoord();
   for (Int_t i = 0; i < fNdimensions; ++i)
      coord[i] = GetAxis(i)->FindBin(x[i]);
   cc->UpdateCoord();

   return GetBinIndexForCurrentBin(allocate);
}


////////////////////////////////////////////////////////////////////////////////
/// Get the bin index for the n dimensional tuple addressed by "name",
/// allocate one if it doesn't exist yet and "allocate" is true.

Long64_t THnSparse::GetBin(const char* name[], Bool_t allocate /* = kTRUE */)
{
   THnSparseCompactBinCoord* cc = GetCompactCoord();
   Int_t *coord = cc->GetCoord();
   for (Int_t i = 0; i < fNdimensions; ++i)
      coord[i] = GetAxis(i)->FindBin(name[i]);
   cc->UpdateCoord();

   return GetBinIndexForCurrentBin(allocate);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the bin index for the n dimensional coordinates coord,
/// allocate one if it doesn't exist yet and "allocate" is true.

Long64_t THnSparse::GetBin(const Int_t* coord, Bool_t allocate /*= kTRUE*/)
{
   GetCompactCoord()->SetCoord(coord);
   return GetBinIndexForCurrentBin(allocate);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the content of the filled bin number "idx".
/// If coord is non-null, it will contain the bin's coordinates for each axis
/// that correspond to the bin.

Double_t THnSparse::GetBinContent(Long64_t idx, Int_t* coord /* = 0 */) const
{
   if (idx >= 0) {
      THnSparseArrayChunk* chunk = GetChunk(idx / fChunkSize);
      idx %= fChunkSize;
      if (chunk && chunk->fContent->GetSize() > idx) {
         if (coord) {
            THnSparseCompactBinCoord* cc = GetCompactCoord();
            Int_t sizeCompact = cc->GetBufferSize();
            cc->SetCoordFromBuffer(chunk->fCoordinates + idx * sizeCompact,
                                                  coord);

         }
         return chunk->fContent->GetAt(idx);
      }
   }
   if (coord)
      memset(coord, -1, sizeof(Int_t) * fNdimensions);
   return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Get square of the error of bin addressed by linidx as
/// \f$\sum weight^{2}\f$
/// If errors are not enabled (via Sumw2() or CalculateErrors())
/// return contents.

Double_t THnSparse::GetBinError2(Long64_t linidx) const {
   if (!GetCalculateErrors())
      return GetBinContent(linidx);

   if (linidx < 0) return 0.;
   THnSparseArrayChunk* chunk = GetChunk(linidx / fChunkSize);
   linidx %= fChunkSize;
   if (!chunk || chunk->fContent->GetSize() < linidx)
      return 0.;

   return chunk->fSumw2->GetAt(linidx);
}


////////////////////////////////////////////////////////////////////////////////
/// Return the index for fCurrentBinIndex.
/// If it doesn't exist then return -1, or allocate a new bin if allocate is set

Long64_t THnSparse::GetBinIndexForCurrentBin(Bool_t allocate)
{
   THnSparseCompactBinCoord* cc = GetCompactCoord();
   ULong64_t hash = cc->GetHash();
   if (fBinContent.GetSize() && !fBins.GetSize())
      FillExMap();
   Long64_t linidx = (Long64_t) fBins.GetValue(hash);
   while (linidx) {
      // fBins stores index + 1!
      THnSparseArrayChunk* chunk = GetChunk((linidx - 1)/ fChunkSize);
      if (chunk->Matches((linidx - 1) % fChunkSize, cc->GetBuffer()))
         return linidx - 1; // we store idx+1, 0 is "TExMap: not found"

      Long64_t nextlinidx = fBinsContinued.GetValue(linidx);
      if (!nextlinidx) break;

      linidx = nextlinidx;
   }
   if (!allocate) return -1;

   ++fFilledBins;

   // allocate bin in chunk
   THnSparseArrayChunk *chunk = (THnSparseArrayChunk*) fBinContent.Last();
   Long64_t newidx = chunk ? ((Long64_t) chunk->GetEntries()) : -1;
   if (!chunk || newidx == (Long64_t)fChunkSize) {
      chunk = AddChunk();
      newidx = 0;
   }
   chunk->AddBin(newidx, cc->GetBuffer());

   // store translation between hash and bin
   newidx += (fBinContent.GetEntriesFast() - 1) * fChunkSize;
   if (!linidx) {
      // fBins didn't find it
      if (2 * GetNbins() > fBins.Capacity())
         fBins.Expand(3 * GetNbins());
      fBins.Add(hash, newidx + 1);
   } else {
      // fBins contains one, but it's the wrong one;
      // add entry to fBinsContinued.
      fBinsContinued.Add(linidx, newidx + 1);
   }
   return newidx;
}

////////////////////////////////////////////////////////////////////////////////
/// Return THnSparseCompactBinCoord object.

THnSparseCompactBinCoord* THnSparse::GetCompactCoord() const
{
   if (!fCompactCoord) {
      Int_t *bins = new Int_t[fNdimensions];
      for (Int_t d = 0; d < fNdimensions; ++d)
         bins[d] = GetAxis(d)->GetNbins();
      const_cast<THnSparse*>(this)->fCompactCoord
         = new THnSparseCompactBinCoord(fNdimensions, bins);
      delete [] bins;
   }
   return fCompactCoord;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of filled bins over all bins

Double_t THnSparse::GetSparseFractionBins() const {
   Double_t nbinsTotal = 1.;
   for (Int_t d = 0; d < fNdimensions; ++d)
      nbinsTotal *= GetAxis(d)->GetNbins() + 2;
   return fFilledBins / nbinsTotal;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the amount of used memory over memory that would be used by a
/// non-sparse n-dimensional histogram. The value is approximate.

Double_t THnSparse::GetSparseFractionMem() const {
   Int_t arrayElementSize = 0;
   if (fFilledBins) {
      TClass* clArray = GetChunk(0)->fContent->IsA();
      TDataMember* dm = clArray ? clArray->GetDataMember("fArray") : 0;
      arrayElementSize = dm ? dm->GetDataType()->Size() : 0;
   }
   if (!arrayElementSize) {
      Warning("GetSparseFractionMem", "Cannot determine type of elements!");
      return -1.;
   }

   Double_t sizePerChunkElement = arrayElementSize + GetCompactCoord()->GetBufferSize();
   if (fFilledBins && GetChunk(0)->fSumw2)
      sizePerChunkElement += sizeof(Double_t); /* fSumw2 */

   Double_t size = 0.;
   size += fBinContent.GetEntries() * (GetChunkSize() * sizePerChunkElement + sizeof(THnSparseArrayChunk));
   size += + 3 * sizeof(Long64_t) * fBins.GetSize() /* TExMap */;

   Double_t nbinsTotal = 1.;
   for (Int_t d = 0; d < fNdimensions; ++d)
      nbinsTotal *= GetAxis(d)->GetNbins() + 2;

   return size / nbinsTotal / arrayElementSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an iterator over all filled bins of a THnSparse.
/// Use THnIter instead.

ROOT::Internal::THnBaseBinIter* THnSparse::CreateIter(Bool_t respectAxisRange) const
{
   return new THnSparseBinIter(respectAxisRange, this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set content of bin with index "bin" to "v"

void THnSparse::SetBinContent(Long64_t bin, Double_t v)
{
   THnSparseArrayChunk* chunk = GetChunk(bin / fChunkSize);
   chunk->fContent->SetAt(v, bin % fChunkSize);
   ++fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Set error of bin with index "bin" to "e", enable errors if needed

void THnSparse::SetBinError2(Long64_t bin, Double_t e2)
{
   THnSparseArrayChunk* chunk = GetChunk(bin / fChunkSize);
   if (!chunk->fSumw2 ) {
      // if fSumw2 is zero GetCalculateErrors should return false
      if (GetCalculateErrors()) {
         Error("SetBinError", "GetCalculateErrors() logic error!");
      }
      Sumw2(); // enable error calculation
   }

   chunk->fSumw2->SetAt(e2, bin % fChunkSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Add "e" to error of bin with index "bin", enable errors if needed

void THnSparse::AddBinError2(Long64_t bin, Double_t e2)
{
   THnSparseArrayChunk* chunk = GetChunk(bin / fChunkSize);
   if (!chunk->fSumw2 ) {
      // if fSumw2 is zero GetCalculateErrors should return false
      if (GetCalculateErrors()) {
         Error("SetBinError", "GetCalculateErrors() logic error!");
      }
      Sumw2(); // enable error calculation
   }

   (*chunk->fSumw2)[bin % fChunkSize] += e2;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable calculation of errors

void THnSparse::Sumw2()
{
   if (GetCalculateErrors()) return;

   fTsumw2 = 0.;
   TIter iChunk(&fBinContent);
   THnSparseArrayChunk* chunk = 0;
   while ((chunk = (THnSparseArrayChunk*) iChunk()))
      chunk->Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the histogram

void THnSparse::Reset(Option_t *option /*= ""*/)
{
   fFilledBins = 0;
   fBins.Delete();
   fBinsContinued.Clear();
   fBinContent.Delete();
   ResetBase(option);
}

