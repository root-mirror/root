// @(#)root/hist:$Id$
// Author: Axel Naumann (2007-09-11)

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THnSparse
#define ROOT_THnSparse

/*************************************************************************

 THnSparse: histogramming multi-dimensional, sparse distributions in
 a memory-efficient way.

*************************************************************************/


#include "THnBase.h"
#include "TExMap.h"
#include "THnSparse_Internal.h"

// needed only for template instantiations of THnSparseT:
#include "TArrayF.h"
#include "TArrayL.h"
#include "TArrayI.h"
#include "TArrayS.h"
#include "TArrayC.h"

class THnSparseCompactBinCoord;

class THnSparse: public THnBase {
 private:
   Int_t      fChunkSize;                   ///<  Number of entries for each chunk
   Long64_t   fFilledBins;                  ///<  Number of filled bins
   TObjArray  fBinContent;                  ///<  Array of THnSparseArrayChunk
   TExMap     fBins;                        ///<! Filled bins
   TExMap     fBinsContinued;               ///<! Filled bins for non-unique hashes, containing pairs of (bin index 0, bin index 1)
   THnSparseCompactBinCoord *fCompactCoord; ///<! Compact coordinate

   THnSparse(const THnSparse&); // Not implemented
   THnSparse& operator=(const THnSparse&); // Not implemented

 protected:

   THnSparse();
   THnSparse(const char* name, const char* title, Int_t dim,
             const Int_t* nbins, const Double_t* xmin, const Double_t* xmax,
             Int_t chunksize);
   THnSparseCompactBinCoord* GetCompactCoord() const;
   THnSparseArrayChunk* GetChunk(Int_t idx) const {
      return (THnSparseArrayChunk*) fBinContent[idx]; }

   THnSparseArrayChunk* AddChunk();
   void Reserve(Long64_t nbins);
   void FillExMap();
   virtual TArray* GenerateArray() const = 0;
   Long64_t GetBinIndexForCurrentBin(Bool_t allocate);

   /// Increment the bin content of "bin" by "w",
   /// return the bin index.
   void FillBin(Long64_t bin, Double_t w) {
      THnSparseArrayChunk* chunk = GetChunk(bin / fChunkSize);
      chunk->AddBinContent(bin % fChunkSize, w);
      FillBinBase(w);
   }
   void InitStorage(Int_t* nbins, Int_t chunkSize);

 public:
   virtual ~THnSparse();

   static THnSparse* CreateSparse(const char* name, const char* title,
                                  const TH1* h1, Int_t chunkSize = 1024 * 16) {
      return (THnSparse*) CreateHnAny(name, title, h1, kTRUE /*sparse*/,
                                       chunkSize);
   }
   static THnSparse* CreateSparse(const char* name, const char* title,
                                  const THnBase* hn, Int_t chunkSize = 1024 * 16) {
      return (THnSparse*) CreateHnAny(name, title, hn, kTRUE /*sparse*/,
                                       chunkSize);
   }

   Int_t GetChunkSize() const { return fChunkSize; }
   Int_t GetNChunks() const { return fBinContent.GetEntriesFast(); }

   ROOT::Internal::THnBaseBinIter* CreateIter(Bool_t respectAxisRange) const;

   Long64_t GetNbins() const { return fFilledBins; }
   void SetFilledBins(Long64_t nbins) { fFilledBins = nbins; }

   Long64_t GetBin(const Int_t* idx) const { return const_cast<THnSparse*>(this)->GetBin(idx, kFALSE); }
   Long64_t GetBin(const Double_t* x) const { return const_cast<THnSparse*>(this)->GetBin(x, kFALSE); }
   Long64_t GetBin(const char* name[]) const { return const_cast<THnSparse*>(this)->GetBin(name, kFALSE); }
   Long64_t GetBin(const Int_t* idx, Bool_t allocate = kTRUE);
   Long64_t GetBin(const Double_t* x, Bool_t allocate = kTRUE);
   Long64_t GetBin(const char* name[], Bool_t allocate = kTRUE);

   /// Forwards to THnBase::SetBinContent().
   /// Non-virtual, CINT-compatible replacement of a using declaration.
   void SetBinContent(const Int_t* idx, Double_t v) {
      THnBase::SetBinContent(idx, v);
   }
   void SetBinContent(Long64_t bin, Double_t v);
   void SetBinError2(Long64_t bin, Double_t e2);

   /// Forwards to THnBase::SetBinContent().
   /// Non-virtual, CINT-compatible replacement of a using declaration.
   void AddBinContent(const Int_t* idx, Double_t v = 1.) {
      THnBase::AddBinContent(idx, v);
   }
   void AddBinContent(Long64_t bin, Double_t v = 1.);
   void AddBinError2(Long64_t bin, Double_t e2);

   /// Forwards to THnBase::GetBinContent() overload.
   /// Non-virtual, CINT-compatible replacement of a using declaration.
   Double_t GetBinContent(const Int_t *idx) const {

      return THnBase::GetBinContent(idx);
   }
   Double_t GetBinContent(Long64_t bin, Int_t* idx = 0) const;
   Double_t GetBinError2(Long64_t linidx) const;

   Double_t GetSparseFractionBins() const;
   Double_t GetSparseFractionMem() const;

   /// Forwards to THnBase::Projection().
   /// Non-virtual, as a CINT-compatible replacement of a using
   /// declaration.
   TH1D*      Projection(Int_t xDim, Option_t* option = "") const{
      return THnBase::Projection(xDim, option);
   }

   /// Forwards to THnBase::Projection().
   /// Non-virtual, as a CINT-compatible replacement of a using
   /// declaration.
   TH2D*      Projection(Int_t yDim, Int_t xDim,
                         Option_t* option = "") const {
      return THnBase::Projection(yDim, xDim, option);
   }

   /// Forwards to THnBase::Projection().
   /// Non-virtual, as a CINT-compatible replacement of a using
   /// declaration.
   TH3D*      Projection(Int_t xDim, Int_t yDim, Int_t zDim,
                         Option_t* option = "") const {
      return THnBase::Projection(xDim, yDim, zDim, option);
   }

   THnSparse* Projection(Int_t ndim, const Int_t* dim,
                         Option_t* option = "") const {
      return (THnSparse*) ProjectionND(ndim, dim, option);
   }

   THnSparse* Rebin(Int_t group) const {
      return (THnSparse*) RebinBase(group);
   }
   THnSparse* Rebin(const Int_t* group) const {
      return (THnSparse*) RebinBase(group);
   }

   void Reset(Option_t* option = "");
   void Sumw2();

   ClassDef(THnSparse, 3); // Interfaces of sparse n-dimensional histogram
};



//______________________________________________________________________________
/** \class THnSparseT
 Templated implementation of the abstract base THnSparse.
 All functionality and the interfaces to be used are in THnSparse!

 THnSparse does not know how to store any bin content itself. Instead, this
 is delegated to the derived, templated class: the template parameter decides
 what the format for the bin content is. In fact it even defines the array
 itself; possible implementations probably derive from TArray.

 Typedefs exist for template parameters with ROOT's generic types:

 Templated name        |    Typedef    |    Bin content type
 ----------------------|---------------|--------------------
 THnSparseT<TArrayC>   |  THnSparseC   |  Char_t
 THnSparseT<TArrayS>   |  THnSparseS   |  Short_t
 THnSparseT<TArrayI>   |  THnSparseI   |  Int_t
 THnSparseT<TArrayL>   |  THnSparseL   |  Long_t
 THnSparseT<TArrayL64> |  THnSparseL64 |  Long64_t
 THnSparseT<TArrayF>   |  THnSparseF   |  Float_t
 THnSparseT<TArrayD>   |  THnSparseD   |  Double_t

 We recommend to use THnSparseC wherever possible, and to map its value space
 of 256 possible values to e.g. float values outside the class. This saves an
 enormous amount of memory. Only if more than 256 values need to be
 distinguished should e.g. THnSparseS or even THnSparseF be chosen.

 Implementation detail: the derived, templated class is kept extremely small
 on purpose. That way the (templated thus inlined) uses of this class will
 only create a small amount of machine code, in contrast to e.g. STL.
*/


template <class CONT>
class THnSparseT: public THnSparse {
 public:
   THnSparseT() {}
   THnSparseT(const char* name, const char* title, Int_t dim,
              const Int_t* nbins, const Double_t* xmin = 0,
              const Double_t* xmax = 0, Int_t chunksize = 1024 * 16):
      THnSparse(name, title, dim, nbins, xmin, xmax, chunksize) {}

   TArray* GenerateArray() const { return new CONT(GetChunkSize()); }
 private:
   ClassDef(THnSparseT, 1); // Sparse n-dimensional histogram with templated content
};

typedef THnSparseT<TArrayD> THnSparseD;
typedef THnSparseT<TArrayF> THnSparseF;
typedef THnSparseT<TArrayL> THnSparseL;
typedef THnSparseT<TArrayL64> THnSparseL64;
typedef THnSparseT<TArrayI> THnSparseI;
typedef THnSparseT<TArrayS> THnSparseS;
typedef THnSparseT<TArrayC> THnSparseC;


#endif //  ROOT_THnSparse
