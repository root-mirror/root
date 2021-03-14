// @(#)root/hist:$Id$
// Author: Axel Naumann (2011-12-20)

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THnBase
#define ROOT_THnBase

/*************************************************************************

 THnBase: Common base class for n-dimensional histogramming.
 Defines interfaces and algorithms.

*************************************************************************/


#include "TNamed.h"
#include "TMath.h"
#include "TFitResultPtr.h"
#include "TObjArray.h"
#include "TArrayD.h"

class TAxis;
class TH1;
class TH1D;
class TH2D;
class TH3D;
class TF1;
class THnIter;

namespace ROOT {
namespace Internal {
   class THnBaseBinIter;
}
}

class THnBase: public TNamed {
protected:
   Int_t      fNdimensions;  // number of dimensions
   TObjArray  fAxes;         // axes of the histogram
   TObjArray  fBrowsables;   //! browser-helpers for each axis
   Double_t   fEntries;      // number of entries, spread over chunks
   Double_t   fTsumw;        // total sum of weights
   Double_t   fTsumw2;       // total sum of weights squared; -1 if no errors are calculated
   TArrayD    fTsumwx;       // total sum of weight*X for each dimension
   TArrayD    fTsumwx2;      // total sum of weight*X*X for each dimension
   std::vector<Double_t>  fIntegral;     //! array with bin weight sums
   enum {
      kNoInt,
      kValidInt,
      kInvalidInt
   } fIntegralStatus;        //! status of integral

 protected:
   THnBase():
      fNdimensions(0), fEntries(0),
      fTsumw(0), fTsumw2(-1.), fIntegral(), fIntegralStatus(kNoInt)
   {}

   THnBase(const char* name, const char* title, Int_t dim,
           const Int_t* nbins, const Double_t* xmin, const Double_t* xmax);

   void UpdateXStat(const Double_t *x, Double_t w = 1.) {
      if (GetCalculateErrors()) {
         for (Int_t d = 0; d < fNdimensions; ++d) {
            const Double_t xd = x[d];
            fTsumwx[d]  += w * xd;
            fTsumwx2[d] += w * xd * xd;
         }
      }
   }

   /// Increment the statistics due to filled weight "w",
   void FillBinBase(Double_t w) {
      fEntries += 1;
      if (GetCalculateErrors()) {
         fTsumw += w;
         fTsumw2 += w*w;
      }
      fIntegralStatus = kInvalidInt;
   }

   virtual void InitStorage(Int_t* nbins, Int_t chunkSize) = 0;
   void Init(const char* name, const char* title,
             const TObjArray* axes, Bool_t keepTargetAxis,
             Int_t chunkSize = 1024 * 16);
   THnBase* CloneEmpty(const char* name, const char* title,
                       const TObjArray* axes, Bool_t keepTargetAxis) const;
   virtual void Reserve(Long64_t /*nbins*/) {}
   virtual void SetFilledBins(Long64_t /*nbins*/) {};

   Bool_t CheckConsistency(const THnBase *h, const char *tag) const;
   TH1* CreateHist(const char* name, const char* title,
                   const TObjArray* axes, Bool_t keepTargetAxis) const;
   TObject* ProjectionAny(Int_t ndim, const Int_t* dim,
                          Bool_t wantNDim, Option_t* option = "") const;
   Bool_t PrintBin(Long64_t idx, Int_t* coord, Option_t* options) const;
   void AddInternal(const THnBase* h, Double_t c, Bool_t rebinned);
   THnBase* RebinBase(Int_t group) const;
   THnBase* RebinBase(const Int_t* group) const;
   void ResetBase(Option_t *option= "");

   static THnBase* CreateHnAny(const char* name, const char* title,
                               const TH1* h1, Bool_t sparse,
                               Int_t chunkSize = 1024 * 16);
   static THnBase* CreateHnAny(const char* name, const char* title,
                               const THnBase* hn, Bool_t sparse,
                               Int_t chunkSize = 1024 * 16);

 public:
   virtual ~THnBase();

   TObjArray* GetListOfAxes() { return &fAxes; }
   const TObjArray* GetListOfAxes() const { return &fAxes; }
   TAxis* GetAxis(Int_t dim) const { return (TAxis*)fAxes[dim]; }

   TFitResultPtr Fit(TF1 *f1 ,Option_t *option = "", Option_t *goption = "");
   TList* GetListOfFunctions() { return 0; }

   virtual ROOT::Internal::THnBaseBinIter* CreateIter(Bool_t respectAxisRange) const = 0;

   virtual Long64_t GetNbins() const = 0;
   Double_t GetEntries() const { return fEntries; }
   Double_t GetWeightSum() const { return fTsumw; }
   Int_t    GetNdimensions() const { return fNdimensions; }
   Bool_t   GetCalculateErrors() const { return fTsumw2 >= 0.; }

   /// Calculate errors (or not if "calc" == kFALSE)
   void     CalculateErrors(Bool_t calc = kTRUE) {
      if (calc) Sumw2();
      else fTsumw2 = -1.;
   }

   Long64_t Fill(const Double_t *x, Double_t w = 1.) {
      UpdateXStat(x, w);
      Long64_t bin = GetBin(x, kTRUE /*alloc*/);
      FillBin(bin, w);
      return bin;
   }
   Long64_t Fill(const char* name[], Double_t w = 1.) {
      Long64_t bin = GetBin(name, kTRUE /*alloc*/);
      FillBin(bin, w);
      return bin;
   }

   // alternate Fill semantics for RDF
   // firstval argument is needed to disambiguate template overloads
   template<typename... MoreTypes> 
   Long64_t Fill(Double_t firstval, MoreTypes ... morevals) {
      const std::array<double, 1 + sizeof...(morevals)> x =  {{ firstval, morevals... }};
      if (x.size() == GetNdimensions()) {
         // without weight
         return Fill(x.data()); 
      }
      else if (x.size() == (GetNdimensions() + 1)) {
         // with weight
         return Fill(x.data(), x.back());
      }
      
      return -1;
   }

   virtual void FillBin(Long64_t bin, Double_t w) = 0;

   void SetBinEdges(Int_t idim, const Double_t* bins);
   Bool_t IsInRange(Int_t *coord) const;
   Double_t GetBinError(const Int_t *idx) const { return GetBinError(GetBin(idx)); }
   Double_t GetBinError(Long64_t linidx) const { return TMath::Sqrt(GetBinError2(linidx)); }
   void SetBinError(const Int_t* idx, Double_t e) { SetBinError(GetBin(idx), e); }
   void SetBinError(Long64_t bin, Double_t e) { SetBinError2(bin, e*e); }
   void AddBinContent(const Int_t* x, Double_t v = 1.) { AddBinContent(GetBin(x), v); }
   void SetEntries(Double_t entries) { fEntries = entries; }
   void SetTitle(const char *title);

   Double_t GetBinContent(const Int_t *idx) const { return GetBinContent(GetBin(idx)); } // intentionally non-virtual
   virtual Double_t GetBinContent(Long64_t bin, Int_t* idx = 0) const = 0;
   virtual Double_t GetBinError2(Long64_t linidx) const = 0;
   virtual Long64_t GetBin(const Int_t* idx) const = 0;
   virtual Long64_t GetBin(const Double_t* x) const = 0;
   virtual Long64_t GetBin(const char* name[]) const = 0;
   virtual Long64_t GetBin(const Int_t* idx, Bool_t /*allocate*/ = kTRUE) = 0;
   virtual Long64_t GetBin(const Double_t* x, Bool_t /*allocate*/ = kTRUE) = 0;
   virtual Long64_t GetBin(const char* name[], Bool_t /*allocate*/ = kTRUE) = 0;

   void SetBinContent(const Int_t* idx, Double_t v) { SetBinContent(GetBin(idx), v); } // intentionally non-virtual
   virtual void SetBinContent(Long64_t bin, Double_t v) = 0;
   virtual void SetBinError2(Long64_t bin, Double_t e2) = 0;
   virtual void AddBinError2(Long64_t bin, Double_t e2) = 0;
   virtual void AddBinContent(Long64_t bin, Double_t v = 1.) = 0;

   Double_t GetSumw() const  { return fTsumw; }
   Double_t GetSumw2() const { return fTsumw2; }
   Double_t GetSumwx(Int_t dim) const  { return fTsumwx[dim]; }
   Double_t GetSumwx2(Int_t dim) const { return fTsumwx2[dim]; }

   /// Project all bins into a 1-dimensional histogram,
   /// keeping only axis "xDim".
   /// If "option" contains:
   ///  - "E" errors will be calculated.
   ///  - "A" ranges of the taget axes will be ignored.
   ///  - "O" original axis range of the taget axes will be
   ///    kept, but only bins inside the selected range
   ///    will be filled.
   TH1D*    Projection(Int_t xDim, Option_t* option = "") const {
      return (TH1D*) ProjectionAny(1, &xDim, false, option);
   }

   /// Project all bins into a 2-dimensional histogram,
   /// keeping only axes "xDim" and "yDim".
   ///
   /// WARNING: just like TH3::Project3D("yx") and TTree::Draw("y:x"),
   /// Projection(y,x) uses the first argument to define the y-axis and the
   /// second for the x-axis!
   ///
   /// If "option" contains "E" errors will be calculated.
   ///                      "A" ranges of the taget axes will be ignored.
   TH2D*    Projection(Int_t yDim, Int_t xDim, Option_t* option = "") const {
      const Int_t dim[2] = {xDim, yDim};
      return (TH2D*) ProjectionAny(2, dim, false, option);
   }

   /// Project all bins into a 3-dimensional histogram,
   /// keeping only axes "xDim", "yDim", and "zDim".
   /// If "option" contains:
   ///  - "E" errors will be calculated.
   ///  - "A" ranges of the taget axes will be ignored.
   ///  - "O" original axis range of the taget axes will be
   ///    kept, but only bins inside the selected range
   ///    will be filled.
   TH3D*    Projection(Int_t xDim, Int_t yDim, Int_t zDim, Option_t* option = "") const {
      const Int_t dim[3] = {xDim, yDim, zDim};
      return (TH3D*) ProjectionAny(3, dim, false, option);
   }

   THnBase* ProjectionND(Int_t ndim, const Int_t* dim,
                         Option_t* option = "") const {
      return (THnBase*)ProjectionAny(ndim, dim, kTRUE /*wantNDim*/, option);
   }

   Long64_t   Merge(TCollection* list);

   void Scale(Double_t c);
   void Add(const THnBase* h, Double_t c=1.);
   void Add(const TH1* hist, Double_t c=1.);
   void Multiply(const THnBase* h);
   void Multiply(TF1* f, Double_t c = 1.);
   void Divide(const THnBase* h);
   void Divide(const THnBase* h1, const THnBase* h2, Double_t c1 = 1., Double_t c2 = 1., Option_t* option="");
   void RebinnedAdd(const THnBase* h, Double_t c=1.);

   virtual void Reset(Option_t* option = "") = 0;
   virtual void Sumw2() = 0;

   Double_t ComputeIntegral();
   void GetRandom(Double_t *rand, Bool_t subBinRandom = kTRUE);

   void Print(Option_t* option = "") const;
   void PrintEntries(Long64_t from = 0, Long64_t howmany = -1, Option_t* options = 0) const;
   void PrintBin(Int_t* coord, Option_t* options) const {
      PrintBin(-1, coord, options);
   }
   void PrintBin(Long64_t idx, Option_t* options) const;

   void Browse(TBrowser *b);
   Bool_t IsFolder() const { return kTRUE; }

   //void Draw(Option_t* option = "");

   ClassDef(THnBase, 1); // Common base for n-dimensional histogram

   friend class THnIter;
};

namespace ROOT {
namespace Internal {
   // Helper class for browsing THnBase objects
   class THnBaseBrowsable: public TNamed {
   public:
      THnBaseBrowsable(THnBase* hist, Int_t axis);
      ~THnBaseBrowsable();
      void Browse(TBrowser *b);
      Bool_t IsFolder() const { return kFALSE; }

   private:
      THnBase* fHist; // Original histogram
      Int_t    fAxis; // Axis to visualize
      TH1*     fProj; // Projection result
      ClassDef(THnBaseBrowsable, 0); // Browser-helper for THnBase
   };

   // Base class for iterating over THnBase bins
   class THnBaseBinIter {
   public:
      THnBaseBinIter(Bool_t respectAxisRange):
         fRespectAxisRange(respectAxisRange), fHaveSkippedBin(kFALSE) {}
      virtual ~THnBaseBinIter();
      Bool_t HaveSkippedBin() const { return fHaveSkippedBin; }
      Bool_t RespectsAxisRange() const { return fRespectAxisRange; }

      virtual Int_t GetCoord(Int_t dim) const = 0;
      virtual Long64_t Next(Int_t* coord = 0) = 0;

   protected:
      Bool_t fRespectAxisRange;
      Bool_t fHaveSkippedBin;
   };
}
}

class THnIter: public TObject {
public:
   THnIter(const THnBase* hist, Bool_t respectAxisRange = kFALSE):
      fIter(hist->CreateIter(respectAxisRange)) {}
   virtual ~THnIter();

   /// Return the next bin's index.
   /// If provided, set coord to that bin's coordinates (bin indexes).
   /// I.e. coord must point to Int_t[hist->GetNdimensions()]
   /// Returns -1 when all bins have been visited.
   Long64_t Next(Int_t* coord = 0) {
      return fIter->Next(coord);
   }

   Int_t GetCoord(Int_t dim) const { return fIter->GetCoord(dim); }
   Bool_t HaveSkippedBin() const { return fIter->HaveSkippedBin(); }
   Bool_t RespectsAxisRange() const { return fIter->RespectsAxisRange(); }

private:
   ROOT::Internal::THnBaseBinIter* fIter;
   ClassDef(THnIter, 0); //Iterator over bins of a THnBase.
};

#endif //  ROOT_THnBase
