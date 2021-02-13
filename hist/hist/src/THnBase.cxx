// @(#)root/hist:$Id$
// Author: Axel Naumann (2011-12-20)

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THnBase.h"

#include "TAxis.h"
#include "TBrowser.h"
#include "TError.h"
#include "TClass.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVirtualPad.h"

#include "HFitInterface.h"
#include "Fit/DataRange.h"
#include "Fit/SparseData.h"
#include "Math/MinimizerOptions.h"
#include "Math/WrappedMultiTF1.h"


/** \class THnBase
    \ingroup Hist
Multidimensional histogram base.
Defines common functionality and interfaces for THn, THnSparse.
*/

ClassImp(THnBase);

////////////////////////////////////////////////////////////////////////////////
/// Construct a THnBase with "dim" dimensions,
/// "nbins" holds the number of bins for each dimension;
/// "xmin" and "xmax" the minimal and maximal value for each dimension.
/// The arrays "xmin" and "xmax" can be NULL; in that case SetBinEdges()
/// must be called for each dimension.

THnBase::THnBase(const char* name, const char* title, Int_t dim,
                 const Int_t* nbins, const Double_t* xmin, const Double_t* xmax):
TNamed(name, title), fNdimensions(dim), fAxes(dim), fBrowsables(dim),
fEntries(0), fTsumw(0), fTsumw2(-1.), fTsumwx(dim), fTsumwx2(dim),
fIntegral(0), fIntegralStatus(kNoInt)
{
   for (Int_t i = 0; i < fNdimensions; ++i) {
      TAxis* axis = new TAxis(nbins[i], xmin ? xmin[i] : 0., xmax ? xmax[i] : 1.);
      axis->SetName(TString::Format("axis%d", i));
      fAxes.AddAtAndExpand(axis, i);
   }
   SetTitle(title);
   fAxes.SetOwner();
}

////////////////////////////////////////////////////////////////////////////////
/// Destruct a THnBase

THnBase::~THnBase() {
   if (fIntegralStatus != kNoInt) delete [] fIntegral;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a new THnBase object that is of the same type as *this,
/// but with dimensions and bins given by axes.
/// If keepTargetAxis is true, the axes will keep their original xmin / xmax,
/// else they will be restricted to the range selected (first / last).

THnBase* THnBase::CloneEmpty(const char* name, const char* title,
                             const TObjArray* axes, Bool_t keepTargetAxis) const
{
   THnBase* ret = (THnBase*)IsA()->New();
   Int_t chunkSize = 1024 * 16;
   if (InheritsFrom(THnSparse::Class())) {
      chunkSize = ((const THnSparse*)this)->GetChunkSize();
   }
   ret->Init(name, title, axes, keepTargetAxis, chunkSize);
   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Initialize axes and name.

void THnBase::Init(const char* name, const char* title,
                   const TObjArray* axes, Bool_t keepTargetAxis,
                   Int_t chunkSize /*= 1024 * 16*/)
{
   SetNameTitle(name, title);

   TIter iAxis(axes);
   const TAxis* axis = 0;
   Int_t pos = 0;
   Int_t *nbins = new Int_t[axes->GetEntriesFast()];
   while ((axis = (TAxis*)iAxis())) {
      TAxis* reqaxis = new TAxis(*axis);
      if (!keepTargetAxis && axis->TestBit(TAxis::kAxisRange)) {
         Int_t binFirst = axis->GetFirst();
         // The lowest egde of the underflow is meaningless.
         if (binFirst == 0)
            binFirst = 1;
         Int_t binLast = axis->GetLast();
         // The overflow edge is implicit.
         if (binLast > axis->GetNbins())
            binLast = axis->GetNbins();
         Int_t nBins = binLast - binFirst + 1;
         if (axis->GetXbins()->GetSize()) {
            // non-uniform bins:
            reqaxis->Set(nBins, axis->GetXbins()->GetArray() + binFirst - 1);
         } else {
            // uniform bins:
            reqaxis->Set(nBins, axis->GetBinLowEdge(binFirst), axis->GetBinUpEdge(binLast));
         }
         reqaxis->ResetBit(TAxis::kAxisRange);
      }

      nbins[pos] = reqaxis->GetNbins();
      fAxes.AddAtAndExpand(new TAxis(*reqaxis), pos++);
   }
   fAxes.SetOwner();

   fNdimensions = axes->GetEntriesFast();
   InitStorage(nbins, chunkSize);
   delete [] nbins;
}


////////////////////////////////////////////////////////////////////////////////
/// Create an empty histogram with name and title with a given
/// set of axes. Create a TH1D/TH2D/TH3D, depending on the number
/// of elements in axes.

TH1* THnBase::CreateHist(const char* name, const char* title,
                         const TObjArray* axes,
                         Bool_t keepTargetAxis ) const {
   const int ndim = axes->GetSize();

   TH1* hist = 0;
   // create hist with dummy axes, we fix them later.
   if (ndim == 1)
      hist = new TH1D(name, title, 1, 0., 1.);
   else if (ndim == 2)
      hist = new TH2D(name, title, 1, 0., 1., 1, 0., 1.);
   else if (ndim == 3)
      hist = new TH3D(name, title, 1, 0., 1., 1, 0., 1., 1, 0., 1.);
   else {
      Error("CreateHist", "Cannot create histogram %s with %d dimensions!", name, ndim);
      return 0;
   }

   TAxis* hax[3] = {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()};
   for (Int_t d = 0; d < ndim; ++d) {
      TAxis* reqaxis = (TAxis*)(*axes)[d];
      hax[d]->SetTitle(reqaxis->GetTitle());
      if (!keepTargetAxis && reqaxis->TestBit(TAxis::kAxisRange)) {
         // axis cannot extend to underflow/overflows (fix ROOT-8781)
         Int_t binFirst = std::max(reqaxis->GetFirst(),1);
         Int_t binLast = std::min(reqaxis->GetLast(), reqaxis->GetNbins() );
         Int_t nBins = binLast - binFirst + 1;
         if (reqaxis->GetXbins()->GetSize()) {
            // non-uniform bins:
            hax[d]->Set(nBins, reqaxis->GetXbins()->GetArray() + binFirst - 1);
         } else {
            // uniform bins:
            hax[d]->Set(nBins, reqaxis->GetBinLowEdge(binFirst), reqaxis->GetBinUpEdge(binLast));
         }
      } else {
         if (reqaxis->GetXbins()->GetSize()) {
            // non-uniform bins:
            hax[d]->Set(reqaxis->GetNbins(), reqaxis->GetXbins()->GetArray());
         } else {
            // uniform bins:
            hax[d]->Set(reqaxis->GetNbins(), reqaxis->GetXmin(), reqaxis->GetXmax());
         }
      }
   }

   hist->Rebuild();

   return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a THn / THnSparse object from a histogram deriving from TH1.

THnBase* THnBase::CreateHnAny(const char* name, const char* title,
                              const TH1* h, Bool_t sparse, Int_t chunkSize)
{
   // Get the dimension of the TH1
   int ndim = h->GetDimension();

   // Axis properties
   int nbins[3] = {0,0,0};
   double minRange[3] = {0.,0.,0.};
   double maxRange[3] = {0.,0.,0.};
   const TAxis* axis[3] = { h->GetXaxis(), h->GetYaxis(), h->GetZaxis() };
   for (int i = 0; i < ndim; ++i) {
      nbins[i]    = axis[i]->GetNbins();
      minRange[i] = axis[i]->GetXmin();
      maxRange[i] = axis[i]->GetXmax();
   }

   // Create the corresponding THnSparse, depending on the storage
   // type of the TH1. The class name will be "TH??\0" where the first
   // ? is 1,2 or 3 and the second ? indicates the storage as C, S,
   // I, F or D.
   THnBase* s = 0;
   const char* cname( h->ClassName() );
   if (cname[0] == 'T' && cname[1] == 'H'
       && cname[2] >= '1' && cname[2] <= '3' && cname[4] == 0) {

#define R__THNBCASE(TAG)                                                \
if (sparse) {                                                     \
s = new _NAME2_(THnSparse,TAG)(name, title, ndim, nbins,       \
minRange, maxRange, chunkSize); \
} else {                                                          \
s = new _NAME2_(THn,TAG)(name, title, ndim, nbins,             \
minRange, maxRange);                  \
}                                                                 \
break;

      switch (cname[3]) {
         case 'F': R__THNBCASE(F);
         case 'D': R__THNBCASE(D);
         case 'I': R__THNBCASE(I);
         case 'S': R__THNBCASE(S);
         case 'C': R__THNBCASE(C);
      }
#undef R__THNBCASE
   }
   if (!s) {
      ::Warning("THnSparse::CreateHnAny", "Unknown Type of Histogram");
      return 0;
   }

   for (int i = 0; i < ndim; ++i) {
      s->GetAxis(i)->SetTitle(axis[i]->GetTitle());
   }

   // Get the array to know the number of entries of the TH1
   const TArray *array = dynamic_cast<const TArray*>(h);
   if (!array) {
      ::Warning("THnSparse::CreateHnAny", "Unknown Type of Histogram");
      return 0;
   }

   s->Add(h);
   return s;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a THnSparse (if "sparse") or THn  from "hn", possibly
/// converting THn <-> THnSparse.

THnBase* THnBase::CreateHnAny(const char* name, const char* title,
                              const THnBase* hn, Bool_t sparse,
                              Int_t chunkSize /*= 1024 * 16*/)
{
   TClass* type = 0;
   if (hn->InheritsFrom(THnSparse::Class())) {
      if (sparse) type = hn->IsA();
      else {
         char bintype;
         if (hn->InheritsFrom(THnSparseD::Class())) bintype = 'D';
         else if (hn->InheritsFrom(THnSparseF::Class())) bintype = 'F';
         else if (hn->InheritsFrom(THnSparseL::Class())) bintype = 'L';
         else if (hn->InheritsFrom(THnSparseI::Class())) bintype = 'I';
         else if (hn->InheritsFrom(THnSparseS::Class())) bintype = 'S';
         else if (hn->InheritsFrom(THnSparseC::Class())) bintype = 'C';
         else {
            hn->Error("CreateHnAny", "Type %s not implemented; please inform the ROOT team!",
                      hn->IsA()->GetName());
            return 0;
         }
         type = TClass::GetClass(TString::Format("THn%c", bintype));
      }
   } else if (hn->InheritsFrom(THn::Class())) {
      if (!sparse) type = hn->IsA();
      else {
         char bintype = 0;
         if (hn->InheritsFrom(THnD::Class())) bintype = 'D';
         else if (hn->InheritsFrom(THnF::Class())) bintype = 'F';
         else if (hn->InheritsFrom(THnC::Class())) bintype = 'C';
         else if (hn->InheritsFrom(THnS::Class())) bintype = 'S';
         else if (hn->InheritsFrom(THnI::Class())) bintype = 'I';
         else if (hn->InheritsFrom(THnL::Class())) bintype = 'L';
         else if (hn->InheritsFrom(THnL64::Class())) {
            hn->Error("CreateHnAny", "Type THnSparse with Long64_t bins is not available!");
            return 0;
         }
         if (bintype) {
            type = TClass::GetClass(TString::Format("THnSparse%c", bintype));
         }
      }
   } else {
      hn->Error("CreateHnAny", "Unhandled type %s, not deriving from THn nor THnSparse!",
                hn->IsA()->GetName());
      return 0;
   }
   if (!type) {
      hn->Error("CreateHnAny", "Unhandled type %s, please inform the ROOT team!",
                hn->IsA()->GetName());
      return 0;
   }

   THnBase* ret = (THnBase*)type->New();
   ret->Init(name, title, hn->GetListOfAxes(),
             kFALSE /*keepTargetAxes*/, chunkSize);

   ret->Add(hn);
   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the THnBase with the bins of hist that have content
/// or error != 0.

void THnBase::Add(const TH1* hist, Double_t c /*=1.*/)
{
   Long64_t nbins = hist->GetNcells();
   int x[3] = {0,0,0};
   for (int i = 0; i < nbins; ++i) {
      double value = hist->GetBinContent(i);
      double error = hist->GetBinError(i);
      if (!value && !error) continue;
      hist->GetBinXYZ(i, x[0], x[1], x[2]);
      SetBinContent(x, value * c);
      SetBinError(x, error * c);
   }
}

////////////////////////////////////////////////////////////////////////////////
///   Fit a THnSparse with function f
///
///   since the data is sparse by default a likelihood fit is performed
///   merging all the regions with empty bins for better performance efficiency
///
///  Since the THnSparse is not drawn no graphics options are passed
///  Here is the list of possible options
///
///                = "I"  Use integral of function in bin instead of value at bin center
///                = "X"  Use chi2 method (default is log-likelihood method)
///                = "U"  Use a User specified fitting algorithm (via SetFCN)
///                = "Q"  Quiet mode (minimum printing)
///                = "V"  Verbose mode (default is between Q and V)
///                = "E"  Perform better Errors estimation using Minos technique
///                = "B"  Use this option when you want to fix one or more parameters
///                       and the fitting function is like "gaus", "expo", "poln", "landau".
///                = "M"  More. Improve fit results
///                = "R"  Use the Range specified in the function range

TFitResultPtr THnBase::Fit(TF1 *f ,Option_t *option ,Option_t *goption)
{

   Foption_t fitOption;

   if (!TH1::FitOptionsMake(option,fitOption)) return 0;

   // The function used to fit cannot be stored in a THnSparse. It
   // cannot be drawn either. Perhaps in the future.
   fitOption.Nostore = true;
   // Use likelihood fit if not specified
   if (!fitOption.Chi2) fitOption.Like = true;
   // create range and minimizer options with default values
   ROOT::Fit::DataRange range(GetNdimensions());
   for ( int i = 0; i < GetNdimensions(); ++i ) {
      TAxis *axis = GetAxis(i);
      range.AddRange(i, axis->GetXmin(), axis->GetXmax());
   }
   ROOT::Math::MinimizerOptions minOption;

   return ROOT::Fit::FitObject(this, f , fitOption , minOption, goption, range);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an n-dimensional random tuple based on the histogrammed
/// distribution. If subBinRandom, the returned tuple will be additionally
/// randomly distributed within the randomized bin, using a flat
/// distribution.

void THnBase::GetRandom(Double_t *rand, Bool_t subBinRandom /* = kTRUE */)
{
   // check whether the integral array is valid
   if (fIntegralStatus != kValidInt)
      ComputeIntegral();

   // generate a random bin
   Double_t p = gRandom->Rndm();
   Long64_t idx = TMath::BinarySearch(GetNbins() + 1, fIntegral, p);
   const Int_t nStaticBins = 40;
   Int_t bin[nStaticBins];
   Int_t* pBin = bin;
   if (GetNdimensions() > nStaticBins) {
      pBin = new Int_t[GetNdimensions()];
   }
   GetBinContent(idx, pBin);

   // convert bin coordinates to real values
   for (Int_t i = 0; i < fNdimensions; i++) {
      rand[i] = GetAxis(i)->GetBinCenter(pBin[i]);

      // randomize the vector within a bin
      if (subBinRandom)
         rand[i] += (gRandom->Rndm() - 0.5) * GetAxis(i)->GetBinWidth(pBin[i]);
   }
   if (pBin != bin) {
      delete [] pBin;
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether bin coord is in range, as defined by TAxis::SetRange().

Bool_t THnBase::IsInRange(Int_t *coord) const
{
   Int_t min = 0;
   Int_t max = 0;
   for (Int_t i = 0; i < fNdimensions; ++i) {
      TAxis *axis = GetAxis(i);
      if (!axis->TestBit(TAxis::kAxisRange)) continue;
      min = axis->GetFirst();
      max = axis->GetLast();
      if (coord[i] < min || coord[i] > max)
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Project all bins into a ndim-dimensional THn / THnSparse (whatever
/// *this is) or if (ndim < 4 and !wantNDim) a TH1/2/3 histogram,
/// keeping only axes in dim (specifying ndim dimensions).
/// If "option" contains:
///  - "E" errors will be calculated.
///  - "A" ranges of the target axes will be ignored.
///  - "O" original axis range of the target axes will be
///    kept, but only bins inside the selected range
///    will be filled.

TObject* THnBase::ProjectionAny(Int_t ndim, const Int_t* dim,
                                Bool_t wantNDim,
                                Option_t* option /*= ""*/) const
{
   TString name(GetName());
   name +="_proj";

   for (Int_t d = 0; d < ndim; ++d) {
      name += "_";
      name += dim[d];
   }

   TString title(GetTitle());
   Ssiz_t posInsert = title.First(';');
   if (posInsert == kNPOS) {
      title += " projection ";
      for (Int_t d = 0; d < ndim; ++d)
         title += GetAxis(dim[d])->GetTitle();
   } else {
      for (Int_t d = ndim - 1; d >= 0; --d) {
         title.Insert(posInsert, GetAxis(d)->GetTitle());
         if (d)
            title.Insert(posInsert, ", ");
      }
      title.Insert(posInsert, " projection ");
   }

   TObjArray newaxes(ndim);
   for (Int_t d = 0; d < ndim; ++d) {
      newaxes.AddAt(GetAxis(dim[d]),d);
   }

   THnBase* hn = 0;
   TH1* hist = 0;
   TObject* ret = 0;

   Bool_t* hadRange = 0;
   Bool_t ignoreTargetRange = (option && (strchr(option, 'A') || strchr(option, 'a')));
   Bool_t keepTargetAxis = ignoreTargetRange || (option && (strchr(option, 'O') || strchr(option, 'o')));
   if (ignoreTargetRange) {
      hadRange = new Bool_t[ndim];
      for (Int_t d = 0; d < ndim; ++d){
         TAxis *axis = GetAxis(dim[d]);
         hadRange[d] = axis->TestBit(TAxis::kAxisRange);
         axis->SetBit(TAxis::kAxisRange, kFALSE);
      }
   }

   if (wantNDim)
      ret = hn = CloneEmpty(name, title, &newaxes, keepTargetAxis);
   else
      ret = hist = CreateHist(name, title, &newaxes, keepTargetAxis);

   if (keepTargetAxis) {
      // make the whole axes visible, i.e. unset the range
      if (wantNDim) {
         for (Int_t d = 0; d < ndim; ++d) {
            hn->GetAxis(d)->SetRange(0, 0);
         }
      } else {
         hist->GetXaxis()->SetRange(0, 0);
         hist->GetYaxis()->SetRange(0, 0);
         hist->GetZaxis()->SetRange(0, 0);
      }
   }

   Bool_t haveErrors = GetCalculateErrors();
   Bool_t wantErrors = haveErrors || (option && (strchr(option, 'E') || strchr(option, 'e')));

   Int_t* bins  = new Int_t[ndim];
   Long64_t myLinBin = 0;

   THnIter iter(this, kTRUE /*use axis range*/);

   while ((myLinBin = iter.Next()) >= 0) {
      Double_t v = GetBinContent(myLinBin);

      for (Int_t d = 0; d < ndim; ++d) {
         bins[d] = iter.GetCoord(dim[d]);
         if (!keepTargetAxis && GetAxis(dim[d])->TestBit(TAxis::kAxisRange)) {
            Int_t binOffset = GetAxis(dim[d])->GetFirst();
            // Don't subtract even more if underflow is alreday included:
            if (binOffset > 0) --binOffset;
            bins[d] -= binOffset;
         }
      }

      Long64_t targetLinBin = -1;
      if (!wantNDim) {
         if (ndim == 1) targetLinBin = bins[0];
         else if (ndim == 2) targetLinBin = hist->GetBin(bins[0], bins[1]);
         else if (ndim == 3) targetLinBin = hist->GetBin(bins[0], bins[1], bins[2]);
      } else {
         targetLinBin = hn->GetBin(bins, kTRUE /*allocate*/);
      }

      if (wantErrors) {
         Double_t err2 = 0.;
         if (haveErrors) {
            err2 = GetBinError2(myLinBin);
         } else {
            err2 = v;
         }
         if (wantNDim) {
            hn->AddBinError2(targetLinBin, err2);
         } else {
            Double_t preverr = hist->GetBinError(targetLinBin);
            hist->SetBinError(targetLinBin, TMath::Sqrt(preverr * preverr + err2));
         }
      }

      // only _after_ error calculation, or sqrt(v) is taken into account!
      if (wantNDim)
         hn->AddBinContent(targetLinBin, v);
      else
         hist->AddBinContent(targetLinBin, v);
   }

   delete [] bins;

   if (wantNDim) {
      hn->SetEntries(fEntries);
   } else {
      if (!iter.HaveSkippedBin()) {
         hist->SetEntries(fEntries);
      } else {
         // re-compute the entries
         // in case of error calculation (i.e. when Sumw2() is set)
         // use the effective entries for the entries
         // since this  is the only way to estimate them
         hist->ResetStats();
         Double_t entries = hist->GetEffectiveEntries();
         if (!wantErrors) {
            // to avoid numerical rounding
            entries = TMath::Floor(entries + 0.5);
         }
         hist->SetEntries(entries);
      }
   }

   if (hadRange) {
      // reset kAxisRange bit:
      for (Int_t d = 0; d < ndim; ++d)
         GetAxis(dim[d])->SetBit(TAxis::kAxisRange, hadRange[d]);

      delete [] hadRange;
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Scale contents and errors of this histogram by c:
/// this = this * c
/// It does not modify the histogram's number of entries.

void THnBase::Scale(Double_t c)
{

   Double_t nEntries = GetEntries();
   // Scale the contents & errors
   Bool_t haveErrors = GetCalculateErrors();
   Long64_t i = 0;
   THnIter iter(this);
   while ((i = iter.Next()) >= 0) {
      // Get the content of the bin from the current histogram
      Double_t v = GetBinContent(i);
      SetBinContent(i, c * v);
      if (haveErrors) {
         Double_t err2 = GetBinError2(i);
         SetBinError2(i, c * c * err2);
      }
   }
   SetEntries(nEntries);
}

////////////////////////////////////////////////////////////////////////////////
/// Add() implementation for both rebinned histograms and those with identical
/// binning. See THnBase::Add().

void THnBase::AddInternal(const THnBase* h, Double_t c, Bool_t rebinned)
{
   if (fNdimensions != h->GetNdimensions()) {
      Warning("RebinnedAdd", "Different number of dimensions, cannot carry out operation on the histograms");
      return;
   }

   // Trigger error calculation if h has it
   if (!GetCalculateErrors() && h->GetCalculateErrors())
      Sumw2();
   Bool_t haveErrors = GetCalculateErrors();

   Double_t* x = 0;
   if (rebinned) {
      x = new Double_t[fNdimensions];
   }
   Int_t* coord = new Int_t[fNdimensions];

   // Expand the exmap if needed, to reduce collisions
   Long64_t numTargetBins = GetNbins() + h->GetNbins();
   Reserve(numTargetBins);

   Long64_t i = 0;
   THnIter iter(h);
   // Add to this whatever is found inside the other histogram
   while ((i = iter.Next(coord)) >= 0) {
      // Get the content of the bin from the second histogram
      Double_t v = h->GetBinContent(i);

      Long64_t mybinidx = -1;
      if (rebinned) {
         // Get the bin center given a coord
         for (Int_t j = 0; j < fNdimensions; ++j)
            x[j] = h->GetAxis(j)->GetBinCenter(coord[j]);

         mybinidx = GetBin(x, kTRUE /* allocate*/);
      } else {
         mybinidx = GetBin(coord, kTRUE /*allocate*/);
      }

      if (haveErrors) {
         Double_t err2 = h->GetBinError2(i) * c * c;
         AddBinError2(mybinidx, err2);
      }
      // only _after_ error calculation, or sqrt(v) is taken into account!
      AddBinContent(mybinidx, c * v);
   }

   delete [] coord;
   delete [] x;

   Double_t nEntries = GetEntries() + c * h->GetEntries();
   SetEntries(nEntries);
}

////////////////////////////////////////////////////////////////////////////////
/// Add contents of h scaled by c to this histogram:
/// this = this + c * h
/// Note that if h has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.

void THnBase::Add(const THnBase* h, Double_t c)
{
   // Check consistency of the input
   if (!CheckConsistency(h, "Add")) return;

   AddInternal(h, c, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Add contents of h scaled by c to this histogram:
/// this = this + c * h
/// Note that if h has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
/// In contrast to Add(), RebinnedAdd() does not require consistent binning of
/// this and h; instead, each bin's center is used to determine the target bin.

void THnBase::RebinnedAdd(const THnBase* h, Double_t c)
{
   AddInternal(h, c, kTRUE);
}


////////////////////////////////////////////////////////////////////////////////
/// Merge this with a list of THnBase's. All THnBase's provided
/// in the list must have the same bin layout!

Long64_t THnBase::Merge(TCollection* list)
{
   if (!list) return 0;
   if (list->IsEmpty()) return (Long64_t)GetEntries();

   Long64_t sumNbins = GetNbins();
   TIter iter(list);
   const TObject* addMeObj = 0;
   while ((addMeObj = iter())) {
      const THnBase* addMe = dynamic_cast<const THnBase*>(addMeObj);
      if (addMe) {
         sumNbins += addMe->GetNbins();
      }
   }
   Reserve(sumNbins);

   iter.Reset();
   while ((addMeObj = iter())) {
      const THnBase* addMe = dynamic_cast<const THnBase*>(addMeObj);
      if (!addMe)
         Error("Merge", "Object named %s is not THnBase! Skipping it.",
               addMeObj->GetName());
      else
         Add(addMe);
   }
   return (Long64_t)GetEntries();
}


////////////////////////////////////////////////////////////////////////////////
/// Multiply this histogram by histogram h
/// this = this * h
/// Note that if h has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.

void THnBase::Multiply(const THnBase* h)
{
   // Check consistency of the input
   if(!CheckConsistency(h, "Multiply"))return;

   // Trigger error calculation if h has it
   Bool_t wantErrors = kFALSE;
   if (GetCalculateErrors() || h->GetCalculateErrors())
      wantErrors = kTRUE;

   if (wantErrors) Sumw2();

   Double_t nEntries = GetEntries();
   // Now multiply the contents: in this case we have the intersection of the sets of bins
   Int_t* coord = new Int_t[fNdimensions];
   Long64_t i = 0;
   THnIter iter(this);
   // Add to this whatever is found inside the other histogram
   while ((i = iter.Next(coord)) >= 0) {
      // Get the content of the bin from the current histogram
      Double_t v1 = GetBinContent(i);
      // Now look at the bin with the same coordinates in h
      Long64_t idxh = h->GetBin(coord);
      Double_t v2 = 0.;
      if (idxh >= 0) v2 = h->GetBinContent(idxh);
      SetBinContent(i, v1 * v2);
      if (wantErrors) {
         Double_t err1 = GetBinError(i) * v2;
         Double_t err2 = 0.;
         if (idxh >= 0) err2 = h->GetBinError(idxh) * v1;
         SetBinError(i, TMath::Sqrt((err2 * err2 + err1 * err1)));
      }
   }
   SetEntries(nEntries);

   delete [] coord;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this*c*f1
/// if errors are defined, errors are also recalculated.
///
/// Only bins inside the function range are recomputed.
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after
/// calling Multiply()

void THnBase::Multiply(TF1* f, Double_t c)
{
   Int_t* coord = new Int_t[fNdimensions];
   Double_t* x = new Double_t[fNdimensions];

   Bool_t wantErrors = GetCalculateErrors();
   if (wantErrors) Sumw2();

   Long64_t i = 0;
   THnIter iter(this);
   // Add to this whatever is found inside the other histogram
   while ((i = iter.Next(coord)) >= 0) {
      Double_t value = GetBinContent(i);

      // Get the bin coordinates given an index array
      for (Int_t j = 0; j < fNdimensions; ++j)
         x[j] = GetAxis(j)->GetBinCenter(coord[j]);

      if (!f->IsInside(x))
         continue;
      TF1::RejectPoint(kFALSE);

      // Evaluate function at points
      Double_t fvalue = f->EvalPar(x, NULL);

      SetBinContent(i, c * fvalue * value);
      if (wantErrors) {
         Double_t error = GetBinError(i);
         SetBinError(i, c * fvalue * error);
      }
   }

   delete [] x;
   delete [] coord;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this histogram by h
/// this = this/(h)
/// Note that if h has Sumw2 set, Sumw2 is automatically called for
/// this if not already set.
/// The resulting errors are calculated assuming uncorrelated content.

void THnBase::Divide(const THnBase *h)
{
   // Check consistency of the input
   if (!CheckConsistency(h, "Divide"))return;

   // Trigger error calculation if h has it
   Bool_t wantErrors=GetCalculateErrors();
   if (!GetCalculateErrors() && h->GetCalculateErrors())
      wantErrors=kTRUE;

   // Remember original histogram statistics
   Double_t nEntries = fEntries;

   if (wantErrors) Sumw2();
   Bool_t didWarn = kFALSE;

   // Now divide the contents: also in this case we have the intersection of the sets of bins
   Int_t* coord = new Int_t[fNdimensions];
   Long64_t i = 0;
   THnIter iter(this);
   // Add to this whatever is found inside the other histogram
   while ((i = iter.Next(coord)) >= 0) {
      // Get the content of the bin from the first histogram
      Double_t v1 = GetBinContent(i);
      // Now look at the bin with the same coordinates in h
      Long64_t hbin = h->GetBin(coord);
      Double_t v2 = h->GetBinContent(hbin);
      if (!v2) {
         v1 = 0.;
         v2 = 1.;
         if (!didWarn) {
            Warning("Divide(h)", "Histogram h has empty bins - division by zero! Setting bin to 0.");
            didWarn = kTRUE;
         }
      }
      SetBinContent(i, v1 / v2);
      if (wantErrors) {
         Double_t err1 = GetBinError(i) * v2;
         Double_t err2 = h->GetBinError(hbin) * v1;
         Double_t b22 = v2 * v2;
         Double_t err = (err1 * err1 + err2 * err2) / (b22 * b22);
         SetBinError2(i, err);
      }
   }
   delete [] coord;
   SetEntries(nEntries);
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this histogram by multiplication of h1 by h2
/// this = (c1*h1)/(c2*h2)
/// Note that if h1 or h2 have Sumw2 set, Sumw2 is automatically called for
/// this if not already set.
/// The resulting errors are calculated assuming uncorrelated content.
/// However, if option ="B" is specified, Binomial errors are computed.
/// In this case c1 and c2 do not make real sense and they are ignored.

void THnBase::Divide(const THnBase *h1, const THnBase *h2, Double_t c1, Double_t c2, Option_t *option)
{

   TString opt = option;
   opt.ToLower();
   Bool_t binomial = kFALSE;
   if (opt.Contains("b")) binomial = kTRUE;

   // Check consistency of the input
   if (!CheckConsistency(h1, "Divide") || !CheckConsistency(h2, "Divide"))return;
   if (!c2) {
      Error("Divide","Coefficient of dividing histogram cannot be zero");
      return;
   }

   Reset();

   // Trigger error calculation if h1 or h2 have it
   if (!GetCalculateErrors() && (h1->GetCalculateErrors()|| h2->GetCalculateErrors() != 0))
      Sumw2();

   // Count filled bins
   Long64_t nFilledBins=0;

   // Now divide the contents: we have the intersection of the sets of bins

   Int_t* coord = new Int_t[fNdimensions];
   memset(coord, 0, sizeof(Int_t) * fNdimensions);
   Bool_t didWarn = kFALSE;

   Long64_t i = 0;
   THnIter iter(h1);
   // Add to this whatever is found inside the other histogram
   while ((i = iter.Next(coord)) >= 0) {
      // Get the content of the bin from the first histogram
      Double_t v1 = h1->GetBinContent(i);
      // Now look at the bin with the same coordinates in h2
      Long64_t h2bin = h2->GetBin(coord);
      Double_t v2 = h2->GetBinContent(h2bin);
      if (!v2) {
         v1 = 0.;
         v2 = 1.;
         if (!didWarn) {
            Warning("Divide(h1, h2)", "Histogram h2 has empty bins - division by zero! Setting bin to 0.");
            didWarn = kTRUE;
         }
      }
      nFilledBins++;
      Long64_t myBin = GetBin(coord);
      SetBinContent(myBin, c1 * v1 / c2 / v2);
      if(GetCalculateErrors()){
         Double_t err1 = h1->GetBinError(i);
         Double_t err2 = h2->GetBinError(h2bin);
         Double_t errSq = 0.;
         if (binomial) {
            if (v1 != v2) {
               Double_t w = v1 / v2;
               err2 *= w;
               errSq = TMath::Abs( ( (1. - 2.*w) * err1 * err1 + err2 * err2 ) / (v2 * v2) );
            }
         } else {
            c1 *= c1;
            c2 *= c2;
            Double_t b22 = v2 * v2 * c2;
            err1 *= v2;
            err2 *= v1;
            errSq = c1 * c2 * (err1 * err1 + err2 * err2) / (b22 * b22);
         }
         SetBinError2(myBin, errSq);
      }
   }

   delete [] coord;
   SetFilledBins(nFilledBins);

   // Set as entries in the result histogram the entries in the numerator
   SetEntries(h1->GetEntries());
}

////////////////////////////////////////////////////////////////////////////////
/// Consistency check on (some of) the parameters of two histograms (for operations).

Bool_t THnBase::CheckConsistency(const THnBase *h, const char *tag) const
{
   if (fNdimensions != h->GetNdimensions()) {
      Warning(tag, "Different number of dimensions, cannot carry out operation on the histograms");
      return kFALSE;
   }
   for (Int_t dim = 0; dim < fNdimensions; dim++){
      if (GetAxis(dim)->GetNbins() != h->GetAxis(dim)->GetNbins()) {
         Warning(tag, "Different number of bins on axis %i, cannot carry out operation on the histograms", dim);
         return kFALSE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the axis # of bins and bin limits on dimension idim

void THnBase::SetBinEdges(Int_t idim, const Double_t* bins)
{
   TAxis* axis = (TAxis*) fAxes[idim];
   axis->Set(axis->GetNbins(), bins);
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the title.
///
/// If title is in the form "stringt;string0;string1;string2 ..."
/// the histogram title is set to stringt, the title of axis0 to string0,
/// of axis1 to string1, of axis2 to string2, etc, just like it is done
/// for TH1/TH2/TH3.
/// To insert the character ";" in one of the titles, one should use "#;"
/// or "#semicolon".

void THnBase::SetTitle(const char *title)
{
   fTitle = title;
   fTitle.ReplaceAll("#;",2,"#semicolon",10);

   Int_t endHistTitle = fTitle.First(';');
   if (endHistTitle >= 0) {
      // title contains a ';' so parse the axis titles
      Int_t posTitle = endHistTitle + 1;
      Int_t lenTitle = fTitle.Length();
      Int_t dim = 0;
      while (posTitle > 0 && posTitle < lenTitle && dim < fNdimensions){
         Int_t endTitle = fTitle.Index(";", posTitle);
         TString axisTitle = fTitle(posTitle, endTitle - posTitle);
         axisTitle.ReplaceAll("#semicolon", 10, ";", 1);
         GetAxis(dim)->SetTitle(axisTitle);
         dim++;
         if (endTitle > 0)
            posTitle = endTitle + 1;
         else
            posTitle = -1;
      }
      // Remove axis titles from histogram title
      fTitle.Remove(endHistTitle, lenTitle - endHistTitle);
   }

   fTitle.ReplaceAll("#semicolon", 10, ";", 1);

}

////////////////////////////////////////////////////////////////////////////////
/// Combine the content of "group" neighboring bins into
/// a new bin and return the resulting THnBase.
/// For group=2 and a 3 dimensional histogram, all "blocks"
/// of 2*2*2 bins will be put into a bin.

THnBase* THnBase::RebinBase(Int_t group) const
{
   Int_t* ngroup = new Int_t[GetNdimensions()];
   for (Int_t d = 0; d < GetNdimensions(); ++d)
      ngroup[d] = group;
   THnBase* ret = RebinBase(ngroup);
   delete [] ngroup;
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Combine the content of "group" neighboring bins for each dimension
/// into a new bin and return the resulting THnBase.
/// For group={2,1,1} and a 3 dimensional histogram, pairs of x-bins
/// will be grouped.

THnBase* THnBase::RebinBase(const Int_t* group) const
{
   Int_t ndim = GetNdimensions();
   TString name(GetName());
   for (Int_t d = 0; d < ndim; ++d)
      name += Form("_%d", group[d]);


   TString title(GetTitle());
   Ssiz_t posInsert = title.First(';');
   if (posInsert == kNPOS) {
      title += " rebin ";
      for (Int_t d = 0; d < ndim; ++d)
         title += Form("{%d}", group[d]);
   } else {
      for (Int_t d = ndim - 1; d >= 0; --d)
         title.Insert(posInsert, Form("{%d}", group[d]));
      title.Insert(posInsert, " rebin ");
   }

   TObjArray newaxes(ndim);
   newaxes.SetOwner();
   for (Int_t d = 0; d < ndim; ++d) {
      newaxes.AddAt(new TAxis(*GetAxis(d) ),d);
      if (group[d] > 1) {
         TAxis* newaxis = (TAxis*) newaxes.At(d);
         Int_t newbins = (newaxis->GetNbins() + group[d] - 1) / group[d];
         if (newaxis->GetXbins() && newaxis->GetXbins()->GetSize()) {
            // variable bins
            Double_t *edges = new Double_t[newbins + 1];
            for (Int_t i = 0; i < newbins + 1; ++i)
               if (group[d] * i <= newaxis->GetNbins())
                  edges[i] = newaxis->GetXbins()->At(group[d] * i);
               else edges[i] = newaxis->GetXmax();
            newaxis->Set(newbins, edges);
            delete [] edges;
         } else {
            newaxis->Set(newbins, newaxis->GetXmin(), newaxis->GetXmax());
         }
      }
   }

   THnBase* h = CloneEmpty(name.Data(), title.Data(), &newaxes, kTRUE);
   Bool_t haveErrors = GetCalculateErrors();
   Bool_t wantErrors = haveErrors;

   Int_t* bins  = new Int_t[ndim];
   Int_t* coord = new Int_t[fNdimensions];

   Long64_t i = 0;
   THnIter iter(this);
   while ((i = iter.Next(coord)) >= 0) {
      Double_t v = GetBinContent(i);
      for (Int_t d = 0; d < ndim; ++d) {
         bins[d] = TMath::CeilNint( (double) coord[d]/group[d] );
      }
      Long64_t idxh = h->GetBin(bins, kTRUE /*allocate*/);

      if (wantErrors) {
         // wantErrors == haveErrors, thus:
         h->AddBinError2(idxh, GetBinError2(i));
      }

      // only _after_ error calculation, or sqrt(v) is taken into account!
      h->AddBinContent(idxh, v);
   }

   delete [] bins;
   delete [] coord;

   h->SetEntries(fEntries);

   return h;

}

////////////////////////////////////////////////////////////////////////////////
/// Clear the histogram

void THnBase::ResetBase(Option_t * /*option = ""*/)
{
   fEntries = 0.;
   fTsumw = 0.;
   fTsumw2 = -1.;
   if (fIntegralStatus != kNoInt) {
      delete [] fIntegral;
      fIntegralStatus = kNoInt;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the integral of the histogram

Double_t THnBase::ComputeIntegral()
{
   // delete old integral
   if (fIntegralStatus != kNoInt) {
      delete [] fIntegral;
      fIntegralStatus = kNoInt;
   }

   // check number of bins
   if (GetNbins() == 0) {
      Error("ComputeIntegral", "The histogram must have at least one bin.");
      return 0.;
   }

   // allocate integral array
   fIntegral = new Double_t [GetNbins() + 1];
   fIntegral[0] = 0.;

   // fill integral array with contents of regular bins (non over/underflow)
   Int_t* coord = new Int_t[fNdimensions];
   Long64_t i = 0;
   THnIter iter(this);
   while ((i = iter.Next(coord)) >= 0) {
      Double_t v = GetBinContent(i);

      // check whether the bin is regular
      bool regularBin = true;
      for (Int_t dim = 0; dim < fNdimensions; dim++) {
         if (coord[dim] < 1 || coord[dim] > GetAxis(dim)->GetNbins()) {
            regularBin = false;
            break;
         }
      }

      // if outlayer, count it with zero weight
      if (!regularBin) v = 0.;

      fIntegral[i + 1] = fIntegral[i] + v;
   }
   delete [] coord;

   // check sum of weights
   if (fIntegral[GetNbins()] == 0.) {
      Error("ComputeIntegral", "No hits in regular bins (non over/underflow).");
      delete [] fIntegral;
      return 0.;
   }

   // normalize the integral array
   for (Long64_t j = 0; j <= GetNbins(); ++j)
      fIntegral[j] = fIntegral[j] / fIntegral[GetNbins()];

   // set status to valid
   fIntegralStatus = kValidInt;
   return fIntegral[GetNbins()];
}

////////////////////////////////////////////////////////////////////////////////
/// Print bin with linex index "idx".
/// For valid options see PrintBin(Long64_t idx, Int_t* bin, Option_t* options).

void THnBase::PrintBin(Long64_t idx, Option_t* options) const
{
   Int_t* coord = new Int_t[fNdimensions];
   PrintBin(idx, coord, options);
   delete [] coord;
}

////////////////////////////////////////////////////////////////////////////////
/// Print one bin. If "idx" is != -1 use that to determine the bin,
/// otherwise (if "idx" == -1) use the coordinate in "bin".
/// If "options" contains:
///  - '0': only print bins with an error or content != 0
/// Return whether the bin was printed (depends on options)

Bool_t THnBase::PrintBin(Long64_t idx, Int_t* bin, Option_t* options) const
{
   Double_t v = -42;
   if (idx == -1) {
      idx = GetBin(bin);
      v = GetBinContent(idx);
   } else {
      v = GetBinContent(idx, bin);
   }

   Double_t err = 0.;
   if (GetCalculateErrors()) {
      if (idx != -1) {
         err = GetBinError(idx);
      }
   }

   if (v == 0. && err == 0. && options && strchr(options, '0')) {
      // suppress zeros, and we have one.
      return kFALSE;
   }

   TString coord;
   for (Int_t dim = 0; dim < fNdimensions; ++dim) {
      coord += bin[dim];
      coord += ',';
   }
   coord.Remove(coord.Length() - 1);

   if (GetCalculateErrors()) {
      Printf("Bin at (%s) = %g (+/- %g)", coord.Data(), v, err);
   } else {
      Printf("Bin at (%s) = %g", coord.Data(), v);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print "howmany" entries starting at "from". If "howmany" is -1, print all.
/// If "options" contains:
///  - 'x': print in the order of axis bins, i.e. (0,0,...,0), (0,0,...,1),...
///  - '0': only print bins with content != 0

void THnBase::PrintEntries(Long64_t from /*=0*/, Long64_t howmany /*=-1*/,
                           Option_t* options /*=0*/) const
{
   if (from < 0) from = 0;
   if (howmany == -1) howmany = GetNbins();

   Int_t* bin = new Int_t[fNdimensions];

   if (options && (strchr(options, 'x') || strchr(options, 'X'))) {
      Int_t* nbins = new Int_t[fNdimensions];
      for (Int_t dim = fNdimensions - 1; dim >= 0; --dim) {
         nbins[dim] = GetAxis(dim)->GetNbins();
         bin[dim] = from % nbins[dim];
         from /= nbins[dim];
      }

      for (Long64_t i = 0; i < howmany; ++i) {
         if (!PrintBin(-1, bin, options))
            ++howmany;
         // Advance to next bin:
         ++bin[fNdimensions - 1];
         for (Int_t dim = fNdimensions - 1; dim >= 0; --dim) {
            if (bin[dim] >= nbins[dim]) {
               bin[dim] = 0;
               if (dim > 0) {
                  ++bin[dim - 1];
               } else {
                  howmany = -1; // aka "global break"
               }
            }
         }
      }
      delete [] nbins;
   } else {
      for (Long64_t i = from; i < from + howmany; ++i) {
         if (!PrintBin(i, bin, options))
            ++howmany;
      }
   }
   delete [] bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a THnBase. If "option" contains:
///  - 'a': print axis details
///  - 'm': print memory usage
///  - 's': print statistics
///  - 'c': print its content, too (this can generate a LOT of output!)
/// Other options are forwarded to PrintEntries().

void THnBase::Print(Option_t* options) const
{
   Bool_t optAxis    = options && (strchr(options, 'A') || (strchr(options, 'a')));
   Bool_t optMem     = options && (strchr(options, 'M') || (strchr(options, 'm')));
   Bool_t optStat    = options && (strchr(options, 'S') || (strchr(options, 's')));
   Bool_t optContent = options && (strchr(options, 'C') || (strchr(options, 'c')));

   Printf("%s (*0x%lx): \"%s\" \"%s\"", IsA()->GetName(), (unsigned long)this, GetName(), GetTitle());
   Printf("  %d dimensions, %g entries in %lld filled bins", GetNdimensions(), GetEntries(), GetNbins());

   if (optAxis) {
      for (Int_t dim = 0; dim < fNdimensions; ++dim) {
         TAxis* axis = GetAxis(dim);
         Printf("    axis %d \"%s\": %d bins (%g..%g), %s bin sizes", dim,
                axis->GetTitle(), axis->GetNbins(), axis->GetXmin(), axis->GetXmax(),
                (axis->GetXbins() ? "variable" : "fixed"));
      }
   }

   if (optStat) {
      Printf("  %s error calculation", (GetCalculateErrors() ? "with" : "without"));
      if (GetCalculateErrors()) {
         Printf("    Sum(w)=%g, Sum(w^2)=%g", GetSumw(), GetSumw2());
         for (Int_t dim = 0; dim < fNdimensions; ++dim) {
            Printf("    axis %d: Sum(w*x)=%g, Sum(w*x^2)=%g", dim, GetSumwx(dim), GetSumwx2(dim));
         }
      }
   }

   if (optMem && InheritsFrom(THnSparse::Class())) {
      const THnSparse* hsparse = dynamic_cast<const THnSparse*>(this);
      Printf("  coordinates stored in %d chunks of %d entries\n    %g of bins filled using %g of memory compared to an array",
             hsparse->GetNChunks(), hsparse->GetChunkSize(),
             hsparse->GetSparseFractionBins(), hsparse->GetSparseFractionMem());
   }

   if (optContent) {
      Printf("  BIN CONTENT:");
      PrintEntries(0, -1, options);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Browse a THnSparse: create an entry (ROOT::THnSparseBrowsable) for each
/// dimension.

void THnBase::Browse(TBrowser *b)
{
   if (fBrowsables.IsEmpty()) {
      for (Int_t dim = 0; dim < fNdimensions; ++dim) {
         fBrowsables.AddAtAndExpand(new ROOT::Internal::THnBaseBrowsable(this, dim), dim);
      }
      fBrowsables.SetOwner();
   }

   for (Int_t dim = 0; dim < fNdimensions; ++dim) {
      b->Add(fBrowsables[dim]);
   }
}




/** \class ROOT::Internal::THnBaseBinIter
  Iterator over THnBase bins (internal implementation).
*/

/// Destruct a bin iterator.

ROOT::Internal::THnBaseBinIter::~THnBaseBinIter() {
   // Not much to do, but pin vtable
}



/** \class THnIter
   Iterator over THnBase bins
*/

ClassImp(THnIter);

THnIter::~THnIter() {
   // Destruct a bin iterator.
   delete fIter;
}


/** \class ROOT::Internal::THnBaseBrowsable
   TBrowser helper for THnBase.
*/


ClassImp(ROOT::Internal::THnBaseBrowsable);

////////////////////////////////////////////////////////////////////////////////
/// Construct a THnBaseBrowsable.

ROOT::Internal::THnBaseBrowsable::THnBaseBrowsable(THnBase* hist, Int_t axis):
fHist(hist), fAxis(axis), fProj(0)
{
   TString axisName = hist->GetAxis(axis)->GetName();
   if (axisName.IsNull()) {
      axisName = TString::Format("axis%d", axis);
   }

   SetNameTitle(axisName,
                TString::Format("Projection on %s of %s", axisName.Data(),
                                hist->IsA()->GetName()).Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Destruct a THnBaseBrowsable.

ROOT::Internal::THnBaseBrowsable::~THnBaseBrowsable()
{
   delete fProj;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse an axis of a THnBase, i.e. draw its projection.

void ROOT::Internal::THnBaseBrowsable::Browse(TBrowser* b)
{
   if (!fProj) {
      fProj = fHist->Projection(fAxis);
   }
   fProj->Draw(b ? b->GetDrawOption() : "");
   gPad->Update();
}

