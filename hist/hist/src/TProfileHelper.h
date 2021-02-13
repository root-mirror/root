// @(#)root/hist:$Id$
// Author: David Gonzalez Maline   18/01/2008

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProfileHelper
#define ROOT_TProfileHelper


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProfileHelper                                                       //
//                                                                      //
// Profile helper class.                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TError.h"
#include "THashList.h"
#include "TMath.h"
#include "TH1Merger.h"

class TProfileHelper {

public:
   template <typename T>
   static Bool_t Add(T* p, const TH1 *h1,  const TH1 *h2, Double_t c1, Double_t c2=1);

   template <typename T>
   static void BuildArray(T* p);

   template <typename T>
   static Double_t GetBinEffectiveEntries(T* p, Int_t bin);

   template <typename T>
   static Long64_t Merge(T* p, TCollection *list);

   template <typename T>
   static T* ExtendAxis(T* p, Double_t x, TAxis *axis);

   template <typename T>
   static void Scale(T* p, Double_t c1, Option_t * option);

   template <typename T>
   static void Sumw2(T* p, Bool_t flag );

   template <typename T>
   static void LabelsDeflate(T* p, Option_t *);

   template <typename T>
   static void LabelsInflate(T* p, Option_t *);

   template <typename T>
   static Double_t GetBinError(T* p, Int_t bin);

   template <typename T>
   static void SetBinEntries(T* p, Int_t bin, Double_t w);

   template <typename T>
   static void SetErrorOption(T* p, Option_t * opt);
};

template <typename T>
Bool_t TProfileHelper::Add(T* p, const TH1 *h1,  const TH1 *h2, Double_t c1, Double_t c2)
{
   // Performs the operation: this = c1*h1 + c2*h2

   T *p1 = (T*)h1;
   T *p2 = (T*)h2;

   // delete buffer if it is there since it will become invalid
   if (p->fBuffer) p->BufferEmpty(1);

// Check profile compatibility
   Int_t nx = p->GetNbinsX();
   Int_t ny = p->GetNbinsY();
   Int_t nz = p->GetNbinsZ();

   if ( nx != p1->GetNbinsX() ||  nx != p2->GetNbinsX() ||
        ny != p1->GetNbinsY() ||  ny != p2->GetNbinsY() ||
        nz != p1->GetNbinsZ() ||  nz != p2->GetNbinsZ() ) {
      Error("TProfileHelper::Add","Attempt to add profiles with different number of bins");
      return kFALSE;
   }

// Add statistics
   Double_t ac1 = TMath::Abs(c1);
   Double_t ac2 = TMath::Abs(c2);
   p->fEntries = ac1*p1->GetEntries() + ac2*p2->GetEntries();
   Double_t s0[TH1::kNstat], s1[TH1::kNstat], s2[TH1::kNstat];
   Int_t i;
   for (i=0;i<TH1::kNstat;i++) {s0[i] = s1[i] = s2[i] = 0;}
   p->GetStats(s0);
   p1->GetStats(s1);
   p2->GetStats(s2);
   for (i=0;i<TH1::kNstat;i++) {
      if (i == 1) s0[i] = c1*c1*s1[i] + c2*c2*s2[i];
      else        s0[i] = ac1*s1[i] + ac2*s2[i];
   }
   p->PutStats(s0);

   // create sumw2 per bin if not set
   if (p->fBinSumw2.fN == 0 && (p1->fBinSumw2.fN != 0 || p2->fBinSumw2.fN != 0)) p->Sumw2();

   // Make the loop over the bins to calculate the Addition
   Double_t *cu1 = p1->GetW();    Double_t *cu2 = p2->GetW();
   Double_t *er1 = p1->GetW2();   Double_t *er2 = p2->GetW2();
   Double_t *en1 = p1->GetB();    Double_t *en2 = p2->GetB();
   Double_t *ew1 = p1->GetB2();   Double_t *ew2 = p2->GetB2();
   // if p1 has not the sum of weight squared/bin stored use just the sum of weights
   if (ew1 == nullptr) ew1 = en1;
   if (ew2 == nullptr) ew2 = en2;
   for (Int_t bin = 0; bin < p->fN; bin++) {
      p->fArray[bin]             = c1*cu1[bin] + c2*cu2[bin];
      p->fSumw2.fArray[bin]      = ac1*er1[bin] + ac2*er2[bin];
      p->fBinEntries.fArray[bin] = ac1*en1[bin] + ac2*en2[bin];
      if (p->fBinSumw2.fN ) p->fBinSumw2.fArray[bin]  = ac1*ac1*ew1[bin] + ac2*ac2*ew2[bin];
   }
   return kTRUE;
}

template <typename T>
void TProfileHelper::BuildArray(T* p) {
   // Build the extra profile data structure in addition to the histograms
   // this are:    array of bin entries:  fBinEntries
   //              array of sum of profiled observable value - squared
   //              stored in TH1::fSumw2
   //              array of some of weight squared (optional) in TProfile::fBinSumw2
   p->fBinEntries.Set(p->fNcells);
   p->fSumw2.Set(p->fNcells);
   if (TH1::GetDefaultSumw2() || p->fBinSumw2.fN > 0 ) p->fBinSumw2.Set(p->fNcells);
}


template <typename T>
Double_t TProfileHelper::GetBinEffectiveEntries(T* p, Int_t bin)
{
//            Return bin effective entries for a weighted filled Profile histogram.
//            In case of an unweighted profile, it is equivalent to the number of entries per bin
//            The effective entries is defined as the square of the sum of the weights divided by the
//            sum of the weights square.
//            TProfile::Sumw2() must be called before filling the profile with weights.
//            Only by calling this method the  sum of the square of the weights per bin is stored.
//

   if (p->fBuffer) p->BufferEmpty();

   if (bin < 0 || bin >= p->fNcells) return 0;
   double sumOfWeights = p->fBinEntries.fArray[bin];
   if ( p->fBinSumw2.fN == 0 || p->fBinSumw2.fN != p->fNcells) {
      // this can happen  when reading an old file
      p->fBinSumw2.Set(0);
      return sumOfWeights;
   }
   double sumOfWeightsSquare = p->fBinSumw2.fArray[bin];
   return ( sumOfWeightsSquare > 0 ?  sumOfWeights * sumOfWeights /   sumOfWeightsSquare : 0 );
}

template <typename T>
Long64_t TProfileHelper::Merge(T* p, TCollection *li) {
   //Merge all histograms in the collection in this histogram.
   //This function computes the min/max for the axes,
   //compute a new number of bins, if necessary,
   //add bin contents, errors and statistics.
   //If overflows are present and limits are different the function will fail.
   //The function returns the total number of entries in the result histogram
   //if the merge is successful, -1 otherwise.
   //
   //IMPORTANT remark. The 2 axis x and y may have different number
   //of bins and different limits, BUT the largest bin width must be
   //a multiple of the smallest bin width and the upper limit must also
   //be a multiple of the bin width.

   if (!li) return 0;
   if (li->IsEmpty()) return (Int_t) p->GetEntries();

   TList inlist;
   inlist.AddAll(li);

   // use TH1Merger class
   TH1Merger merger(*p, *li, "");
   Bool_t ret = merger();

   return (ret) ? p->GetEntries() : -1;

#ifdef OLD_PROFILE_MERGE

   TAxis newXAxis;
   TAxis newYAxis;
   TAxis newZAxis;
   Bool_t initialLimitsFound = kFALSE;
   Bool_t allSameLimits = kTRUE;
   Bool_t allHaveLimits = kTRUE;
//   Bool_t firstNonEmptyHist = kTRUE;

   TIter next(&inlist);
   T* h = p;

   do {

      // skip empty histograms
      // if (h->fTsumw == 0 && h->GetEntries() == 0) continue;

      Bool_t hasLimits = h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax();
      allHaveLimits = allHaveLimits && hasLimits;

      if (hasLimits) {
         h->BufferEmpty();

#ifdef LATER
         // this is done in case the first histograms are empty and
         // the histogram have different limits
         if (firstNonEmptyHist ) {
            // set axis limits in the case the first histogram was empty
            if (h != p ) {
               if (!p->SameLimitsAndNBins(p->fXaxis, *(h->GetXaxis())) )
                  p->fXaxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
               if (!p->SameLimitsAndNBins(p->fYaxis, *(h->GetYaxis())) )
                  p->fYaxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
               if (!p->SameLimitsAndNBins(p->fZaxis, *(h->GetZaxis())) )
                  p->fZaxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),h->GetZaxis()->GetXmax());
            }
            firstNonEmptyHist = kFALSE;
         }
#endif

         // this is executed the first time an histogram with limits is found
         // to set some initial values on the new axis
         if (!initialLimitsFound) {
            initialLimitsFound = kTRUE;
            newXAxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
                     h->GetXaxis()->GetXmax());
            if ( p->GetDimension() >= 2 )
            newYAxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
                     h->GetYaxis()->GetXmax());
            if ( p->GetDimension() >= 3 )
            newZAxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),
                     h->GetZaxis()->GetXmax());
         }
         else {
            // check first if histograms have same bins
            if (!p->SameLimitsAndNBins(newXAxis, *(h->GetXaxis())) ||
                !p->SameLimitsAndNBins(newYAxis, *(h->GetYaxis())) ||
                !p->SameLimitsAndNBins(newZAxis, *(h->GetZaxis())) ) {

               allSameLimits = kFALSE;
               // recompute in this case the optimal limits
               // The condition to works is that the histogram have same bin with
               // and one common bin edge

               if (!p->RecomputeAxisLimits(newXAxis, *(h->GetXaxis()))) {
                  Error("TProfileHelper::Merge", "Cannot merge profiles %d dim - limits are inconsistent:\n "
                     "first: (%d, %f, %f), second: (%d, %f, %f)",p->GetDimension(),
                        newXAxis.GetNbins(), newXAxis.GetXmin(), newXAxis.GetXmax(),
                        h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
                        h->GetXaxis()->GetXmax());
                  return -1;
               }
               if (p->GetDimension() >= 2 && !p->RecomputeAxisLimits(newYAxis, *(h->GetYaxis()))) {
                  Error("TProfileHelper::Merge", "Cannot merge profiles %d dim - limits are inconsistent:\n "
                        "first: (%d, %f, %f), second: (%d, %f, %f)",p->GetDimension(),
                        newYAxis.GetNbins(), newYAxis.GetXmin(), newYAxis.GetXmax(),
                        h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
                        h->GetYaxis()->GetXmax());
                  return -1;
               }
               if (p->GetDimension() >= 3 && !p->RecomputeAxisLimits(newZAxis, *(h->GetZaxis()))) {
                  Error("TProfileHelper::Merge", "Cannot merge profiles %d dim - limits are inconsistent:\n "
                        "first: (%d, %f, %f), second: (%d, %f, %f)",p->GetDimension(),
                        newZAxis.GetNbins(), newZAxis.GetXmin(), newZAxis.GetXmax(),
                        h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),
                        h->GetZaxis()->GetXmax());
                  return -1;
               }
            }
         }
      }
   }  while ( ( h = dynamic_cast<T*> ( next() ) ) != NULL );
   if (!h && (*next) ) {
      Error("TProfileHelper::Merge","Attempt to merge object of class: %s to a %s",
            (*next)->ClassName(),p->ClassName());
      return -1;
   }

   next.Reset();

   // In the case of histogram with different limits
   // newX(Y)Axis will now have the new found limits
   // but one needs first to clone this histogram to perform the merge
   // The clone is not needed when all histograms have the same limits
   T * hclone = 0;
   if (!allSameLimits) {
      // We don't want to add the clone to gDirectory,
      // so remove our kMustCleanup bit temporarily
      Bool_t mustCleanup = p->TestBit(kMustCleanup);
      if (mustCleanup) p->ResetBit(kMustCleanup);
      hclone = (T*)p->IsA()->New();
      R__ASSERT(hclone);
      hclone->SetDirectory(0);
      p->Copy(*hclone);
      if (mustCleanup) p->SetBit(kMustCleanup);
      p->BufferEmpty(1);         // To remove buffer.
      p->Reset();                // BufferEmpty sets limits so we can't use it later.
      p->SetEntries(0);
      inlist.AddFirst(hclone);
   }

   if (!allSameLimits && initialLimitsFound) {
      Int_t b[] = { newXAxis.GetNbins(), newYAxis.GetNbins(), newZAxis.GetNbins() };
      Double_t v[] = { newXAxis.GetXmin(), newXAxis.GetXmax(),
                       newYAxis.GetXmin(), newYAxis.GetXmax(),
                       newZAxis.GetXmin(), newZAxis.GetXmax() };
      p->SetBins(b, v);
   }

   if (!allHaveLimits) {
      // fill this histogram with all the data from buffers of histograms without limits
      while ( (h = dynamic_cast<T*> (next()) ) ) {
         if (h->GetXaxis()->GetXmin() >= h->GetXaxis()->GetXmax() && h->fBuffer) {
             // no limits
            Int_t nbentries = (Int_t)h->fBuffer[0];
            Double_t v[5];
            for (Int_t i = 0; i < nbentries; i++)
               if ( p->GetDimension() == 3 ) {
                  v[0] = h->fBuffer[5*i + 2];
                  v[1] = h->fBuffer[5*i + 3];
                  v[2] = h->fBuffer[5*i + 4];
                  v[3] = h->fBuffer[5*i + 5];
                  v[4] = h->fBuffer[5*i + 1];
                  p->Fill(v);
               } else if ( p->GetDimension() == 2 ) {
                  v[0] = h->fBuffer[4*i + 2];
                  v[1] = h->fBuffer[4*i + 3];
                  v[2] = h->fBuffer[4*i + 4];
                  v[3] = h->fBuffer[4*i + 1];
                  v[4] = 0;
                  p->Fill(v);
               }
               else if ( p->GetDimension() == 1 ) {
                  v[0] = h->fBuffer[3*i + 2];
                  v[1] = h->fBuffer[3*i + 3];
                  v[2] = h->fBuffer[3*i + 1];
                  v[3] = v[4] = 0;
                  p->Fill(v);
               }
         }
      }
      if (!initialLimitsFound) {
         if (hclone) {
            inlist.Remove(hclone);
            delete hclone;
         }
         return (Int_t) p->GetEntries();  // all histograms have been processed
      }
      next.Reset();
   }

   //merge bin contents and errors
   Double_t stats[TH1::kNstat], totstats[TH1::kNstat];
   for (Int_t i=0;i<TH1::kNstat;i++) {totstats[i] = stats[i] = 0;}
   p->GetStats(totstats);
   Double_t nentries = p->GetEntries();
   Bool_t canExtend = p->CanExtendAllAxes();
   p->SetCanExtend(TH1::kNoAxis); // reset, otherwise setting the under/overflow will extend the axis

   while ( (h=static_cast<T*>(next())) ) {
      // process only if the histogram has limits; otherwise it was processed before

      if (h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax()) {
         // import statistics
         h->GetStats(stats);
         for (Int_t i = 0; i < TH1::kNstat; i++)
            totstats[i] += stats[i];
         nentries += h->GetEntries();

         for ( Int_t hbin = 0; hbin < h->fN; ++hbin ) {
            Int_t pbin = hbin;
            if (!allSameLimits) {
               // histogram have different limits:
               // find global bin number in p given the x,y,z axis bin numbers in h
               // in case of non equal axes
               // we can use FindBin on p axes because SetCanExtend(TH1::kNoAxis) has been called
               if ( h->GetW()[hbin] != 0 && (h->IsBinUnderflow(hbin) || h->IsBinOverflow(hbin)) ) {
                  // reject cases where underflow/overflow are there and bin content is not zero
                  Error("TProfileHelper::Merge", "Cannot merge profiles - they have"
                        " different limits and underflows/overflows are present."
                        " The initial profile is now broken!");
                  return -1;
               }
               Int_t hbinx, hbiny, hbinz;
               h->GetBinXYZ(hbin, hbinx, hbiny, hbinz);

               pbin = p->GetBin( p->fXaxis.FindBin( h->GetXaxis()->GetBinCenter(hbinx) ),
                                 p->fYaxis.FindBin( h->GetYaxis()->GetBinCenter(hbiny) ),
                                 p->fZaxis.FindBin( h->GetZaxis()->GetBinCenter(hbinz) ) );
            }


            p->fArray[pbin]             += h->GetW()[hbin];
            p->fSumw2.fArray[pbin]      += h->GetW2()[hbin];
            p->fBinEntries.fArray[pbin] += h->GetB()[hbin];
            if (p->fBinSumw2.fN) {
               if ( h->GetB2() ) p->fBinSumw2.fArray[pbin] += h->GetB2()[hbin];
               else p->fBinSumw2.fArray[pbin] += h->GetB()[hbin];
            }
         }
      }
   }
   if (canExtend) p->SetCanExtend(TH1::kAllAxes);

   //copy merged stats
   p->PutStats(totstats);
   p->SetEntries(nentries);
   if (hclone) {
      inlist.Remove(hclone);
      delete hclone;
   }
   return (Long64_t)nentries;
#endif
}

template <typename T>
T* TProfileHelper::ExtendAxis(T* p, Double_t x, TAxis *axis)
{
// Profile histogram is resized along axis such that x is in the axis range.
// The new axis limits are recomputed by doubling iteratively
// the current axis range until the specified value x is within the limits.
// The algorithm makes a copy of the histogram, then loops on all bins
// of the old histogram to fill the extended histogram.
// Takes into account errors (Sumw2) if any.
// The axis must be extendable before invoking this function.
// Ex: h->GetXaxis()->SetCanExtend(kTRUE)


   if (!axis->CanExtend()) return 0;
   if (axis->GetXmin() >= axis->GetXmax()) return 0;
   if (axis->GetNbins() <= 0) return 0;
   if (TMath::IsNaN(x)) { // x may be a NaN
      return 0;
   }

   Double_t xmin, xmax;
   if (!p->FindNewAxisLimits(axis, x, xmin, xmax))
      return 0;

   //save a copy of this histogram
   T* hold = (T*)p->IsA()->New();
   R__ASSERT(hold);
   hold->SetDirectory(0);
   p->Copy(*hold);
   //set new axis limits but keep same number of bins
   axis->SetLimits(xmin,xmax);
   if (p->fBinSumw2.fN) hold->Sumw2();

   // total bins (inclusing underflow /overflow)
   Int_t  nx = p->fXaxis.GetNbins() + 2;
   Int_t  ny = (p->GetDimension() > 1) ? p->fYaxis.GetNbins() + 2 : 1;
   Int_t  nz = (p->GetDimension() > 2) ? p->fZaxis.GetNbins() + 2 : 1;

   Int_t iaxis = 0;
   if (axis == p->GetXaxis()) iaxis = 1;
   if (axis == p->GetYaxis()) iaxis = 2;
   if (axis == p->GetZaxis()) iaxis = 3;
   Bool_t firstw = kTRUE;

   //now loop on all bins and refill
   p->Reset("ICE"); //reset only Integral, contents and Errors

   // need to consider also underflow/overflow in the non-extending axes
   Double_t xc,yc,zc;
   Int_t ix, iy, iz, binx, biny, binz;
   for (binz=0;binz< nz;binz++) {
      zc  = hold->GetZaxis()->GetBinCenter(binz);
      iz  = p->fZaxis.FindFixBin(zc);
      for (biny=0;biny<ny;biny++) {
         yc  = hold->GetYaxis()->GetBinCenter(biny);
         iy  = p->fYaxis.FindFixBin(yc);
         for (binx=0;binx<nx;binx++) {
            xc = hold->GetXaxis()->GetBinCenter(binx);
            ix  = p->fXaxis.FindFixBin(xc);
            Int_t sourceBin = hold->GetBin(binx,biny,binz);
            // skip empty bins
            if (hold->fBinEntries.fArray[sourceBin] == 0) continue;
            if (hold->IsBinUnderflow(sourceBin, iaxis) || hold->IsBinOverflow(sourceBin, iaxis)) {
               if (firstw) {
                  Warning("ExtendAxis",
                          "Histogram %s has underflow or overflow in the %s that is extendable"
                          " their content will be lost",p->GetName(),axis->GetName());
                  firstw = kFALSE;
               }
               continue;
            }
            Int_t destinationBin = p->GetBin(ix,iy,iz);
            p->AddBinContent(destinationBin, hold->fArray[sourceBin]);
            p->fBinEntries.fArray[destinationBin] += hold->fBinEntries.fArray[sourceBin];
            p->fSumw2.fArray[destinationBin] += hold->fSumw2.fArray[sourceBin];
            if (p->fBinSumw2.fN) p->fBinSumw2.fArray[destinationBin] += hold->fBinSumw2.fArray[sourceBin];
         }
      }
   }
   return hold;
}

template <typename T>
void TProfileHelper::Scale(T* p, Double_t c1, Option_t *)
{
   Double_t ac1 = TMath::Abs(c1);

   // Make the loop over the bins to calculate the Addition
   Int_t bin;
   Double_t *cu1 = p->GetW();
   Double_t *er1 = p->GetW2();
   Double_t *en1 = p->GetB();
   for (bin=0;bin<p->fN;bin++) {
      p->fArray[bin]             = c1*cu1[bin];
      p->fSumw2.fArray[bin]      = ac1*ac1*er1[bin];
      p->fBinEntries.fArray[bin] = en1[bin];
   }
}

template <typename T>
void TProfileHelper::Sumw2(T* p, Bool_t flag)
{
   // Create/Delete structure to store sum of squares of weights per bin  *-*-*-*-*-*-*-*
   //   This is needed to compute  the correct statistical quantities
   //    of a profile filled with weights
   //
   //
   //  This function is automatically called when the histogram is created
   //  if the static function TH1::SetDefaultSumw2 has been called before.

   if (!flag) {
      // clear array if existing or do nothing
      if (p->fBinSumw2.fN > 0 ) p->fBinSumw2.Set(0);
      return;
   }

   if ( p->fBinSumw2.fN == p->fNcells) {
      if (!p->fgDefaultSumw2)
         Warning("Sumw2","Sum of squares of profile bin weights structure already created");
      return;
   }

   p->fBinSumw2.Set(p->fNcells);

   // by default fill with the sum of weights which are stored in fBinEntries
   for (Int_t bin=0; bin<p->fNcells; bin++) {
      p->fBinSumw2.fArray[bin] = p->fBinEntries.fArray[bin];
   }
}

template <typename T>
void TProfileHelper::LabelsDeflate(T* p, Option_t *ax)
{
   // Reduce the number of bins for this axis to the number of bins having a label.
   // Works only for the given axis passed in the option
   // The method will remove only the extra bins existing after the last "labeled" bin.
   // Note that if there are "un-labeled" bins present between "labeled" bins they will not be removed


   TAxis *axis = p->GetXaxis();
   if (ax[0] == 'y' || ax[0] == 'Y') axis = p->GetYaxis();
   if (ax[0] == 'z' || ax[0] == 'Z') axis = p->GetZaxis();
   if (!axis) {
      Error("TProfileHelper::LabelsDeflate","Invalid axis option %s",ax);
      return;
   }
   if (!axis->GetLabels()) return;

   // find bin with last labels
   // bin number is object ID in list of labels
   // therefore max bin number is number of bins of the deflated histograms
   TIter next(axis->GetLabels());
   TObject *obj;
   Int_t nbins = 0;
   while ((obj = next())) {
      Int_t ibin = obj->GetUniqueID();
      if (ibin > nbins) nbins = ibin;
   }
   if (nbins < 1) nbins = 1;

   // do nothing in case it was the last bin
   if (nbins==axis->GetNbins()) return;

   T *hold = (T*)p->IsA()->New();;
   hold->SetDirectory(0);
   p->Copy(*hold);


   Double_t xmin = axis->GetXmin();
   Double_t xmax = axis->GetBinUpEdge(nbins);
   axis->SetRange(0,0);
   // set the new bins and range
   axis->Set(nbins,xmin,xmax);
   p->SetBinsLength(-1); // reset the number of cells
   p->fBinEntries.Set(p->fN);
   p->fSumw2.Set(p->fN);
   if (p->fBinSumw2.fN)  p->fBinSumw2.Set(p->fN);

   // reset the content
   p->Reset("ICE");

   //now loop on all bins and refill
   Int_t bin,binx,biny,binz;
   for (bin =0; bin < hold->fN; ++bin)
   {
      hold->GetBinXYZ(bin, binx, biny, binz);
      Int_t ibin = p->GetBin(binx, biny, binz);
      p->fArray[ibin] += hold->fArray[bin];
      p->fBinEntries.fArray[ibin] += hold->fBinEntries.fArray[bin];
      p->fSumw2.fArray[ibin] += hold->fSumw2.fArray[bin];
      if (p->fBinSumw2.fN) p->fBinSumw2.fArray[ibin] += hold->fBinSumw2.fArray[bin];
   }

   delete hold;
}

template <typename T>
void TProfileHelper::LabelsInflate(T* p, Option_t *ax)
{
// Double the number of bins for axis.
// Refill histogram
// This function is called by TAxis::FindBin(const char *label)
// Works only for the given axis

   if (gDebug) Info("LabelsInflate","Inflate label for axis %s of profile %s",ax,p->GetName());

   Int_t iaxis = p->AxisChoice(ax);
   TAxis *axis = 0;
   if (iaxis == 1) axis = p->GetXaxis();
   if (iaxis == 2) axis = p->GetYaxis();
   if (iaxis == 3) axis = p->GetZaxis();
   if (!axis) return;
   // TAxis *axis = p->GetXaxis();
   // if (ax[0] == 'y' || ax[0] == 'Y') axis = p->GetYaxis();

   T *hold = (T*)p->IsA()->New();;
   hold->SetDirectory(0);
   p->Copy(*hold);


//   Int_t  nbxold = p->fXaxis.GetNbins();
//   Int_t  nbyold = p->fYaxis.GetNbins();
   Int_t  nbins  = axis->GetNbins();
   Double_t xmin = axis->GetXmin();
   Double_t xmax = axis->GetXmax();
   xmax = xmin + 2*(xmax-xmin);
   axis->SetRange(0,0);
   // double the number of bins
   nbins *= 2;
   axis->Set(nbins,xmin,xmax);
   // reset the array of content according to the axis
   p->SetBinsLength(-1);
   Int_t ncells = p->fN;
   p->fBinEntries.Set(ncells);
   p->fSumw2.Set(ncells);
   if (p->fBinSumw2.fN)  p->fBinSumw2.Set(ncells);

   p->Reset("ICE");  // reset content and error

   //now loop on all old bins and refill excluding underflow/overflow in
   // the axis that has the bin doubled
   Int_t binx, biny, binz = 0;
   for (Int_t ibin =0; ibin < hold->fNcells; ibin++) {
      // get the binx,y,z values . The x-y-z (axis) bin values will stay the same between new-old after the expanding
      hold->GetBinXYZ(ibin,binx,biny,binz);
      Int_t bin = p->GetBin(binx,biny,binz);

      // underflow and overflow will be cleaned up because their meaning has been altered
      if (hold->IsBinUnderflow(ibin,iaxis) || hold->IsBinOverflow(ibin,iaxis)) {
         if (gDebug && hold->fBinEntries.fArray[ibin] > 0) Info("LabelsInflate","Content for underflow/overflow of bin (%d,%d,%d) will be lost",binx,biny,binz);
         continue;
      }
      else {
         p->fArray[bin] = hold->fArray[ibin];
         p->fBinEntries.fArray[bin] = hold->fBinEntries.fArray[ibin];
         p->fSumw2.fArray[bin] = hold->fSumw2.fArray[ibin];
         if (p->fBinSumw2.fN) p->fBinSumw2.fArray[bin] = hold->fBinSumw2.fArray[ibin];
         if (gDebug) Info("LabelsInflate","Copy Content from bin (%d,%d,%d) from %d in %d (%f,%f)",binx,biny,binz, ibin, bin, hold->fArray[ibin],hold->fBinEntries.fArray[ibin] );
      }
   }
   delete hold;
}

template <typename T>
void TProfileHelper::SetErrorOption(T* p, Option_t * option) {
   // set the profile option
   TString opt = option;
   opt.ToLower();
   p->fErrorMode = kERRORMEAN;
   if (opt.Contains("s")) p->fErrorMode = kERRORSPREAD;
   if (opt.Contains("i")) p->fErrorMode = kERRORSPREADI;
   if (opt.Contains("g")) p->fErrorMode = kERRORSPREADG;
}

template <typename T>
Double_t TProfileHelper::GetBinError(T* p, Int_t bin)
{
   // compute bin error of profile histograms

   if (p->fBuffer) p->BufferEmpty();

   if (bin < 0 || bin >= p->fNcells) return 0;
   Double_t cont = p->fArray[bin];                  // sum of bin w *y
   Double_t sum  = p->fBinEntries.fArray[bin];      // sum of bin weights
   Double_t err2 = p->fSumw2.fArray[bin];           // sum of bin w * y^2
   Double_t neff = p->GetBinEffectiveEntries(bin);  // (sum of w)^2 / (sum of w^2)
   if (sum == 0) return 0;      // for empty bins
   // case the values y are gaussian distributed y +/- sigma and w = 1/sigma^2
   if (p->fErrorMode == kERRORSPREADG) {
      return 1./TMath::Sqrt(sum);
   }
   // compute variance in y (eprim2) and standard deviation in y (eprim)
   Double_t contsum = cont/sum;
   Double_t eprim2  = TMath::Abs(err2/sum - contsum*contsum);
   Double_t eprim   = TMath::Sqrt(eprim2);

   if (p->fErrorMode == kERRORSPREADI) {
      if (eprim != 0) return eprim/TMath::Sqrt(neff);
      // in case content y is an integer (so each my has an error +/- 1/sqrt(12)
      // when the std(y) is zero
      return 1/TMath::Sqrt(12*neff);
   }

   // if approximate compute the sums (of w, wy and wy2) using all the bins
   //  when the variance in y is zero
   Double_t test = 1;
   if (err2 != 0 && neff < 5) test = eprim2*sum/err2;
   //Int_t cellLimit = (p->GetDimension() == 3)?1000404:10404;
   if (p->fgApproximate && (test < 1.e-4 || eprim2 <= 0)) {
      Double_t stats[TH1::kNstat] = {0};
      p->GetStats(stats);
      Double_t ssum = stats[0];
      // for 1D profile
      int index = 4;  // index in the stats array for 1D
      if (p->GetDimension() == 2) index = 7;   // for 2D
      if (p->GetDimension() == 3) index = 11;   // for 3D
      Double_t scont = stats[index];
      Double_t serr2 = stats[index+1];

      // compute mean and variance in y
      Double_t scontsum = scont/ssum;                                  // global mean
      Double_t seprim2  = TMath::Abs(serr2/ssum - scontsum*scontsum);  // global variance
      eprim           = 2*TMath::Sqrt(seprim2);                        // global std (why factor of 2 ??)
      sum = ssum;
   }
   sum = TMath::Abs(sum);

   // case option "S" return standard deviation in y
   if (p->fErrorMode == kERRORSPREAD) return eprim;

   // default case : fErrorMode = kERRORMEAN
   // return standard error on the mean of y
   //if (neff == 0) std::cerr << "NEFF = 0 for bin " << bin << "   " << eprim << "  " << neff << "  " << std::endl;
   return eprim/TMath::Sqrt(neff);

}


template <typename T>
void TProfileHelper::SetBinEntries(T* p, Int_t bin, Double_t w) {
//    Set the number of entries in bin for the profile
//    In case the profile stores the sum of weight squares - set the sum of weight square to the number entries
//    (i.e. assume the entries have been filled all with a weight == 1 )

   if (bin < 0 || bin >= p->fNcells) return;
   p->fBinEntries.fArray[bin] = w;
   if (p->fBinSumw2.fN) p->fBinSumw2.fArray[bin] = w;

}

#endif
