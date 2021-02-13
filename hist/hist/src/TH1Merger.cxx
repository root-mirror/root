// Helper clas implementing some of the TH1 functionality

#include "TH1Merger.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TError.h"
#include "THashList.h"
#include "TClass.h"
#include <iostream>
#include <limits>
#include <utility>

#define PRINTRANGE(a, b, bn)                                                                                          \
   Printf(" base: %f %f %d, %s: %f %f %d", a->GetXmin(), a->GetXmax(), a->GetNbins(), bn, b->GetXmin(), b->GetXmax(), \
          b->GetNbins());

Bool_t TH1Merger::AxesHaveLimits(const TH1 * h) {
   Bool_t hasLimits = h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax();
   if (h->GetDimension() > 1) hasLimits &=  h->GetYaxis()->GetXmin() < h->GetYaxis()->GetXmax();
   if (h->GetDimension() > 2) hasLimits &=  h->GetZaxis()->GetXmin() < h->GetZaxis()->GetXmax();
   return hasLimits;
}

/// Function performing the actual merge
Bool_t TH1Merger::operator() () {


   EMergerType type = ExamineHistograms();

   if (gDebug) Info("Merge","Histogram Merge type is %d and new axis flag is %d",(int) type,(int) fNewAxisFlag);

   if (type == kNotCompatible) return kFALSE;

   if (type == kAllSameAxes)
      return SameAxesMerge();

   if (type == kAllLabel)
      return LabelMerge();

   if (type == kAllNoLimits)
      return BufferMerge();

   if (type == kAutoP2HaveLimits || (type == kAutoP2NeedLimits && AutoP2BufferMerge()))
      return AutoP2Merge();

   // this is the mixed case - more complicated
   if (type == kHasNewLimits) {
      // we need to define some new axes
      DefineNewAxes();
      // we might need to merge some histogram using the buffer
      Bool_t ret =  BufferMerge();
      // if ret is true the merge is completed and we can exit
      if (ret) return kTRUE;
      // in the other cases then we merge using FindBin
      return DifferentAxesMerge();
   }
   Error("TH1Merger","Unknown type of Merge for histogram %s",fH0->GetName());
   return kFALSE;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Determine final boundaries and number of bins for histograms created in power-of-2
/// autobin mode.
///
/// Return kTRUE if compatible, updating fNewXaxis accordingly; return kFALSE if something
/// wrong.
///
/// The histograms are not merge-compatible if
///
///       1. have different variable-size bins
///       2. larger bin size is not an integer multiple of the smaller one
///       3. the final estimated range is smalle then the bin size
///

Bool_t TH1Merger::AutoP2BuildAxes(TH1 *h)
{
   // They must be both defined
   if (!h) {
      Error("AutoP2BuildAxes", "undefined histogram: %p", h);
      return kFALSE;
   }

   // They must be created in power-of-2 autobin mode
   if (!h->TestBit(TH1::kAutoBinPTwo)) {
      Error("AutoP2BuildAxes", "not in autobin-power-of-2 mode!");
      return kFALSE;
   }

   // Point to axes
   TAxis *a0 = &fNewXAxis, *a1 = h->GetXaxis();

   // This is for future merging of detached ranges (only possible if no over/underflows)
   Bool_t canextend = (h->GetBinContent(0) > 0 || h->GetBinContent(a1->GetNbins() + 1) > 0) ? kFALSE : kTRUE;

   // The first time we just copy the boundaries and bins
   if (a0->GetFirst() == a0->GetLast()) {
      a0->Set(a1->GetNbins(), a1->GetXmin(), a1->GetXmax());
      // This is for future merging of detached ranges (only possible if no over/underflows)
      a0->SetCanExtend(canextend);
      return kTRUE;
   }

   // Bin sizes must be in integer ratio
   Double_t bwmax = (a0->GetXmax() - a0->GetXmin()) / a0->GetNbins();
   Double_t bwmin = (a1->GetXmax() - a1->GetXmin()) / a1->GetNbins();
   Bool_t b0 = kTRUE;
   if (bwmin > bwmax) {
      std::swap(bwmax, bwmin);
      b0 = kFALSE;
   }
   if (!(bwmin > 0.)) {
      PRINTRANGE(a0, a1, h->GetName());
      Error("AutoP2BuildAxes", "minimal bin width negative or null: %f", bwmin);
      return kFALSE;
   }

   Double_t rt;
   Double_t re = std::modf(bwmax / bwmin, &rt);
   if (re > std::numeric_limits<Double_t>::epsilon()) {
      PRINTRANGE(a0, a1, h->GetName());
      Error("AutoP2BuildAxes", "bin widths not in integer ratio: %f", re);
      return kFALSE;
   }

   // Range of the merged histogram, taking into account overlaps
   Bool_t domax = kFALSE;
   Double_t xmax, xmin;
   if (a0->GetXmin() < a1->GetXmin()) {
      if (a0->GetXmax() < a1->GetXmin()) {
         if (!a0->CanExtend() || !canextend) {
            PRINTRANGE(a0, a1, h->GetName());
            Error("AutoP2BuildAxes", "ranges are disconnected and under/overflows: cannot merge");
            return kFALSE;
         }
         xmax = a1->GetXmax();
         xmin = a0->GetXmin();
         domax = b0 ? kTRUE : kFALSE;
      } else {
         if (a0->GetXmax() >= a1->GetXmax()) {
            xmax = a1->GetXmax();
            xmin = a1->GetXmin();
            domax = !b0 ? kTRUE : kFALSE;
         } else {
            xmax = a0->GetXmax();
            xmin = a1->GetXmin();
            domax = !b0 ? kTRUE : kFALSE;
         }
      }
   } else {
      if (a1->GetXmax() < a0->GetXmin()) {
         if (!a0->CanExtend() || !canextend) {
            PRINTRANGE(a0, a1, h->GetName());
            Error("AutoP2BuildAxes", "ranges are disconnected and under/overflows: cannot merge");
            return kFALSE;
         }
         xmax = a0->GetXmax();
         xmin = a1->GetXmin();
         domax = !b0 ? kTRUE : kFALSE;
      } else {
         if (a1->GetXmax() >= a0->GetXmax()) {
            xmax = a0->GetXmax();
            xmin = a0->GetXmin();
            domax = b0 ? kTRUE : kFALSE;
         } else {
            xmax = a1->GetXmax();
            xmin = a0->GetXmin();
            domax = b0 ? kTRUE : kFALSE;
         }
      }
   }
   Double_t range = xmax - xmin;

   re = std::modf(range / bwmax, &rt);
   if (rt < 1.) {
      PRINTRANGE(a0, a1, h->GetName());
      Error("MergeCompatibleHistograms", "range smaller than bin width: %f %f %f", range, bwmax, rt);
      return kFALSE;
   }
   if (re > std::numeric_limits<Double_t>::epsilon()) {
      if (domax) {
         xmax -= bwmax * re;
      } else {
         xmin += bwmax * re;
      }
   }
   // Number of bins
   Int_t nb = (Int_t)rt;

   // Set the result
   a0->Set(nb, xmin, xmax);

   // This is for future merging of detached ranges (only possible if no over/underflows)
   if (!a0->CanExtend())
      a0->SetCanExtend(canextend);

   // Done
   return kTRUE;
}

/**
   Examine the list of histograms to find out which type of Merge we need to do
   Pass the input list containing the histogram to merge and h0 which is the initial histogram
   on which all the histogram of the list will be merged
   This are the possible cases:
    - 1. All histogram have the same axis (allSameLimits = true)
    - 2. Histogram have different axis but compatible  (allSameLimits = false) and sameLimitsX,Y,Z specifies which axis
         has different limits
    - 3. Histogram do not have limits (so the Buffer is used)  allHaveLimits = false
    - 3b. One histogram has limits the other not : allHaveLimits = false AND initialLimitsFound = true
    - 4. Histogram Have labels  = allHaveLabels = true


*/
TH1Merger::EMergerType TH1Merger::ExamineHistograms() {



   Bool_t initialLimitsFound = kFALSE;
   Bool_t allHaveLabels = kTRUE;  // assume all histo have labels and check later
   UInt_t labelAxisType = TH1::kNoAxis;    // type of axes that have label
   Bool_t allHaveLimits = kTRUE;
   Bool_t allSameLimits = kTRUE;
   Bool_t sameLimitsX = kTRUE;
   Bool_t sameLimitsY = kTRUE;
   Bool_t sameLimitsZ = kTRUE;
   Bool_t foundLabelHist = kFALSE;
   Bool_t haveWeights = kFALSE;

   Bool_t isAutoP2 = kFALSE;

   // TAxis newXAxis;
   // TAxis newYAxis;
   // TAxis newZAxis;

   TIter next(&fInputList);
   TH1 * h = fH0;  // start with fH0

   int dimension = fH0->GetDimension();

   isAutoP2 = fH0->TestBit(TH1::kAutoBinPTwo) ? kTRUE : kFALSE;

   // if the option alphanumeric merge is set
   // we assume we do not have labels
   if (fNoLabelMerge)  allHaveLabels = kFALSE;

   // start looping on the histograms

   do  {

      // check first histogram compatibility
      if (h != fH0) {
         if (h->GetDimension() != dimension) {
            Error("Merge", "Cannot merge histogram - dimensions are different\n "
                  "%s has dim=%d and %s has dim=%d",fH0->GetName(),dimension,h->GetName(),h->GetDimension());
            return kNotCompatible;
         }
      }

      // check if one of the histogram is weighted
      haveWeights |= h->GetSumw2N() != 0;

      // do not skip anymore empty histograms
      // since are used to set the limits
      Bool_t hasLimits = TH1Merger::AxesHaveLimits(h);
      allHaveLimits = allHaveLimits && hasLimits;
      allSameLimits &= allHaveLimits;

      if (isAutoP2 && !h->TestBit(TH1::kAutoBinPTwo)) {
         Error("Merge", "Cannot merge histogram - some are in autobin-power-of-2 mode, but not %s!", h->GetName());
         return kNotCompatible;
      }
      if (!isAutoP2 && h->TestBit(TH1::kAutoBinPTwo)) {
         Error("Merge", "Cannot merge histogram - %s is in autobin-power-of-2 mode, but not the previous ones",
               h->GetName());
         return kNotCompatible;
      }

      if (hasLimits) {
         h->BufferEmpty();

//          // this is done in case the first histograms are empty and
//          // the histogram have different limits
// #ifdef LATER
//          if (firstHistWithLimits ) {
//             // set axis limits in the case the first histogram did not have limits
//             if (h != this && !SameLimitsAndNBins( fXaxis, *h->GetXaxis()) ) {
//               if (h->GetXaxis()->GetXbins()->GetSize() != 0) fXaxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXbins()->GetArray());
//               else                                           fXaxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
//             }
//             firstHistWithLimits = kFALSE;
//          }
// #endif

         // this is executed the first time an histogram with limits is found
         // to set some initial values on the new axis
         if (!initialLimitsFound) {
            initialLimitsFound = kTRUE;
            if (h->GetXaxis()->GetXbins()->GetSize() != 0)
               fNewXAxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXbins()->GetArray());
            else
               fNewXAxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            if (dimension > 1) {
               if (h->GetYaxis()->GetXbins()->GetSize() != 0)
                  fNewYAxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXbins()->GetArray());
               else
                  fNewYAxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
            }
            if (dimension > 2) {
               if (h->GetZaxis()->GetXbins()->GetSize() != 0)
                  fNewZAxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXbins()->GetArray());
               else
                  fNewZAxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax());

            }
         }
         else {
            // check first if histograms have same bins in X
            if (!TH1::SameLimitsAndNBins(fNewXAxis, *(h->GetXaxis())) ) {
               sameLimitsX = kFALSE;
               // recompute the limits in this case the optimal limits
               // The condition to works is that the histogram have same bin with
               // and one common bin edge
               if (!TH1::RecomputeAxisLimits(fNewXAxis, *(h->GetXaxis()))) {
                  Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
                        "first: %s (%d, %f, %f), second: %s (%d, %f, %f)", fH0->GetName(),
                        fNewXAxis.GetNbins(), fNewXAxis.GetXmin(), fNewXAxis.GetXmax(),
                        h->GetName(),h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
                        h->GetXaxis()->GetXmax());
                  return kNotCompatible;
               }
            }
            // check first if histograms have same bins in Y
            if (dimension > 1 && !TH1::SameLimitsAndNBins(fNewYAxis, *(h->GetYaxis()))) {
             sameLimitsY = kFALSE;
             // recompute in this case the optimal limits
             // The condition to works is that the histogram have same bin with
             // and one common bin edge
             if (!TH1::RecomputeAxisLimits(fNewYAxis, *(h->GetYaxis()))) {
               Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
                     "first: %s (%d, %f, %f), second: %s (%d, %f, %f)", fH0->GetName(),
                     fNewYAxis.GetNbins(), fNewYAxis.GetXmin(), fNewYAxis.GetXmax(),
                     h->GetName(), h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
                     h->GetYaxis()->GetXmax());
               return kNotCompatible;
             }
           }
           if(dimension > 2 && !fH0->SameLimitsAndNBins(fNewZAxis, *(h->GetZaxis()))) {
             sameLimitsZ = kFALSE;
             if (!TH1::RecomputeAxisLimits(fNewZAxis, *(h->GetZaxis()))) {
               Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
                     "first: %s (%d, %f, %f), second: %s (%d, %f, %f)", fH0->GetName(),
                     fNewZAxis.GetNbins(), fNewZAxis.GetXmin(), fNewZAxis.GetXmax(),
                     h->GetName(),h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),
                     h->GetZaxis()->GetXmax());
               return kNotCompatible;
             }
           }
           allSameLimits = sameLimitsX && sameLimitsY && sameLimitsZ;


         }
      }
      Bool_t histoIsEmpty = h->IsEmpty();
      // std::cout << "considering histo " << h->GetName() << "  labels - " << allHaveLabels << " is empty "
      //           <<  histoIsEmpty << std::endl;

      // if histogram is empty it does not matter if it has label or not
      if (allHaveLabels && !histoIsEmpty) {
         THashList* hlabelsX = h->GetXaxis()->GetLabels();
         THashList* hlabelsY = (dimension > 1) ? h->GetYaxis()->GetLabels() : nullptr;
         THashList* hlabelsZ = (dimension > 2) ? h->GetZaxis()->GetLabels() : nullptr;
         Bool_t haveOneLabelX = hlabelsX != nullptr;
         Bool_t haveOneLabelY = hlabelsY != nullptr;
         Bool_t haveOneLabelZ = hlabelsZ != nullptr;
         Bool_t haveOneLabel = haveOneLabelX || haveOneLabelY || haveOneLabelZ;
         // do here to print message only one time
         if (foundLabelHist && allHaveLabels && !haveOneLabel) {
            Warning("Merge","Not all histograms have labels. I will ignore labels,"
            " falling back to bin numbering mode.");
         }

         allHaveLabels &= (haveOneLabel);

         if (haveOneLabel) {
            foundLabelHist = kTRUE;
            UInt_t type = 0;
            if (haveOneLabelX) type |= TH1::kXaxis;
            if (haveOneLabelY) type |= TH1::kYaxis;
            if (haveOneLabelZ) type |= TH1::kZaxis;
            if (labelAxisType == TH1::kNoAxis) labelAxisType = type;
            // check if all histogram have consistent label axis
            // this means that there is at least one axis where boith histogram have labels
            Bool_t consistentLabels = (type & labelAxisType) != TH1::kNoAxis;
            allHaveLabels &= consistentLabels;
            if (!consistentLabels)
               Warning("TH1Merger::ExamineHistogram","Histogram %s has inconsistent labels: %d is not consistent with  %d",
                     h->GetName(), (int) type, (int) labelAxisType );
            if (gDebug && consistentLabels)
               Info("TH1Merger::ExamineHistogram","Histogram %s has consistent labels",h->GetName() );
         }

         // Check compatibility of axis that have labels with axis that can be extended
         UInt_t extendAxisType = TH1::kNoAxis;
         if (fH0->GetXaxis()->CanExtend()) extendAxisType |= TH1::kXaxis;
         if (dimension > 1 && fH0->GetYaxis()->CanExtend()) extendAxisType |= TH1::kYaxis;
         if (dimension > 2 && fH0->GetZaxis()->CanExtend()) extendAxisType |= TH1::kZaxis;
         // it is sufficient to have a consistent label axis that can be extended
         Bool_t labelAxisCanBeExtended = ((extendAxisType & labelAxisType) != TH1::kNoAxis);
         // If histograms have labels but corresponding axes cannot be extended use bin center mode
         if (allHaveLabels && !labelAxisCanBeExtended)   {
            // special case for this histogram when is empty
            // and axis cannot be extended (because it is the default)
            if ( fH0->IsEmpty()  ) {
                if (gDebug)
                   Info("TH1Merger::ExamineHistogram","Histogram %s to be merged is empty and we are merging with %s that has labels. Force the axis to be extended",fH0->GetName(),h->GetName());
               fH0->SetCanExtend( labelAxisType );
            }
            else { // histogram is not empty
               if (gDebug)
                  Info("TH1Merger::ExamineHistogram","Histogram %s to be merged has labels but corresponding axis cannot be extended - using bin numeric mode to merge. Call TH1::SetCanExtend(TH1::kAllAxes) if want to merge using label mode",fH0->GetName());
               allHaveLabels = kFALSE;
            }
         }
         // we don;t need to check anymore for case of non=empty histograms with some labels.
         // If we have some labels set ans axis is not extendable the LabelsMerge function handles
         // that case correctly
#if 0
         if (allHaveLabels ) {
            // count number of bins with non-null content
            Int_t non_zero_bins = 0;
            // loop on axis that have labels. Support this only for 1D histogram
            if (hlabelsX && dimension == 1 && !h->GetXaxis()->CanExtend()) {
               Int_t nbins = h->GetXaxis()->GetNbins();
               if (nbins > hlabelsX->GetEntries() ) {
                  for (Int_t i = 1; i <= nbins; i++) {
                     if (h->RetrieveBinContent(i) != 0 || (h->fSumw2.fN && h->GetBinError(i) != 0) ) {
                        non_zero_bins++;
                     }
                  }
                  if (non_zero_bins > hlabelsX->GetEntries() ) {
                     Warning("TH1Merger::ExamineHistograms","Histogram %s contains non-empty bins without labels - falling back to bin numbering mode",h->GetName() );
                     allHaveLabels = kFALSE;
                  }
               }
            }
            // for multidimensional case check that labels size is less than axis size

            if (dimension > 1 ) {
               if (hlabelsX && !h->GetXaxis()->CanExtend() && hlabelsX->GetEntries() < h->GetXaxis()->GetNbins()-1 ) { // use -1 because one bin without label is like a dummy label
                  Warning("TH1Merger::ExamineHistograms","Histogram %s has the X axis containing more than one bin without labels - falling back to bin numbering mode",h->GetName() );
                  allHaveLabels = kFALSE;
               }
               if (hlabelsY && !h->GetYaxis()->CanExtend() && hlabelsY->GetEntries() < h->GetYaxis()->GetNbins()-1 ) { // use -1 because one bin without label is like a dummy label
                  Warning("TH1Merger::ExamineHistograms","Histogram %s has the Y axis containing more than one bin without labels - falling back to bin numbering mode",h->GetName() );
                  allHaveLabels = kFALSE;
               }
               if (hlabelsZ && !h->GetZaxis()->CanExtend() && hlabelsZ->GetEntries() < h->GetZaxis()->GetNbins()-1 ) { // use -1 because one bin without label is like a dummy label
                  Warning("TH1Merger::ExamineHistograms","Histogram %s has the Z axis containing more than one bin without labels - falling back to bin numbering mode",h->GetName() );
                  allHaveLabels = kFALSE;
               }
            }
         }
#endif
      }
      if (gDebug)
         Info("TH1Merger::ExamineHistogram","Examine histogram %s - labels %d - same limits %d - axis found %d",h->GetName(),allHaveLabels,allSameLimits,initialLimitsFound );

   }    while ( ( h = dynamic_cast<TH1*> ( next() ) ) != NULL );

   if (!h && (*next) ) {
      Error("Merge","Attempt to merge object of class: %s to a %s",
            (*next)->ClassName(),fH0->ClassName());
      return kNotCompatible;
   }

   // in case of weighted histogram set Sumw2() on fH0 is is not weighted
   if (haveWeights && fH0->GetSumw2N() == 0)
      fH0->Sumw2();

   // AutoP2
   if (isAutoP2) {
      if (allHaveLimits)
         return kAutoP2HaveLimits;
      return kAutoP2NeedLimits;
   }

   // return the type of merge
   if (allHaveLabels) return kAllLabel;
   if (allSameLimits) return kAllSameAxes;
   if (!initialLimitsFound) {
      R__ASSERT(!allHaveLimits);
      // case where no limits are found and the buffer is used
      return kAllNoLimits;
   }
   // remaining case should be the mixed one. Some histogram have limits some not
   fNewAxisFlag = 0;
   if (!sameLimitsX) fNewAxisFlag |= TH1::kXaxis;
   if (!sameLimitsY) fNewAxisFlag |= TH1::kYaxis;
   if (!sameLimitsZ) fNewAxisFlag |= TH1::kZaxis;

   // we need to set the flag also in case this histogram has no limits
   // we need to set explicitly the flag to re-define  a new axis
   if (fH0->GetXaxis()->GetXmin() >= fH0->GetXaxis()->GetXmax()) fNewAxisFlag |= TH1::kXaxis;
   if (dimension > 1 && fH0->GetYaxis()->GetXmin() >= fH0->GetYaxis()->GetXmax()) fNewAxisFlag |= TH1::kYaxis;
   if (dimension > 2 && fH0->GetZaxis()->GetXmin() >= fH0->GetZaxis()->GetXmax()) fNewAxisFlag |= TH1::kZaxis;


   return kHasNewLimits;

}

/**
   Function to define new histogram axis when merging
   It is call only in case of merging with different axis or with the
   buffer  (kHasNewLimits)
*/

void TH1Merger::DefineNewAxes() {

   // first we need to create a copy of the histogram in case is not empty

   if (!fH0->IsEmpty() ) {
      Bool_t mustCleanup = fH0->TestBit(kMustCleanup);
      if (mustCleanup) fH0->ResetBit(kMustCleanup);
      fHClone = (TH1*)fH0->IsA()->New();
      fHClone->SetDirectory(0);
      fH0->Copy(*fHClone);
      if (mustCleanup) fH0->SetBit(kMustCleanup);
      fH0->BufferEmpty(1);         // To remove buffer.
      fH0->Reset();                // BufferEmpty sets limits so we can't use it later.
      fH0->SetEntries(0);
      fInputList.AddFirst(fHClone);

   }

   bool newLimitsX = (fNewAxisFlag & TH1::kXaxis);
   bool newLimitsY = (fNewAxisFlag & TH1::kYaxis);
   bool newLimitsZ = (fNewAxisFlag & TH1::kZaxis);
   if (newLimitsX) {
      fH0->fXaxis.SetRange(0,0);
      if (fNewXAxis.GetXbins()->GetSize() != 0)
          fH0->fXaxis.Set(fNewXAxis.GetNbins(),fNewXAxis.GetXbins()->GetArray());
       else
          fH0->fXaxis.Set(fNewXAxis.GetNbins(),fNewXAxis.GetXmin(), fNewXAxis.GetXmax());
   }
   if (newLimitsY) {
      fH0->fYaxis.SetRange(0,0);
      if (fNewYAxis.GetXbins()->GetSize() != 0)
          fH0->fYaxis.Set(fNewYAxis.GetNbins(),fNewYAxis.GetXbins()->GetArray());
       else
          fH0->fYaxis.Set(fNewYAxis.GetNbins(),fNewYAxis.GetXmin(), fNewYAxis.GetXmax());
   }
   if (newLimitsZ) {
      fH0->fZaxis.SetRange(0,0);
      if (fNewZAxis.GetXbins()->GetSize() != 0)
          fH0->fZaxis.Set(fNewZAxis.GetNbins(),fNewZAxis.GetXbins()->GetArray());
       else
          fH0->fZaxis.Set(fNewZAxis.GetNbins(),fNewZAxis.GetXmin(), fNewZAxis.GetXmax());
   }

   // we need to recompute fNcells and set the array size (as in TH1::SetBins)
   fH0->fNcells = fH0->fXaxis.GetNbins()+2;
   if (fH0->fDimension > 1) fH0->fNcells *= fH0->fYaxis.GetNbins()+2;
   if (fH0->fDimension > 2) fH0->fNcells *= fH0->fZaxis.GetNbins()+2;
   fH0->SetBinsLength(fH0->fNcells);
   if (fH0->fSumw2.fN) fH0->fSumw2.Set(fH0->fNcells);
   // set dummy Y and Z axis for lower dim histogras
   if (fH0->fDimension < 3)  fH0->fZaxis.Set(1,0,1);
   if (fH0->fDimension < 2)  fH0->fYaxis.Set(1,0,1);

   if (gDebug) {
      if (newLimitsX) Info("DefineNewAxis","A new X axis has been defined Nbins=%d , [%f,%f]", fH0->fXaxis.GetNbins(),
                           fH0->fXaxis.GetXmin(), fH0->fXaxis.GetXmax() );
      if (newLimitsY) Info("DefineNewAxis","A new Y axis has been defined Nbins=%d , [%f,%f]", fH0->fYaxis.GetNbins(),
                           fH0->fYaxis.GetXmin(), fH0->fYaxis.GetXmax() );
      if (newLimitsZ) Info("DefineNewAxis","A new Z axis has been defined Nbins=%d , [%f,%f]", fH0->fZaxis.GetNbins(),
                           fH0->fZaxis.GetXmin(), fH0->fZaxis.GetXmax() );
   }

   return;

}

void TH1Merger::CopyBuffer(TH1 *hsrc, TH1 *hdes)
{
   // Check inputs
   //if (!hsrc || !hsrc->fBuffer || !hdes || !hdes->fBuffer) {
   if (!hsrc || !hsrc->fBuffer || !hdes ) {
      void *p1 = hsrc ? hsrc->fBuffer : 0;
      //void *p2 = hdes ? hdes->fBuffer : 0;
      //Warning("TH1Merger::CopyMerge", "invalid inputs: %p, %p, %p, %p -> do nothing", hsrc, hdes, p1, p2);
      Warning("TH1Merger::CopyMerge", "invalid inputs: %p, %p, %p, -> do nothing", hsrc, hdes, p1);
   }

   // Entries from buffers have to be filled one by one
   // because FillN doesn't resize histograms.
   Int_t nbentries = (Int_t)hsrc->fBuffer[0];
   if (hdes->fDimension == 1) {
      for (Int_t i = 0; i < nbentries; i++)
         hdes->Fill(hsrc->fBuffer[2 * i + 2], hsrc->fBuffer[2 * i + 1]);
   }
   if (hdes->fDimension == 2) {
      auto h2 = dynamic_cast<TH2 *>(hdes);
      R__ASSERT(h2);
      for (Int_t i = 0; i < nbentries; i++)
         h2->Fill(hsrc->fBuffer[3 * i + 2], hsrc->fBuffer[3 * i + 3], hsrc->fBuffer[3 * i + 1]);
   }
   if (hdes->fDimension == 3) {
      auto h3 = dynamic_cast<TH3 *>(hdes);
      R__ASSERT(h3);
      for (Int_t i = 0; i < nbentries; i++)
         h3->Fill(hsrc->fBuffer[4 * i + 2], hsrc->fBuffer[4 * i + 3], hsrc->fBuffer[4 * i + 4],
                  hsrc->fBuffer[4 * i + 1]);
   }
}

Bool_t TH1Merger::AutoP2BufferMerge()
{

   TH1 *href = 0, *hist = 0;
   TIter nextref(&fInputList);
   if (TH1Merger::AxesHaveLimits(fH0)) {
      href = fH0;
   } else {
      while ((hist = (TH1 *)nextref()) && !href) {
         if (TH1Merger::AxesHaveLimits(hist))
            href = hist;
      }
   }
   Bool_t resetfH0 = kFALSE;
   if (!href) {
      // Merge all histograms to fH0 and do a final projection
      href = fH0;
   } else {
      if (href != fH0) {
         // Temporary add fH0 to the list for buffer merging
         fInputList.Add(fH0);
         resetfH0 = kTRUE;
      }
   }
   TIter next(&fInputList);
   while ((hist = (TH1 *)next())) {
      if (!TH1Merger::AxesHaveLimits(hist) && hist->fBuffer) {
         if (gDebug)
            Info("AutoP2BufferMerge", "merging buffer of %s into %s", hist->GetName(), href->GetName());
         CopyBuffer(hist, href);
         fInputList.Remove(hist);
      }
   }
   // Final projection
   if (href->fBuffer)
      href->BufferEmpty(1);
   // Reset fH0, if already added, to avoid double counting
   if (resetfH0)
      fH0->Reset("ICES");
   // Done, all histos have been processed
   return kTRUE;
}

Bool_t TH1Merger::AutoP2Merge()
{

   Double_t stats[TH1::kNstat], totstats[TH1::kNstat];
   for (Int_t i = 0; i < TH1::kNstat; i++) {
      totstats[i] = stats[i] = 0;
   }

   TIter next(&fInputList);
   TH1 *hist = 0;
   // Calculate boundaries and bins
   Double_t xmin = 0., xmax = 0.;
   if (!(fH0->IsEmpty())) {
      hist = fH0;
   } else {
      while ((hist = (TH1 *)next())) {
         if (!hist->IsEmpty())
            break;
      }
   }

   if (!hist) {
      if (gDebug)
         Info("TH1Merger::AutoP2Merge", "all histograms look empty!");
      return kFALSE;
   }

   // Start building the axes from the reference histogram
   if (!AutoP2BuildAxes(hist)) {
      Error("TH1Merger::AutoP2Merge", "cannot create axes from %s", hist->GetName());
      return kFALSE;
   }
   TH1 *h = 0;
   while ((h = (TH1 *)next())) {
      if (!AutoP2BuildAxes(h)) {
         Error("TH1Merger::AutoP2Merge", "cannot merge histogram %s: not merge compatible", h->GetName());
         return kFALSE;
      }
   }
   xmin = fNewXAxis.GetXmin();
   xmax = fNewXAxis.GetXmax();
   Int_t nbins = fNewXAxis.GetNbins();

   // Prepare stats
   fH0->GetStats(totstats);
   // Clone fH0 and add it to the list
   if (!fH0->IsEmpty())
      fInputList.Add(fH0->Clone());

   // reset fH0
   fH0->Reset("ICES");
   // Set the new boundaries
   fH0->SetBins(nbins, xmin, xmax);

   next.Reset();
   Double_t nentries = 0.;
   while ((hist = (TH1 *)next())) {
      // process only if the histogram has limits; otherwise it was processed before
      // in the case of an existing buffer (see if statement just before)

      if (gDebug)
         Info("TH1Merger::AutoP2Merge", "merging histogram %s into %s (entries: %f)", hist->GetName(), fH0->GetName(),
              hist->GetEntries());

      // skip empty histograms
      if (hist->IsEmpty())
         continue;

      // import statistics
      hist->GetStats(stats);
      for (Int_t i = 0; i < TH1::kNstat; i++)
         totstats[i] += stats[i];
      nentries += hist->GetEntries();

      // Int_t nx = hist->GetXaxis()->GetNbins();
      // loop on bins of the histogram and do the merge
      for (Int_t ibin = 0; ibin < hist->fNcells; ibin++) {

         Double_t cu = hist->RetrieveBinContent(ibin);
         Double_t e1sq = TMath::Abs(cu);
         if (fH0->fSumw2.fN)
            e1sq = hist->GetBinErrorSqUnchecked(ibin);

         Double_t xu = hist->GetBinCenter(ibin);
         Int_t jbin = fH0->FindBin(xu);

         fH0->AddBinContent(jbin, cu);
         if (fH0->fSumw2.fN)
            fH0->fSumw2.fArray[jbin] += e1sq;
      }
   }
   // copy merged stats
   fH0->PutStats(totstats);
   fH0->SetEntries(nentries);

   return kTRUE;
}

Bool_t TH1Merger::BufferMerge()
{

   TIter next(&fInputList);
   while (TH1* hist = (TH1*)next()) {
      // support also case where some histogram have limits and some have the buffer
      if ( !TH1Merger::AxesHaveLimits(hist) && hist->fBuffer  ) {

         if (gDebug)
            Info("TH1Merger::BufferMerge","Merging histogram %s into %s",hist->GetName(), fH0->GetName() );
         CopyBuffer(hist, fH0);
         fInputList.Remove(hist);
      }
   }
   // return true if the merge is completed
   if (fInputList.GetSize() == 0) {
      // all histo have been merged
      return kTRUE;
   }
   // we need to reset the buffer in case of merging later on
   // is this really needed ???
   if (fH0->fBuffer) fH0->BufferEmpty(1);

   return kFALSE;
}

Bool_t TH1Merger::SameAxesMerge() {


   Double_t stats[TH1::kNstat], totstats[TH1::kNstat];
   for (Int_t i=0;i<TH1::kNstat;i++) {
      totstats[i] = stats[i] = 0;
   }
   fH0->GetStats(totstats);
   Double_t nentries = fH0->GetEntries();

   TIter next(&fInputList);
   while (TH1* hist=(TH1*)next()) {
      // process only if the histogram has limits; otherwise it was processed before
      // in the case of an existing buffer (see if statement just before)

      if (gDebug)
         Info("TH1Merger::SameAxesMerge","Merging histogram %s into %s",hist->GetName(), fH0->GetName() );

      // skip empty histograms
      if (hist->IsEmpty()) continue;

      // import statistics
      hist->GetStats(stats);
      for (Int_t i=0; i<TH1::kNstat; i++)
         totstats[i] += stats[i];
      nentries += hist->GetEntries();

         //Int_t nx = hist->GetXaxis()->GetNbins();
         // loop on bins of the histogram and do the merge
      for (Int_t ibin = 0; ibin < hist->fNcells; ibin++) {
         MergeBin(hist, ibin, ibin);
      }
   }
   //copy merged stats
   fH0->PutStats(totstats);
   fH0->SetEntries(nentries);

   return kTRUE;
}


/**
   Merged histogram when axis can be different.
   Histograms are merged looking at bin center positions

 */
Bool_t TH1Merger::DifferentAxesMerge() {

   Double_t stats[TH1::kNstat], totstats[TH1::kNstat];
   for (Int_t i=0;i<TH1::kNstat;i++) {totstats[i] = stats[i] = 0;}
   fH0->GetStats(totstats);
   Double_t nentries = fH0->GetEntries();

   TIter next(&fInputList);
   while (TH1* hist=(TH1*)next()) {

      if (gDebug)
         Info("TH1Merger::DifferentAxesMerge","Merging histogram %s into %s",hist->GetName(), fH0->GetName() );

      // skip empty histograms
      if (hist->IsEmpty()) continue;

      // import statistics
      hist->GetStats(stats);
      for (Int_t i=0;i<TH1::kNstat;i++)
         totstats[i] += stats[i];
      nentries += hist->GetEntries();

      // loop on bins of the histogram and do the merge
      for (Int_t ibin = 0; ibin < hist->fNcells; ibin++) {

         // if bin is empty we can skip it
         if (IsBinEmpty(hist,ibin)) continue;

         Int_t binx,biny,binz;
         hist->GetBinXYZ(ibin, binx, biny, binz);

         // case of underflow/overflows in the histogram being merged
         if (binx <= 0 || binx >= hist->GetNbinsX() + 1) {
            if (fH0->fXaxis.CanExtend() || ( hist->fXaxis.GetBinCenter(binx) > fH0->fXaxis.GetXmin() && hist->fXaxis.GetBinCenter(binx) < fH0->fXaxis.GetXmax()) ) {
               Error("TH1Merger::DifferentAxesMerge", "Cannot merge histograms - the histograms %s can extend the X axis or have"
                     " different limits and underflows/overflows are present in the histogram %s.",fH0->GetName(),hist->GetName());
                  return kFALSE;
            }
         }
         if (biny <= 0 || biny >= hist->GetNbinsY() + 1) {
            if (fH0->fYaxis.CanExtend() || ( hist->fYaxis.GetBinCenter(biny) > fH0->fYaxis.GetXmin() && hist->fYaxis.GetBinCenter(biny) < fH0->fYaxis.GetXmax()) ) {
               Error("TH1Merger::DifferentAxesMerge", "Cannot merge histograms - the histograms %s can extend the Y axis or have"
                     " different limits and underflows/overflows are present in the histogram %s.",fH0->GetName(),hist->GetName());
                  return kFALSE;
            }
         }
         if (binz <= 0 || binz >= hist->GetNbinsZ() + 1) {
            if (fH0->fZaxis.CanExtend() || ( hist->fZaxis.GetBinCenter(binz) > fH0->fZaxis.GetXmin() && hist->fZaxis.GetBinCenter(binz) < fH0->fZaxis.GetXmax()) ) {
               Error("TH1Merger::DifferentAxesMerge", "Cannot merge histograms - the histograms %s can extend the Z axis or have"
                     " different limits and underflows/overflows are present in the histogram %s.",fH0->GetName(),hist->GetName());
                  return kFALSE;
            }
         }

         Int_t ix = 0;
         Int_t iy = 0;
         Int_t iz = 0;

         // we can extend eventually the axis if histogram is capable of doing it
         // by using FindBin
          ix = fH0->fXaxis.FindBin(hist->GetXaxis()->GetBinCenter(binx));
          if (fH0->fDimension > 1)
             iy = fH0->fYaxis.FindBin(hist->GetYaxis()->GetBinCenter(biny));
          if (fH0->fDimension > 2)
             iz = fH0->fZaxis.FindBin(hist->GetZaxis()->GetBinCenter(binz));

         Int_t ib = fH0->GetBin(ix,iy,iz);
         if (ib < 0 || ib > fH0->fNcells) {
            Fatal("TH1Merger::LabelMerge","Fatal error merging histogram %s - bin number is %d and array size is %d",
                  fH0->GetName(), ib,fH0->fNcells);
         }

         MergeBin(hist, ibin, ib);

      }
   }
   //copy merged stats
   fH0->PutStats(totstats);
   fH0->SetEntries(nentries);

   return kTRUE;
}

/**
   Find a duplicate labels in an axis label list
*/
Bool_t TH1Merger::HasDuplicateLabels(const THashList * labels) {

   if (!labels) return kFALSE;

   for (const auto * obj: *labels) {
      auto objList = labels->GetListForObject(obj);
      //objList->ls();
      if (objList->GetSize() > 1 ) {
         // check here if in the list we have duplicates
         std::unordered_set<std::string> s;
         for ( const auto * o: *objList) {
            auto ret = s.insert(std::string(o->GetName() ));
            if (!ret.second) return kTRUE;
         }
      }
   }
   return kFALSE;
}

/**
 Check if histogram has duplicate labels
 Return an integer with bit set correponding
  on the axis that has duplicate labels
  e.g. duplicate labels on x axis : return 1
       duplicate labels on x and z axis : return 5

*/
Int_t TH1Merger::CheckForDuplicateLabels(const TH1 * hist) {

   R__ASSERT(hist != nullptr);

   auto labelsX = hist->GetXaxis()->GetLabels();
   auto labelsY = hist->GetYaxis()->GetLabels();
   auto labelsZ = hist->GetZaxis()->GetLabels();

   Int_t res = 0;
   if (HasDuplicateLabels(labelsX) ) {
      Warning("TH1Merger::CheckForDuplicateLabels","Histogram %s has duplicate labels in the x axis. "
              "Bin contents will be merged in a single bin",hist->GetName());
      res |= 1;
   }
   if (HasDuplicateLabels(labelsY) ) {
      Warning("TH1Merger::CheckForDuplicateLabels","Histogram %s has duplicate labels in the y axis. "
              "Bin contents will be merged in a single bin",hist->GetName());
      res |= 2;
   }
   if (HasDuplicateLabels(labelsZ) ) {
      Warning("TH1Merger::CheckForDuplicateLabels","Histogram %s has duplicate labels in the z axis. "
              "Bin contents will be merged in a single bin",hist->GetName());
      res |= 4;
   }
   return res;
}

/**
   Merge histograms with labels
*/
Bool_t TH1Merger::LabelMerge() {

   Double_t stats[TH1::kNstat], totstats[TH1::kNstat];
   for (Int_t i=0;i<TH1::kNstat;i++) {totstats[i] = stats[i] = 0;}
   fH0->GetStats(totstats);
   Double_t nentries = fH0->GetEntries();

   // check for duplicate labels
   if (!fNoCheck && nentries > 0) CheckForDuplicateLabels(fH0);

   TIter next(&fInputList);
   while (TH1* hist=(TH1*)next()) {

      if (gDebug)
         Info("TH1Merger::LabelMerge","Merging histogram %s into %s",hist->GetName(), fH0->GetName() );

      // skip empty histograms
      if (hist->IsEmpty()) continue;

      // import statistics
      hist->GetStats(stats);
      for (Int_t i=0;i<TH1::kNstat;i++)
         totstats[i] += stats[i];
      nentries += hist->GetEntries();

      auto labelsX = hist->GetXaxis()->GetLabels();
      auto labelsY = hist->GetYaxis()->GetLabels();
      auto labelsZ = hist->GetZaxis()->GetLabels();
      R__ASSERT(!( labelsX == nullptr  && labelsY == nullptr && labelsZ == nullptr));

      Bool_t mergeLabelsX = labelsX && fH0->fXaxis.CanExtend() && hist->fXaxis.CanExtend();
      Bool_t mergeLabelsY = labelsY && fH0->fYaxis.CanExtend() && hist->fYaxis.CanExtend();
      Bool_t mergeLabelsZ = labelsZ && fH0->fZaxis.CanExtend() && hist->fZaxis.CanExtend();

      if (gDebug) {
         if (mergeLabelsX)
            Info("TH1Merger::LabelMerge","Merging X axis in label mode");
         else
            Info("TH1Merger::LabelMerge","Merging X axis in numeric mode");
         if (mergeLabelsY)
            Info("TH1Merger::LabelMerge","Merging Y axis in label mode");
         else if (hist->GetDimension() > 1)
            Info("TH1Merger::LabelMerge","Merging Y axis in numeric mode");
         if (mergeLabelsZ)
            Info("TH1Merger::LabelMerge","Merging Z axis in label mode" );
         else if (hist->GetDimension() > 2)
            Info("TH1Merger::LabelMerge","Merging Z axis in numeric mode");
      }

      // check if histogram has duplicate labels
      if (!fNoCheck && hist->GetEntries() > 0) CheckForDuplicateLabels(hist);

      // loop on bins of the histogram and do the merge
      if (gDebug) {
         // print bins original histogram
         std::cout << "Bins of original histograms\n";
         for (int ix = 1; ix <= fH0->GetXaxis()->GetNbins(); ++ix) {
            for (int iy = 1; iy <= fH0->GetYaxis()->GetNbins(); ++iy) {
               for (int iz = 1; iz <= fH0->GetZaxis()->GetNbins(); ++iz) {
                  int i = fH0->GetBin(ix,iy,iz);
                  std::cout << "bin" << ix << "," << iy << "," << iz
                     << "  : " << fH0->RetrieveBinContent(i) /* << " , " << fH0->fBinEntries.fArray[i] */ << std::endl;
               }
            }
         }
      }
      for (Int_t ibin = 0; ibin < hist->fNcells; ibin++) {

         // if bin is empty we can skip it
         if (IsBinEmpty(hist,ibin)) continue;

         Int_t binx,biny,binz;
         hist->GetBinXYZ(ibin, binx, biny, binz);

         // here only in the case of bins with labels
         const char * labelX = 0;
         const char * labelY = 0;
         const char * labelZ = 0;
         labelX=hist->GetXaxis()->GetBinLabel(binx);
         if (fH0->fDimension > 1) labelY = hist->GetYaxis()->GetBinLabel(biny);
         if (fH0->fDimension > 2) labelZ = hist->GetYaxis()->GetBinLabel(binz);
         // do we need to support case when there are bins with labels and bins without them ??
         // this case should have been detected before when examining the histograms


         Int_t ix = -1;
         Int_t iy = (fH0->fDimension > 1) ? -1 : 0;
         Int_t iz = (fH0->fDimension > 2) ? -1 : 0;

         // special case for underflow/overflows which have normally empty labels
         if (binx == 0 && TString(labelX) == "" ) ix = 0;
         if (binx == hist->fXaxis.GetNbins() +1 && TString(labelX) == "" ) ix = fH0->fXaxis.GetNbins() +1;
         if (fH0->fDimension > 1 ) {
            if (biny == 0 && TString(labelY) == "" ) iy = 0;
            if (biny == hist->fYaxis.GetNbins() +1 && TString(labelY) == "" ) iy = fH0->fYaxis.GetNbins() +1;
         }
         if (fH0->fDimension > 2 ) {
            if (binz == 0 && TString(labelZ) == "" ) iz = 0;
            if (binz == hist->fZaxis.GetNbins() +1 && TString(labelZ) == "" ) iz = fH0->fZaxis.GetNbins() +1;
         }



         // find corresponding case (in case bin is not overflow)
         // and see if for that axis we need to merge using labels or bin numbers
         if (ix == -1) {
            if (mergeLabelsX) {
               // std::cout << "find bin for label " << labelX << "  " << fH0->GetBinContent(1,1,1) << " nbins "
               // << fH0->GetXaxis()->GetNbins() << std::endl;
               ix = fH0->fXaxis.FindBin(labelX);
               // std::cout << "bin for label " << ix << "  " << fH0->GetBinContent(1,1,1) << " nbins "
               // << fH0->GetXaxis()->GetNbins() << std::endl;
            }
            else
               ix = FindFixBinNumber(binx, hist->fXaxis, fH0->fXaxis);
         }

         if (iy == -1 && fH0->fDimension> 1 ) { // check on dim should not be needed
            if (mergeLabelsY)
               iy= fH0->fYaxis.FindBin(labelY);
            else
               iy = FindFixBinNumber(biny, hist->fYaxis, fH0->fYaxis);
         }
         if (iz == -1 && fH0->fDimension> 2)  {
            if (mergeLabelsZ)
               iz= fH0->fZaxis.FindBin(labelZ);
            else
               iz = FindFixBinNumber(binz, hist->fZaxis, fH0->fZaxis);
         }

         if (gDebug)
            Info("TH1Merge::LabelMerge","Merge bin [%d,%d,%d] with label [%s,%s,%s] into bin [%d,%d,%d]",
                 binx,biny,binz,labelX,labelY,labelZ,ix,iy,iz);


         Int_t ib = fH0->GetBin(ix,iy,iz);
         if (ib < 0 || ib >= fH0->fNcells) {
            Fatal("TH1Merger::LabelMerge","Fatal error merging histogram %s - bin number is %d and array size is %d",
                  fH0->GetName(), ib,fH0->fNcells);
         }

         MergeBin(hist, ibin, ib);
      }
   }
   //copy merged stats
   fH0->PutStats(totstats);
   fH0->SetEntries(nentries);

   return kTRUE;
}

/// helper function for merging

Bool_t TH1Merger::IsBinEmpty(const TH1 * hist, Int_t ibin) {
   Double_t cu = hist->RetrieveBinContent(ibin);
   Double_t e1sq = (hist->fSumw2.fN) ?  hist->GetBinErrorSqUnchecked(ibin) : cu;
   return cu == 0 && e1sq == 0;
}

// merge input bin (ibin) of histograms hist ibin into current bin cbin of this histogram
void TH1Merger::MergeBin(const TH1 *hist, Int_t ibin, Int_t cbin)
{
   if (!fIsProfileMerge) {
      Double_t cu = hist->RetrieveBinContent(ibin);
      fH0->AddBinContent(cbin, cu);
      if (fH0->fSumw2.fN) {
         Double_t e1sq = (hist->fSumw2.fN) ? hist->GetBinErrorSqUnchecked(ibin) : cu;
         fH0->fSumw2.fArray[cbin] += e1sq;
      }
   } else {
      if (fIsProfile1D)
         MergeProfileBin(static_cast<const TProfile *> (hist), ibin, cbin);
      else if (fIsProfile2D)
         MergeProfileBin(static_cast<const TProfile2D *> (hist), ibin, cbin);
      else if (fIsProfile3D)
         MergeProfileBin(static_cast<const TProfile3D *> (hist), ibin, cbin);
   }
   return;
}

// merge profile input bin (ibin) of histograms hist ibin into current bin cbin of this histogram
template<class TProfileType>
void TH1Merger::MergeProfileBin(const TProfileType *h, Int_t hbin, Int_t pbin)
{
   TProfileType *p = static_cast<TProfileType *>(fH0);
   p->fArray[pbin] += h->fArray[hbin];
   p->fSumw2.fArray[pbin] += h->fSumw2.fArray[hbin];
   p->fBinEntries.fArray[pbin] += h->fBinEntries.fArray[hbin];
   if (p->fBinSumw2.fN) {
      if (h->fBinSumw2.fN)
         p->fBinSumw2.fArray[pbin] += h->fBinSumw2.fArray[hbin];
      else
         p->fBinSumw2.fArray[pbin] += h->fArray[hbin];
   }
   if (gDebug)
      Info("TH1Merge::MergeProfileBin", "Merge bin %d of profile %s with content %f in bin %d - result is %f", hbin,
           h->GetName(), h->fArray[hbin], pbin, p->fArray[pbin]);
}