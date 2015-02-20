// @(#)root/hist:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "Riostream.h"
#include "TAxis.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TError.h"
#include "THashList.h"
#include "TH1.h"
#include "TObjString.h"
#include "TDatime.h"
#include "TTimeStamp.h"
#include "TROOT.h"
#include "TClass.h"
#include "TMath.h"
#include <time.h>
#include <cassert>

ClassImp(TAxis)

//______________________________________________________________________________
//
// This class manages histogram axis. It is referenced by TH1 and TGraph.
// To make a graphical representation of an histogram axis, this class
// references the TGaxis class.
//
// TAxis supports axis with fixed or variable bin sizes.
// Labels may be associated to individual bins.
//
//    see examples of various axis representations drawn by class TGaxis.
//

//______________________________________________________________________________
TAxis::TAxis(): TNamed(), TAttAxis()
{
   // Default constructor.

   fNbins   = 1;
   fXmin    = 0;
   fXmax    = 1;
   fFirst   = 0;
   fLast    = 0;
   fParent  = 0;
   fLabels  = 0;
   fBits2   = 0;
   fTimeDisplay = 0;
}

//______________________________________________________________________________
TAxis::TAxis(Int_t nbins,Double_t xlow,Double_t xup): TNamed(), TAttAxis()
{
   // Axis constructor for axis with fix bin size

   fParent  = 0;
   fLabels  = 0;
   Set(nbins,xlow,xup);
}

//______________________________________________________________________________
TAxis::TAxis(Int_t nbins,const Double_t *xbins): TNamed(), TAttAxis()
{
   // Axis constructor for variable bin size

   fParent  = 0;
   fLabels  = 0;
   Set(nbins,xbins);
}

//______________________________________________________________________________
TAxis::~TAxis()
{
   // Destructor.

   if (fLabels) {
      fLabels->Delete();
      delete fLabels;
      fLabels = 0;
   }
}

//______________________________________________________________________________
TAxis::TAxis(const TAxis &axis) : TNamed(axis), TAttAxis(axis), fLabels(0)
{
   // Copy constructor.

   axis.Copy(*this);
}

//______________________________________________________________________________
TAxis& TAxis::operator=(const TAxis &orig)
{
   // Assignment operator.

   orig.Copy( *this );
   return *this;
}


//______________________________________________________________________________
const char *TAxis::ChooseTimeFormat(Double_t axislength)
{
   // Choose a reasonable time format from the coordinates in the active pad
   // and the number of divisions in this axis
   // If orientation = "X", the horizontal axis of the pad will be used for ref.
   // If orientation = "Y", the vertical axis of the pad will be used for ref.

   const char *formatstr = nullptr;
   Int_t reasformat = 0;
   Int_t ndiv,nx1,nx2,n;
   Double_t awidth;
   Double_t length;

   if (!axislength) {
      length = gPad->GetUxmax() - gPad->GetUxmin();
   } else {
      length = axislength;
   }

   ndiv = GetNdivisions();
   if (ndiv > 1000) {
      nx2   = ndiv/100;
      nx1   = TMath::Max(1, ndiv%100);
      ndiv = 100*nx2 + Int_t(Double_t(nx1)*gPad->GetAbsWNDC());
   }
   ndiv = TMath::Abs(ndiv);
   n = ndiv - (ndiv/100)*100;
   awidth = length/n;

//  width in seconds ?
   if (awidth>=.5) {
      reasformat = 1;
//  width in minutes ?
      if (awidth>=30) {
         awidth /= 60; reasformat = 2;
//   width in hours ?
         if (awidth>=30) {
            awidth /=60; reasformat = 3;
//   width in days ?
            if (awidth>=12) {
               awidth /= 24; reasformat = 4;
//   width in months ?
               if (awidth>=15.218425) {
                  awidth /= 30.43685; reasformat = 5;
//   width in years ?
                  if (awidth>=6) {
                     awidth /= 12; reasformat = 6;
                     if (awidth>=2) {
                        awidth /= 12; reasformat = 7;
                     }
                  }
               }
            }
         }
      }
   }
//   set reasonable format
   switch (reasformat) {
      case 0:
        formatstr = "%S";
        break;
      case 1:
        formatstr = "%Mm%S";
        break;
      case 2:
        formatstr = "%Hh%M";
        break;
      case 3:
        formatstr = "%d-%Hh";
        break;
      case 4:
        formatstr = "%d/%m";
        break;
      case 5:
        formatstr = "%d/%m/%y";
        break;
      case 6:
        formatstr = "%d/%m/%y";
        break;
      case 7:
        formatstr = "%m/%y";
        break;
   }
   return formatstr;
}

//______________________________________________________________________________
void TAxis::Copy(TObject &obj) const
{
   // Copy axis structure to another axis

   TNamed::Copy(obj);
   TAttAxis::Copy(((TAxis&)obj));
   TAxis &axis( ((TAxis&)obj) );
   axis.fNbins  = fNbins;
   axis.fXmin   = fXmin;
   axis.fXmax   = fXmax;
   axis.fFirst  = fFirst;
   axis.fLast   = fLast;
   axis.fBits2  = fBits2;
   fXbins.Copy(axis.fXbins);
   axis.fTimeFormat   = fTimeFormat;
   axis.fTimeDisplay  = fTimeDisplay;
   axis.fParent       = fParent;
   if (axis.fLabels) {
      axis.fLabels->Delete();
      delete axis.fLabels;
      axis.fLabels = 0;
   }
   if (fLabels) {
      //Properly handle case where not all bins have labels
      TIter next(fLabels);
      TObjString *label;
      if(! axis.fLabels) {
         axis.fLabels = new THashList(axis.fNbins, 3);
      }
      while( (label=(TObjString*)next()) ) {
         TObjString *copyLabel = new TObjString(*label);
         axis.fLabels->Add(copyLabel);
         copyLabel->SetUniqueID(label->GetUniqueID());
      }
   }
}

//______________________________________________________________________________
Int_t TAxis::DistancetoPrimitive(Int_t, Int_t)
{
   // Compute distance from point px,py to an axis

   return 9999;
}

//______________________________________________________________________________
void TAxis::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   // Execute action corresponding to one event
   //
   //  This member function is called when an axis is clicked with the locator
   //
   //  The axis range is set between the position where the mouse is pressed
   //  and the position where it is released.
   //  If the mouse position is outside the current axis range when it is released
   //  the axis is unzoomed with the corresponding proportions.
   //  Note that the mouse does not need to be in the pad or even canvas
   //  when it is released.

   if (!gPad) return;
   gPad->ExecuteEventAxis(event,px,py,this);
}

//______________________________________________________________________________
Int_t TAxis::FindBin(Double_t x)
{
   // Find bin number corresponding to abscissa x. NOTE: this method does not work with alphanumeric bins !!!
   //
   // If x is underflow or overflow, attempt to extend the axis if TAxis::kCanExtend is true.
   // Otherwise, return 0 or fNbins+1.

   Int_t bin;
   // NOTE: This should not be allowed for Alphanumeric histograms,
   // but it is heavily used (legacy) in the TTreePlayer to fill alphanumeric histograms.
   // but in case of alphanumeric do-not extend the axis. It makes no sense 
   if (IsAlphanumeric() && gDebug) Info("FindBin","Numeric query on alphanumeric axis - Sorting the bins or extending the axes / rebinning can alter the correspondence between the label and the bin interval.");
   if (x < fXmin) {              //*-* underflow
      bin = 0;
      if (fParent == 0) return bin;
      if (!CanExtend() || IsAlphanumeric() ) return bin;   
      ((TH1*)fParent)->ExtendAxis(x,this);
      return FindFixBin(x);
   } else  if ( !(x < fXmax)) {     //*-* overflow  (note the way to catch NaN)
      bin = fNbins+1;
      if (fParent == 0) return bin;
      if (!CanExtend() || IsAlphanumeric() ) return bin;
      ((TH1*)fParent)->ExtendAxis(x,this);
      return FindFixBin(x);
   } else {
      if (!fXbins.fN) {        //*-* fix bins
         bin = 1 + int (fNbins*(x-fXmin)/(fXmax-fXmin) );
      } else {                  //*-* variable bin sizes
         //for (bin =1; x >= fXbins.fArray[bin]; bin++);
         bin = 1 + TMath::BinarySearch(fXbins.fN,fXbins.fArray,x);
      }
   }
   return bin;
}

//______________________________________________________________________________
Int_t TAxis::FindBin(const char *label)
{
   // Find bin number with label.
   // If the List of labels does not exist create it and make the axis alphanumeric
   // If one wants just to add a single label- just call TAxis::SetBinLabel
   // If label is not in the list of labels do the following depending on the
   // bit TAxis::kCanExtend; of the axis.
   //   - if the bit is set add the new label and if the number of labels exceeds
   //      the number of bins, double the number of bins via TH1::LabelsInflate
   //   - if the bit is not set and the histogram has labels in each bin
   //        set the bit automatically and consider the histogram as alphanumeric
   //     if histogram has only some bins with labels then the histogram is not
   //     consider alphanumeric and return -1
   //
   // -1 is returned only when the Axis has no parent histogram

   //create list of labels if it does not exist yet
   if (!fLabels) {
      if (!fParent) return -1;
      fLabels = new THashList(fNbins,3);
      // we set the axis alphanumeric
      // when list of labels does not exist
      // do we want to do this also when histogram is not empty ?????
      if (CanBeAlphanumeric() ) { 
         SetCanExtend(kTRUE);
         SetAlphanumeric(kTRUE);
         if (fXmax <= fXmin) {
            //L.M. Dec 2010 in case of no min and max specified use 0 ->NBINS
            fXmin = 0;
            fXmax = fNbins;
         }
      }
   }

   // search for label in the existing list and return it if it exists
   TObjString *obj = (TObjString*)fLabels->FindObject(label);
   if (obj) return (Int_t)obj->GetUniqueID();

   // if labels is not in the list and we have already labels
   if (!IsAlphanumeric()) {
      // if bins without labels exist or if the axis cannot be set to alphanumeric
      if (HasBinWithoutLabel() || !CanBeAlphanumeric() ) {         
         Info("FindBin","Label %s is not in the list and the axis is not alphanumeric - ignore it",label);
         return -1;
      }
      else {
         Info("FindBin","Label %s not in the list. It will be added to the histogram",label);
         SetCanExtend(kTRUE);
         SetAlphanumeric(kTRUE);
      }
   }

   //Not yet in the list. Can we extend the axis ?
   assert ( CanExtend() && IsAlphanumeric() );
   // {
   //    if (gDebug>0)
   //       Info("FindBin","Label %s is not in the list and the axis cannot be extended - the entry will be added in the underflow bin",label);
   //    return 0;
   // }

   Int_t n = fLabels->GetEntries();

   //may be we have to resize the histogram (doubling number of channels)
   if (n >= fNbins) ((TH1*)fParent)->LabelsInflate(GetName());

   //add new label to the list: assign bin number
   obj = new TObjString(label);
   fLabels->Add(obj);
   obj->SetUniqueID(n+1);
   return n+1;
}

//______________________________________________________________________________
Int_t TAxis::FindFixBin(const char *label) const
{
   // Find bin number with label.
   // If the List of labels does not exist or the label doe not exist just return -1 .
   // Do not attempt to modify the axis. This is different than FindBin

   //create list of labels if it does not exist yet
   if (!fLabels) return -1; 
 
   // search for label in the existing list and return it if it exists
   TObjString *obj = (TObjString*)fLabels->FindObject(label);
   if (obj) return (Int_t)obj->GetUniqueID();
   return -1;
}   


//______________________________________________________________________________
Int_t TAxis::FindFixBin(Double_t x) const
{
   // Find bin number corresponding to abscissa x
   //
   // Identical to TAxis::FindBin except that if x is an underflow/overflow
   // no attempt is made to extend the axis.

   Int_t bin;
   if (x < fXmin) {              //*-* underflow
      bin = 0;
   } else  if ( !(x < fXmax)) {     //*-* overflow  (note the way to catch NaN
      bin = fNbins+1;
   } else {
      if (!fXbins.fN) {        //*-* fix bins
         bin = 1 + int (fNbins*(x-fXmin)/(fXmax-fXmin) );
      } else {                  //*-* variable bin sizes
//         for (bin =1; x >= fXbins.fArray[bin]; bin++);
         bin = 1 + TMath::BinarySearch(fXbins.fN,fXbins.fArray,x);
      }
   }
   return bin;
}

//______________________________________________________________________________
const char *TAxis::GetBinLabel(Int_t bin) const
{
   // Return label for bin

   if (!fLabels) return "";
   if (bin <= 0 || bin > fNbins) return "";
   TIter next(fLabels);
   TObjString *obj;
   while ((obj=(TObjString*)next())) {
      Int_t binid = (Int_t)obj->GetUniqueID();
      if (binid == bin) return obj->GetName();
   }
   return "";
}

//______________________________________________________________________________
Int_t TAxis::GetFirst() const
{
   //             return first bin on the axis
   //       ie 1 if no range defined
   //       NOTE: in some cases a zero is returned (see TAxis::SetRange)

   if (!TestBit(kAxisRange)) return 1;
   return fFirst;
}

//______________________________________________________________________________
Int_t TAxis::GetLast() const
{
   //             return last bin on the axis
   //       ie fNbins if no range defined
   //       NOTE: in some cases a zero is returned (see TAxis::SetRange)

   if (!TestBit(kAxisRange)) return fNbins;
   return fLast;
}

//______________________________________________________________________________
Double_t TAxis::GetBinCenter(Int_t bin) const
{
   // Return center of bin

   Double_t binwidth;
   if (!fXbins.fN || bin<1 || bin>fNbins) {
      binwidth = (fXmax - fXmin) / Double_t(fNbins);
      return fXmin + (bin-1) * binwidth + 0.5*binwidth;
   } else {
      binwidth = fXbins.fArray[bin] - fXbins.fArray[bin-1];
      return fXbins.fArray[bin-1] + 0.5*binwidth;
   }
}

//______________________________________________________________________________
Double_t TAxis::GetBinCenterLog(Int_t bin) const
{
   // Return center of bin in log
   // With a log-equidistant binning for a bin with low and up edges, the mean is :
   // 0.5*(ln low + ln up) i.e. sqrt(low*up) in logx (e.g. sqrt(10^0*10^2) = 10).
   //Imagine a bin with low=1 and up=100 :
   // - the center in lin is (100-1)/2=50.5
   // - the center in log would be sqrt(1*100)=10 (!=log(50.5))
   // NB: if the low edge of the bin is negative, the function returns the bin center
   //     as computed by TAxis::GetBinCenter

   Double_t low,up;
   if (!fXbins.fN || bin<1 || bin>fNbins) {
      Double_t binwidth = (fXmax - fXmin) / Double_t(fNbins);
      low = fXmin + (bin-1) * binwidth;
      up  = low+binwidth;
   } else {
      low = fXbins.fArray[bin-1];
      up  = fXbins.fArray[bin];
   }
   if (low <=0 ) return GetBinCenter(bin);
   return TMath::Sqrt(low*up);
}
//______________________________________________________________________________
Double_t TAxis::GetBinLowEdge(Int_t bin) const
{
   // Return low edge of bin

   if (fXbins.fN && bin > 0 && bin <=fNbins) return fXbins.fArray[bin-1];
   Double_t binwidth = (fXmax - fXmin) / Double_t(fNbins);
   return fXmin + (bin-1) * binwidth;
}

//______________________________________________________________________________
Double_t TAxis::GetBinUpEdge(Int_t bin) const
{
   // Return up edge of bin

   if (!fXbins.fN || bin < 1 || bin>fNbins) {
      Double_t binwidth = (fXmax - fXmin) / Double_t(fNbins);
      return fXmin + bin*binwidth;
   }
   return fXbins.fArray[bin];
}

//______________________________________________________________________________
Double_t TAxis::GetBinWidth(Int_t bin) const
{
   // Return bin width

   if (fNbins <= 0) return 0;
   if (fXbins.fN <= 0)  return (fXmax - fXmin) / Double_t(fNbins);
   if (bin >fNbins) bin = fNbins;
   if (bin <1 ) bin = 1;
   return fXbins.fArray[bin] - fXbins.fArray[bin-1];
}


//______________________________________________________________________________
void TAxis::GetCenter(Double_t *center) const
{
   // Return an array with the center of all bins

   Int_t bin;
   for (bin=1; bin<=fNbins; bin++) *(center + bin-1) = GetBinCenter(bin);
}

//______________________________________________________________________________
void TAxis::GetLowEdge(Double_t *edge) const
{
   // Return an array with the lod edge of all bins

   Int_t bin;
   for (bin=1; bin<=fNbins; bin++) *(edge + bin-1) = GetBinLowEdge(bin);
}

//______________________________________________________________________________
const char *TAxis::GetTimeFormatOnly() const
{
   // Return *only* the time format from the string fTimeFormat

   static TString timeformat;
   Int_t idF = fTimeFormat.Index("%F");
   if (idF>=0) {
      timeformat = fTimeFormat(0,idF);
   } else {
      timeformat = fTimeFormat;
   }
   return timeformat.Data();
}

//______________________________________________________________________________
const char *TAxis::GetTicks() const
{
   // Return the ticks option (see SetTicks)

   if (TestBit(kTickPlus) && TestBit(kTickMinus)) return "+-";
   if (TestBit(kTickMinus)) return "-";
   return "+";
}

//______________________________________________________________________________
Bool_t TAxis::HasBinWithoutLabel() const
{
   // this helper function checks if there is a bin without a label
   // if all bins have labels, the axis can / will become alphanumeric
   return fLabels->GetSize() != fNbins;
}

//______________________________________________________________________________
void TAxis::LabelsOption(Option_t *option)
{
   //  Set option(s) to draw axis with labels
   //  option = "a" sort by alphabetic order
   //         = ">" sort by decreasing values
   //         = "<" sort by increasing values
   //         = "h" draw labels horizonthal
   //         = "v" draw labels vertical
   //         = "u" draw labels up (end of label right adjusted)
   //         = "d" draw labels down (start of label left adjusted)

   if (!fLabels) {
      Warning("Sort","Cannot sort. No labels");
      return;
   }
   TH1 *h = (TH1*)GetParent();
   if (!h) {
      Error("Sort","Axis has no parent");
      return;
   }

   h->LabelsOption(option,GetName());
}

//______________________________________________________________________________
void TAxis::ImportAttributes(const TAxis *axis)
{
   // Copy axis attributes to this

   SetTitle(axis->GetTitle());
   SetNdivisions(axis->GetNdivisions());
   SetAxisColor(axis->GetAxisColor());
   SetLabelColor(axis->GetLabelColor());
   SetLabelFont(axis->GetLabelFont());
   SetLabelOffset(axis->GetLabelOffset());
   SetLabelSize(axis->GetLabelSize());
   SetTickLength(axis->GetTickLength());
   SetTitleOffset(axis->GetTitleOffset());
   SetTitleSize(axis->GetTitleSize());
   SetTitleColor(axis->GetTitleColor());
   SetTitleFont(axis->GetTitleFont());
   SetBit(TAxis::kCenterTitle,   axis->TestBit(TAxis::kCenterTitle));
   SetBit(TAxis::kCenterLabels,  axis->TestBit(TAxis::kCenterLabels));
   SetBit(TAxis::kRotateTitle,   axis->TestBit(TAxis::kRotateTitle));
   SetBit(TAxis::kNoExponent,    axis->TestBit(TAxis::kNoExponent));
   SetBit(TAxis::kTickPlus,      axis->TestBit(TAxis::kTickPlus));
   SetBit(TAxis::kTickMinus,     axis->TestBit(TAxis::kTickMinus));
   SetBit(TAxis::kMoreLogLabels, axis->TestBit(TAxis::kMoreLogLabels));
   SetBit(TAxis::kDecimals,      axis->TestBit(TAxis::kDecimals));
   SetTimeFormat(axis->GetTimeFormat());
}



//______________________________________________________________________________
void TAxis::SaveAttributes(std::ostream &out, const char *name, const char *subname)
{
    // Save axis attributes as C++ statement(s) on output stream out

   char quote = '"';
   if (strlen(GetTitle())) {
      TString t(GetTitle());
      t.ReplaceAll("\\","\\\\");
      out<<"   "<<name<<subname<<"->SetTitle("<<quote<<t.Data()<<quote<<");"<<std::endl;
   }
   if (fTimeDisplay) {
      out<<"   "<<name<<subname<<"->SetTimeDisplay(1);"<<std::endl;
      out<<"   "<<name<<subname<<"->SetTimeFormat("<<quote<<GetTimeFormat()<<quote<<");"<<std::endl;
   }
   if (fLabels) {
      TIter next(fLabels);
      TObjString *obj;
      while ((obj=(TObjString*)next())) {
         out<<"   "<<name<<subname<<"->SetBinLabel("<<obj->GetUniqueID()<<","<<quote<<obj->GetName()<<quote<<");"<<std::endl;
      }
   }

   if (fFirst || fLast) {
      out<<"   "<<name<<subname<<"->SetRange("<<fFirst<<","<<fLast<<");"<<std::endl;
   }

   if (TestBit(kLabelsHori)) {
      out<<"   "<<name<<subname<<"->SetBit(TAxis::kLabelsHori);"<<std::endl;
   }

   if (TestBit(kLabelsVert)) {
      out<<"   "<<name<<subname<<"->SetBit(TAxis::kLabelsVert);"<<std::endl;
   }

   if (TestBit(kLabelsDown)) {
      out<<"   "<<name<<subname<<"->SetBit(TAxis::kLabelsDown);"<<std::endl;
   }

   if (TestBit(kLabelsUp)) {
      out<<"   "<<name<<subname<<"->SetBit(TAxis::kLabelsUp);"<<std::endl;
   }

   if (TestBit(kCenterLabels)) {
      out<<"   "<<name<<subname<<"->CenterLabels(true);"<<std::endl;
   }

   if (TestBit(kCenterTitle)) {
      out<<"   "<<name<<subname<<"->CenterTitle(true);"<<std::endl;
   }

   if (TestBit(kRotateTitle)) {
      out<<"   "<<name<<subname<<"->RotateTitle(true);"<<std::endl;
   }

   if (TestBit(kDecimals)) {
      out<<"   "<<name<<subname<<"->SetDecimals();"<<std::endl;
   }

   if (TestBit(kMoreLogLabels)) {
      out<<"   "<<name<<subname<<"->SetMoreLogLabels();"<<std::endl;
   }

   if (TestBit(kNoExponent)) {
      out<<"   "<<name<<subname<<"->SetNoExponent();"<<std::endl;
   }

   TAttAxis::SaveAttributes(out,name,subname);
}

//______________________________________________________________________________
void TAxis::Set(Int_t nbins, Double_t xlow, Double_t xup)
{
   // Initialize axis with fix bins

   fNbins   = nbins;
   fXmin    = xlow;
   fXmax    = xup;
   if (!fParent) SetDefaults();
   if (fXbins.fN > 0) fXbins.Set(0);
}

//______________________________________________________________________________
void TAxis::Set(Int_t nbins, const Float_t *xbins)
{
   // Initialize axis with variable bins

   Int_t bin;
   fNbins  = nbins;
   fXbins.Set(fNbins+1);
   for (bin=0; bin<= fNbins; bin++)
      fXbins.fArray[bin] = xbins[bin];
   for (bin=1; bin<= fNbins; bin++)
      if (fXbins.fArray[bin] < fXbins.fArray[bin-1])
         Error("TAxis::Set", "bins must be in increasing order");
   fXmin      = fXbins.fArray[0];
   fXmax      = fXbins.fArray[fNbins];
   if (!fParent) SetDefaults();
}

//______________________________________________________________________________
void TAxis::Set(Int_t nbins, const Double_t *xbins)
{
   // Initialize axis with variable bins

   Int_t bin;
   fNbins  = nbins;
   fXbins.Set(fNbins+1);
   for (bin=0; bin<= fNbins; bin++)
      fXbins.fArray[bin] = xbins[bin];
   for (bin=1; bin<= fNbins; bin++)
      if (fXbins.fArray[bin] < fXbins.fArray[bin-1])
         Error("TAxis::Set", "bins must be in increasing order");
   fXmin      = fXbins.fArray[0];
   fXmax      = fXbins.fArray[fNbins];
   if (!fParent) SetDefaults();
}

//______________________________________________________________________________
void TAxis::SetAlphanumeric(Bool_t alphanumeric)
{
   if (alphanumeric) fBits2 |= kAlphanumeric;
   else fBits2 &= ~kAlphanumeric;

   // clear underflow and overflow (in an alphanumeric situation they do not make sense)
   // NOTE: using AddBinContent instead of SetBinContent in order to not change
   //  the number of entries
   //((TH1 *)fParent)->ClearUnderflowAndOverflow();
   // L.M. 26.1.15 Keep underflow and overflows (see ROOT-7034)
   if (gDebug && fParent) {
      TH1 * h = dynamic_cast<TH1*>( fParent);
      if (!h) return;
      double s[TH1::kNstat];
      h->GetStats(s);
      if (s[0] != 0. && gDebug > 0)
         Info("SetAlphanumeric","Histogram %s is set alphanumeric but has non-zero content",GetName());
   }
}


//______________________________________________________________________________
void TAxis::SetDefaults()
{
   // Set axis default values (from TStyle)

   fFirst   = 0;
   fLast    = 0;
   fBits2   = 0;
   char name[2];
   strlcpy(name,GetName(),2);
   name[1] = 0;
   TAttAxis::ResetAttAxis(name);
   fTimeDisplay = 0;
   SetTimeFormat();
}

//______________________________________________________________________________
void TAxis::SetBinLabel(Int_t bin, const char *label)
{
   // Set label for bin
   // If no label list exists, it is created. If all the bins have labels, the
   // axis becomes alphanumeric and extendable.
   // New labels will not be added with the Fill method but will end-up in the
   // underflow bin. See documentation of TAxis::FindBin(const char*)

   if (!fLabels) fLabels = new THashList(fNbins,3);

   if (bin <= 0 || bin > fNbins) {
      Error("SetBinLabel","Illegal bin number: %d",bin);
      return;
   }

   // Check whether this bin already has a label.
   TIter next(fLabels);
   TObjString *obj;
   while ((obj=(TObjString*)next())) {
      if ( obj->GetUniqueID()==(UInt_t)bin ) {
         // It does. Overwrite it.
         obj->SetString(label);
         // LM need to rehash the labels list (see ROOT-5025)
         fLabels->Rehash(fLabels->GetSize() );
         return;
      }
   }
   // It doesn't. Add this new label.
   obj = new TObjString(label);
   fLabels->Add(obj);
   obj->SetUniqueID((UInt_t)bin);

   // check for Alphanumeric case (labels for each bin)
   if (CanBeAlphanumeric() && fLabels->GetSize() == fNbins) {      
      SetAlphanumeric(kTRUE);
      SetCanExtend(kTRUE);
   }
}


//______________________________________________________________________________
void TAxis::SetRange(Int_t first, Int_t last)
{
   //  Set the viewing range for the axis from bin first to last
   //  To set a range using the axis coordinates, use TAxis::SetRangeUser.

   //  If first == last == 0 or if last < first or if the range specified does
   //  not intersect at all with the maximum available range [0, fNbins + 1],
   //  then the range is reset by removing the bit TAxis::kAxisRange. In this
   //  case the functions TAxis::GetFirst() and TAxis::GetLast() will return 1
   //  and fNbins.

   //  If the range specified partially intersects [0, fNbins + 1], then the
   //  intersection range is set. For instance, if first == -2 and last == fNbins,
   //  then the set range is [0, fNbins] (fFirst = 0 and fLast = fNbins).
   //
   //  NOTE: for historical reasons, SetRange(0,0) resets the range even though Bin 0 is
   //       technically reserved for the underflow; in order to set the range of the axis
   //       so that it only includes the underflow, use SetRange(a,0), where a < 0

   Int_t nCells = fNbins + 1; // bins + overflow

   // special reset range cases
   if (last < first || (first < 0 && last < 0) ||
         (first > nCells && last > nCells) || (first == 0 && last == 0)
   ) {
      fFirst = 1;
      fLast = fNbins;
      SetBit(kAxisRange, 0);
   } else {
      fFirst = std::max(first, 0);
      fLast = std::min(last, nCells);
      SetBit(kAxisRange, 1);
   }

}


//______________________________________________________________________________
void TAxis::SetRangeUser(Double_t ufirst, Double_t ulast)
{
   //  Set the viewing range for the axis from ufirst to ulast (in user coordinates)
   //  To set a range using the axis bin numbers, use TAxis::SetRange.

   if (!strstr(GetName(),"xaxis")) {
      TH1 *hobj = (TH1*)GetParent();
      if (hobj &&
          ((hobj->GetDimension() == 2 && strstr(GetName(),"zaxis"))
           || (hobj->GetDimension() == 1 && strstr(GetName(),"yaxis")))) {
         hobj->SetMinimum(ufirst);
         hobj->SetMaximum(ulast);
         return;
      }
   }
   Int_t ifirst = FindFixBin(ufirst);
   Int_t ilast = FindFixBin(ulast);
   // fixes for numerical error and for https://savannah.cern.ch/bugs/index.php?99777
   if (GetBinUpEdge(ifirst) <= ufirst ) ifirst += 1;
   if (GetBinLowEdge(ilast) >= ulast ) ilast -= 1;
   SetRange(ifirst, ilast);
}

//______________________________________________________________________________
void TAxis::SetTicks(Option_t *option)
{
   //  set ticks orientation
   //  option = "+"  ticks drawn on the "positive side" (default)
   //  option = "-"  ticks drawn on the "negative side"
   //  option = "+-" ticks drawn on both sides

   ResetBit(kTickPlus);
   ResetBit(kTickMinus);
   if (strchr(option,'+')) SetBit(kTickPlus);
   if (strchr(option,'-')) SetBit(kTickMinus);
}

//______________________________________________________________________________
void TAxis::SetTimeFormat(const char *tformat)
{
   // Change the format used for time plotting
   //
   //  The format string for date and time use the same options as the one used
   //  in the standard strftime C function, i.e. :
   //    for date :
   //      %a abbreviated weekday name
   //      %b abbreviated month name
   //      %d day of the month (01-31)
   //      %m month (01-12)
   //      %y year without century
   //
   //    for time :
   //      %H hour (24-hour clock)
   //      %I hour (12-hour clock)
   //      %p local equivalent of AM or PM
   //      %M minute (00-59)
   //      %S seconds (00-61)
   //      %% %
   //
   //    This function allows also to define the time offset. It is done via %F
   //    which should be appended at the end of the format string. The time
   //    offset has the following format: 'yyyy-mm-dd hh:mm:ss'
   //    Example:
   //
   //          h = new TH1F("Test","h",3000,0.,200000.);
   //          h->GetXaxis()->SetTimeDisplay(1);
   //          h->GetXaxis()->SetTimeFormat("%d\/%m\/%y%F2000-02-28 13:00:01");
   //
   //    This defines the time format being "dd/mm/yy" and the time offset as the
   //    February 28th 2003 at 13:00:01
   //
   //    If %F is not specified, the time offset used will be the one defined by:
   //    gStyle->SetTimeOffset. For example like that:
   //
   //          TDatime da(2003,02,28,12,00,00);
   //          gStyle->SetTimeOffset(da.Convert());

   TString timeformat = tformat;

   if (timeformat.Index("%F")>=0 || timeformat.IsNull()) {
      fTimeFormat = timeformat;
      return;
   }

   Int_t idF = fTimeFormat.Index("%F");
   if (idF>=0) {
      Int_t lnF = fTimeFormat.Length();
      TString stringtimeoffset = fTimeFormat(idF,lnF);
      fTimeFormat = tformat;
      fTimeFormat.Append(stringtimeoffset);
   } else {
      fTimeFormat = tformat;
      SetTimeOffset(gStyle->GetTimeOffset());
   }
}


//______________________________________________________________________________
void TAxis::SetTimeOffset(Double_t toffset, Option_t *option)
{
   // Change the time offset
   // If option = "gmt", set display mode to GMT.

   TString opt = option;
   opt.ToLower();

   char tmp[20];
   time_t timeoff;
   struct tm* utctis;
   Int_t idF = fTimeFormat.Index("%F");
   if (idF>=0) fTimeFormat.Remove(idF);
   fTimeFormat.Append("%F");

   timeoff = (time_t)((Long_t)(toffset));
   // offset is always saved in GMT to allow file transport
   // to different time zones
   utctis = gmtime(&timeoff);

   strftime(tmp,20,"%Y-%m-%d %H:%M:%S",utctis);
   fTimeFormat.Append(tmp);

   // append the decimal part of the time offset
   Double_t ds = toffset-(Int_t)toffset;
   snprintf(tmp,20,"s%g",ds);
   fTimeFormat.Append(tmp);

   // add GMT/local option
   if (opt.Contains("gmt")) fTimeFormat.Append(" GMT");
}


//______________________________________________________________________________
void TAxis::Streamer(TBuffer &R__b)
{
   // Stream an object of class TAxis.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 5) {
         R__b.ReadClassBuffer(TAxis::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      TAttAxis::Streamer(R__b);
      R__b >> fNbins;
      if (R__v < 5) {
         Float_t xmin,xmax;
         R__b >> xmin; fXmin = xmin;
         R__b >> xmax; fXmax = xmax;
         Float_t *xbins = 0;
         Int_t n = R__b.ReadArray(xbins);
         fXbins.Set(n);
         for (Int_t i=0;i<n;i++) fXbins.fArray[i] = xbins[i];
         delete [] xbins;
      } else {
         R__b >> fXmin;
         R__b >> fXmax;
         fXbins.Streamer(R__b);
      }
      if (R__v > 2) {
         R__b >> fFirst;
         R__b >> fLast;
          // following lines required to repair for a bug in Root version 1.03
         if (fFirst < 0 || fFirst > fNbins) fFirst = 0;
         if (fLast  < 0 || fLast  > fNbins) fLast  = 0;
         if (fLast  < fFirst) { fFirst = 0; fLast = 0;}
         if (fFirst ==0 && fLast == 0) SetBit(kAxisRange,0);
      }
      if (R__v > 3) {
         R__b >> fTimeDisplay;
         fTimeFormat.Streamer(R__b);
      } else {
         SetTimeFormat();
      }
      R__b.CheckByteCount(R__s, R__c, TAxis::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TAxis::Class(),this);
   }
}

//______________________________________________________________________________
void TAxis::UnZoom()
{
   // Reset first & last bin to the full range

   if (!gPad) return;
   gPad->SetView();

   //unzoom object owning this axis
   SetRange(0,0);
   TH1 *hobj1 = (TH1*)GetParent();
   if (!strstr(GetName(),"xaxis")) {
      if (!hobj1) return;
      if (hobj1->GetDimension() == 2) {
         if (strstr(GetName(),"zaxis")) {
            hobj1->SetMinimum();
            hobj1->SetMaximum();
            hobj1->ResetBit(TH1::kIsZoomed);
         }
         return;
      }
      if (strcmp(hobj1->GetName(),"hframe") == 0 ) {
         hobj1->SetMinimum(fXmin);
         hobj1->SetMaximum(fXmax);
      } else {
         if (fXmin==hobj1->GetMinimum() && fXmax==hobj1->GetMaximum()) {
            hobj1->SetMinimum(fXmin);
            hobj1->SetMaximum(fXmax);
         } else {
            hobj1->SetMinimum();
            hobj1->SetMaximum();
         }
         hobj1->ResetBit(TH1::kIsZoomed);
      }
   }
   //must unzoom all histograms in the pad
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   while ((obj= next())) {
      if (!obj->InheritsFrom(TH1::Class())) continue;
      TH1 *hobj = (TH1*)obj;
      if (hobj == hobj1) continue;
      if (!strstr(GetName(),"xaxis")) {
         if (hobj->GetDimension() == 2) {
            if (strstr(GetName(),"zaxis")) {
               hobj->SetMinimum();
               hobj->SetMaximum();
               hobj->ResetBit(TH1::kIsZoomed);
            } else {
               hobj->GetYaxis()->SetRange(0,0);
            }
            return;
         }
         if (strcmp(hobj->GetName(),"hframe") == 0 ) {
            hobj->SetMinimum(fXmin);
            hobj->SetMaximum(fXmax);
         } else {
            hobj->SetMinimum();
            hobj->SetMaximum();
            hobj->ResetBit(TH1::kIsZoomed);
         }
      } else {
         hobj->GetXaxis()->SetRange(0,0);
      }
   }
}

//______________________________________________________________________________
void TAxis::ZoomOut(Double_t factor, Double_t offset)
{
   // Zoom out by a factor of 'factor' (default =2)
   //   uses previous zoom factor by default
   // Keep center defined by 'offset' fixed
   //   ie. -1 at left of current range, 0 in center, +1 at right


   if (factor <= 0) factor = 2;
   Double_t center = (GetFirst()*(1-offset) + GetLast()*(1+offset))/2.;
   Int_t first = int(TMath::Floor(center+(GetFirst()-center)*factor + 0.4999999));
   Int_t last  = int(TMath::Floor(center+(GetLast() -center)*factor + 0.5000001));
   if (first==GetFirst() && last==GetLast()) { first--; last++; }
   SetRange(first,last);
}
