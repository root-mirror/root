// @(#)root/treeviewer:$Id$
// Author: Bastien Dalla Piazza  02/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TParallelCoordVar.h"
#include "TParallelCoord.h"
#include "TParallelCoordRange.h"

#include <iostream>

#include "TLatex.h"
#include "TLine.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TH1.h"
#include "TStyle.h"
#include "TBox.h"
#include "TCollection.h"
#include "TList.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "snprintf.h"

ClassImp(TParallelCoordVar);

/** \class TParallelCoordVar

TParallelCoord axes. Class containing a variable for the TParallelCoord.

Options can be defined each axis separately using the right mouse click. These
options can be applied to every axes using the editor.

  - Axis width: If set to 0, the axis is simply a line. If higher, a color
    histogram is drawn on the axis.
  - Axis histogram height: If not 0, a usual bar histogram is drawn on the plot.

The order in which the variables are drawn is essential to see the clusters. The
axes can be dragged to change their position. A zoom is also available. The
logarithm scale is also available by right clicking on the axis.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TParallelCoordVar::TParallelCoordVar()
   :TNamed(), TAttLine(), TAttFill(kOrange+9,0)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
///Destructor.

TParallelCoordVar::~TParallelCoordVar()
{
   if (fHistogram) delete fHistogram;
   if (fRanges){
      TIter next(fRanges);
      TParallelCoordRange* range;
      while((range = (TParallelCoordRange*)next())) fParallel->CleanUpSelections(range);
      fRanges->Delete();
      delete fRanges;
   }
   if (fVal) delete [] fVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor. By default, the title and the name are the expression
/// given to TTree::Draw. The name can be changed by the user (the label on the
/// plot) but not the title.

TParallelCoordVar::TParallelCoordVar(Double_t *val, const char* title, Int_t id, TParallelCoord* parallel)
   :TNamed(title,title), TAttLine(1,1,1), TAttFill(kOrange+9,3001)
{
   Init();
   fId = id;
   fParallel      = parallel;
   fRanges        = new TList();
   fNentries       = fParallel->GetNentries();

   fVal = new Double_t[fParallel->GetNentries()];

   for(Long64_t ui = 0;ui<fParallel->GetNentries();++ui) fVal[ui]=val[ui];

   GetMinMaxMean();
   GetHistogram();
   GetQuantiles();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a range to the current selection on the axis.

void TParallelCoordVar::AddRange(TParallelCoordRange* range)
{
   if (!range) {
      TParallelCoordSelect *select = fParallel->GetCurrentSelection();
      if (select) {
         range = new TParallelCoordRange(this,0,0,select);
         fRanges->Add(range);
         range->GetSelection()->Add(range);
      } else {
         Error("AddRange","You must create a selection before adding ranges.");
      }
   } else {
      fRanges->Add(range);
      range->GetSelection()->Add(range);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete variables.

void TParallelCoordVar::DeleteVariable()
{
   fParallel->RemoveVariable(this);
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the distance from the axis.

Int_t TParallelCoordVar::DistancetoPrimitive(Int_t px, Int_t py)
{
   if(!gPad) return 9999;
   Double_t xx = gPad->AbsPixeltoX(px);
   Double_t yy = gPad->AbsPixeltoY(py);
   Double_t dist = 9999;

   if (fX1==fX2) {
      if (yy>fY1+0.01 && yy<fY2-0.01) dist = TMath::Abs(xx - fX1);
   } else {
      if (xx>fX1+0.01 && xx<fX2-0.01) dist = TMath::Abs(yy - fY1);
   }
   if (dist<=0.005) return 0;
   else             return 9999;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the axis.

void TParallelCoordVar::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the entry is within the range(s) of "select". Two ranges on a
/// single axis are conjugated as a "or": to be selected, the entry must be in
/// one of the ranges.

Bool_t TParallelCoordVar::Eval(Long64_t evtidx, TParallelCoordSelect *select)
{
   if (fRanges->GetSize() > 0){
      TIter next(fRanges);
      Bool_t inarange = kFALSE;
      Bool_t noOwnedRange = kTRUE;
      TParallelCoordRange *range;
      while ((range = (TParallelCoordRange*)next())){
         if(select->Contains(range)) {
            noOwnedRange = kFALSE;
            if(range->IsIn(fVal[evtidx])) inarange = kTRUE;
         }
      }
      if (noOwnedRange) return kTRUE;
      else return inarange;
   }
   else return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the corresponding entry.

void TParallelCoordVar::ExecuteEvent(Int_t entry, Int_t px, Int_t py)
{
   if (!gPad) return;
   if (!gPad->IsEditable() && entry!=kMouseEnter) return;

   static Int_t pxold, pyold;
   static Int_t zoom; // -1:nothing zoomed or translated, 0:translating the axis, 1:zooming
   static Int_t pzoomold;
   static Bool_t first = kTRUE;

   Int_t px1,px2,py1,py2,n=-1;
   px1 = gPad->XtoAbsPixel(fX1);
   px2 = gPad->XtoAbsPixel(fX2);
   py1 = gPad->YtoAbsPixel(fY1);
   py2 = gPad->YtoAbsPixel(fY2);

   if(fX1 == fX2) {
      if(gPad->AbsPixeltoX(px)-fX1 > 0) gPad->SetCursor(kArrowVer);
      else                              gPad->SetCursor(kArrowHor);
   } else {
      if(gPad->AbsPixeltoY(py)-fY1 > 0) gPad->SetCursor(kArrowHor);
      else                              gPad->SetCursor(kArrowVer);
   }

   gVirtualX->SetLineColor(-1);
   switch (entry) {
      case kButton1Down:
         if (fX1==fX2){
            ((TCanvas*)gPad)->Selected(gPad,fParallel,1);
            if(gPad->AbsPixeltoX(px)-fX1 > 0){
               zoom = 1;
               gVirtualX->DrawLine(gPad->XtoAbsPixel(fX1-0.05),py,gPad->XtoAbsPixel(fX1+0.05),py);
               first = kTRUE;
               pzoomold = py;
            } else {
               zoom = 0;
               gVirtualX->DrawLine(px,py1,px,py2);
            }
         } else {
            if(gPad->AbsPixeltoY(py)-fY1 > 0){
               zoom = 1;
               gVirtualX->DrawLine(px,gPad->YtoAbsPixel(fY1-0.05),px,gPad->YtoAbsPixel(fY1+0.05));
               first=kTRUE;
               pzoomold = px;
            } else {
               zoom = 0;
               gVirtualX->DrawLine(px1,py,px2,py);
            }
         }
         pxold = px;
         pyold = py;
         break;
      case kButton1Up: {
         Double_t xx = gPad->AbsPixeltoX(px);
         Double_t yy = gPad->AbsPixeltoY(py);
         TFrame *frame = gPad->GetFrame();
         if (fX1==fX2) {
            if(zoom == 0){
               Double_t axisSpace = (frame->GetX2() - frame->GetX1())/(fParallel->GetNvar() - 1);
               Double_t pos = (xx - frame->GetX1())/axisSpace;
               if (pos < 0) n = -1;
               else         n = (Int_t)pos;
            } else {
               Double_t min = GetValuefromXY(xx,yy);
               Double_t max = GetValuefromXY(xx,gPad->AbsPixeltoY(pzoomold));
               if(TMath::Abs(min-max) < 0.00001) return;       // Avoid zooming if the axis is just clicked.
               if (fParallel->TestBit(TParallelCoord::kGlobalScale)) {
                  if (min>max) {
                     Double_t mem = min;
                     min = max; max = mem;
                  }
                  fParallel->SetGlobalMin(min);
                  fParallel->SetGlobalMax(max);
               } else {
                  SetCurrentLimits(min,max);
               }
            }
         } else {
            if(zoom == 0) {
               Double_t axisSpace = (frame->GetY2() - frame->GetY1())/(fParallel->GetNvar() - 1);
               Double_t pos = (yy-frame->GetY1())/axisSpace;
               if (pos < 0) n= -1;
               else         n = (Int_t)pos;
            } else {
               Double_t min = GetValuefromXY(xx,yy);
               Double_t max = GetValuefromXY(gPad->AbsPixeltoX(pzoomold),yy);
               SetCurrentLimits(min,max);
            }
         }
         if(zoom == 0){
            if (n>=0 && (UInt_t)n>=fParallel->GetNvar()) --n;
            else if (n<fParallel->GetVarList()->IndexOf(this)) ++n;
            fParallel->GetVarList()->Remove(this);
            fParallel->GetVarList()->AddAt(this,n);
         }
         gPad->Modified();
         break;
      }
      case kMouseMotion:
         pxold=px;
         pyold=py;
         break;
      case kButton1Motion:
         if(fX1==fX2){
            if(zoom==0){
               gPad->SetCursor(kArrowHor);
               gVirtualX->DrawLine(pxold,py1,pxold,py2);
               gVirtualX->DrawLine(px,py1,px,py2);
            } else if(zoom==1) {
               gPad->SetCursor(kArrowVer);
               if(!first) gVirtualX->DrawLine(gPad->XtoAbsPixel(fX1-0.05),pyold,gPad->XtoAbsPixel(fX1+0.05),pyold);
               gVirtualX->DrawLine(gPad->XtoAbsPixel(fX1-0.05),py,gPad->XtoAbsPixel(fX1+0.05),py);
               first = kFALSE;
            }
         } else {
            if(zoom==0){
               gPad->SetCursor(kArrowVer);
               gVirtualX->DrawLine(px1,pyold,px2,pyold);
               gVirtualX->DrawLine(px1,py,px2,py);
            } else if(zoom==1){
               gPad->SetCursor(kArrowHor);
               if(!first) gVirtualX->DrawLine(pxold,gPad->YtoAbsPixel(fY1-0.05),pxold,gPad->YtoAbsPixel(fY1+0.05));
               gVirtualX->DrawLine(px,gPad->YtoAbsPixel(fY1-0.05),px,gPad->YtoAbsPixel(fY1+0.05));
               first = kFALSE;
            }
         }
         pxold = px;
         pyold = py;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the position of the variable on the graph for the n'th entry.

void TParallelCoordVar::GetEntryXY(Long64_t n, Double_t & x, Double_t & y)
{
   if(fX1==fX2){
      x = fX1;
      if (fMinCurrent != fMaxCurrent) {
         if (TestBit(kLogScale)) y = fY1 + (fY2 - fY1) *
                                    (TMath::Log10(fVal[n]/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
         else                    y = fY1 + (fY2 - fY1) *
                                    (fVal[n] - fMinCurrent) / (fMaxCurrent - fMinCurrent);
      } else {
         y = fY1 + 0.5*(fY2-fY1);
      }
   } else {
      y = fY1;
      if (fMinCurrent != fMaxCurrent) {
         if (TestBit(kLogScale)) x = fX1 + (fX2 - fX1) *
                                    (TMath::Log10(fVal[n]/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
         else                    x = fX1 + (fX2 - fX1) *
                                    (fVal[n] - fMinCurrent) / (fMaxCurrent - fMinCurrent);
      } else {
         x = fX1 + 0.5*(fX2-fX1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the entry weight: The weight of an entry for a given variable
/// is the bin content of the histogram bin the entry is going through.

Int_t TParallelCoordVar::GetEntryWeight(Long64_t evtidx)
{
   Int_t bin = 1 + (Int_t)((fVal[evtidx] - fMinCurrent)/((fMaxCurrent-fMinCurrent)/fNbins));
   return (Int_t)fHistogram->GetBinContent(bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Create or recreate the histogram.

TH1F* TParallelCoordVar::GetHistogram()
{
   if (fHistogram) delete fHistogram;
   fHistogram = NULL;
   fHistogram = new TH1F("hpa", "hpa", fNbins, fMinCurrent, fMaxCurrent+0.0001*(fMaxCurrent-fMinCurrent));
   fHistogram->SetDirectory(0);
   Long64_t first = fParallel->GetCurrentFirst();
   Long64_t nentries = fParallel->GetCurrentN();
   for(Long64_t li=first; li<first+nentries;++li) {
      if(fVal[li] >= fMinCurrent && fVal[li] <= fMaxCurrent) fHistogram->Fill(fVal[li]);
   }
   return fHistogram;
}

////////////////////////////////////////////////////////////////////////////////
/// Get mean, min and max of those variable.

void TParallelCoordVar::GetMinMaxMean()
{
   Double_t min,max,ave = 0;
   min = DBL_MAX;
   max = -DBL_MAX;
   Long64_t first,nentries;
   first = fParallel->GetCurrentFirst();
   nentries = fParallel->GetCurrentN();
   for(Long64_t li=first; li<first+nentries;++li){
      if(fVal[li]<min) min = fVal[li];
      if(fVal[li]>max) max = fVal[li];
      ave+=fVal[li];
   }

   fMean = ave/((Double_t)nentries);
   fMinCurrent = fMinInit = min;
   fMaxCurrent = fMaxInit = max;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns info about this axis.

char* TParallelCoordVar::GetObjectInfo(Int_t px, Int_t py) const
{
   static char info[128];
   info[0] = 0;

   if (!gPad) return info;
   Double_t xx = gPad->AbsPixeltoX(px);
   Double_t yy = gPad->AbsPixeltoY(py);
   if (fX1 == fX2) {
      if (yy<fY1) {
         snprintf(info,128,"%s = %f", GetTitle(), fMinCurrent);
      } else if (yy>fY2) {
         snprintf(info,128,"%s = %f", GetTitle(), fMaxCurrent);
      } else {
         Double_t axislength = fY2-fY1;
         Double_t pos = (yy-fY1)/axislength;
         snprintf(info,128,"%s = %f", GetTitle(), fMinCurrent + pos*(fMaxCurrent-fMinCurrent));
      }
   } else {
      if (xx<fX1) {
         snprintf(info,128,"%s = %f", GetTitle(), fMinCurrent);
      } else if(xx>fX2) {
         snprintf(info,128,"%s = %f", GetTitle(), fMaxCurrent);
      } else {
         Double_t axislength = fX2-fX1;
         Double_t pos = (xx-fX1)/axislength;
         snprintf(info,128,"%s = %f", GetTitle(), pos*(fMaxCurrent-fMinCurrent));
      }
   }
   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the box plot values (quantiles).

void TParallelCoordVar::GetQuantiles()
{
   Double_t *quantiles = new Double_t[3];
   quantiles[0]=0.; quantiles[1]=0.; quantiles[2] = 0.;
   Double_t *prob = new Double_t[3];
   prob[0]=0.25; prob[1]=0.5; prob[2] = 0.75;
   Long64_t first = fParallel->GetCurrentFirst();
   Long64_t nentries = fParallel->GetCurrentN();
   if (!TestBit(kLogScale) && first==0 && nentries==fNentries) TMath::Quantiles(fNentries,3,fVal,quantiles,prob,kFALSE);
   else {
      Double_t* val = new Double_t[nentries];
      Int_t selected = 0;
      if(fMinInit<=0) {
         for (Long64_t n=first;n<first+nentries;++n) {
            if (fVal[n] >= fMinCurrent) {
               if (TestBit(kLogScale)) val[selected] = TMath::Log10(fVal[n]);
               else                    val[selected] = fVal[n];
               ++selected;
            }
         }
      } else {
         for (Long64_t n=first;n<first+nentries;++n) {
            if (TestBit(kLogScale)) val[selected] = TMath::Log10(fVal[n]);
            else                    val[selected] = fVal[n];
            ++selected;
         }
      }
      TMath::Quantiles(selected,3,val,quantiles,prob,kFALSE);
      delete [] val;
   }
   fQua1 = quantiles[0];
   fMed = quantiles[1];
   fQua3 = quantiles[2];
   delete [] quantiles;
   delete [] prob;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value corresponding to the position.

Double_t TParallelCoordVar::GetValuefromXY(Double_t x,Double_t y)
{
   Double_t pos;
   if (fMinCurrent == fMaxCurrent) return fMinCurrent;
   if (fX1 == fX2) {
      if (y<=fY1)      pos = fMinCurrent;
      else if (y>=fY2) pos = fMaxCurrent;
      else             pos = fMinCurrent + ((y-fY1)/(fY2-fY1))*(fMaxCurrent-fMinCurrent);
   } else {
      if (x<=fX1)      pos = fMinCurrent;
      else if (x>=fX2) pos = fMaxCurrent;
      else             pos = fMinCurrent + ((x-fX1)/(fX2-fX1))*(fMaxCurrent-fMinCurrent);
   }
   return pos;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a position corresponding to the value on the axis.

void TParallelCoordVar::GetXYfromValue(Double_t value, Double_t & x, Double_t & y)
{
   if(value < fMinCurrent || value > fMaxCurrent) return;

   if (fX1==fX2) {
      x = fX1;
      if (fMinCurrent != fMaxCurrent) {
         if (TestBit(kLogScale)) y = fY1 + (fY2 - fY1) *
                                    (TMath::Log10(value/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
         else                    y = fY1 + (fY2 - fY1) *
                                    (value - fMinCurrent) / (fMaxCurrent - fMinCurrent);
      } else {
         y = fY1 + 0.5*(fY2-fY1);
      }
   } else {
      y = fY1;
      if (fMinCurrent != fMaxCurrent) {
         if (TestBit(kLogScale)) x = fX1 + (fX2 - fX1) *
                                    (TMath::Log10(value/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
         else                    x = fX1 + (fX2 - fX1) *
                                    (value - fMinCurrent) / (fMaxCurrent - fMinCurrent);
      } else {
         x = fX1 + 0.5*(fX2-fX1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialise the TParallelVar variables.

void TParallelCoordVar::Init()
{
   fX1         = 0;
   fX2         = 0;
   fY1         = 0;
   fY2         = 0;
   fId         = 0;
   fVal        = NULL;
   fMean       = 0;
   fMinInit    = 0;
   fMinCurrent = 0;
   fMaxInit    = 0;
   fMaxCurrent = 0;
   fMed        = 0;
   fQua1       = 0;
   fQua3       = 0;
   fNentries   = 0;
   fParallel   = NULL;
   fHistogram  = NULL;
   fNbins      = 100;
   fHistoLW    = 2;
   fHistoHeight     = 0.5;
   fRanges     = NULL;
   SetBit(kLogScale,kFALSE);
   SetBit(kShowBox,kFALSE);
   SetBit(kShowBarHisto,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the axis.

void TParallelCoordVar::Paint(Option_t* option)
{
   TIter next(fRanges);
   TParallelCoordRange* range;
   while ((range = (TParallelCoordRange*)next())) range->Paint(option);

   PaintHistogram();
   if (TestBit(kShowBox)) PaintBoxPlot();
   PaintLabels();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the boxes in the case of a candle chart.

void TParallelCoordVar::PaintBoxPlot()
{
   TLine *line = new TLine();
   line->SetLineColor(GetLineColor());
   line->SetLineWidth(1);
   TBox *box = new TBox();
   box->SetLineWidth(1);
   box->SetLineColor(GetLineColor());
   box->SetLineStyle(1);
   box->SetFillStyle(0);

   TFrame* frame = gPad->GetFrame();

   Double_t boxSize;
   if (fParallel->GetNvar() > 1) {
      if (fX1==fX2) boxSize = fHistoHeight*((frame->GetY2()-frame->GetY1())/(fParallel->GetNvar()-1));
      else          boxSize = fHistoHeight*((frame->GetX2()-frame->GetX1())/(fParallel->GetNvar()-1));
      if (boxSize >= 0.03) boxSize = 0.03;
   }
   else boxSize = 0.03;

   Double_t qua1,med,qua3,max,min;
   Double_t a,b,maxinit,mininit;
   if (TestBit(kLogScale)) {
      a = TMath::Log10(fMinCurrent);
      b = TMath::Log10(fMaxCurrent/fMinCurrent);
      if(fMinInit > 0) mininit = TMath::Log10(fMinInit);
      else             mininit = TMath::Log10(fMinCurrent);
      maxinit = TMath::Log10(fMaxInit);
   } else {
      a = fMinCurrent;
      b = fMaxCurrent-fMinCurrent;
      mininit = fMinInit;
      maxinit = fMaxInit;
   }
   if(fX1==fX2) {
      qua1 = fY1 + ((fQua1-a)/b)*(fY2-fY1);
      qua3 = fY1 + ((fQua3-a)/b)*(fY2-fY1);
      med  = fY1 + ((fMed-a)/b)*(fY2-fY1);
      max  = fY1 + ((maxinit-a)/b)*(fY2-fY1);
      min  = fY1 + ((mininit-a)/b)*(fY2-fY1);
   } else {
      qua1 = fX1 + ((fQua1-a)/b)*(fX2-fX1);
      qua3 = fX1 + ((fQua3-a)/b)*(fX2-fX1);
      med  = fX1 + ((fMed-a)/b)*(fX2-fX1);
      max  = fX1 + ((maxinit-a)/b)*(fX2-fX1);
      min  = fX1 + ((mininit-a)/b)*(fX2-fX1);
   }

   // min and max lines.
   if (fX1==fX2) {
      line->PaintLine(fX1-boxSize,min,fX1+boxSize,min);
      line->PaintLine(fX2-boxSize,max,fX2+boxSize,max);
   } else {
      line->PaintLine(min,fY1-boxSize,min,fY1+boxSize);
      line->PaintLine(max,fY2-boxSize,max,fY2+boxSize);
   }

   // lines from min and max to the box.
   line->SetLineStyle(7);
   if (fX1==fX2) {
      if (min<frame->GetY1()) min = frame->GetY1();
      if (max>frame->GetY2()) max = frame->GetY2();
      line->PaintLine(fX1,min,fX1,qua1);
      line->PaintLine(fX1,qua3,fX1,max);
   } else {
      if (min<frame->GetX1()) min = frame->GetX1();
      if (max>frame->GetX2()) max = frame->GetX2();
      line->PaintLine(min,fY1,qua1,fY2);
      line->PaintLine(qua3,fY1,max,fY2);
   }

   // Box
   if(fX1==fX2) box->PaintBox(fX1-boxSize,qua1,fX1+boxSize,qua3);
   else box->PaintBox(qua1,fY1-boxSize,qua3,fY1+boxSize);

   // Median line
   line->SetLineStyle(1);
   if(fX1==fX2) line->PaintLine(fX1-boxSize,med,fX1+boxSize,med);
   else line->PaintLine(med,fY1-boxSize,med,fY1+boxSize);

   // Paint average
   if (!TestBit(kLogScale) || (TestBit(kLogScale) && fMean > 0)) {
      Double_t mean;
      if (TestBit(kLogScale)) mean = TMath::Log10(fMean);
      else mean = fMean;
      TMarker *mark = NULL;
      if(fX1==fX2) mark = new TMarker(fX1,fY1 + ((mean-a)/b)*(fY2-fY1),24);
      else         mark = new TMarker(fX1 + ((mean-a)/b)*(fX2-fX1),fY1,24);
      mark->Paint();
      delete mark;
   }

   delete line;
   delete box;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the histogram on the axis.

void TParallelCoordVar::PaintHistogram()
{
   Int_t i;

   TFrame *frame = gPad->GetFrame();

   if (!fHistogram) GetHistogram();

   // Paint the axis body.
   if (fHistoHeight!=0 && TestBit(kShowBarHisto)) {
      // Paint the axis body using bar chart.
      TBox *b = new TBox();
      b->SetFillStyle(GetFillStyle());
      b->SetFillColor(GetFillColor());
      b->SetLineStyle(1);
      b->SetLineColor(GetFillColor());
      b->SetLineWidth(1);
      Double_t hmin = fHistogram->GetMinimum();
      Double_t hmax = fHistogram->GetMaximum();
      if (fX1 == fX2) {
         // Vertical case.
         Double_t dy = (fY2-fY1)/fNbins;
         Double_t dv = (fMaxCurrent - fMinCurrent)/fNbins;
         Double_t v = fMinCurrent;
         Double_t y1 = fY1,x2,y2;
         for (i=1; i<=fNbins; i++) {
            x2 = fX1+((fHistogram->GetBinContent(i)-hmin)/(hmax-hmin))*fHistoHeight*
                 ((frame->GetX2()-frame->GetX1())/(fParallel->GetNvar()-1));
            if(TestBit(kLogScale)) y2 = fY1 + (fY2-fY1)*(TMath::Log10((v+dv)/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
            else y2=y1+dy;
            b->PaintBox(fX1,y1,x2,y2,"l");
            y1=y2;
            v += dv;
         }
      } else {
         // Horizontal case.
         Double_t dx = (fX2-fX1)/fNbins;
         Double_t dv = (fMaxCurrent - fMinCurrent)/fNbins;
         Double_t v = fMinCurrent;
         Double_t x1 = fX1,x2,y2;
         for (i=1; i<=fNbins; i++) {
            y2 = fY1+((fHistogram->GetBinContent(i)-hmin)/(hmax-hmin))*fHistoHeight*((frame->GetY2()-frame->GetY1())/(fParallel->GetNvar()-1));
            if(TestBit(kLogScale)) x2 = fX1 + (fX2-fX1)*(TMath::Log10((v+dv)/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
            else x2=x1+dx;
            b->PaintBox(x1,fY1,x2,y2,"l");
            x1=x2;
            v+=dv;
         }
      }
      delete b;
   }
   if (fHistoLW==0 && !TestBit(kShowBox)) {
      // Paint the axis body as a simple line.
      TLine* l = new TLine(fX1,fY1,fX2,fY2);
      l->SetLineWidth(GetLineWidth());
      l->SetLineColor(GetLineColor());
      l->SetLineStyle(GetLineColor());
      l->Paint();
      delete l;
   } else if (fHistoLW!=0){
      // Paint the axis body using the color palette.
      TLine *lb = new TLine();
      lb->SetLineWidth(fHistoLW);
      Double_t hmin = fHistogram->GetMinimum();
      Double_t hmax = fHistogram->GetMaximum();
      Int_t theColor;
      Int_t ncolors = gStyle->GetNumberOfColors();
      if (fX1 == fX2) {
         // Vertical case.
         Double_t dy = (fY2-fY1)/fNbins;
         Double_t y1 = fY1,y2;
         Double_t dv = (fMaxCurrent - fMinCurrent)/fNbins;
         Double_t v = fMinCurrent;
         for (i=1; i<=fNbins; i++) {
            theColor = (Int_t)( ((fHistogram->GetBinContent(i)-hmin)/(hmax-hmin))*(ncolors-1) );
            if(TestBit(kLogScale)) y2 = fY1 + (fY2-fY1)*(TMath::Log10((v+dv)/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
            else y2=y1+dy;
            lb->SetLineColor(gStyle->GetColorPalette(theColor));
            lb->PaintLine(fX1,y1,fX1,y2);
            y1=y2;
            v+=dv;
         }
      } else {
         // Horizontal case.
         Double_t dx = (fX2-fX1)/fNbins;
         Double_t dv = (fMaxCurrent - fMinCurrent)/fNbins;
         Double_t v = fMinCurrent;
         Double_t x1 = fX1,x2;
         for (i=1; i<=fNbins; i++) {
            theColor = (Int_t)( ((fHistogram->GetBinContent(i)-hmin)/(hmax-hmin))*(ncolors-1) );
            lb->SetLineColor(gStyle->GetColorPalette(theColor));
            if(TestBit(kLogScale)) x2 = fX1 + (fX2-fX1)*(TMath::Log10((v+dv)/fMinCurrent)) / (TMath::Log10(fMaxCurrent/fMinCurrent));
            else x2=x1+dx;
            lb->PaintLine(x1,fY1,x2,fY1);
            x1=x2;
            v+=dv;
         }
      }
      delete lb;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the axis labels and titles.

void TParallelCoordVar::PaintLabels()
{
   TLatex* t = new TLatex();
   TFrame *frame = gPad->GetFrame();
   t->SetTextSize(0.03);
   if (fX1==fX2) {
      t->SetText(fX1,frame->GetY1() - 0.04 - t->GetTextSize(),GetName());
      Double_t tlength = t->GetXsize();
      if (fX1-0.5*tlength<0.01) {
         t->SetTextAlign(11);
         t->SetText(0.01, frame->GetY1() - 0.04 - t->GetTextSize(), GetName());
         t->Paint();
      } else if (fX1+0.5*tlength > 0.99) {
         t->SetTextAlign(31);
         t->SetText(0.99,frame->GetY1() - 0.04 - t->GetTextSize(),GetName());
         t->Paint();
      } else {
         t->SetTextAlign(21);
         t->PaintLatex(fX1,frame->GetY1() - 0.04 - t->GetTextSize(),0,0.03,GetName());
      }
      if (!fParallel->TestBit(TParallelCoord::kCandleChart)) {
         t->SetTextAlign(21);
         t->PaintLatex(fX1,frame->GetY2() + 0.005,0,0.025,Form("%g",fMaxCurrent));
         t->SetTextAlign(23);
         t->PaintLatex(fX1,frame->GetY1() - 0.005,0,0.025,Form("%g",fMinCurrent));
      }
   } else {
      t->SetText(fX1-0.04,fY1+0.02,GetName());
      t->SetTextSize(0.03);
      Double_t tlength = t->GetXsize();
      if (fX1-0.04-tlength<0.01) {
         t->SetTextAlign(12);
         t->SetText(0.01,fY1+0.02,GetName());
         t->Paint();
      } else {
         t->SetTextAlign(32);
         t->PaintLatex(fX1-0.04,fY1+0.02,0,0.03,GetName());
      }
      if (!fParallel->TestBit(TParallelCoord::kCandleChart)) {
         t->SetTextAlign(12);
         t->PaintLatex(0.01,fY1-0.02,0,0.025,Form("%g",fMinCurrent));
         t->SetTextAlign(32);
         t->PaintLatex(0.99,fY1-0.02,0,0.025,Form("%g",fMaxCurrent));
      }
   }
   delete t;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the axis main data.

void TParallelCoordVar::Print(Option_t* /*option*/) const
{
   printf("**************variable #%d**************\n",fParallel->GetVarList()->IndexOf(this));
   printf("at x1=%f, y1=%f, x2=%f, y2=%f.\n",fX1,fY1,fX2,fY2);
   printf("min = %f, Q1 = %f, Med = %f, Q3 = %f, Max = %f\n", fMinInit, fQua1, fMed, fQua3, fMaxInit);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the TParallelCoordVar as a macro. Can be used only in the context
/// of TParallelCoord::SavePrimitive (pointer "TParallelCoord* para" is
/// defined in TParallelCoord::SavePrimitive) with the option "pcalled".

void TParallelCoordVar::SavePrimitive(std::ostream & out, Option_t* options)
{
   TString opt = options;
   if (opt.Contains("pcalled")) {
      out<<"   var->SetBit(TParallelCoordVar::kLogScale,"<<TestBit(kLogScale)<<");"<<std::endl;
      out<<"   var->SetBit(TParallelCoordVar::kShowBox,"<<TestBit(kShowBox)<<");"<<std::endl;
      out<<"   var->SetBit(TParallelCoordVar::kShowBarHisto,"<<TestBit(kShowBarHisto)<<");"<<std::endl;
      out<<"   var->SetHistogramBinning("<<fNbins<<");"<<std::endl;
      out<<"   var->SetHistogramLineWidth("<<fHistoLW<<");"<<std::endl;
      out<<"   var->SetInitMin("<<fMinInit<<");"<<std::endl;
      out<<"   var->SetInitMax("<<fMaxInit<<");"<<std::endl;
      out<<"   var->SetHistogramHeight("<<fHistoHeight<<");"<<std::endl;
      out<<"   var->GetMinMaxMean();"<<std::endl;
      out<<"   var->GetHistogram();"<<std::endl;
      out<<"   var->SetFillStyle("<<GetFillStyle()<<");"<<std::endl;
      out<<"   var->SetFillColor("<<GetFillColor()<<");"<<std::endl;
      out<<"   var->SetLineColor("<<GetLineColor()<<");"<<std::endl;
      out<<"   var->SetLineWidth("<<GetLineWidth()<<");"<<std::endl;
      out<<"   var->SetLineStyle("<<GetLineStyle()<<");"<<std::endl;
      if (TestBit(kShowBox)) out<<"   var->GetQuantiles();"<<std::endl;
      TIter next(fRanges);
      TParallelCoordRange* range;
      Int_t i = 1;
      while ((range = (TParallelCoordRange*)next())) {
         out<<"   //***************************************"<<std::endl;
         out<<"   // Create the "<<i<<"th range owned by the axis \""<<GetTitle()<<"\"."<<std::endl;
         out<<"   TParallelCoordSelect* sel = para->GetSelection(\""<<range->GetSelection()->GetTitle()<<"\");"<<std::endl;
         out<<"   TParallelCoordRange* newrange = new TParallelCoordRange(var,"<<range->GetMin()<<","<<range->GetMax()<<",sel);"<<std::endl;
         out<<"   var->AddRange(newrange);"<<std::endl;
         out<<"   sel->Add(newrange);"<<std::endl;
         ++i;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the axis to display a candle.

void TParallelCoordVar::SetBoxPlot(Bool_t box)
{
   SetBit(kShowBox,box);
   if (box) SetHistogramHeight(0.5);
   else {
      SetHistogramHeight(0);
      SetHistogramLineWidth(2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the histogram binning.

void TParallelCoordVar::SetHistogramBinning(Int_t n)
{
   if (n < 0 || n == fNbins) return;
   fNbins = n;
   GetHistogram();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the height of the bar histogram.

void TParallelCoordVar::SetHistogramHeight(Double_t h)
{
   fHistoHeight = h;
   if (!fParallel->TestBit(TParallelCoord::kCandleChart)){
      if(h!=0) SetBit(kShowBarHisto,kTRUE);
      else SetBit(kShowBarHisto,kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current minimum of the axis.

void TParallelCoordVar::SetCurrentMin(Double_t min)
{
   fMinCurrent = min;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current maximum of the axis.

void TParallelCoordVar::SetCurrentMax(Double_t max)
{
   fMaxCurrent = max;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the limits within which one the entries must be painted.

void TParallelCoordVar::SetCurrentLimits(Double_t min, Double_t max)
{
   if (min>max) {
      Double_t mem = min;
      min = max;
      max = mem;
   }
   if(TestBit(kLogScale) && max<=0) return;
   if(TestBit(kLogScale) && min<=0) min = 0.00001*max;
   fMinCurrent = min;
   fMaxCurrent = max;

   delete fHistogram;
   fHistogram = NULL;
   GetHistogram();

   if (fParallel->TestBit(TParallelCoord::kGlobalScale)) {
      fParallel->SetGlobalMin(min);
      fParallel->SetGlobalMax(max);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If true, the pad is updated while the motion of a dragged range.

void TParallelCoordVar::SetLiveRangesUpdate(Bool_t on)
{
   TIter next(fRanges);
   TParallelCoordRange* range;
   while ((range = (TParallelCoordRange*)next())) range->SetBit(TParallelCoordRange::kLiveUpdate,on);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the axis in log scale.

void TParallelCoordVar::SetLogScale(Bool_t log)
{
   if (log == TestBit (kLogScale)) return;
   if (fMaxInit < 0)             SetBit(kLogScale,kFALSE);
   else if (log) {
      if (fMaxCurrent < 0 ) fMaxCurrent = fMaxInit;
      if (fMinCurrent < 0 ) fMinCurrent = 0.00001*fMaxCurrent;
      SetBit(kLogScale,kTRUE);
      SetCurrentMin(fMinCurrent);
      SetCurrentMax(fMaxCurrent);
   } else {
      SetBit(kLogScale,kFALSE);
      SetCurrentMin(fMinInit);
      SetCurrentMax(fMaxInit);
   }
   GetQuantiles();
   GetHistogram();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the variable values.

void TParallelCoordVar::SetValues(Long64_t length, Double_t* val)
{
   if (fVal) delete [] fVal;
   fVal = new Double_t[length];
   fNentries = length;
   for (Long64_t li = 0; li < length; ++li) fVal[li] = val[li];
   GetMinMaxMean();
   GetHistogram();
   if (TestBit(kShowBox)) GetQuantiles();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the X position of the axis in the case of a vertical axis.
/// and rotate the axis if it was horizontal.

void TParallelCoordVar::SetX(Double_t x, Bool_t gl)
{
   TFrame *frame = gPad->GetFrame();
   if (!gl) {
      fY1 = frame->GetY1();
      fY2 = frame->GetY2();
   } else {
      Double_t gmin = fParallel->GetGlobalMin();
      Double_t gmax = fParallel->GetGlobalMax();
      fY1 = frame->GetY1() + ((fMinCurrent-gmin)/(gmax-gmin))*(frame->GetY2()-frame->GetY1());
      fY2 = frame->GetY1() + ((fMaxCurrent-gmin)/(gmax-gmin))*(frame->GetY2()-frame->GetY1());
   }
   fX1 = fX2 = x;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the Y position of the axis in the case of a horizontal axis.
/// and rotate the axis if it was vertical.

void TParallelCoordVar::SetY(Double_t y, Bool_t gl)
{
   TFrame *frame = gPad->GetFrame();
   if (!gl) {
      fX1 = frame->GetX1();
      fX2 = frame->GetX2();
   } else {
      Double_t gmin = fParallel->GetGlobalMin();
      Double_t gmax = fParallel->GetGlobalMax();
      fX1 = frame->GetX1() + ((fMinCurrent-gmin)/(gmax-gmin))*(frame->GetX2()-frame->GetX1());
      fX2 = frame->GetX1() + ((fMaxCurrent-gmin)/(gmax-gmin))*(frame->GetX2()-frame->GetX1());
   }
   fY1 = fY2 = y;
}
