// @(#)root/treeviewer:$Id$
// Author: Bastien Dalla Piazza  02/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TParallelCoordRange.h"
#include "TParallelCoord.h"
#include "TParallelCoordVar.h"

#include "TPolyLine.h"
#include "TList.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TPoint.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TString.h"

#include <iostream>

ClassImp(TParallelCoordRange);

/** \class TParallelCoordRange
A TParallelCoordRange is a range used for parallel coordinates plots.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TParallelCoordRange::TParallelCoordRange()
   :TNamed("Range","Range"), TAttLine(), fSize(0.01)
{
   fMin = 0;
   fMax = 0;
   fVar = NULL;
   fSelect = NULL;
   SetBit(kShowOnPad,kTRUE);
   SetBit(kLiveUpdate,kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TParallelCoordRange::~TParallelCoordRange()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TParallelCoordRange::TParallelCoordRange(TParallelCoordVar *var, Double_t min, Double_t max, TParallelCoordSelect *sel)
   :TNamed("Range","Range"), TAttLine(1,1,1), fSize(0.01)
{
   if(min == max) {
      min = var->GetCurrentMin();
      max = var->GetCurrentMax();
   }
   fMin = min;
   fMax = max;

   fVar = var;
   fSelect = NULL;

   if (!sel) {
      TParallelCoordSelect* s = var->GetParallel()->GetCurrentSelection();
      if (s) fSelect = s;
      else return;
   } else {
      fSelect = sel;
   }

   SetLineColor(fSelect->GetLineColor());

   SetBit(kShowOnPad,kTRUE);
   SetBit(kLiveUpdate,var->GetParallel()->TestBit(TParallelCoord::kLiveUpdate));
}

////////////////////////////////////////////////////////////////////////////////
/// Make the selection which owns the range to be drawn on top of the others.

void TParallelCoordRange::BringOnTop()
{
   TList *list = fVar->GetParallel()->GetSelectList();
   list->Remove(fSelect);
   list->AddLast(fSelect);
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the range.

void TParallelCoordRange::Delete(const Option_t* /*options*/)
{
   fVar->GetRanges()->Remove(this);
   fVar->GetParallel()->CleanUpSelections(this);
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the distance to the primitive.

Int_t TParallelCoordRange::DistancetoPrimitive(Int_t px, Int_t py)
{
   if(TestBit(kShowOnPad)){
      Double_t xx,yy,thisx=0,thisy=0;
      xx = gPad->AbsPixeltoX(px);
      yy = gPad->AbsPixeltoY(py);
      fVar->GetXYfromValue(fMin,thisx,thisy);
      Int_t dist = 9999;
      if(fVar->GetVert()){
         if(xx > thisx-2*fSize && xx < thisx && yy > thisy-fSize && yy<thisy+fSize) dist = 0;
         fVar->GetXYfromValue(fMax,thisx,thisy);
         if(xx > thisx-2*fSize && xx < thisx && yy > thisy-fSize && yy<thisy+fSize) dist = 0;
      } else {
         if(yy > thisy-2*fSize && yy < thisy && xx > thisx-fSize && xx<thisx+fSize) dist = 0;
         fVar->GetXYfromValue(fMax,thisx,thisy);
         if(yy > thisy-2*fSize && yy < thisy && xx > thisx-fSize && xx<thisx+fSize) dist = 0;
      }
      return dist;
   } else return 9999;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a TParallelCoordRange.

void TParallelCoordRange::Draw(Option_t* options)
{
   AppendPad(options);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the entry.

void TParallelCoordRange::ExecuteEvent(Int_t entry, Int_t px, Int_t py)
{
   if (!gPad) return;
   if (!gPad->IsEditable() && entry!=kMouseEnter) return;

   Bool_t vert = fVar->GetVert();
   static Int_t pxold, pyold;
   static Int_t mindragged = -1; //-1:nothing dragged, 0:max dragged, 1:mindragged, 2:both dragged;
   Int_t plx1,plx2,ply1,ply2;

   Double_t xx,yy,txxmin=0,txxmax=0,tyymin=0,tyymax=0;
   TFrame *frame = gPad->GetFrame();
   xx = gPad->AbsPixeltoX(px);
   yy = gPad->AbsPixeltoY(py);
   fVar->GetXYfromValue(fMin,txxmin,tyymin);
   fVar->GetXYfromValue(fMax,txxmax,tyymax);
   if (vert) {
      plx1 = gPad->XtoAbsPixel(txxmin-2*fSize);
      plx2 = gPad->XtoAbsPixel(txxmax-2*fSize);
      ply1 = gPad->YtoAbsPixel(tyymin+fSize);
      ply2 = gPad->YtoAbsPixel(tyymax-fSize);
   } else {
      plx1 = gPad->XtoAbsPixel(txxmin+fSize);
      plx2 = gPad->XtoAbsPixel(txxmax-fSize);
      ply1 = gPad->YtoAbsPixel(tyymin-2*fSize);
      ply2 = gPad->YtoAbsPixel(tyymax-2*fSize);
   }

   gPad->SetCursor(kPointer);
   gVirtualX->SetLineColor(-1);
   gVirtualX->SetLineWidth(1);
   TPoint *p = NULL;
   switch (entry) {
      case kButton1Down:
         fVar->GetParallel()->SetCurrentSelection(fSelect);
         ((TCanvas*)gPad)->Selected(gPad,fVar->GetParallel(),1);
         if ((vert && yy<tyymax-fSize) || (!vert && xx < txxmax-fSize)) {     //checks if the min slider is clicked.
            mindragged = 1;
            p = GetSliderPoints(fMin);
            gVirtualX->DrawPolyLine(5,p);
            delete [] p;
         } else {
            mindragged = 0;
            p = GetSliderPoints(fMax);
            gVirtualX->DrawPolyLine(5,p);
            delete [] p;
         }
         gVirtualX->DrawLine(plx1,ply1,plx2,ply2);
         break;
      case kButton1Up: {
         Double_t min = fMin, max= fMax;
         if (mindragged == 1) min = fVar->GetValuefromXY(xx,yy);
         if (mindragged == 0) max = fVar->GetValuefromXY(xx,yy);
         if(fMin!=min || fMax != max) {
            if (min>max) {
               Double_t mem = min;
               min = max;
               max = mem;
            }
            fMin = min;
            fMax = max;
            gPad->Modified();
         }
         mindragged = -1;
         break;
      }
      case kMouseMotion:
         pxold = px;
         pyold = py;
         break;
      case kButton1Motion:
         if((vert && yy > frame->GetY1() && yy < frame->GetY2()) ||
            (!vert && xx > frame->GetX1() && xx < frame->GetX2())){
            if (vert) p = GetSliderPoints(pyold);
            else      p = GetSliderPoints(pxold);
            gVirtualX->DrawPolyLine(5,p);
            delete [] p;
            if (vert) p = GetBindingLinePoints(pyold,mindragged);
            else p = GetBindingLinePoints(pxold,mindragged);
            gVirtualX->DrawPolyLine(2,p);
            delete [] p;
            if (vert) p = GetSliderPoints(py);
            else      p = GetSliderPoints(px);
            gVirtualX->DrawPolyLine(5,p);
            delete [] p;
            if (vert) p = GetBindingLinePoints(py,mindragged);
            else p = GetBindingLinePoints(px,mindragged);
            gVirtualX->DrawPolyLine(2,p);
            delete [] p;
            if (TestBit(kLiveUpdate)){
               Double_t min = fMin, max= fMax;
               if (mindragged == 1) min = fVar->GetValuefromXY(xx,yy);
               if (mindragged == 0) max = fVar->GetValuefromXY(xx,yy);
               if(fMin!=min || fMax != max) {
                  if (min>max) {
                     Double_t mem = min;
                     min = max;
                     max = mem;
                  }
                  fMin = min;
                  fMax = max;
                  gPad->Modified();
                  gPad->Update();
               }
            }
         }
         pxold = px;
         pyold = py;
         break;
      default:
         //std::cout<<"entry: "<<entry<<std::endl;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the points of the line binding the two needles of the range.

TPoint* TParallelCoordRange::GetBindingLinePoints(Int_t pos,Int_t mindragged)
{
   Double_t txx,tyy,txxo,tyyo=0;
   if (fVar->GetVert()){
      txx = fVar->GetX();
      tyy = gPad->AbsPixeltoY(pos);
   } else {
      tyy = fVar->GetY();
      txx = gPad->AbsPixeltoX(pos);
   }
   if (mindragged==1) fVar->GetXYfromValue(fMax,txxo,tyyo);
   else fVar->GetXYfromValue(fMin,txxo,tyyo);

   TPoint *bindline = new TPoint[2];
   if (fVar->GetVert()) {
      if (mindragged==1) {
         bindline[0] = TPoint(gPad->XtoAbsPixel(txx-2*fSize),gPad->YtoAbsPixel(tyy+fSize));
         bindline[1] = TPoint(gPad->XtoAbsPixel(txx-2*fSize),gPad->YtoAbsPixel(tyyo-fSize));
      } else {
         bindline[0] = TPoint(gPad->XtoAbsPixel(txx-2*fSize),gPad->YtoAbsPixel(tyyo+fSize));
         bindline[1] = TPoint(gPad->XtoAbsPixel(txx-2*fSize),gPad->YtoAbsPixel(tyy-fSize));
      }
   } else {
      if (mindragged==1) {
         bindline[0] = TPoint(gPad->XtoAbsPixel(txx+fSize),gPad->YtoAbsPixel(tyy-2*fSize));
         bindline[1] = TPoint(gPad->XtoAbsPixel(txxo-fSize),gPad->YtoAbsPixel(tyy-2*fSize));
      } else {
         bindline[0] = TPoint(gPad->XtoAbsPixel(txxo+fSize),gPad->YtoAbsPixel(tyy-2*fSize));
         bindline[1] = TPoint(gPad->XtoAbsPixel(txx-fSize),gPad->YtoAbsPixel(tyy-2*fSize));
      }
   }
   return bindline;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the points to paint the needles at "value".

TPoint* TParallelCoordRange::GetSliderPoints(Double_t value)
{
   Double_t txx=0,tyy=0;
   fVar->GetXYfromValue(value,txx,tyy);
   Int_t tx[5];
   Int_t ty[5];
   if (fVar->GetVert()) {
      tx[0]=gPad->XtoAbsPixel(txx);
      tx[1]=tx[4]=gPad->XtoAbsPixel(txx-fSize);
      ty[0]=ty[1]=ty[4]=gPad->YtoAbsPixel(tyy);
      tx[2]=tx[3]=gPad->XtoAbsPixel(txx-2*fSize);
      ty[2]=gPad->YtoAbsPixel(tyy+fSize);
      ty[3]=gPad->YtoAbsPixel(tyy-fSize);
   } else {
      ty[0]=gPad->YtoAbsPixel(tyy);
      ty[1]=ty[4]=gPad->YtoAbsPixel(tyy-fSize);
      tx[0]=tx[1]=tx[4]=gPad->XtoAbsPixel(txx);
      ty[2]=ty[3]=gPad->YtoAbsPixel(tyy-2*fSize);
      tx[2]=gPad->XtoAbsPixel(txx-fSize);
      tx[3]=gPad->XtoAbsPixel(txx+fSize);
   }
   TPoint *slider = new TPoint[5];
   for(UInt_t ui=0;ui<5;++ui) slider[ui] = TPoint(tx[ui],ty[ui]);
   return slider;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the points to paint the needle at "pos".

TPoint* TParallelCoordRange::GetSliderPoints(Int_t pos)
{
   Double_t txx,tyy;
   if (fVar->GetVert()){
      txx = fVar->GetX();
      tyy = gPad->AbsPixeltoY(pos);
   } else {
      tyy = fVar->GetY();
      txx = gPad->AbsPixeltoX(pos);
   }

   Int_t tx[5];
   Int_t ty[5];
   if (fVar->GetVert()) {
      tx[0]=gPad->XtoAbsPixel(txx);
      tx[1]=tx[4]=gPad->XtoAbsPixel(txx-fSize);
      ty[0]=ty[1]=ty[4]=gPad->YtoAbsPixel(tyy);
      tx[2]=tx[3]=gPad->XtoAbsPixel(txx-2*fSize);
      ty[2]=gPad->YtoAbsPixel(tyy+fSize);
      ty[3]=gPad->YtoAbsPixel(tyy-fSize);
   } else {
      ty[0]=gPad->YtoAbsPixel(tyy);
      ty[1]=ty[4]=gPad->YtoAbsPixel(tyy-fSize);
      tx[0]=tx[1]=tx[4]=gPad->XtoAbsPixel(txx);
      ty[2]=ty[3]=gPad->YtoAbsPixel(tyy-2*fSize);
      tx[2]=gPad->XtoAbsPixel(txx-fSize);
      tx[3]=gPad->XtoAbsPixel(txx+fSize);
   }
   TPoint *slider = new TPoint[5];
   for(UInt_t ui=0;ui<5;++ui) slider[ui] = TPoint(tx[ui],ty[ui]);
   return slider;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate if the given value is within the range or not.

Bool_t TParallelCoordRange::IsIn(Double_t evtval)
{
   return evtval>=fMin && evtval<=fMax;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint a TParallelCoordRange.

void TParallelCoordRange::Paint(Option_t* /*options*/)
{
   if(TestBit(kShowOnPad)){
      PaintSlider(fMin,kTRUE);
      PaintSlider(fMax,kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint a slider.

void TParallelCoordRange::PaintSlider(Double_t value, Bool_t fill)
{
   SetLineColor(fSelect->GetLineColor());

   TPolyLine *p= new TPolyLine();
   p->SetLineStyle(1);
   p->SetLineColor(1);
   p->SetLineWidth(1);

   Double_t *x = new Double_t[5];
   Double_t *y = new Double_t[5];

   Double_t xx,yy;

   fVar->GetXYfromValue(value,xx,yy);
   if(fVar->GetVert()){
      x[0] = xx; x[1]=x[4]=xx-fSize; x[2]=x[3]=xx-2*fSize;
      y[0]=y[1]=y[4]=yy; y[2] = yy+fSize; y[3] = yy-fSize;
   } else {
      y[0] = yy; y[1]=y[4]=yy-fSize; y[2]=y[3]= yy-2*fSize;
      x[0]=x[1]=x[4]=xx; x[2]=xx-fSize; x[3] = xx+fSize;
   }
   if (fill) {
      p->SetFillStyle(1001);
      p->SetFillColor(0);
      p->PaintPolyLine(4,&x[1],&y[1],"f");
      p->SetFillColor(GetLineColor());
      p->SetFillStyle(3001);
      p->PaintPolyLine(4,&x[1],&y[1],"f");
   }
   p->PaintPolyLine(5,x,y);

   delete p;
   delete [] x;
   delete [] y;
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about the range.

void TParallelCoordRange::Print(Option_t* /*options*/) const
{
   printf("On \"%s\" : min = %f, max = %f\n", fVar->GetTitle(), fMin, fMax);
}

////////////////////////////////////////////////////////////////////////////////
/// Make the selection which owns the range to be drawn under all the others.

void TParallelCoordRange::SendToBack()
{
   TList *list = fVar->GetParallel()->GetSelectList();
   list->Remove(fSelect);
   list->AddFirst(fSelect);
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the selection line color.

void  TParallelCoordRange::SetLineColor(Color_t col)
{
   fSelect->SetLineColor(col);
   TAttLine::SetLineColor(col);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the selection line width.

void  TParallelCoordRange::SetLineWidth(Width_t wid)
{
   fSelect->SetLineWidth(wid);
}


ClassImp(TParallelCoordSelect);

/** \class TParallelCoordSelect
A TParallelCoordSelect is a specialised TList to hold TParallelCoordRanges used
by TParallelCoord.

Selections of specific entries can be defined over the data se using parallel
coordinates. With that representation, a selection is an ensemble of ranges
defined on the axes. Ranges defined on the same axis are conjugated with OR
(an entry must be in one or the other ranges to be selected). Ranges on
different axes are are conjugated with AND (an entry must be in all the ranges
to be selected). Several selections can be defined with different colors. It is
possible to generate an entry list from a given selection and apply it to the
tree using the editor ("Apply to tree" button).
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TParallelCoordSelect::TParallelCoordSelect()
   : TList(), TAttLine(kBlue,1,1)
{
   fTitle = "Selection";
   SetBit(kActivated,kTRUE);
   SetBit(kShowRanges,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TParallelCoordSelect::TParallelCoordSelect(const char* title)
   : TList(), TAttLine(kBlue,1,1)
{
   fTitle = title;
   SetBit(kActivated,kTRUE);
   SetBit(kShowRanges,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TParallelCoordSelect::~TParallelCoordSelect()
{
   TIter next(this);
   TParallelCoordRange* range;
   while ((range = (TParallelCoordRange*)next())) range->GetVar()->GetRanges()->Remove(range);
   TList::Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Activate the selection.

void TParallelCoordSelect::SetActivated(Bool_t on)
{
   TIter next(this);
   TParallelCoordRange* range;
   while ((range = (TParallelCoordRange*)next())) range->SetBit(TParallelCoordRange::kShowOnPad,on);
   SetBit(kActivated,on);
}

////////////////////////////////////////////////////////////////////////////////
/// Show the ranges needles.

void TParallelCoordSelect::SetShowRanges(Bool_t s)
{
   TIter next(this);
   TParallelCoordRange* range;
   while ((range = (TParallelCoordRange*)next())) range->SetBit(TParallelCoordRange::kShowOnPad,s);
   SetBit(kShowRanges,s);
}
