// @(#)root/graf:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TROOT.h"
#include "TBuffer.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TPolyLine.h"

ClassImp(TPolyLine);

/** \class TPolyLine
\ingroup BasicGraphics

Defined by an array on N points in a 2-D space.

One can draw the contour of the polyline or/and its fill area.
Example:
Begin_Macro(source)
{
   Double_t x[5] = {.2,.7,.6,.25,.2};
   Double_t y[5] = {.5,.1,.9,.7,.5};
   TPolyLine *pline = new TPolyLine(5,x,y);
   pline->SetFillColor(38);
   pline->SetLineColor(2);
   pline->SetLineWidth(4);
   pline->Draw("f");
   pline->Draw();
}
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// PolyLine default constructor.

TPolyLine::TPolyLine(): TObject()
{
   fN = 0;
   fX = 0;
   fY = 0;
   fLastPoint = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// PolyLine normal constructor without initialisation.
/// Allocates n points.  The option string is ignored.

TPolyLine::TPolyLine(Int_t n, Option_t *option)
      :TObject(), TAttLine(), TAttFill()
{
   fOption = option;
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      fX = fY = 0;
      return;
   }
   fN = n;
   fX = new Double_t[fN];
   fY = new Double_t[fN];
}

////////////////////////////////////////////////////////////////////////////////
/// PolyLine normal constructor (single precision).
/// Makes n points with (x, y) coordinates from x and y.
/// The option string is ignored.

TPolyLine::TPolyLine(Int_t n, Float_t *x, Float_t *y, Option_t *option)
      :TObject(), TAttLine(), TAttFill()
{
   fOption = option;
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      fX = fY = 0;
      return;
   }
   fN = n;
   fX = new Double_t[fN];
   fY = new Double_t[fN];
   if (!x || !y) return;
   for (Int_t i=0; i<fN;i++) { fX[i] = x[i]; fY[i] = y[i];}
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// PolyLine normal constructor (double precision).
/// Makes n points with (x, y) coordinates from x and y.
/// The option string is ignored.

TPolyLine::TPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option)
      :TObject(), TAttLine(), TAttFill()
{
   fOption = option;
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      fX = fY = 0;
      return;
   }
   fN = n;
   fX = new Double_t[fN];
   fY = new Double_t[fN];
   if (!x || !y) return;
   for (Int_t i=0; i<fN;i++) { fX[i] = x[i]; fY[i] = y[i];}
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TPolyLine& TPolyLine::operator=(const TPolyLine& pl)
{
   if(this!=&pl) {
      pl.Copy(*this);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// PolyLine default destructor.

TPolyLine::~TPolyLine()
{
   if (fX) delete [] fX;
   if (fY) delete [] fY;
}

////////////////////////////////////////////////////////////////////////////////
/// PolyLine copy constructor.

TPolyLine::TPolyLine(const TPolyLine &polyline) : TObject(polyline), TAttLine(polyline), TAttFill(polyline)
{
   fN = 0;
   fX = 0;
   fY = 0;
   fLastPoint = -1;
   ((TPolyLine&)polyline).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this polyline to polyline.

void TPolyLine::Copy(TObject &obj) const
{
   TObject::Copy(obj);
   TAttLine::Copy(((TPolyLine&)obj));
   TAttFill::Copy(((TPolyLine&)obj));
   ((TPolyLine&)obj).fN = fN;
   delete [] ((TPolyLine&)obj).fX;
   delete [] ((TPolyLine&)obj).fY;
   if (fN > 0) {
      ((TPolyLine&)obj).fX = new Double_t[fN];
      ((TPolyLine&)obj).fY = new Double_t[fN];
      for (Int_t i=0; i<fN;i++)  {((TPolyLine&)obj).fX[i] = fX[i]; ((TPolyLine&)obj).fY[i] = fY[i];}
   } else {
      ((TPolyLine&)obj).fX = 0;
      ((TPolyLine&)obj).fY = 0;
   }
   ((TPolyLine&)obj).fOption = fOption;
   ((TPolyLine&)obj).fLastPoint = fLastPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns closest distance in pixels from point (px, py) to a polyline.
///
/// First looks for distances to the points of the polyline.  Stops search
/// and returns if a vertex of the polyline is found to be closer than 10
/// pixels.  Thus the return value may depend on the ordering of points
/// in the polyline.
///
/// Then looks for distances to the lines of the polyline.  There is no
/// arbitrary cutoff; any distance may be found.
///
/// Finally checks whether (px, py) is inside a closed and filled polyline.
/// (Must be EXACTLY closed.  "Filled" means fill color and fill style are
/// both non-zero.) If so, returns zero.
///
/// Returns 9999 if the polyline has no points.

Int_t TPolyLine::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t big = 9999;
   const Int_t kMaxDiff = 10;

   // check if point is near one of the points
   Int_t i, pxp, pyp, d;
   Int_t distance = big;
   if (Size() <= 0) return distance;

   for (i=0;i<Size();i++) {
      pxp = gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      pyp = gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
      d   = TMath::Abs(pxp-px) + TMath::Abs(pyp-py);
      if (d < distance) distance = d;
   }
   if (distance < kMaxDiff) return distance;

   // check if point is near one of the connecting lines
   for (i=0;i<Size()-1;i++) {
      d = DistancetoLine(px, py, gPad->XtoPad(fX[i]), gPad->YtoPad(fY[i]), gPad->XtoPad(fX[i+1]), gPad->YtoPad(fY[i+1]));
      if (d < distance) distance = d;
   }

   // in case of a closed and filled polyline, check if we are inside
   if (fFillColor && fFillStyle && fX[0] == fX[fLastPoint] && fY[0] == fY[fLastPoint]) {
      if (TMath::IsInside(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),fLastPoint+1,fX,fY)) distance = 0;
   }
   return distance;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this polyline with its current attributes.

void TPolyLine::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this polyline with new coordinates.

void TPolyLine::DrawPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option)
{
   TPolyLine *newpolyline = new TPolyLine(n,x,y);
   TAttLine::Copy(*newpolyline);
   TAttFill::Copy(*newpolyline);
   newpolyline->fOption = fOption;
   newpolyline->SetBit(kCanDelete);
   newpolyline->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function is called when a polyline is clicked with the locator
///
///  If Left button clicked on one of the line end points, this point
///     follows the cursor until button is released.
///
///  if Middle button clicked, the line is moved parallel to itself
///     until the button is released.

void TPolyLine::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   Int_t i, d;
   Double_t xmin, xmax, ymin, ymax, dx, dy, dxr, dyr;
   const Int_t kMaxDiff = 10;
   static Bool_t middle;
   static Int_t ipoint, pxp, pyp;
   static Int_t px1,px2,py1,py2;
   static Int_t pxold, pyold, px1old, py1old, px2old, py2old;
   static Int_t dpx, dpy;
   static Int_t *x=0, *y=0;
   Bool_t opaque  = gPad->OpaqueMoving();

   if (!gPad->IsEditable()) return;

   Int_t np = Size();

   switch (event) {

   case kButton1Down:
      gVirtualX->SetLineColor(-1);
      TAttLine::Modify();  //Change line attributes only if necessary
      px1 = gPad->XtoAbsPixel(gPad->GetX1());
      py1 = gPad->YtoAbsPixel(gPad->GetY1());
      px2 = gPad->XtoAbsPixel(gPad->GetX2());
      py2 = gPad->YtoAbsPixel(gPad->GetY2());
      ipoint = -1;


      if (x || y) break;
      x = new Int_t[np+1];
      y = new Int_t[np+1];
      for (i=0;i<np;i++) {
         pxp = gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
         pyp = gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
         if (!opaque) {
            gVirtualX->DrawLine(pxp-4, pyp-4, pxp+4,  pyp-4);
            gVirtualX->DrawLine(pxp+4, pyp-4, pxp+4,  pyp+4);
            gVirtualX->DrawLine(pxp+4, pyp+4, pxp-4,  pyp+4);
            gVirtualX->DrawLine(pxp-4, pyp+4, pxp-4,  pyp-4);
         }
         x[i] = pxp;
         y[i] = pyp;
         d   = TMath::Abs(pxp-px) + TMath::Abs(pyp-py);
         if (d < kMaxDiff) ipoint =i;
      }
      dpx = 0;
      dpy = 0;
      pxold = px;
      pyold = py;
      if (ipoint < 0) return;
      if (ipoint == 0) {
         px1old = 0;
         py1old = 0;
         px2old = gPad->XtoAbsPixel(fX[1]);
         py2old = gPad->YtoAbsPixel(fY[1]);
      } else if (ipoint == fN-1) {
         px1old = gPad->XtoAbsPixel(gPad->XtoPad(fX[fN-2]));
         py1old = gPad->YtoAbsPixel(gPad->YtoPad(fY[fN-2]));
         px2old = 0;
         py2old = 0;
      } else {
         px1old = gPad->XtoAbsPixel(gPad->XtoPad(fX[ipoint-1]));
         py1old = gPad->YtoAbsPixel(gPad->YtoPad(fY[ipoint-1]));
         px2old = gPad->XtoAbsPixel(gPad->XtoPad(fX[ipoint+1]));
         py2old = gPad->YtoAbsPixel(gPad->YtoPad(fY[ipoint+1]));
      }
      pxold = gPad->XtoAbsPixel(gPad->XtoPad(fX[ipoint]));
      pyold = gPad->YtoAbsPixel(gPad->YtoPad(fY[ipoint]));

      break;


   case kMouseMotion:

      middle = kTRUE;
      for (i=0;i<np;i++) {
         pxp = gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
         pyp = gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
         d   = TMath::Abs(pxp-px) + TMath::Abs(pyp-py);
         if (d < kMaxDiff) middle = kFALSE;
      }


   // check if point is close to an axis
      if (middle) gPad->SetCursor(kMove);
      else gPad->SetCursor(kHand);
      break;

   case kButton1Motion:
      if (!opaque) {
         if (middle) {
            for(i=0;i<np-1;i++) {
               gVirtualX->DrawLine(x[i]+dpx, y[i]+dpy, x[i+1]+dpx, y[i+1]+dpy);
               pxp = x[i]+dpx;
               pyp = y[i]+dpy;
               gVirtualX->DrawLine(pxp-4, pyp-4, pxp+4,  pyp-4);
               gVirtualX->DrawLine(pxp+4, pyp-4, pxp+4,  pyp+4);
               gVirtualX->DrawLine(pxp+4, pyp+4, pxp-4,  pyp+4);
               gVirtualX->DrawLine(pxp-4, pyp+4, pxp-4,  pyp-4);
            }
            pxp = x[np-1]+dpx;
            pyp = y[np-1]+dpy;
            gVirtualX->DrawLine(pxp-4, pyp-4, pxp+4,  pyp-4);
            gVirtualX->DrawLine(pxp+4, pyp-4, pxp+4,  pyp+4);
            gVirtualX->DrawLine(pxp+4, pyp+4, pxp-4,  pyp+4);
            gVirtualX->DrawLine(pxp-4, pyp+4, pxp-4,  pyp-4);
            dpx += px - pxold;
            dpy += py - pyold;
            pxold = px;
            pyold = py;
            for(i=0;i<np-1;i++) {
               gVirtualX->DrawLine(x[i]+dpx, y[i]+dpy, x[i+1]+dpx, y[i+1]+dpy);
               pxp = x[i]+dpx;
               pyp = y[i]+dpy;
               gVirtualX->DrawLine(pxp-4, pyp-4, pxp+4,  pyp-4);
               gVirtualX->DrawLine(pxp+4, pyp-4, pxp+4,  pyp+4);
               gVirtualX->DrawLine(pxp+4, pyp+4, pxp-4,  pyp+4);
               gVirtualX->DrawLine(pxp-4, pyp+4, pxp-4,  pyp-4);
            }
            pxp = x[np-1]+dpx;
            pyp = y[np-1]+dpy;
            gVirtualX->DrawLine(pxp-4, pyp-4, pxp+4,  pyp-4);
            gVirtualX->DrawLine(pxp+4, pyp-4, pxp+4,  pyp+4);
            gVirtualX->DrawLine(pxp+4, pyp+4, pxp-4,  pyp+4);
            gVirtualX->DrawLine(pxp-4, pyp+4, pxp-4,  pyp-4);
         } else {
            if (px1old) gVirtualX->DrawLine(px1old, py1old, pxold,  pyold);
            if (px2old) gVirtualX->DrawLine(pxold,  pyold,  px2old, py2old);
            gVirtualX->DrawLine(pxold-4, pyold-4, pxold+4,  pyold-4);
            gVirtualX->DrawLine(pxold+4, pyold-4, pxold+4,  pyold+4);
            gVirtualX->DrawLine(pxold+4, pyold+4, pxold-4,  pyold+4);
            gVirtualX->DrawLine(pxold-4, pyold+4, pxold-4,  pyold-4);
            pxold = px;
            pxold = TMath::Max(pxold, px1);
            pxold = TMath::Min(pxold, px2);
            pyold = py;
            pyold = TMath::Max(pyold, py2);
            pyold = TMath::Min(pyold, py1);
            if (px1old) gVirtualX->DrawLine(px1old, py1old, pxold,  pyold);
            if (px2old) gVirtualX->DrawLine(pxold,  pyold,  px2old, py2old);
            gVirtualX->DrawLine(pxold-4, pyold-4, pxold+4,  pyold-4);
            gVirtualX->DrawLine(pxold+4, pyold-4, pxold+4,  pyold+4);
            gVirtualX->DrawLine(pxold+4, pyold+4, pxold-4,  pyold+4);
            gVirtualX->DrawLine(pxold-4, pyold+4, pxold-4,  pyold-4);
         }
      } else {
         if (middle) {
            for(i=0;i<np-1;i++) {
               pxp = x[i]+dpx;
               pyp = y[i]+dpy;
            }
            pxp = x[np-1]+dpx;
            pyp = y[np-1]+dpy;
            dpx += px - pxold;
            dpy += py - pyold;
            pxold = px;
            pyold = py;
         } else {
            pxold = px;
            pxold = TMath::Max(pxold, px1);
            pxold = TMath::Min(pxold, px2);
            pyold = py;
            pyold = TMath::Max(pyold, py2);
            pyold = TMath::Min(pyold, py1);
         }
         if (x && y) {
            if (middle) {
               for(i=0;i<np;i++) {
                  fX[i] = gPad->PadtoX(gPad->AbsPixeltoX(x[i]+dpx));
                  fY[i] = gPad->PadtoY(gPad->AbsPixeltoY(y[i]+dpy));
               }
            } else {
               fX[ipoint] = gPad->PadtoX(gPad->AbsPixeltoX(pxold));
               fY[ipoint] = gPad->PadtoY(gPad->AbsPixeltoY(pyold));
            }
         }
         gPad->Modified(kTRUE);
      }
      break;

   case kButton1Up:

   // Compute x,y range
      xmin = gPad->GetUxmin();
      xmax = gPad->GetUxmax();
      ymin = gPad->GetUymin();
      ymax = gPad->GetUymax();
      dx   = xmax-xmin;
      dy   = ymax-ymin;
      dxr  = dx/(1 - gPad->GetLeftMargin() - gPad->GetRightMargin());
      dyr  = dy/(1 - gPad->GetBottomMargin() - gPad->GetTopMargin());

   // Range() could change the size of the pad pixmap and therefore should
   // be called before the other paint routines
         gPad->Range(xmin - dxr*gPad->GetLeftMargin(),
                     ymin - dyr*gPad->GetBottomMargin(),
                     xmax + dxr*gPad->GetRightMargin(),
                     ymax + dyr*gPad->GetTopMargin());
         gPad->RangeAxis(xmin, ymin, xmax, ymax);

      if (x && y) {
         if (middle) {
            for(i=0;i<np;i++) {
               fX[i] = gPad->PadtoX(gPad->AbsPixeltoX(x[i]+dpx));
               fY[i] = gPad->PadtoY(gPad->AbsPixeltoY(y[i]+dpy));
            }
         } else {
            fX[ipoint] = gPad->PadtoX(gPad->AbsPixeltoX(pxold));
            fY[ipoint] = gPad->PadtoY(gPad->AbsPixeltoY(pyold));
         }
         delete [] x; x = 0;
         delete [] y; y = 0;
      }
      gPad->Modified(kTRUE);
      gVirtualX->SetLineColor(-1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List this polyline with its attributes.
/// The option string is ignored.

void TPolyLine::ls(Option_t *) const
{
   TROOT::IndentLevel();
   printf("TPolyLine  N=%d\n",fN);
}

////////////////////////////////////////////////////////////////////////////////
/// Merge polylines in the collection in this polyline

Int_t TPolyLine::Merge(TCollection *li)
{
   if (!li) return 0;
   TIter next(li);

   //first loop to count the number of entries
   TPolyLine *pl;
   Int_t npoints = 0;
   while ((pl = (TPolyLine*)next())) {
      if (!pl->InheritsFrom(TPolyLine::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s",pl->ClassName(),this->ClassName());
         return -1;
      }
      npoints += pl->Size();
   }

   //extend this polyline to hold npoints
   if (npoints > 1) SetPoint(npoints-1,0,0);

   //merge all polylines
   next.Reset();
   while ((pl = (TPolyLine*)next())) {
      Int_t np = pl->Size();
      Double_t *x = pl->GetX();
      Double_t *y = pl->GetY();
      for (Int_t i=0;i<np;i++) {
         SetPoint(i,x[i],y[i]);
      }
   }

   return npoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this polyline with its current attributes.

void TPolyLine::Paint(Option_t *option)
{
   if (TestBit(kPolyLineNDC)) {
      if (strlen(option) > 0) PaintPolyLineNDC(fLastPoint+1, fX, fY, option);
      else                    PaintPolyLineNDC(fLastPoint+1, fX, fY, fOption.Data());
   } else {
      if (strlen(option) > 0) PaintPolyLine(fLastPoint+1, fX, fY, option);
      else                    PaintPolyLine(fLastPoint+1, fX, fY, fOption.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this polyline with new coordinates.
///
///  If option = 'f' or 'F' the fill area is drawn.
///  The default is to draw the lines only.

void TPolyLine::PaintPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option)
{
   if (!gPad) return;
   if (n <= 0) return;
   TAttLine::Modify();  //Change line attributes only if necessary
   TAttFill::Modify();  //Change fill area attributes only if necessary
   Double_t *xx = x;
   Double_t *yy = y;
   if (gPad->GetLogx()) {
      xx = new Double_t[n];
      for (Int_t ix=0;ix<n;ix++) xx[ix] = gPad->XtoPad(x[ix]);
   }
   if (gPad->GetLogy()) {
      yy = new Double_t[n];
      for (Int_t iy=0;iy<n;iy++) yy[iy] = gPad->YtoPad(y[iy]);
   }
   if (*option == 'f' || *option == 'F') gPad->PaintFillArea(n,xx,yy,option);
   else                                  gPad->PaintPolyLine(n,xx,yy,option);
   if (x != xx) delete [] xx;
   if (y != yy) delete [] yy;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this polyline with new coordinates in NDC.

void TPolyLine::PaintPolyLineNDC(Int_t n, Double_t *x, Double_t *y, Option_t *option)
{
   TAttLine::Modify();  //Change line attributes only if necessary
   TAttFill::Modify();  //Change fill area attributes only if necessary
   if (*option == 'f' || *option == 'F') gPad->PaintFillAreaNDC(n,x,y,option);
   else                                  gPad->PaintPolyLineNDC(n,x,y,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this polyline with its attributes.
/// The option string is ignored.

void TPolyLine::Print(Option_t *) const
{
   printf("PolyLine  N=%d\n",fN);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TPolyLine::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TPolyLine::Class())) {
      out<<"   ";
   } else {
      out<<"   Double_t *dum = 0;"<<std::endl;
      out<<"   TPolyLine *";
   }
   out<<"pline = new TPolyLine("<<fN<<",dum,dum,"<<quote<<fOption<<quote<<");"<<std::endl;

   SaveFillAttributes(out,"pline",0,1001);
   SaveLineAttributes(out,"pline",1,1,1);

   for (Int_t i=0;i<Size();i++) {
      out<<"   pline->SetPoint("<<i<<","<<fX[i]<<","<<fY[i]<<");"<<std::endl;
   }
   out<<"   pline->Draw("
      <<quote<<option<<quote<<");"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NDC mode on if isNDC = kTRUE, off otherwise

void TPolyLine::SetNDC(Bool_t isNDC)
{
   ResetBit(kPolyLineNDC);
   if (isNDC) SetBit(kPolyLineNDC);
}

////////////////////////////////////////////////////////////////////////////////
/// Set point following LastPoint to x, y.
/// Returns index of the point (new last point).

Int_t TPolyLine::SetNextPoint(Double_t x, Double_t y)
{
   fLastPoint++;
   SetPoint(fLastPoint, x, y);
   return fLastPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point number n to (x, y)
/// If n is greater than the current size, the arrays are automatically
/// extended.

void TPolyLine::SetPoint(Int_t n, Double_t x, Double_t y)
{
   if (n < 0) return;
   if (!fX || !fY || n >= fN) {
      // re-allocate the object
      Int_t newN = TMath::Max(2*fN,n+1);
      Double_t *savex = new Double_t [newN];
      Double_t *savey = new Double_t [newN];
      if (fX && fN){
         memcpy(savex,fX,fN*sizeof(Double_t));
         memset(&savex[fN],0,(newN-fN)*sizeof(Double_t));
         delete [] fX;
      }
      if (fY && fN){
         memcpy(savey,fY,fN*sizeof(Double_t));
         memset(&savey[fN],0,(newN-fN)*sizeof(Double_t));
         delete [] fY;
      }
      fX = savex;
      fY = savey;
      fN = newN;
   }
   fX[n] = x;
   fY[n] = y;
   fLastPoint = TMath::Max(fLastPoint,n);
}

////////////////////////////////////////////////////////////////////////////////
/// Resize this polyline to size n.
/// If n <= 0 the current arrays of points are deleted.
/// If n is greater than the current size, the new points are set to (0, 0)

void TPolyLine::SetPolyLine(Int_t n)
{
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      delete [] fX;
      delete [] fY;
      fX = fY = 0;
      return;
   }
   if (n < fN) {
      fN = n;
      fLastPoint = n - 1;
   } else {
      SetPoint(n-1,0,0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set new values for this polyline (single precision).
///
/// If n <= 0 the current arrays of points are deleted.

void TPolyLine::SetPolyLine(Int_t n, Float_t *x, Float_t *y, Option_t *option)
{
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      delete [] fX;
      delete [] fY;
      fX = fY = 0;
      return;
   }
   fN =n;
   if (fX) delete [] fX;
   if (fY) delete [] fY;
   fX = new Double_t[fN];
   fY = new Double_t[fN];
   for (Int_t i=0; i<fN;i++) {
      if (x) fX[i] = (Double_t)x[i];
      if (y) fY[i] = (Double_t)y[i];
   }
   fOption = option;
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set new values for this polyline (double precision).
///
/// If n <= 0 the current arrays of points are deleted.

void TPolyLine::SetPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option)
{
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      delete [] fX;
      delete [] fY;
      fX = fY = 0;
      return;
   }
   fN =n;
   if (fX) delete [] fX;
   if (fY) delete [] fY;
   fX = new Double_t[fN];
   fY = new Double_t[fN];
   for (Int_t i=0; i<fN;i++) {
      if (x) fX[i] = x[i];
      if (y) fY[i] = y[i];
   }
   fOption = option;
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TPolyLine::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         b.ReadClassBuffer(TPolyLine::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TObject::Streamer(b);
      TAttLine::Streamer(b);
      TAttFill::Streamer(b);
      b >> fN;
      fX = new Double_t[fN];
      fY = new Double_t[fN];
      Float_t *x = new Float_t[fN];
      Float_t *y = new Float_t[fN];
      b.ReadFastArray(x,fN);
      b.ReadFastArray(y,fN);
      for (Int_t i=0;i<fN;i++) {
         fX[i] = x[i];
         fY[i] = y[i];
      }
      fOption.Streamer(b);
      b.CheckByteCount(R__s, R__c, TPolyLine::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TPolyLine::Class(),this);
   }
}
