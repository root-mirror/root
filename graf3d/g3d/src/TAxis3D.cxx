// @(#)root/g3d:$Id$
// Author: Valery Fine(fine@mail.cern.ch)   07/01/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cctype>
#include <cassert>
#include <iostream>

#include "TMath.h"
#include "TList.h"
#include "TAxis3D.h"
#include "TClass.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TView.h"
#include "TVirtualX.h"
#include "TBrowser.h"
#include "TStyle.h"
#include "strlcpy.h"

/** \class TAxis3D
\ingroup g3d

The 3D axis painter class

This class provide up to 3 axes to any 3D ROOT plot and "ZOOM" service.
ExecuteEvent() method does provide zooming and moving a projection
3D object within TPad client area. With Zoom mode on the user can access
TAxis3D context menu and set /change the attributes of axes all together
or separately.

To add the 3D rulers to any 3D view one has to create
an instance of this class and Draw it.

~~~ {.cpp}
 TAxis3D rulers;
 rulers.Draw();
~~~

One can use a static method to create ruler and attach it to the current gPad

~~~ {.cpp}
 TAxis3D::ToggleRulers(); // Brings the 3D axes up
 TAxis3D::ToggleRulers(); // next calls remove the rulers from the TPad etc
~~~

To activate Zoomer one may call

~~~ {.cpp}
 TAxis3D::ToggleZoom();
~~~

each time one needs move or zoom the image. Then the user can:

  - move:
\image html g3d_axis3d_01.png
  - zoom:
\image html g3d_axis3d_02.png

its 3D view with <left-mouse button> press / move.
The "Zoom" deactivates itself just the user release the <left-mouse button>

To change attributes of the rulers attached to the current Pad, one may
query its pointer first:

~~~ {.cpp}
TAxis3D *axis = TAxis3D::GetPadAxis(); // Ask axis pointer
if (axis) {
  TAxis3D::ToggleRulers()     // To pop axes down
  axis->SetLabelColor(kBlue); // Paint the axes labels with blue color
  axis->SetAxisColor(kRed);   // Paint the axes itself with blue color
  TAxis3D::ToggleRulers()     // To pop axes up
}
~~~

The attributes of the created axes are affected by the current style
(see TStyle class ) and Set... methods of this class

For example:
~~~ {.cpp}
  gStyle->SetAxisColor(kYellow,"X");
  gStyle->SetAxisColor(kYellow,"Y");
  gStyle->SetAxisColor(kYellow,"Z");

  gStyle->SetLabelColor(kYellow,"X");
  gStyle->SetLabelColor(kYellow,"Y");
  gStyle->SetLabelColor(kYellow,"Z");

  TAxis3D::ToggleRulers();
  TAxis3D::ToggleRulers();
~~~

will draw all axes and labels with yellow color.
*/

const Char_t *TAxis3D::fgRulerName = "axis3druler";
ClassImp(TAxis3D);

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TAxis3D::TAxis3D() : TNamed(TAxis3D::fgRulerName,"ruler")
{
   fSelected = 0;
   fZoomMode = kFALSE;
   fStickyZoom = kFALSE;
   InitSet();
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TAxis3D::TAxis3D(Option_t *) : TNamed(TAxis3D::fgRulerName,"ruler")
{
   fSelected = 0;
   InitSet();
   fZoomMode = kFALSE;
   fStickyZoom = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy  constructor.

TAxis3D::TAxis3D(const TAxis3D &axis) : TNamed(axis)
{
   ((TAxis3D&)axis).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy axis3d

void TAxis3D::Copy(TObject &obj) const
{
   TNamed::Copy(obj);
   for (Int_t i=0;i<2;i++) fAxis[i].Copy(((TAxis3D&)obj).fAxis[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization.

void TAxis3D::InitSet()
{
   fAxis[0].SetName("xaxis");
   fAxis[1].SetName("yaxis");
   fAxis[2].SetName("zaxis");

   fAxis[0].Set(1,0.,1.);
   fAxis[1].Set(1,0.,1.);
   fAxis[2].Set(1,0.,1.);
   UseCurrentStyle();
}

////////////////////////////////////////////////////////////////////////////////
/// Add all 3 axes to the TBrowser

void TAxis3D::Browse(TBrowser *b)
{
   for (Int_t i=0;i<3;i++) b->Add(&fAxis[i],fAxis[i].GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a line.

Int_t TAxis3D::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t dist = 9;
   for (int i=0;i<3;i++) {
      Int_t axDist = fAxis[i].DistancetoPrimitive(px,py);
      if (dist > axDist) { dist = axDist; fSelected = &fAxis[i]; }
   }
   if (fZoomMode)
      return 0;
   else
      return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when an axis is clicked with the locator

void TAxis3D::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   if (fSelected) fSelected->ExecuteEvent(event,px,py);

   // Execute action corresponding to the mouse event

   static Double_t x0, y0, x1, y1;

   static Int_t pxold, pyold;
   static Int_t px0, py0;
   static Int_t linedrawn;

   if (!fZoomMode) return;

   // something to zoom ?

   gPad->SetCursor(kCross);

   switch (event) {

   case kButton1Down:
      gVirtualX->SetLineColor(-1);
      gPad->TAttLine::Modify();  //Change line attributes only if necessary
      ((TPad *)gPad)->AbsPixeltoXY(px,py,x0,y0);
      px0   = px; py0   = py;
      pxold = px; pyold = py;
      linedrawn = 0;
      break;

   case kButton1Motion:
      if (linedrawn) gVirtualX->DrawBox(px0, py0, pxold, pyold, TVirtualX::kHollow);
      pxold = px;
      pyold = py;
      linedrawn = 1;
      gVirtualX->DrawBox(px0, py0, pxold, pyold, TVirtualX::kHollow);
      break;

   case kButton1Up: {
      Int_t i;
      gPad->SetDoubleBuffer(1);
      gVirtualX->SetDrawMode(TVirtualX::kCopy); // set drawing mode back to normal (copy) mode
      TView *view = gPad->GetView();
      if (!view) break;                       // no 3D view yet

      Double_t min[3],max[3],viewCenter[3],viewCenterNDC[3];

      view->GetRange(min,max);
      for (i =0; i<3;i++) viewCenter[i] = (max[i]+min[i])/2;
      view->WCtoNDC(viewCenter,viewCenterNDC);
      // Define the center
      Double_t center[3],pointNDC[3],size[3],oldSize[3];
      ((TPad *)gPad)->AbsPixeltoXY(px,py,x1,y1);
      pointNDC[0] = (x0+x1)/2; pointNDC[1] = (y0+y1)/2;
      pointNDC[2] = viewCenterNDC[2];
      view->NDCtoWC(pointNDC, center);

      for (i =0; i<3;i++) oldSize[i] = size[i]= (max[i]-min[i])/2;

      // If there was a small motion, move the center only, do not change a scale
      if (TMath::Abs(px-px0)+TMath::Abs(py - py0) > 4 ) {
         Double_t newEdge[3];
         for (i =0; i<3;i++) size[i] = -1;

         pointNDC[0] = x0; pointNDC[1] = y0;

         view->NDCtoWC(pointNDC, newEdge);
         for (i =0; i<3;i++) {
            Double_t newSize = TMath::Abs(newEdge[i]-center[i]);
            if ( newSize/oldSize[i] > 0.002)
               size[i] = TMath::Max(size[i], newSize);
            else
               size[i] = oldSize[i];
         }

         pointNDC[0] = x1; pointNDC[1] = y1;

         view->NDCtoWC(pointNDC, newEdge);
         for (i =0; i<3;i++) {
            Double_t newSize = TMath::Abs(newEdge[i]-center[i]);
            if ( newSize/oldSize[i] > 0.002)
               size[i] = TMath::Max(size[i], newSize);
            else
               size[i] = oldSize[i];
         }
#if 0
         if (fZooms == kMAXZOOMS) fZoom = 0;
         fZooms++;
         memcpy(fZoomMin[fZooms],min,3*sizeof(Float_t));
         memcpy(fZoomMax[fZooms],max,3*sizeof(Float_t));
#endif
      }
      for (i =0; i<3;i++) {
         max[i] = center[i] + size[i];
         min[i] = center[i] - size[i];
      }
      view->SetRange(min,max);

      if(!fStickyZoom)SwitchZoom();
         gPad->Modified(kTRUE);
         gPad->Update();
         break;
      }
      default: break;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Dummy method
///  returns the const char * to "axis3d"

char *TAxis3D::GetObjectInfo(Int_t , Int_t ) const
{
   return (char*)"axis3d";
}

////////////////////////////////////////////////////////////////////////////////
/// Paint axis over 3D view on the TPad

void TAxis3D::Paint(Option_t *)
{
   TGaxis axis;
   PaintAxis(&axis, 90);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the axis for TView object.
///
/// The original idea belongs:
///
/// void THistPainter::PaintLegoAxis(TGaxis *axis, Double_t ang)

void TAxis3D::PaintAxis(TGaxis *axis, Float_t ang)
{
   static Double_t epsil = 0.001;

   Double_t cosa, sina;
   Double_t bmin, bmax;
   Double_t r[24]       /* was [3][8] */;
   Int_t ndiv, i;
   Double_t x1[3], x2[3], y1[3], y2[3], z1[3], z2[3], av[24]  /*  was [3][8] */;
   char chopax[10];
   Int_t ix1, ix2, iy1, iy2, iz1, iz2;
   Double_t rad;

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintAxis", "no TView in current pad");
      return;
   }

   rad  = TMath::ATan(1.) * 4. / 180.;
   cosa = TMath::Cos(ang*rad);
   sina = TMath::Sin(ang*rad);

   view->AxisVertex(ang, av, ix1, ix2, iy1, iy2, iz1, iz2);
   for (i = 1; i <= 8; ++i) {
      r[i*3 - 3] = av[i*3 - 3] + av[i*3 - 2]*cosa;
      r[i*3 - 2] = av[i*3 - 2]*sina;
      r[i*3 - 1] = av[i*3 - 1];
   }

   view->WCtoNDC(&r[ix1*3 - 3], x1);
   view->WCtoNDC(&r[ix2*3 - 3], x2);
   view->WCtoNDC(&r[iy1*3 - 3], y1);
   view->WCtoNDC(&r[iy2*3 - 3], y2);
   view->WCtoNDC(&r[iz1*3 - 3], z1);
   view->WCtoNDC(&r[iz2*3 - 3], z2);

   view->SetAxisNDC(x1, x2, y1, y2, z1, z2);

   Double_t *rmin = view->GetRmin();
   Double_t *rmax = view->GetRmax();

   axis->SetLineWidth(1);

   for (i=0;i<3;i++) {

      // X axis drawing
      Double_t ax[2], ay[2];
      Bool_t logAx = kFALSE;
      memset(chopax,0,sizeof(chopax));
      switch (i) {
         case 0 :
            ax[0] = x1[0]; ax[1] = x2[0];
            ay[0] = x1[1]; ay[1] = x2[1];
            logAx = gPad->GetLogx();
            break;
         case 1 :
            if (TMath::Abs(y1[0] - y2[0]) < epsil)  y2[0] = y1[0];
            ax[0] = y1[0]; ax[1] = y2[0];
            ay[0] = y1[1]; ay[1] = y2[1];
            logAx = gPad->GetLogy();
            break;
         case 2 :
            ax[0] = z1[0]; ax[1] = z2[0];
            ay[0] = z1[1]; ay[1] = z2[1];
            strlcpy(chopax, "SDH+=",10);
            logAx = gPad->GetLogz();
            break;
         default:
            assert(0);
            continue;
      };

      // If the axis is too short - skip it
      if ( ( TMath::Abs(ax[0] - ax[1]) + TMath::Abs(ay[0] - ay[1]))  < epsil  ) continue;

      if (i != 2 ) {
         if (ax[0] > ax[1]) strlcpy(chopax, "SDHV=+",10);
         else               strlcpy(chopax, "SDHV=-",10);
      }

      if (i==1 && (TMath::Abs(z1[0] - z2[0]) + TMath::Abs(z1[1] - z2[1])) < epsil)
         strlcpy(chopax, "SDH+=",10);

      // Initialize the axis options
      if (logAx) {
         strlcat(chopax,"G",10);
         bmin = TMath::Power(10, rmin[i]);
         bmax = TMath::Power(10, rmax[i]);
      } else {
         bmin = rmin[i];
         bmax = rmax[i];
      }

      axis->SetLineColor(  fAxis[i].GetAxisColor());
      axis->SetTextFont(   fAxis[i].GetTitleFont());
      axis->SetTextColor(  fAxis[i].GetTitleColor());
      axis->SetTickSize(   fAxis[i].GetTickLength());
      axis->SetLabelColor( fAxis[i].GetLabelColor());
      axis->SetLabelFont(  fAxis[i].GetLabelFont());
      axis->SetLabelOffset(fAxis[i].GetLabelOffset()+fAxis[i].GetTickLength());
      axis->SetLabelSize(  fAxis[i].GetLabelSize());
      axis->SetTitle(      fAxis[i].GetTitle());
      axis->SetTitleOffset(fAxis[i].GetTitleOffset());
      axis->SetTitleSize(  fAxis[i].GetTitleSize());
      enum { kCenterTitle = BIT(12) }; // to be removed with the last version of ROOT
      axis->SetBit(kCenterTitle, fAxis[i].TestBit(kCenterTitle));

      //*-*-    Initialize the number of divisions. If the
      //*-*-    number of divisions is negative, option 'N' is required.
      ndiv = fAxis[i].GetNdivisions();
      if (ndiv < 0) {
         ndiv = -ndiv;
         chopax[6] = 'N';
      }

      // Option time display is required ?
      if (fAxis[i].GetTimeDisplay()) {
         strlcat(chopax,"t",10);
         if (strlen(fAxis[i].GetTimeFormatOnly()) == 0) {
            axis->SetTimeFormat(fAxis[i].ChooseTimeFormat(bmax-bmin));
         } else {
            axis->SetTimeFormat(fAxis[i].GetTimeFormat());
         }
      }
      axis->SetOption(chopax);
      axis->PaintAxis(ax[0], ay[0], ax[1], ay[1], bmin, bmax, ndiv, chopax);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert "screen pixel" coordinates to some center of 3D WC coordinate
/// if view and gPad present

Double_t *TAxis3D::PixeltoXYZ(Double_t px, Double_t py, Double_t *point3D, TView *view)
{
   Double_t *thisPoint = 0;
   if (!view && gPad) view = gPad->GetView();
   if (view) {
      Double_t x[3] = {px,py,0.5}; // ((TPad *)thisPad)->AbsPixeltoXY(px,py,x[0],x[1]);
      Double_t min[3], max[3];
      view->GetRange(min,max);
      Int_t i;
      for (i =0; i<3;i++) min[i] = (max[i]+min[i])/2;
      view->WCtoNDC(min,max);
      min[0] = x[0]; min[1] = x[1];
      min[2] = max[2];
      view->NDCtoWC(min, x);
      for (i=0;i<3;i++) point3D[i] = x[i];
      thisPoint = point3D;
   }
   return thisPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TAxis3D::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   fAxis[0].SaveAttributes(out,GetName(),"->GetXaxis()");
   fAxis[1].SaveAttributes(out,GetName(),"->GetYaxis()");
   fAxis[2].SaveAttributes(out,GetName(),"->GetZaxis()");
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current attributes by current style.

void TAxis3D::UseCurrentStyle()
{
   if (gStyle->IsReading()) {
      fAxis[0].ResetAttAxis("X");
      fAxis[1].ResetAttAxis("Y");
      fAxis[2].ResetAttAxis("Z");

      fAxis[0].SetTitle("x");
      fAxis[0].SetLabelColor(kRed);   fAxis[0].SetAxisColor(kRed);
      fAxis[1].SetLabelColor(kGreen); fAxis[1].SetAxisColor(kGreen);
      fAxis[2].SetLabelColor(kBlue);  fAxis[2].SetAxisColor(kBlue);
   } else {
      gStyle->SetNdivisions (fAxis[0].GetNdivisions(), "x");
      gStyle->SetAxisColor  (fAxis[0].GetAxisColor(),  "x");
      gStyle->SetLabelColor (fAxis[0].GetLabelColor(), "x");
      gStyle->SetLabelFont  (fAxis[0].GetLabelFont(),  "x");
      gStyle->SetLabelOffset(fAxis[0].GetLabelOffset(),"x");
      gStyle->SetLabelSize  (fAxis[0].GetLabelSize(),  "x");
      gStyle->SetTickLength (fAxis[0].GetTickLength(), "x");
      gStyle->SetTitleOffset(fAxis[0].GetTitleOffset(),"x");
      gStyle->SetTitleSize  (fAxis[0].GetTitleSize(),  "x");
      gStyle->SetTitleColor (fAxis[0].GetTitleColor(), "x");
      gStyle->SetTitleFont  (fAxis[0].GetTitleFont(),  "x");

      gStyle->SetNdivisions (fAxis[1].GetNdivisions(), "y");
      gStyle->SetAxisColor  (fAxis[1].GetAxisColor(),  "y");
      gStyle->SetLabelColor (fAxis[1].GetLabelColor(), "y");
      gStyle->SetLabelFont  (fAxis[1].GetLabelFont(),  "y");
      gStyle->SetLabelOffset(fAxis[1].GetLabelOffset(),"y");
      gStyle->SetLabelSize  (fAxis[1].GetLabelSize(),  "y");
      gStyle->SetTickLength (fAxis[1].GetTickLength(), "y");
      gStyle->SetTitleOffset(fAxis[1].GetTitleOffset(),"y");
      gStyle->SetTitleSize  (fAxis[1].GetTitleSize(),  "y");
      gStyle->SetTitleColor (fAxis[1].GetTitleColor(), "y");
      gStyle->SetTitleFont  (fAxis[1].GetTitleFont(),  "y");

      gStyle->SetNdivisions (fAxis[2].GetNdivisions(), "z");
      gStyle->SetAxisColor  (fAxis[2].GetAxisColor(),  "z");
      gStyle->SetLabelColor (fAxis[2].GetLabelColor(), "z");
      gStyle->SetLabelFont  (fAxis[2].GetLabelFont(),  "z");
      gStyle->SetLabelOffset(fAxis[2].GetLabelOffset(),"z");
      gStyle->SetLabelSize  (fAxis[2].GetLabelSize(),  "z");
      gStyle->SetTickLength (fAxis[2].GetTickLength(), "z");
      gStyle->SetTitleOffset(fAxis[2].GetTitleOffset(),"z");
      gStyle->SetTitleSize  (fAxis[2].GetTitleSize(),  "z");
      gStyle->SetTitleColor (fAxis[2].GetTitleColor(), "z");
      gStyle->SetTitleFont  (fAxis[2].GetTitleFont(),  "z");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the axis index by its name

Int_t TAxis3D::AxisChoice( Option_t *axis) const
{
   char achoice = toupper(axis[0]);
   if (achoice == 'X') return 0;
   if (achoice == 'Y') return 1;
   if (achoice == 'Z') return 2;
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of divisions.

Int_t TAxis3D::GetNdivisions( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetNdivisions();
}

////////////////////////////////////////////////////////////////////////////////
/// Get axis color.

Color_t TAxis3D::GetAxisColor( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetAxisColor();
}

////////////////////////////////////////////////////////////////////////////////
/// Get label color.

Color_t TAxis3D::GetLabelColor( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetLabelColor();
}

////////////////////////////////////////////////////////////////////////////////
/// Get label font.

Style_t TAxis3D::GetLabelFont( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetLabelFont();
}

////////////////////////////////////////////////////////////////////////////////
/// Get label offset.

Float_t TAxis3D::GetLabelOffset( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetLabelOffset();
}

////////////////////////////////////////////////////////////////////////////////
/// Get label size.

Float_t TAxis3D::GetLabelSize( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetLabelSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Get tick mark length.

Float_t TAxis3D::GetTickLength( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetTickLength();
}

////////////////////////////////////////////////////////////////////////////////
/// Get title offset.

Float_t TAxis3D::GetTitleOffset( Option_t *axis) const
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return 0;
   return fAxis[ax].GetTitleOffset();
}

////////////////////////////////////////////////////////////////////////////////

#define AXISCHOICE                \
   Int_t i = AxisChoice(axis);    \
   Int_t nax = 1;                 \
   if (i == -1) { i = 0; nax = 3;}\
   for (Int_t ax=i;ax<nax+i;ax++)

////////////////////////////////////////////////////////////////////////////////
/// Set number of divisions.

void TAxis3D::SetNdivisions(Int_t n, Option_t *axis)
{
   AXISCHOICE {fAxis[ax].SetNdivisions(n);}
}

////////////////////////////////////////////////////////////////////////////////
/// Set axis color.

void TAxis3D::SetAxisColor(Color_t color, Option_t *axis)
{
   AXISCHOICE {fAxis[ax].SetAxisColor(color);}
}

////////////////////////////////////////////////////////////////////////////////
/// Set axis range.

void TAxis3D::SetAxisRange(Double_t xmin, Double_t xmax, Option_t *axis)
{
   Int_t ax = AxisChoice(axis);
   if (ax < 0) return;
   TAxis *theAxis = &fAxis[ax];
   Int_t bin1 = theAxis->FindBin(xmin);
   Int_t bin2 = theAxis->FindBin(xmax);
   theAxis->SetRange(bin1, bin2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set label color.

void TAxis3D::SetLabelColor(Color_t color, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetLabelColor(color); }
}

////////////////////////////////////////////////////////////////////////////////
/// Set label font.

void TAxis3D::SetLabelFont(Style_t font, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetLabelFont(font); }
}

////////////////////////////////////////////////////////////////////////////////
/// Set label offset.

void TAxis3D::SetLabelOffset(Float_t offset, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetLabelOffset(offset); }
}

////////////////////////////////////////////////////////////////////////////////
/// Set label size.

void TAxis3D::SetLabelSize(Float_t size, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetLabelSize(size); }
}

////////////////////////////////////////////////////////////////////////////////
/// Set tick mark length.

void TAxis3D::SetTickLength(Float_t length, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetTickLength(length); }
}

////////////////////////////////////////////////////////////////////////////////
/// Set title offset.

void TAxis3D::SetTitleOffset(Float_t offset, Option_t *axis)
{
   AXISCHOICE { fAxis[ax].SetTitleOffset(offset); }
}
#undef AXISCHOICE

////////////////////////////////////////////////////////////////////////////////
/// Returns the "pad" Axis3D object pointer if any.

TAxis3D *TAxis3D::GetPadAxis(TVirtualPad *pad)
{
   TObject *obj = 0;
   TVirtualPad *thisPad=pad;
   if (!thisPad) thisPad = gPad;
   if (thisPad) {
      // Find axis in the current thisPad
      obj = thisPad->FindObject(TAxis3D::fgRulerName);
      if (!(obj && obj->InheritsFrom(Class()->GetName()))) obj = 0;
   }
   return (TAxis3D *)obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn ON / OFF the "Ruler", TAxis3D object attached
/// to the current pad

TAxis3D *TAxis3D::ToggleRulers(TVirtualPad *pad)
{
   TAxis3D *ax = 0;
   TVirtualPad *thisPad=pad;
   if (!thisPad) thisPad = gPad;
   if (thisPad && thisPad->GetView() ) {
      TAxis3D *a =  GetPadAxis(pad);
      if (a)  delete a;
      else {
         ax = new TAxis3D;
         ax->SetBit(kCanDelete);
         ax->Draw();
      }
      thisPad->Modified();
      thisPad->Update();
   }
   return ax;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn ON / OFF the "Ruler" and "zoom mode" of the TAxis3D object attached
/// to the current pad (if pad = 0; gPad is used "by default")
///
/// User is given a chance to either:
///  1.  move the center of the 3D scene at the cursor position
///  2.  zoom view with mouse "drugging" the bounder rectangle with "left" mouse
///  3.  Change the axis attributes via TContextMenu with "right mouse button click"

TAxis3D *TAxis3D::ToggleZoom(TVirtualPad *pad)
{
   TAxis3D *ax = 0;
   TVirtualPad *thisPad=pad;
   if (!thisPad) thisPad = gPad;
   if (thisPad && thisPad->GetView()) {
      // Find axis in the current thisPad
      TList *l = thisPad->GetListOfPrimitives();
      TObject *o = l->FindObject(TAxis3D::fgRulerName);
      if (o && o->InheritsFrom(Class()->GetName())) { // Find axis
         if (o != l->Last()) { // make sure the TAxis3D is the last object of the Pad.
            l->Remove(o);
            l->AddLast(o);
         }
         ax = (TAxis3D *)o;
      } else { // There is no
         ax = new TAxis3D;
         ax->SetBit(kCanDelete);
         ax->Draw();
      }
      ax->SwitchZoom();
   }
   return ax;
}
