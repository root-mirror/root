// @(#)root/histpainter:$Id$
// Author: Rene Brun   15/11/2002

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TPaletteAxis.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TStyle.h"
#include "TMath.h"
#include "TView.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "snprintf.h"

#include <iostream>

ClassImp(TPaletteAxis);


////////////////////////////////////////////////////////////////////////////////

/*! \class TPaletteAxis
\ingroup Histpainter
\brief The palette painting class.

A `TPaletteAxis` object is used to display the color palette when
drawing 2-d histograms.

The `TPaletteAxis` is automatically created drawn when drawing a 2-D
histogram when the option "Z" is specified.

A `TPaletteAxis` object is added to the histogram list of functions and
can be retrieved doing:

    TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");

then the pointer `palette` can be used to change the palette attributes.

Because the palette is created at painting time only, one must issue a:

    gPad->Update();

before retrieving the palette pointer in order to create the palette. The following
macro gives an example.

Begin_Macro(source)
{
   TCanvas *c1 = new TCanvas("c1","c1",600,400);
   TH2F *h2 = new TH2F("h2","Example of a resized palette ",40,-4,4,40,-20,20);
   Float_t px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      h2->Fill(px,5*py);
   }
   gStyle->SetPalette(1);
   h2->Draw("COLZ");
   gPad->Update();
   TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
   palette->SetY2NDC(0.7);
   return c1;
}
End_Macro

`TPaletteAxis` inherits from `TBox` and `TPave`. The methods
allowing to specify the palette position are inherited from these two classes.

The palette can be interactively moved and resized. The context menu
can be used to set the axis attributes.

It is possible to select a range on the axis to set the min/max in z

As default labels and ticks are drawn by `TGAxis` at equidistant (lin or log)
points as controlled by SetNdivisions.
If option "CJUST" is given labels and ticks are justified at the
color boundaries defined by the contour levels.
In this case no optimization can be done. It is responsiblity of the
user to adjust minimum, maximum of the histogram and/or the contour levels
to get a reasonable look of the plot.
Only overlap of the labels is avoided if too many contour levels are used.

This option is especially useful with user defined contours.
An example is shown here:
Begin_Macro(source)
{
   gStyle->SetOptStat(0);
   TCanvas *c1 = new TCanvas("c1","exa_CJUST",300,10,400,400);
   TH2F *hpxpy = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
   // Fill histograms randomly
   TRandom3 randomNum;
   Float_t px, py;
   for (Int_t i = 0; i < 25000; i++) {
      randomNum.Rannor(px,py);
      hpxpy->Fill(px,py);
   }
   hpxpy->SetMaximum(200);
   Double_t zcontours[5] = {0, 20, 40, 80, 120};
   hpxpy->SetContour(5, zcontours);
   hpxpy->GetZaxis()->SetTickSize(0.01);
   hpxpy->GetZaxis()->SetLabelOffset(0.01);
   gPad->SetRightMargin(0.13);
   hpxpy->SetTitle("User contours, CJUST");
   hpxpy->Draw("COL Z CJUST");
}
End_Macro
*/


////////////////////////////////////////////////////////////////////////////////
/// Palette default constructor.

TPaletteAxis::TPaletteAxis(): TPave()
{
   fH  = 0;
   SetName("");
}


////////////////////////////////////////////////////////////////////////////////
/// Palette normal constructor.

TPaletteAxis::TPaletteAxis(Double_t x1, Double_t y1, Double_t x2, Double_t  y2, TH1 *h)
   : TPave(x1, y1, x2, y2)
{
   fH = h;
   SetName("palette");
   TAxis *zaxis = fH->GetZaxis();
   fAxis.ImportAxisAttributes(zaxis);
   if (gPad->GetView()) SetBit(kHasView);
}


////////////////////////////////////////////////////////////////////////////////
/// Palette destructor.

TPaletteAxis::~TPaletteAxis()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Palette copy constructor.

TPaletteAxis::TPaletteAxis(const TPaletteAxis &palette) : TPave(palette)
{
   ((TPaletteAxis&)palette).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TPaletteAxis& TPaletteAxis::operator=(const TPaletteAxis &orig)
{
   orig.Copy( *this );
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy a palette to a palette.

void TPaletteAxis::Copy(TObject &obj) const
{
   TPave::Copy(obj);
   ((TPaletteAxis&)obj).fH    = fH;
   ((TPaletteAxis&)obj).fName = fName;
}


////////////////////////////////////////////////////////////////////////////////
/// Check if mouse on the axis region.

Int_t TPaletteAxis::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t plxmax = gPad->XtoAbsPixel(fX2);
   Int_t plymin = gPad->YtoAbsPixel(fY1);
   Int_t plymax = gPad->YtoAbsPixel(fY2);
   if (px > plxmax && px < plxmax + 30 && py >= plymax && py <= plymin) return px - plxmax;

   //otherwise check if inside the box
   return TPave::DistancetoPrimitive(px, py);
}


////////////////////////////////////////////////////////////////////////////////
/// Check if mouse on the axis region.

void TPaletteAxis::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   static Int_t kmode = 0;
   Int_t plxmin = gPad->XtoAbsPixel(fX1);
   Int_t plxmax = gPad->XtoAbsPixel(fX2);
   if (kmode != 0 || px <= plxmax) {
      if (event == kButton1Down) kmode = 1;
      TBox::ExecuteEvent(event, px, py);
      if (event == kButton1Up) kmode = 0;
      // In case palette coordinates have been modified, recompute NDC coordinates
      Double_t dpx  = gPad->GetX2() - gPad->GetX1();
      Double_t dpy  = gPad->GetY2() - gPad->GetY1();
      Double_t xp1  = gPad->GetX1();
      Double_t yp1  = gPad->GetY1();
      fX1NDC = (fX1 - xp1) / dpx;
      fY1NDC = (fY1 - yp1) / dpy;
      fX2NDC = (fX2 - xp1) / dpx;
      fY2NDC = (fY2 - yp1) / dpy;
      return;
   }
   gPad->SetCursor(kHand);
   static Double_t ratio1, ratio2;
   static Int_t px1old, py1old, px2old, py2old;
   Double_t temp, xmin, xmax;

   switch (event) {

      case kButton1Down:
         ratio1 = (gPad->AbsPixeltoY(py) - fY1) / (fY2 - fY1);
         py1old = gPad->YtoAbsPixel(fY1 + ratio1 * (fY2 - fY1));
         px1old = plxmin;
         px2old = plxmax;
         py2old = py1old;
         gVirtualX->DrawBox(px1old, py1old, px2old, py2old, TVirtualX::kHollow);
         gVirtualX->SetLineColor(-1);
         // No break !!!

      case kButton1Motion:
         gVirtualX->DrawBox(px1old, py1old, px2old, py2old, TVirtualX::kHollow);
         ratio2 = (gPad->AbsPixeltoY(py) - fY1) / (fY2 - fY1);
         py2old = gPad->YtoAbsPixel(fY1 + ratio2 * (fY2 - fY1));
         gVirtualX->DrawBox(px1old, py1old, px2old, py2old, TVirtualX::kHollow);
         break;

      case kButton1Up:
         if (gROOT->IsEscaped()) {
            gROOT->SetEscape(kFALSE);
            break;
         }

         ratio2 = (gPad->AbsPixeltoY(py) - fY1) / (fY2 - fY1);
         xmin = ratio1;
         xmax = ratio2;
         if (xmin > xmax) {
            temp   = xmin;
            xmin   = xmax;
            xmax   = temp;
            temp   = ratio1;
            ratio1 = ratio2;
            ratio2 = temp;
         }
         if (ratio2 - ratio1 > 0.05) {
            if (fH->GetDimension() == 2) {
               Double_t zmin = fH->GetMinimum();
               Double_t zmax = fH->GetMaximum();
               if (gPad->GetLogz()) {
                  if (zmin <= 0 && zmax > 0) zmin = TMath::Min((Double_t)1,
                                                                  (Double_t)0.001 * zmax);
                  zmin = TMath::Log10(zmin);
                  zmax = TMath::Log10(zmax);
               }
               Double_t newmin = zmin + (zmax - zmin) * ratio1;
               Double_t newmax = zmin + (zmax - zmin) * ratio2;
               if (newmin < zmin)newmin = fH->GetBinContent(fH->GetMinimumBin());
               if (newmax > zmax)newmax = fH->GetBinContent(fH->GetMaximumBin());
               if (gPad->GetLogz()) {
                  newmin = TMath::Exp(2.302585092994 * newmin);
                  newmax = TMath::Exp(2.302585092994 * newmax);
               }
               fH->SetMinimum(newmin);
               fH->SetMaximum(newmax);
               fH->SetBit(TH1::kIsZoomed);
            }
            gPad->Modified(kTRUE);
         }
         gVirtualX->SetLineColor(-1);
         kmode = 0;
         break;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the color index of the bin (i,j).
///
/// This function should be used after an histogram has been plotted with the
/// option COL or COLZ like in the following example:
///
///     h2->Draw("COLZ");
///     gPad->Update();
///     TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
///     Int_t ci = palette->GetBinColor(20,15);
///
/// Then it is possible to retrieve the RGB components in the following way:
///
///     TColor *c = gROOT->GetColor(ci);
///     float x,y,z;
///     c->GetRGB(x,y,z);

Int_t TPaletteAxis::GetBinColor(Int_t i, Int_t j)
{
   Double_t zc = fH->GetBinContent(i, j);
   return GetValueColor(zc);
}


////////////////////////////////////////////////////////////////////////////////
/// Displays the z value corresponding to cursor position py.

char *TPaletteAxis::GetObjectInfo(Int_t /* px */, Int_t py) const
{
   Double_t z;
   static char info[64];

   Double_t zmin = fH->GetMinimum();
   Double_t zmax = fH->GetMaximum();
   Int_t   y1   = gPad->GetWh() - gPad->VtoPixel(fY1NDC);
   Int_t   y2   = gPad->GetWh() - gPad->VtoPixel(fY2NDC);
   Int_t   y    = gPad->GetWh() - py;

   if (gPad->GetLogz()) {
      if (zmin <= 0 && zmax > 0) zmin = TMath::Min((Double_t)1,
                                                      (Double_t)0.001 * zmax);
      Double_t zminl = TMath::Log10(zmin);
      Double_t zmaxl = TMath::Log10(zmax);
      Double_t zl    = (zmaxl - zminl) * ((Double_t)(y - y1) / (Double_t)(y2 - y1)) + zminl;
      z = TMath::Power(10., zl);
   } else {
      z = (zmax - zmin) * ((Double_t)(y - y1) / (Double_t)(y2 - y1)) + zmin;
   }

   snprintf(info, 64, "(z=%g)", z);
   return info;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the color index of the given z value
///
/// This function should be used after an histogram has been plotted with the
/// option COL or COLZ like in the following example:
///
///     h2->Draw("COLZ");
///     gPad->Update();
///     TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
///     Int_t ci = palette->GetValueColor(30.);
///
/// Then it is possible to retrieve the RGB components in the following way:
///
///     TColor *c = gROOT->GetColor(ci);
///     float x,y,z;
///     c->GetRGB(x,y,z);

Int_t TPaletteAxis::GetValueColor(Double_t zc)
{
   Double_t wmin  = fH->GetMinimum();
   Double_t wmax  = fH->GetMaximum();
   Double_t wlmin = wmin;
   Double_t wlmax = wmax;

   if (gPad->GetLogz()) {
      if (wmin <= 0 && wmax > 0) wmin = TMath::Min((Double_t)1,
                                                      (Double_t)0.001 * wmax);
      wlmin = TMath::Log10(wmin);
      wlmax = TMath::Log10(wmax);
   }

   Int_t ncolors = gStyle->GetNumberOfColors();
   Int_t ndivz = fH->GetContour();
   if (ndivz == 0) return 0;
   ndivz = TMath::Abs(ndivz);
   Int_t theColor, color;
   Double_t scale = ndivz / (wlmax - wlmin);

   if (fH->TestBit(TH1::kUserContour) && gPad->GetLogz()) zc = TMath::Log10(zc);
   if (zc < wlmin) zc = wlmin;

   color = Int_t(0.01 + (zc - wlmin) * scale);

   theColor = Int_t((color + 0.99) * Double_t(ncolors) / Double_t(ndivz));
   return gStyle->GetColorPalette(theColor);
}


////////////////////////////////////////////////////////////////////////////////
/// Paint the palette.

void TPaletteAxis::Paint(Option_t *)
{
   ConvertNDCtoPad();

   SetFillStyle(1001);
   Double_t ymin = fY1;
   Double_t ymax = fY2;
   Double_t xmin = fX1;
   Double_t xmax = fX2;
   Double_t wmin = fH->GetMinimum();
   Double_t wmax = fH->GetMaximum();
   Double_t wlmin = wmin;
   Double_t wlmax = wmax;
   Double_t y1, y2, w1, w2, zc;

   if ((wlmax - wlmin) <= 0) {
      Double_t mz = wlmin * 0.1;
      if (mz == 0) mz = 0.1;
      wlmin = wlmin - mz;
      wlmax = wlmax + mz;
      wmin  = wlmin;
      wmax  = wlmax;
   }

   if (gPad->GetLogz()) {
      if (wmin <= 0 && wmax > 0) wmin = TMath::Min((Double_t)1,
                                                      (Double_t)0.001 * wmax);
      wlmin = TMath::Log10(wmin);
      wlmax = TMath::Log10(wmax);
   }
   Double_t ws    = wlmax - wlmin;
   Int_t ncolors = gStyle->GetNumberOfColors();
   Int_t ndivz = fH->GetContour();
   if (ndivz == 0) return;
   ndivz = TMath::Abs(ndivz);
   Int_t theColor, color;
   // import Attributes already here since we might need them for CJUST
   if (fH->GetDimension() == 2) fAxis.ImportAxisAttributes(fH->GetZaxis());
   // case option "CJUST": put labels directly at color boundaries
   TLatex *label = NULL;
   TLine *line = NULL;
   Double_t prevlab = 0;
   TString opt(fH->GetDrawOption());
   if (opt.Contains("CJUST", TString::kIgnoreCase)) {
      label = new TLatex();
      label->SetTextFont(fAxis.GetLabelFont());
      label->SetTextColor(fAxis.GetLabelColor());
      label->SetTextAlign(kHAlignLeft+kVAlignCenter);
      line = new TLine();
      line->SetLineColor(fAxis.GetLineColor());
      line->PaintLine(xmax, ymin, xmax, ymax);
   }
   Double_t scale = ndivz / (wlmax - wlmin);
   for (Int_t i = 0; i < ndivz; i++) {

      zc = fH->GetContourLevel(i);
      if (fH->TestBit(TH1::kUserContour) && gPad->GetLogz())
         zc = TMath::Log10(zc);
      w1 = zc;
      if (w1 < wlmin) w1 = wlmin;

      w2 = wlmax;
      if (i < ndivz - 1) {
         zc = fH->GetContourLevel(i + 1);
         if (fH->TestBit(TH1::kUserContour) && gPad->GetLogz())
            zc = TMath::Log10(zc);
         w2 = zc;
      }

      if (w2 <= wlmin) continue;
      y1 = ymin + (w1 - wlmin) * (ymax - ymin) / ws;
      y2 = ymin + (w2 - wlmin) * (ymax - ymin) / ws;

      if (fH->TestBit(TH1::kUserContour)) {
         color = i;
      } else {
         color = Int_t(0.01 + (w1 - wlmin) * scale);
      }

      theColor = Int_t((color + 0.99) * Double_t(ncolors) / Double_t(ndivz));
      SetFillColor(gStyle->GetColorPalette(theColor));
      TAttFill::Modify();
      gPad->PaintBox(xmin, y1, xmax, y2);
      // case option "CJUST": put labels directly
      if (label) {
         Double_t lof = fAxis.GetLabelOffset()*(gPad->GetUxmax()-gPad->GetUxmin());
         // the following assumes option "S"
         Double_t tlength = fAxis.GetTickSize() * (gPad->GetUxmax()-gPad->GetUxmin());
         Double_t lsize = fAxis.GetLabelSize();
         Double_t lsize_user = lsize*(gPad->GetUymax()-gPad->GetUymin());
         Double_t zlab = fH->GetContourLevel(i);
         if (gPad->GetLogz()&& !fH->TestBit(TH1::kUserContour)) {
            zlab = TMath::Power(10, zlab);
         }
         // make sure labels dont overlap
         if (i == 0 || (y1 - prevlab) > 1.5*lsize_user) {
            label->PaintLatex(xmax + lof, y1, 0, lsize, Form("%g", zlab));
            prevlab = y1;
         }
         line->PaintLine(xmax-tlength, y1, xmax, y1);
         if (i == ndivz-1) {
            // label + tick at top of axis
            if ((y2 - prevlab > 1.5*lsize_user))
               label->PaintLatex(xmax + lof, y2, 0, lsize, Form("%g",fH->GetMaximum()));
            line->PaintLine(xmax-tlength, y2, xmax, y2);
         }
      }
   }

   // Take primary divisions only
   Int_t ndiv = fH->GetZaxis()->GetNdivisions();
   Bool_t isOptimized = ndiv>0;
   Int_t absDiv = abs(ndiv);
   Int_t maxD = absDiv/1000000;
   ndiv = absDiv%100 + maxD*1000000;
   if (!isOptimized) ndiv  = -ndiv;

   char chopt[6] = "S   ";
   chopt[1] = 0;
   strncat(chopt, "+L", 3);
   if (ndiv < 0) {
      ndiv = TMath::Abs(ndiv);
      strncat(chopt, "N", 2);
   }
   if (gPad->GetLogz()) {
      wmin = TMath::Power(10., wlmin);
      wmax = TMath::Power(10., wlmax);
      strncat(chopt, "G", 2);
   }
   if (label) {
   // case option "CJUST", cleanup
      delete label;
      delete line;
   } else {
      // default
      fAxis.PaintAxis(xmax, ymin, xmax, ymax, wmin, wmax, ndiv, chopt);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TPaletteAxis::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   //char quote = '"';
   out << "   " << std::endl;
   if (gROOT->ClassSaved(TPaletteAxis::Class())) {
      out << "   ";
   } else {
      out << "   " << ClassName() << " *";
   }
   if (fOption.Contains("NDC")) {
      out << "palette = new " << ClassName() << "(" << fX1NDC << "," << fY1NDC << "," << fX2NDC << "," << fY2NDC
          << "," << fH->GetName() << ");" << std::endl;
   } else {
      out << "palette = new " << ClassName() << "(" << fX1 << "," << fY1 << "," << fX2 << "," << fY2
          << "," << fH->GetName() << ");" << std::endl;
   }
   out << "   palette->SetLabelColor(" << fAxis.GetLabelColor() << ");" << std::endl;
   out << "   palette->SetLabelFont("  << fAxis.GetLabelFont() << ");" << std::endl;
   out << "   palette->SetLabelOffset(" << fAxis.GetLabelOffset() << ");" << std::endl;
   out << "   palette->SetLabelSize("  << fAxis.GetLabelSize() << ");" << std::endl;
   out << "   palette->SetTitleOffset(" << fAxis.GetTitleOffset() << ");" << std::endl;
   out << "   palette->SetTitleSize("  << fAxis.GetTitleSize() << ");" << std::endl;
   SaveFillAttributes(out, "palette", -1, -1);
   SaveLineAttributes(out, "palette", 1, 1, 1);
}


////////////////////////////////////////////////////////////////////////////////
/// Unzoom the palette

void TPaletteAxis::UnZoom()
{
   TView *view = gPad->GetView();
   if (view) {
      delete view;
      gPad->SetView(0);
   }
   fH->GetZaxis()->SetRange(0, 0);
   if (fH->GetDimension() == 2) {
      fH->SetMinimum();
      fH->SetMaximum();
      fH->ResetBit(TH1::kIsZoomed);
   }
}
