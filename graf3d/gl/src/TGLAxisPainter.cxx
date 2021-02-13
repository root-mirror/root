// @(#)root/gl:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLAxisPainter.h"

#include "TGLRnrCtx.h"
#include "TGLCamera.h"
#include "TGLIncludes.h"
#include "TGLFontManager.h"

#include "TAttAxis.h"
#include "TAxis.h"
#include "TH1.h"
#include "THLimitsFinder.h"

#include "TMath.h"
#include "TPRegexp.h"

/** \class TGLAxisPainterBox
\ingroup opengl
Utility class to paint axis in GL.
*/

ClassImp(TGLAxisPainter);

////////////////////////////////////////////////////////////////////////////////

TGLAxisPainter::TGLAxisPainter():
   fExp(0),
   fMaxDigits(5),
   fDecimals(0),

   fAttAxis(0), fUseAxisColors(kTRUE),

   fFontMode(TGLFont::kTexture),
   fDir(1, 0, 0),
   fTMNDim(1),
   fLabelPixelFontSize(14), fLabel3DFontSize(1.0),
   fTitlePixelFontSize(14), fTitle3DFontSize(1.0),

   fLabelAlignH(TGLFont::kCenterH),
   fLabelAlignV(TGLFont::kCenterV),
   fAllZeroesRE(0)
{
   // Constructor.

   fAllZeroesRE = new TPMERegexp("[-+]?0\\.0*$", "o");
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLAxisPainter::~TGLAxisPainter()
{
   delete fAllZeroesRE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set label align.

void TGLAxisPainter::SetLabelAlign(TGLFont::ETextAlignH_e h, TGLFont::ETextAlignV_e v)
{
   fLabelAlignH = h;
   fLabelAlignV = v;
}

////////////////////////////////////////////////////////////////////////////////
/// Find first and last character of a label.

void TGLAxisPainter::LabelsLimits(const char *label, Int_t &first, Int_t &last) const
{
   last = strlen(label) - 1;
   for (Int_t i = 0; i <= last; i++) {
      if (strchr("1234567890-+.", label[i])) {
         first = i;
         return;
      }
   }
   Error("LabelsLimits", "attempt to draw a blank label");
}

////////////////////////////////////////////////////////////////////////////////
/// Returns formatted text suitable for display of value.

void TGLAxisPainter::FormAxisValue(Double_t  val, TString &s) const
{
   s.Form(fFormat, val);
   s = s.Strip(TString::kLeading);

   if (s == "-." || s == "-0")
   {
      s  = "0";
      return;
   }

   Ssiz_t ld = s.Last('.') + 1;
   if (s.Length() - ld > fDecimals)
      s.Remove(ld + fDecimals);


   if (fDecimals == 0 && s.EndsWith("."))
      s.Remove(s.Length() -1);

   fAllZeroesRE->Substitute(s, "0", kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct print format from given primary bin width.

void TGLAxisPainter::SetTextFormat(Double_t min, Double_t max, Double_t bw1)
{
   Double_t absMax = TMath::Max(TMath::Abs(min), TMath::Abs(max));
   Double_t epsilon = 1e-5;
   Double_t absMaxLog = TMath::Log10(absMax) + epsilon;

   fExp   = 0;
   Int_t if1, if2;
   Double_t xmicros = TMath::Power(10, -fMaxDigits);
   if (bw1 < xmicros && absMaxLog < 0) {
      // First case : bin width less than 0.001
      fExp = (Int_t)absMaxLog;
      if (fExp % 3 == 1) fExp += TMath::Sign(2, fExp);
      if (fExp % 3 == 2) fExp += TMath::Sign(1, fExp);
      if1     = fMaxDigits;
      if2     = fMaxDigits - 2;
   } else {
      // Use x 10 n format. (only powers of 3 allowed)
      Float_t af = (absMax > 1) ? absMaxLog : TMath::Log10(absMax * 0.0001);
      af += epsilon;
      Int_t clog = Int_t(af) + 1;

      if (clog > fMaxDigits) {
         while (1) {
            fExp++;
            absMax    /= 10;
            if (fExp % 3 == 0 && absMax <= TMath::Power(10, fMaxDigits - 1)) break;
         }
      } else if (clog < -fMaxDigits) {
         Double_t rne   = 1 / TMath::Power(10, fMaxDigits - 2);
         while (1) {
            fExp--;
            absMax  *= 10;
            if (fExp % 3 == 0 && absMax >= rne) break;
         }
      }

      Int_t na = 0;
      for (Int_t i = fMaxDigits - 1; i > 0; i--) {
         if (TMath::Abs(absMax) < TMath::Power(10, i)) na = fMaxDigits - i;
      }
      Double_t size =  TMath::Abs(max - min);
      Int_t ndyn = (Int_t)(size / bw1);
      while (ndyn) {
         if (size / ndyn <= 0.999 && na < fMaxDigits - 2) {
            na++;
            ndyn /= 10;
         } else break;
      }
      if2 = na;
      if1 = TMath::Max(clog + na, fMaxDigits) + 1;
   }

   // compose text format
   if (TMath::Min(min, max) < 0)if1 = if1 + 1;
   if1 = TMath::Min(if1, 32);

   // In some cases, if1 and if2 are too small....
   Double_t dwlabel = bw1 * TMath::Power(10, -fExp);
   while (dwlabel < TMath::Power(10, -if2)) {
      if1++;
      if2++;
   }
   if (if1 > 14) if1 = 14;
   if (if2 > 14) if2 = 14;
   if (if2) fFormat.Form("%%%d.%df", if1, if2);
   else     fFormat.Form("%%%d.%df", if1 + 1, 1);

   // get decimal number
   TString chtemp;
   chtemp.Form("%g", dwlabel);
   fDecimals = 0;
   if (chtemp.First('.') != kNPOS)
      fDecimals = chtemp.Length() - chtemp.First('.') - 1;
}

// Utility functions.

////////////////////////////////////////////////////////////////////////////////
/// Render text at the given position. Offset depends of text alignment.

void TGLAxisPainter::RnrText(const TString &txt, const TGLVector3 &p, TGLFont::ETextAlignH_e aH, TGLFont::ETextAlignV_e aV, const TGLFont &font) const
{
   if (fFontMode == TGLFont::kPixmap || fFontMode ==  TGLFont::kBitmap)
   {
     font.Render(txt, p.X(), p.Y(), p.Z(), aH, aV);
   }
   else
   {
      // In case of non pixmap font, size is adjusted to the projected view in order to
      // be visible on zoom out. In other words texture and polygon fonts imitate
      // pixmap font behaviour.
      glPushMatrix();
      glTranslated(p.X(), p.Y(), p.Z());
      Double_t sc = fLabel3DFontSize/fLabelPixelFontSize;
      glScaled(sc, sc, 1);
      font.Render(txt, 0, 0, 0, aH, aV);
      glPopMatrix();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set label font derived from TAttAxis.

void TGLAxisPainter::SetLabelFont(TGLRnrCtx &rnrCtx, const char* fontName, Int_t fontSize, Double_t size3d)
{
   rnrCtx.RegisterFontNoScale(fontSize, fontName, fFontMode, fLabelFont);
   fLabel3DFontSize = size3d;
   fLabelPixelFontSize = fLabelFont.GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Render label reading prepared list ov value-pos pairs.

void TGLAxisPainter::RnrLabels() const
{
   if (fUseAxisColors)
      TGLUtil::Color(fAttAxis->GetLabelColor());

   glPushMatrix();

   Float_t off = fAttAxis->GetLabelOffset() +  fAttAxis->GetTickLength();
   TGLVector3 offVec = fTMOff[0] * off;
   glTranslated(offVec.X(), offVec.Y(), offVec.Z());

   fLabelFont.PreRender();
   Double_t p = 0.;
   TString s;
   for (LabVec_t::const_iterator it = fLabVec.begin(); it != fLabVec.end(); ++it) {
      FormAxisValue((*it).second, s);
      p = (*it).first;
      RnrText(s, fDir*p, fLabelAlignH, fLabelAlignV, fLabelFont);
   }

   fLabelFont.PostRender();
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Set title font derived from TAttAxis.

void TGLAxisPainter::SetTitleFont(TGLRnrCtx &rnrCtx, const char* fontName,
                                  Int_t fontSize, Double_t size3d)
{
   rnrCtx.RegisterFontNoScale(fontSize, fontName, fFontMode, fTitleFont);
   fTitlePixelFontSize = fTitleFont.GetSize();
   fTitle3DFontSize = size3d;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw title at given position.

void TGLAxisPainter::RnrTitle(const TString &txt, TGLVector3 &pos , TGLFont::ETextAlignH_e aH, TGLFont::ETextAlignV_e aV) const
{
   if (fUseAxisColors)
      TGLUtil::Color(fAttAxis->GetTitleColor());

   TString title = (fExp) ? Form("%s [10^%d]", txt.Data(), fExp) : txt;
   fTitleFont.PreRender();
   RnrText(title, pos, aH, aV, fTitleFont);
   fTitleFont.PostRender();
}

////////////////////////////////////////////////////////////////////////////////
/// Render axis main line and tick-marks.

void TGLAxisPainter::RnrLines() const
{
   if (fUseAxisColors)
      TGLUtil::Color(fAttAxis->GetAxisColor());

   TGLUtil::LineWidth(1);
   glBegin(GL_LINES);

   // Main line.
   //
   Float_t min = fTMVec.front().first;
   Float_t max = fTMVec.back().first;
   TGLVector3 start = fDir * min;
   TGLVector3 end   = fDir * max;
   glVertex3dv(start.Arr());
   glVertex3dv(end.Arr());

   // Tick-marks.
   // Support three possible directions and two orders.
   //
   Float_t tmsOrderFirst  = fAttAxis->GetTickLength();
   Float_t tmsOrderSecond = tmsOrderFirst * 0.5;
   TGLVector3 pos;
   TMVec_t::const_iterator it = fTMVec.begin();
   Int_t nt =  fTMVec.size()-1;
   ++it;
   for (Int_t t = 1; t < nt; ++t, ++it) {
      pos = fDir * ((*it).first);
      for (Int_t dim = 0; dim < fTMNDim; dim++) {
         glVertex3dv(pos.Arr());
         if ((*it).second)
            glVertex3dv((pos + fTMOff[dim]*tmsOrderSecond).Arr());
         else
            glVertex3dv((pos + fTMOff[dim]*tmsOrderFirst).Arr());
      }
   }
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
/// GL render TAxis.

void TGLAxisPainter::PaintAxis(TGLRnrCtx &rnrCtx, TAxis* ax)
{
   fAttAxis = ax;
   Double_t min = ax->GetXmin();
   Double_t max = ax->GetXmax();
   if (min == max)
   {
      Error("TGLAxisPainter::PaintAxis", "axis without range");
      return;
   }

   //___________________________________________________________________________
   // Fill labels value-pos and tick-marks position-length.

   Int_t n1a = TMath::FloorNint(fAttAxis->GetNdivisions() / 100);
   Int_t n2a = fAttAxis->GetNdivisions() - n1a * 100;
   Int_t bn1, bn2;
   Double_t bw1, bw2;                   // primary , secondary bin width
   Double_t bl1=0, bh1=0, bl2=0, bh2=0; // bin low, high values

   // Read limits from users range
   THLimitsFinder::Optimize(min, max,       n1a, bl1, bh1, bn1, bw1);
   THLimitsFinder::Optimize(bl1, bl1 + bw1, n2a, bl2, bh2, bn2, bw2);

   //___________________________________________________________________________

   // Get TM. First and last values are reserved for axis range
   //
   fTMVec.clear();
   fLabVec.clear();

   fTMVec.push_back(TM_t(min, -1));

   Double_t v1 = bl1;
   Double_t v2 = 0;
   for (Int_t t1 = 0; t1 <= bn1; t1++)
   {
      fTMVec.push_back(TM_t(v1, 0));
      fLabVec.push_back(Lab_t(v1, v1));
      v2 = v1 + bw2;
      for (Int_t t2 = 1; t2 < bn2; t2++)
      {
         if (v2 > max) break;
         fTMVec.push_back(TM_t(v2, 1));
         v2 += bw2;
      }
      v1 += bw1;
   }

   // complete low edges for 1.st order TM
   v2 = bl1 -bw2;
   while (v2 > min) {
      fTMVec.push_back(TM_t(v2, 1));
      v2 -= bw2;
   }

   fTMVec.push_back(TM_t(max, -1));

   //___________________________________________________________________________
   // Get labels. In this case trivial one-one mapping.

   Double_t p = bl1;
   fLabVec.clear();
   SetTextFormat(min, max, bw1);
   for (Int_t i = 0; i <= bn1; i++) {
      fLabVec.push_back(Lab_t(p, p));
      p += bw1;
   }

   //___________________________________________________________________________
   // Set font.

   // First projected axis length needed if use relative font size.
   const char* labFontName   = TGLFontManager::GetFontNameFromId(fAttAxis->GetLabelFont());
   const char* titleFontName = TGLFontManager::GetFontNameFromId(fAttAxis->GetTitleFont());

   // pixel font size is set externally for pixmap and bitmap fonts
   // for texture and polygon fonts font size is set here, to get font resolution
   if (fFontMode == TGLFont::kPolygon || fFontMode == TGLFont::kTexture)
   {
      GLdouble mm[16], pm[16];
      GLint    vp[4];
      glGetDoublev(GL_MODELVIEW_MATRIX,  mm);
      glGetDoublev(GL_PROJECTION_MATRIX, pm);
      glGetIntegerv(GL_VIEWPORT, vp);

      GLdouble dn[3], up[3];
      gluProject(fDir.X()*min, fDir.Y()*min, fDir.Z()*min, mm, pm, vp, &dn[0], &dn[1], &dn[2]);
      gluProject(fDir.X()*max, fDir.Y()*max, fDir.Z()*max, mm, pm, vp, &up[0], &up[1], &up[2]);
      Double_t len = TMath::Sqrt((up[0] - dn[0]) * (up[0] - dn[0]) +
                                 (up[1] - dn[1]) * (up[1] - dn[1]) +
                                 (up[2] - dn[2]) * (up[2] - dn[2]));

      fLabelPixelFontSize = TMath::Nint(len*fAttAxis->GetLabelSize());
      fTitlePixelFontSize = TMath::Nint(len*fAttAxis->GetTitleSize());
   }

   SetLabelFont(rnrCtx, labFontName,   fLabelPixelFontSize, (max - min)*fAttAxis->GetLabelSize());
   SetTitleFont(rnrCtx, titleFontName, fTitlePixelFontSize, (max - min)*fAttAxis->GetTitleSize());

   //___________________________________________________________________________
   // Draw.

   if (!fUseAxisColors)
      TGLUtil::Color(rnrCtx.ColorSet().Markup());

   glDisable(GL_LIGHTING);
   RnrLines();
   RnrLabels();

   if (ax->GetTitle())
      RnrTitle(ax->GetTitle(), fTitlePos, fLabelAlignH, fLabelAlignV);
}


/** \class TGLAxisPainterBox
\ingroup opengl
Painter class for axes encompassing a 3D box.
*/

ClassImp(TGLAxisPainterBox);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLAxisPainterBox::TGLAxisPainterBox() :
   TGLAxisPainter()
{
   fAxis[0] = fAxis[1] = fAxis[2] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLAxisPainterBox::~TGLAxisPainterBox()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Get position of axes and titles from projected corners.

void TGLAxisPainterBox::SetAxis3DTitlePos(TGLRnrCtx &rnrCtx)
{
   Double_t x0 =  fAxis[0]->GetXmin();
   Double_t x1 =  fAxis[0]->GetXmax();

   Double_t y0 =  fAxis[1]->GetXmin();
   Double_t y1 =  fAxis[1]->GetXmax();

   Double_t z0 =  fAxis[2]->GetXmin();
   Double_t z1 =  fAxis[2]->GetXmax();

   // project corner points
   const GLdouble *pm = rnrCtx.RefCamera().RefLastNoPickProjM().CArr();
   GLdouble mm[16];
   GLint    vp[4];
   glGetDoublev(GL_MODELVIEW_MATRIX,  mm);
   glGetIntegerv(GL_VIEWPORT, vp);
   GLdouble projX[4], projY[4], projZ[4];
   GLdouble cornerX[4];
   GLdouble cornerY[4];
   cornerX[0] = x0; cornerY[0] = y0;
   cornerX[1] = x1; cornerY[1] = y0;
   cornerX[2] = x1; cornerY[2] = y1;
   cornerX[3] = x0; cornerY[3] = y1;
   gluProject(cornerX[0], cornerY[0], z0, mm, pm, vp, &projX[0], &projY[0], &projZ[0]);
   gluProject(cornerX[1], cornerY[1], z0, mm, pm, vp, &projX[1], &projY[1], &projZ[1]);
   gluProject(cornerX[2], cornerY[2], z0, mm, pm, vp, &projX[2], &projY[2], &projZ[2]);
   gluProject(cornerX[3], cornerY[3], z0, mm, pm, vp, &projX[3], &projY[3], &projZ[3]);


   // Z axis location (left most corner)
   //
   Int_t idxLeft = 0;
   Float_t xt = projX[0];
   for (Int_t i = 1; i < 4; ++i) {
      if (projX[i] < xt) {
         xt  = projX[i];
         idxLeft = i;
      }
   }
   fAxisTitlePos[2].Set(cornerX[idxLeft], cornerY[idxLeft], z1);


   // XY axis location (closest to eye) first
   //
   Float_t zt = 1.f;
   Float_t zMin = 0.f;
   Int_t idxFront = 0;
   for (Int_t i = 0; i < 4; ++i) {
      if (projZ[i] < zt) {
         zt  = projZ[i];
         idxFront = i;
      }
      if (projZ[i] > zMin) zMin = projZ[i];
   }
   Int_t xyIdx = idxFront;
   if (zMin - zt < 1e-2) xyIdx = 0; // avoid flipping in front view


   switch (xyIdx) {
      case 0:
         fAxisTitlePos[0].Set(x1, y0, z0);
         fAxisTitlePos[1].Set(x0, y1, z0);
         break;
      case 1:
         fAxisTitlePos[0].Set(x1, y0, z0);
         fAxisTitlePos[1].Set(x0, y1, z0);
         break;
      case 2:
         fAxisTitlePos[0].Set(x0, y1, z0);
         fAxisTitlePos[1].Set(x1, y0, z0);
         break;
      case 3:
         fAxisTitlePos[0].Set(x1, y1, z0);
         fAxisTitlePos[1].Set(x0, y0, z0);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XYZ axis with bitmap font.

void TGLAxisPainterBox::DrawAxis3D(TGLRnrCtx &rnrCtx)
{
   // set font size first depending on size of projected axis
   TGLMatrix mm;
   GLdouble pm[16];
   GLint    vp[4];
   glGetDoublev(GL_MODELVIEW_MATRIX, mm.Arr());
   glGetDoublev(GL_PROJECTION_MATRIX, pm);
   glGetIntegerv(GL_VIEWPORT, vp);

   // determine bitmap font size from length of projected vertical
   GLdouble dn[3];
   GLdouble up[3];
   gluProject(fAxisTitlePos[2].X(), fAxisTitlePos[2].Y(), fAxis[2]->GetXmin(), mm.Arr(), pm, vp, &dn[0], &dn[1], &dn[2]);
   gluProject(fAxisTitlePos[2].X(), fAxisTitlePos[2].Y(), fAxis[2]->GetXmax(), mm.Arr(), pm, vp, &up[0], &up[1], &up[2]);
   Double_t len = TMath::Sqrt((up[0] - dn[0]) * (up[0] - dn[0]) +
                              (up[1] - dn[1]) * (up[1] - dn[1]) +
                              (up[2] - dn[2]) * (up[2] - dn[2]));
   SetLabelPixelFontSize(TMath::CeilNint(len*fAxis[2]->GetLabelSize()));
   SetTitlePixelFontSize(TMath::CeilNint(len*fAxis[2]->GetTitleSize()));


   // Z axis
   //
   // tick-mark vector = 10 pixels left
   fAxis[2]->SetTickLength(1.); // leave this relative factor neutral
   TGLVertex3 worldRef(fAxisTitlePos[2].X(), fAxisTitlePos[2].Y(), fAxisTitlePos[2].Z());
   RefTMOff(0) = rnrCtx.RefCamera().ViewportDeltaToWorld(worldRef, -10, 0, &mm);
   SetTMNDim(1);
   RefDir().Set(0., 0., 1.);
   SetLabelAlign(TGLFont::kRight,  TGLFont::kBottom);
   glPushMatrix();
   glTranslatef(fAxisTitlePos[2].X(), fAxisTitlePos[2].Y(), 0);
   RefTitlePos().Set(RefTMOff(0).X(), RefTMOff(0).Y(),fAxisTitlePos[2].Z());
   PaintAxis(rnrCtx, fAxis[2]);
   glPopMatrix();

   // XY Axis
   //
   SetTMNDim(2);
   RefTMOff(1).Set(0, 0, fAxis[2]->GetXmin()- fAxis[2]->GetXmax());
   SetLabelAlign(TGLFont::kCenterH, TGLFont::kBottom);
   // X
   glPushMatrix();
   RefDir().Set(1, 0, 0);
   Float_t yOff = fAxis[0]->GetXmax() - fAxis[0]->GetXmin();
   yOff *= 0.5f;
   if (fAxisTitlePos[0].Y() < fAxis[1]->GetXmax()) yOff = -yOff;
   RefTMOff(0).Set(0, yOff, 0);
   glTranslatef(0, fAxisTitlePos[0].Y(), fAxisTitlePos[0].Z());
   RefTitlePos().Set(fAxisTitlePos[0].X(), yOff*1.5*fAxis[0]->GetTickLength(), 0);
   PaintAxis(rnrCtx, fAxis[0]);
   glPopMatrix();

   // Y
   glPushMatrix();
   RefDir().Set(0, 1, 0);
   Float_t xOff = fAxis[1]->GetXmax() - fAxis[1]->GetXmin();
   if (fAxisTitlePos[1].X() < fAxis[0]->GetXmax()) xOff = -xOff;
   RefTMOff(0).Set(xOff, 0, 0);
   glTranslatef(fAxisTitlePos[1].X(), 0, fAxisTitlePos[1].Z());
   RefTitlePos().Set(xOff*1.5*fAxis[1]->GetTickLength(), fAxisTitlePos[1].Y(), 0);
   PaintAxis(rnrCtx, fAxis[1]);
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLAxisPainterBox::PlotStandard(      TGLRnrCtx      &rnrCtx,
                                     TH1                  *histo,
                                     const TGLBoundingBox &bbox)
{
   fAxis[0] = histo->GetXaxis();
   fAxis[1] = histo->GetYaxis();
   fAxis[2] = histo->GetZaxis();
   // fAxis[2]->SetTitle("Z");
   // fAxis[2]->SetLabelSize(0.04);
   // fAxis[2]->SetTitleSize(0.05);

   Double_t sx = (bbox.XMax() - bbox.XMin()) / (fAxis[0]->GetXmax() - fAxis[0]->GetXmin());
   Double_t sy = (bbox.YMax() - bbox.YMin()) / (fAxis[1]->GetXmax() - fAxis[1]->GetXmin());
   Double_t sz = (bbox.ZMax() - bbox.ZMin()) / (fAxis[2]->GetXmax() - fAxis[2]->GetXmin());

   // draw
   glPushMatrix();
   glScaled(sx, sy, sz);
   SetAxis3DTitlePos(rnrCtx);
   DrawAxis3D(rnrCtx);
   glPopMatrix();
}
