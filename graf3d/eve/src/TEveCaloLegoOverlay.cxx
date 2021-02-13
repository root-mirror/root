// @(#)root/eve:$Id$
// Author: Alja Mrak-Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveCaloLegoOverlay.h"
#include "TEveCaloLegoGL.h"

#include "TAxis.h"
#include "THLimitsFinder.h"

#include "TGLRnrCtx.h"
#include "TGLIncludes.h"
#include "TGLSelectRecord.h"
#include "TGLUtil.h"
#include "TGLViewerBase.h"
#include "TGLCamera.h"
#include "TGLAxisPainter.h"
#include "TGLFontManager.h"

#include "TEveCalo.h"
#include "TEveCaloData.h"
#include "TEveRGBAPalette.h"

#include <KeySymbols.h>


/** \class TEveCaloLegoOverlay
\ingroup TEve
GL-overlay control GUI for TEveCaloLego.
*/

ClassImp(TEveCaloLegoOverlay);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveCaloLegoOverlay::TEveCaloLegoOverlay() :
   TGLCameraOverlay(),

   fCalo(0),

   fShowScales(kTRUE),
   fScaleColor(-1), fScaleTransparency(0),
   fScaleCoordX(0.85), fScaleCoordY(0.65),
   fScaleW(0), fScaleH(0),
   fCellX(-1), fCellY(-1),

   fFrameColor(-1), fFrameLineTransp(70), fFrameBgTransp(90),

   fMouseX(0),  fMouseY(0),
   fInDrag(kFALSE),

   fHeaderSelected(kFALSE),

   fPlaneAxis(0), fAxisPlaneColor(kGray),
   fShowPlane(kFALSE),

   fMenuW(0.08),
   fButtonW(0.5),
   fShowSlider(kFALSE),
   fSliderH(0.6),
   fSliderPosY(0.15),
   fSliderVal(0),

   fActiveID(-1), fActiveCol(kRed-4)
{
   fPlaneAxis = new TAxis();
}

/******************************************************************************/
// Virtual event handlers from TGLOverlayElement
/******************************************************************************/

Bool_t TEveCaloLegoOverlay::SetSliderVal(Event_t* event, TGLRnrCtx &rnrCtx)
{
   // Set height of horizontal plane in the calorimeter.

   TGLRect& wprt = rnrCtx.RefCamera().RefViewport();
   fSliderVal = (1 -event->fY*1./wprt.Height() -fSliderPosY)/fSliderH;

   if (fSliderVal < 0 )
      fSliderVal = 0;
   else if (fSliderVal > 1)
      fSliderVal = 1;

   fCalo->SetHPlaneVal(fSliderVal);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle overlay event.
/// Return TRUE if event was handled.

Bool_t TEveCaloLegoOverlay::Handle(TGLRnrCtx          & rnrCtx,
                                   TGLOvlSelectRecord & selRec,
                                   Event_t            * event)
{
   if (selRec.GetN() < 2) return kFALSE;


   if (rnrCtx.RefCamera().IsOrthographic())
   {
      switch (event->fType)
      {      case kButtonPress:
         {
            fMouseX = event->fX;
            fMouseY = event->fY;
            fInDrag = kTRUE;
            return kTRUE;
         }
         case kButtonRelease:
         {
            fInDrag = kFALSE;
            return kTRUE;
         }
         case kMotionNotify:
         {
            if (fInDrag)
            {
               const TGLRect& vp = rnrCtx.RefCamera().RefViewport();
               fScaleCoordX += (Float_t)(event->fX - fMouseX) / vp.Width();
               fScaleCoordY -= (Float_t)(event->fY - fMouseY) / vp.Height();
               fMouseX = event->fX;
               fMouseY = event->fY;
               // Make sure we don't go offscreen (use fDraw variables set in draw)
               if (fScaleCoordX < 0)
                  fScaleCoordX = 0;
               else if (fScaleCoordX + fScaleW > 1.0f)
                  fScaleCoordX = 1.0f - fScaleW;
               if (fScaleCoordY < 0)
                  fScaleCoordY = 0;
               else if (fScaleCoordY + fScaleH > 1.0f)
                  fScaleCoordY = 1.0f - fScaleH;
            }
            return kTRUE;
         }
         default:
            break;
      }
   }

   else
   {
      switch (event->fType)
      {
         case kMotionNotify:
         {
            Int_t item = selRec.GetN() < 2 ? -1 : (Int_t)selRec.GetItem(1);
            if (fActiveID != item) {
               fActiveID = item;
               return kTRUE;
            } else {
               if (fActiveID == 2 && event->fState == 256)
                  return SetSliderVal(event, rnrCtx);

               return kFALSE;
            }
            break;
         }
         case kButtonPress:
         {
            if (event->fCode != kButton1) {
               return kFALSE;
            }
            switch (selRec.GetItem(1))
            {
               case 1:
                  fShowSlider = !fShowSlider;
                  fCalo->SetDrawHPlane(fShowSlider);
                  break;
               case 2:
                  return SetSliderVal(event, rnrCtx);
               case 3:
                  fHeaderSelected = !fHeaderSelected;
               default:
                  break;
            }
         }
         default:
            break;
      }
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has entered overlay area.

Bool_t TEveCaloLegoOverlay::MouseEnter(TGLOvlSelectRecord& /*rec*/)
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has left overlay area.

void TEveCaloLegoOverlay::MouseLeave()
{
   fActiveID = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color and transparency of scales.

void TEveCaloLegoOverlay::SetScaleColorTransparency(Color_t colIdx, Char_t transp)
{
   fScaleColor = colIdx;
   fScaleTransparency = transp;
}

////////////////////////////////////////////////////////////////////////////////
/// Set scale coordinates in range [0,1].

void TEveCaloLegoOverlay::SetScalePosition(Double_t x, Double_t y)
{
   fScaleCoordX = x;
   fScaleCoordY = y;
}

////////////////////////////////////////////////////////////////////////////////
/// Set frame attributes.

void TEveCaloLegoOverlay:: SetFrameAttribs(Color_t frameColor, Char_t lineTransp, Char_t bgTransp)
{
   fFrameColor = frameColor;
   fFrameLineTransp = lineTransp;
   fFrameBgTransp = bgTransp;
}

//==============================================================================
void TEveCaloLegoOverlay::RenderHeader(TGLRnrCtx& rnrCtx)
{
   // Render text on top right corner of the screen.

   TGLRect &vp = rnrCtx.GetCamera()->RefViewport();

   TGLFont font;
   Int_t fs = TMath::Max(TMath::Nint(vp.Height()*0.035), 12);
   rnrCtx.RegisterFontNoScale(fs, "arial", TGLFont::kPixmap, font);
   font.PreRender();
   Float_t off = fs*0.2;
   Float_t bb[6];
   font.BBox(fHeaderTxt.Data(), bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
   Float_t x = vp.Width()  -bb[3] -off;
   Float_t y = vp.Height() -bb[4] -off;
   if (rnrCtx.Selection())
   {
      glPushName(0);
      glLoadName(3);
      glBegin(GL_QUADS);
      glVertex2f(x/vp.Width(), y/ vp.Height());
      glVertex2f(1,  y/ vp.Height());
      glVertex2f(1, 1);
      glVertex2f(x/vp.Width(), 1);
      glEnd();
      glPopName();
   }
   else
   {
      TGLUtil::Color(fHeaderSelected ? fActiveCol : fCalo->GetFontColor());
      glRasterPos2i(0, 0);
      glBitmap(0, 0, 0, 0, x, y, 0);
      font.Render(fHeaderTxt.Data());
   }
   font.PostRender();
}

////////////////////////////////////////////////////////////////////////////////
/// Render menu for plane-value and the plane if marked.

void TEveCaloLegoOverlay::RenderPlaneInterface(TGLRnrCtx &rnrCtx)
{
   glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_LINE_BIT | GL_POINT_BIT);
   glEnable(GL_POINT_SMOOTH);
   glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(0.1, 1);
   glDisable(GL_CULL_FACE);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   // move to the center of menu
   Double_t maxVal = fCalo->GetMaxVal();

   // button
   glPushMatrix();
   glTranslatef(1 -fMenuW, (1-fButtonW )*fMenuW*0.8, 0);

   glPushName(0);
   glLoadName(1);
   Float_t a=0.6;
   (fActiveID == 1) ? TGLUtil::Color(fActiveCol):TGLUtil::Color4f(0, 1, 0, a);
   Float_t bw = fButtonW*fMenuW*0.5;
   Float_t bwt = bw*0.8;
   Float_t bh = fButtonW*fMenuW;
   glBegin(GL_QUADS);
   glVertex2f(-bw, 0);
   glVertex2f( bw, 0);
   glVertex2f( bwt, bh);
   glVertex2f(-bwt, bh);
   glEnd();


   TGLUtil::LineWidth(1);
   glBegin(GL_LINES);
   TGLUtil::Color(4);
   glVertex2f(0, 0); glVertex2f(0, bh);
   glVertex2f((bw+bwt)*0.5, bh*0.5); glVertex2f(-(bw+bwt)*0.5, bh*0.5);
   glEnd();

   TGLUtil::LineWidth(2);
   glBegin(GL_LINE_LOOP);
   glVertex2f(-bw, 0);
   glVertex2f( bw, 0);
   glVertex2f( bwt, bh);
   glVertex2f(-bwt, bh);
   glEnd();
   TGLUtil::LineWidth(1);

   glTranslatef(0, fSliderPosY, 0.5);

   if (fShowSlider)
   {
      // event handler
      if (rnrCtx.Selection())
      {
         glLoadName(2);
         Float_t w = fButtonW*fMenuW*0.5f;
         glBegin(GL_QUADS);
         glVertex2f(-w, 0);
         glVertex2f( w, 0);
         glVertex2f( w, fSliderH);
         glVertex2f(-w, fSliderH);
         glEnd();
      }

      // slider axis
      fAxisPainter->SetLabelPixelFontSize(TMath::CeilNint(rnrCtx.GetCamera()->RefViewport().Height()*GetAttAxis()->GetLabelSize()));
      fAxisPainter->RefDir().Set(0, 1, 0);
      fAxisPainter->RefTMOff(0).Set(1, 0, 0);
      fAxisPainter->SetLabelAlign(TGLFont::kLeft, TGLFont::kCenterV);
      fPlaneAxis->SetRangeUser(0, maxVal);
      fPlaneAxis->SetLimits(0, maxVal);
      fPlaneAxis->SetNdivisions(710);
      fPlaneAxis->SetTickLength(0.02*maxVal);
      fPlaneAxis->SetLabelOffset(0.02*maxVal);
      fPlaneAxis->SetLabelSize(0.05);

      glPushMatrix();
      glScalef(fSliderH/(maxVal), fSliderH/maxVal, 1.);
      fAxisPainter->PaintAxis(rnrCtx, fPlaneAxis);
      glPopMatrix();

      // marker
      TGLUtil::Color((fActiveID == 2) ? fActiveCol : 3);
      TGLUtil::PointSize(8);
      glBegin(GL_POINTS);
      glVertex3f(0, fSliderVal*fSliderH, -0.1);
      glEnd();
   }

   glPopName();
   glPopMatrix();
   glPopAttrib();
}

/******************************************************************************/
void TEveCaloLegoOverlay::RenderLogaritmicScales(TGLRnrCtx& rnrCtx)
{
   // Draw slider of calo 2D in mode TEveCalo:fValSize.

   TGLRect &vp = rnrCtx.GetCamera()->RefViewport();

   Double_t maxVal = fCalo->GetMaxVal();
   Int_t    maxe   = TMath::CeilNint(TMath::Log10(maxVal+1)); // max round exponent
   Double_t sqv    = TMath::Power(10, maxe)+1; // starting max square value
   Double_t fc     = TMath::Log10(sqv)/TMath::Log10(fCalo->GetMaxVal()+1);
   Double_t cellX = fCellX*fc;
   Double_t cellY = fCellY*fc;

   Double_t scaleStepY = 0.1; // step is 10% of screen
   Double_t scaleStepX = scaleStepY*vp.Height()/vp.Width(); // step is 10% of screen

   Double_t frameOff = 0.01;

   // define max starting exponent not to take more than scalStepY height
   while(cellY > scaleStepY)
   {
      fc = TMath::Log10(TMath::Power(10, maxe-1)+1)/TMath::Log10(TMath::Power(10, maxe)+1);
      maxe --;
      cellX *= fc;
      cellY *= fc;
   }

   sqv =  TMath::Power(10, maxe)+1;
   glPushMatrix();
   glTranslatef(fScaleCoordX + 0.5*scaleStepX + frameOff, fScaleCoordY + 0.5*scaleStepY + frameOff, 0); // translate to lower left corner

   glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_LINE_BIT | GL_POINT_BIT);
   glEnable(GL_BLEND);
   glDisable(GL_CULL_FACE);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(0.1, 1);

   glPushName(0);
   glLoadName(1);

   // draw cells
   Color_t color = fScaleColor > -1 ? fScaleColor : rnrCtx.ColorSet().Markup().GetColorIndex();
   TGLUtil::ColorTransparency(color, fScaleTransparency);

   Float_t pos, dx, dy;
   glBegin(GL_QUADS);
   Int_t ne = 3; // max number of columns
   for (Int_t i=0; i < ne; ++i)
   {
      Float_t valFac = TMath::Log10(TMath::Power(10, maxe-i)+1)/TMath::Log10(sqv);
      dx = 0.5* cellX * valFac;
      dy = 0.5* cellY * valFac;
      pos = i* scaleStepY;
      glVertex2f( - dx, pos - dy);
      glVertex2f( - dx, pos + dy);
      glVertex2f( + dx, pos + dy);
      glVertex2f( + dx, pos - dy);
   }
   glEnd();

   // draw points in case square is below pixels
   glBegin(GL_POINTS);
   for (Int_t i=0; i < ne; ++i)
      glVertex2f(0, i* scaleStepY);
   glEnd();

   // draw numbers
   TGLFont fontB;
   Int_t fsb = TMath::Max(TMath::Nint(vp.Height()*0.03), 12);
   rnrCtx.RegisterFontNoScale(fsb, "arial", TGLFont::kPixmap, fontB);
   TGLFont fontE;
   Int_t fsE = TMath::Max(TMath::Nint(vp.Height()*0.01), 8);
   rnrCtx.RegisterFontNoScale(fsE, "arial", TGLFont::kPixmap, fontE);

   Float_t llx, lly, llz, urx, ury, urz;
   fontB.BBox("10", llx, lly, llz, urx, ury, urz);
   Float_t expX = urx/vp.Width();
   Float_t expY = (ury-lly)*0.5/vp.Height();
   Float_t expOff = 1;
   fontB.PreRender();
   fontE.PreRender();
   glPushMatrix();
   glTranslatef(0.5*scaleStepX, 0, 0.1);
   for (Int_t i = 0; i < ne; ++i)
   {
      if (i == maxe)
      {
         fontB.Render("1", 0, i*scaleStepY, 0, TGLFont::kLeft, TGLFont::kCenterV);
      }
      else if ( i == (maxe -1))
      {
         fontB.Render("10", 0, i*scaleStepY, 0, TGLFont::kLeft, TGLFont::kCenterV);
      }
      else
      {
         fontB.Render("10", 0, i*scaleStepY, 0, TGLFont::kLeft, TGLFont::kCenterV);
         fontB.BBox(Form("%d",  maxe-i), llx, lly, llz, urx, ury, urz);
         if (expOff >  urx/vp.Width()) expOff = urx/vp.Width();
         fontE.Render(Form("%d",  maxe-i), expX , i*scaleStepY+expY, 0, TGLFont::kLeft, TGLFont::kCenterV);
      }
   }
   glPopMatrix();
   fontB.PostRender();
   fontE.PostRender();
   if (expOff < 1)  expX += expOff;
   glPopMatrix();

   // draw frame
   {
      fScaleW = scaleStepX + expX+ frameOff*2;
      fScaleH = scaleStepY * ne + frameOff*2;
      Double_t x0 = fScaleCoordX;
      Double_t x1 = x0 + fScaleW;
      Double_t y0 = fScaleCoordY;
      Double_t y1 = y0 + fScaleH;
      Double_t zf = +0.2;

      color = fFrameColor > -1 ?  fFrameColor : rnrCtx.ColorSet().Markup().GetColorIndex();
      TGLUtil::ColorTransparency(color, fFrameLineTransp);

      glBegin(GL_LINE_LOOP);
      glVertex3f(x0, y0, zf); glVertex3f(x1, y0, zf);
      glVertex3f(x1, y1, zf); glVertex3f(x0, y1, zf);
      glEnd();

      TGLUtil::ColorTransparency(color, fFrameBgTransp);
      glBegin(GL_QUADS);
      glVertex2f(x0, y0); glVertex2f(x1, y0);
      glVertex2f(x1, y1); glVertex2f(x0, y1);
      glEnd();
   }
   glPopName();

   glPopAttrib();
} // end draw scales

////////////////////////////////////////////////////////////////////////////////

void TEveCaloLegoOverlay::RenderPaletteScales(TGLRnrCtx& rnrCtx)
{
   // Draw slider of calo 2D in mode TEveCalo:fValColor.

   glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_LINE_BIT);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(0.1, 1);

   TGLRect& vp = rnrCtx.RefCamera().RefViewport();
   Double_t maxVal = fCalo->GetMaxVal();
   Int_t    bn = 0;
   Double_t bw = 0;         // bin with first second order
   Double_t bl = 0, bh = 0; // bin low, high first
   THLimitsFinder::Optimize(0, maxVal, 10, bl, bh, bn, bw);
   bn = TMath::CeilNint(maxVal/bw) + 1;

   fScaleH = 0.25; // relative height of the scale
   fScaleW = fScaleH*1.5/(bn*vp.Aspect());
   Float_t h = 0.5 * bw  ;
   Float_t w = h * 1.5/ vp.Aspect();

   glPushMatrix();
   glTranslatef(fScaleCoordX + fScaleW*0.5, fScaleCoordY + fScaleH/bn*0.5, 0); // translate to lower left corner
   glScalef(fScaleH/(bn*bw), fScaleH/(bn*bw), 1.);

   glPushName(0);
   glLoadName(1);
   TGLAxisPainter::LabVec_t &labVec = fAxisPainter->RefLabVec();
   labVec.clear();
   Float_t val = 0;
   for (Int_t l= 0; l<bn; l++) {
      labVec.push_back( TGLAxisPainter::Lab_t(val, val));
      val += bw;
   }

   TGLUtil::Color(rnrCtx.ColorSet().Markup().GetColorIndex());
   fAxisPainter->RefDir().Set(0, 1, 0);
   Int_t fs = TMath::CeilNint(rnrCtx.GetCamera()->RefViewport().Height()*0.02);
   fAxisPainter->SetLabelFont(rnrCtx, "arial", fs);
   fAxisPainter->SetTextFormat(0, maxVal, bw);
   fAxisPainter->SetLabelAlign(TGLFont::kCenterH, TGLFont::kCenterV);
   TAttAxis att;
   fAxisPainter->SetAttAxis(&att);
   fAxisPainter->RnrLabels();

   UChar_t c[4];
   Float_t y;
   Double_t zf = +0.2;
   glBegin(GL_QUADS);
   for (TGLAxisPainter::LabVec_t::iterator it = labVec.begin(); it != labVec.end(); ++it)
   {
      fCalo->GetPalette()->ColorFromValue((Int_t)((*it).first), c);
      glColor4ub( c[0], c[1], c[2], c[3]);

      y = (*it).second;
      glVertex3f( -w, y - h, zf); glVertex3f( +w, y - h, zf);
      glVertex3f( +w, y + h, zf); glVertex3f( -w, y + h, zf);
   }
   glEnd();

   TGLUtil::Color(rnrCtx.ColorSet().Markup().GetColorIndex());
   glBegin(GL_LINE_LOOP);
   for (TGLAxisPainter::LabVec_t::iterator it = labVec.begin(); it != labVec.end(); ++it)
   {
      y = (*it).second;
      glVertex3f( -w, y - h, zf); glVertex3f( +w, y - h, zf);
      glVertex3f( +w, y + h, zf); glVertex3f( -w, y + h, zf);
   }
   glEnd();

   glPopName();
   glPopMatrix();
   glPopAttrib();
}

////////////////////////////////////////////////////////////////////////////////

void TEveCaloLegoOverlay::Render(TGLRnrCtx& rnrCtx)
{
   // Draw calorimeter scale info and plane interface.

   if ( fCalo == 0 || fCalo->GetData()->Empty()) return;

   Float_t old_depth_range[2];
   glGetFloatv(GL_DEPTH_RANGE, old_depth_range);
   glDepthRange(0, 0.001);

   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   if (rnrCtx.Selection())
   {
      TGLRect rect(*rnrCtx.GetPickRectangle());
      rnrCtx.GetCamera()->WindowToViewport(rect);
      gluPickMatrix(rect.X(), rect.Y(), rect.Width(), rect.Height(),
                    (Int_t*) rnrCtx.GetCamera()->RefViewport().CArr());
   }
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   glTranslatef(-1, -1, 0);
   glScalef(2, 2, 1);


   TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE);
   TGLCamera& cam = rnrCtx.RefCamera();
   Bool_t drawOverlayAxis = kTRUE;

   if (cam.IsOrthographic())
   {
      // in 2D need pixel cell dimension
      // project lego eta-phi boundaries
      TGLVector3 rng(fCalo->GetEtaRng(), fCalo->GetPhiRng(), 0);
      TGLVertex3 p;
      TGLVector3 res = cam.WorldDeltaToViewport(p, rng);

      TEveCaloLegoGL* lgl = dynamic_cast<TEveCaloLegoGL*>(rnrCtx.RefViewer().FindLogicalInScenes(fCalo));
      if (fShowScales && lgl)
      {

         // get smallest bin
         Double_t sq = 1e4;
         if (lgl->fBinStep == 1)
         {
            TEveCaloData::CellData_t cellData;
            for ( TEveCaloData::vCellId_t::iterator i = fCalo->fCellList.begin(); i != fCalo->fCellList.end(); ++i)
            {
               fCalo->fData->GetCellData(*i, cellData);
               if (sq > cellData.EtaDelta()) sq = cellData.EtaDelta();
               if (sq > cellData.PhiDelta()) sq = cellData.PhiDelta();
            }
         }
         else
         {
            TAxis* a;
            Int_t nb;
            a = fCalo->GetData()->GetEtaBins();
            nb = a->GetNbins();
            for (Int_t i=1 ; i<=nb; i++)
            {
               if (sq > a->GetBinWidth(i)) sq = a->GetBinWidth(i);
            }

            a = fCalo->GetData()->GetPhiBins();
            nb = a->GetNbins();
            for (Int_t i=1 ; i<=nb; i++)
            {
               if (sq > a->GetBinWidth(i)) sq = a->GetBinWidth(i);
            }

            sq *= lgl->fBinStep;
         }
         fCellX = (res.X()*sq)/(fCalo->GetEtaRng()*1.*cam.RefViewport().Width());
         fCellY = (res.Y()*sq)/(fCalo->GetPhiRng()*1.*cam.RefViewport().Height());
         // printf("bin width %f cells size %f %f\n", sq, fCellX, fCellY);
         if (fCalo->Get2DMode() == TEveCaloLego::kValSize)
            RenderLogaritmicScales(rnrCtx);
         else if (fCalo->GetPalette())
            RenderPaletteScales(rnrCtx);
      }

      // draw camera overlay if projected lego bbox to large
      SetFrustum(cam);
      if (   fCalo->GetEtaMin() > fFrustum[0] && fCalo->GetEtaMax() < fFrustum[2]
          && fCalo->GetPhiMin() > fFrustum[1] && fCalo->GetPhiMax() < fFrustum[3])
            drawOverlayAxis = kFALSE;
   }

   if (cam.IsPerspective() && fShowPlane)
   {
      RenderPlaneInterface(rnrCtx);
   }

   // draw info text on top right corner
   if (fHeaderTxt.Length())
   {
      RenderHeader(rnrCtx);
   }

   glPopMatrix();
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
   glMatrixMode(GL_MODELVIEW);

   glDepthRange(old_depth_range[0], old_depth_range[1]);

   if (drawOverlayAxis) TGLCameraOverlay::Render(rnrCtx);
}
