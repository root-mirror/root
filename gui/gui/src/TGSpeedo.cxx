// @(#)root/gui:$Id: TGSpeedo.cxx
// Author: Bertrand Bellenot   26/10/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGSpeedo                                                             //
//                                                                      //
// TGSpeedo is a widget looking like a speedometer, with a needle,      //
// a counter and a small odometer window.                               //
//                                                                      //
//Begin_Html
/*
<img src="gif/speedometer.gif">
*/
//End_Html                                                              //
//                                                                      //
// Three thresholds are configurable, with their glowing color          //
// A peak mark can be enabled, allowing to keep track of the highest    //
// value displayed. The mark can be reset by right-clicking on the      //
// widget.                                                              //
// Two signals are available:                                           //
//    OdoClicked(): when user click on the small odometer window        //
//    LedClicked(): when user click on the small led near the counter   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include "TGResourcePool.h"
#include "TImage.h"
#include "TMath.h"
#include "TVirtualX.h"
#include "snprintf.h"

#include "TGSpeedo.h"


ClassImp(TGSpeedo);

////////////////////////////////////////////////////////////////////////////////
/// TGSpeedo widget constructor.

TGSpeedo::TGSpeedo(const TGWindow *p, int id)
   : TGFrame(p, 1, 1), TGWidget (id), fImage(0), fImage2(0), fBase(0)
{
   fAngleMin = -133.5;
   fAngleMax =  133.5;
   fAngle    = -133.5;
   fScaleMin = 0.0;
   fScaleMax = 100.0;
   fValue    = 0.0;
   fCounter  = 0;
   fPeakMark = kFALSE;
   fMeanMark = kFALSE;
   fPeakVal  = 0.0;
   fMeanVal  = 0.0;
   fTextFS = fCounterFS = kNone;
   fThreshold[0] = fThreshold[1] = fThreshold[2] = 0.0;
   fThresholdColor[0] = kGreen;
   fThresholdColor[1] = kOrange;
   fThresholdColor[2] = kRed;
   fThresholdActive = kFALSE;
   fPicName  = "speedo.gif";
   fImage = TImage::Open(fPicName);
   if (!fImage || !fImage->IsValid())
      Error("TGSpeedo::Build", "%s not found", fPicName.Data());
   Build();
   AddInput(kButtonPressMask | kButtonReleaseMask);
}

////////////////////////////////////////////////////////////////////////////////
/// TGSpeedo widget constructor.

TGSpeedo::TGSpeedo(const TGWindow *p, Float_t smin, Float_t smax,
                   const char *lbl1, const char *lbl2, const char *dsp1,
                   const char *dsp2, int id)
   : TGFrame(p, 1, 1), TGWidget (id), fImage(0), fImage2(0), fBase(0)
{
   fAngleMin = -133.5;
   fAngleMax =  133.5;
   fAngle    = -133.5;
   fScaleMin = smin;
   fScaleMax = smax;
   fValue    = smin;
   fCounter  = 0;
   fLabel1   = lbl1;
   fLabel2   = lbl2;
   fDisplay1 = dsp1;
   fDisplay2 = dsp2;
   fPeakMark = kFALSE;
   fMeanMark = kFALSE;
   fPeakVal  = 0.0;
   fMeanVal  = 0.0;
   fTextFS = fCounterFS = kNone;
   fThreshold[0] = fThreshold[1] = fThreshold[2] = 0.0;
   fThresholdColor[0] = kGreen;
   fThresholdColor[1] = kOrange;
   fThresholdColor[2] = kRed;
   fThresholdActive = kFALSE;
   fPicName  = "speedo.gif";
   fImage = TImage::Open(fPicName);
   if (!fImage || !fImage->IsValid())
      Error("TGSpeedo::Build", "%s not found", fPicName.Data());
   Build();
   AddInput(kButtonPressMask | kButtonReleaseMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Build TGSpeedo widget.

void TGSpeedo::Build()
{
   TString sc;
   Float_t step, mark[5];
   TString ar = "arialbd.ttf";
   Int_t i, nexe, offset;

   const TGFont *counterFont = fClient->GetFont("-*-helvetica-bold-r-*-*-12-*-*-*-*-*-*-*");
   if (!counterFont) return;
   fCounterFS = counterFont->GetFontStruct();

   const TGFont *textFont = fClient->GetFont("-*-helvetica-bold-r-*-*-8-*-*-*-*-*-*-*");
   if (!textFont) return;
   fTextFS = textFont->GetFontStruct();

   const TGFont *labelFont = fClient->GetFont("-*-helvetica-bold-r-*-*-14-*-*-*-*-*-*-*");
   if (!labelFont) return;
   FontStruct_t labelFS = labelFont->GetFontStruct();

   if (fImage && fImage->IsValid()) {
      fBase = fClient->GetPicturePool()->GetPicture(gSystem->BaseName(fPicName.Data()),
              fImage->GetPixmap(), fImage->GetMask());
      // center of the image
      Float_t xc = (Float_t)(fBase ? (fBase->GetWidth() + 1) / 2 : 96.0);
      Float_t yc = (Float_t)(fBase ? (fBase->GetHeight() + 1) / 2 : 96.0);

      // compute scale ticks steps
      step = (fScaleMax - fScaleMin) / 4.0;
      mark[0] = fScaleMin;
      mark[4] = fScaleMax;
      for (i=1; i<4; i++) {
         mark[i] = mark[i-1] + step;
      }
      // format tick labels
      if (fScaleMax >= 1000.0) {
         nexe = 0;
         while (1) {
            nexe++;
            for (i=0; i<5; i++) {
               mark[i] /= 10.0;
            }
            // coverity[loop_condition]: ignore - false positive
            if (mark[4] < 1000.0) break;
         }
         // draw multiplier
         fImage->DrawText((Int_t)xc - 11, (Int_t)yc + 15, "x10", 12, "#ffffff", ar);
         sc.Form("%d", nexe);
         fImage->DrawText((Int_t)xc + 11, (Int_t)yc + 13, sc.Data(), 10, "#ffffff", ar);
      }
      else if (fScaleMax < 100.0) {
         nexe = 0;
         while (1) {
            nexe--;
            for (i=0; i<5; i++) {
               mark[i] *= 10.0;
            }
            // coverity[loop_condition]: ignore - false positive
            if (mark[4] > 99.9 ) break;
         }
         // draw multiplier
         fImage->DrawText((Int_t)xc - 11, (Int_t)yc + 15, "x10", 12, "#ffffff", ar);
         sc.Form("%d", nexe);
         fImage->DrawText((Int_t)xc + 11, (Int_t)yc + 13, sc.Data(), 10, "#ffffff", ar);
      }
      // Format and draw scale tickmarks
      sc.Form("%d",(Int_t)mark[0]);
      fImage->DrawText((Int_t)xc - 51, (Int_t)yc + 30, sc.Data(), 14, "#ffffff", ar);
      sc.Form("%d",(Int_t)mark[1]);
      fImage->DrawText((Int_t)xc - 59, (Int_t)yc - 29, sc.Data(), 14, "#ffffff", ar);
      sc.Form("%d",(Int_t)mark[2]);
      offset = gVirtualX->TextWidth(labelFS, sc.Data(), sc.Length()) / 2;
      fImage->DrawText((Int_t)xc - offset, (Int_t)yc - 65, sc.Data(), 14, "#ffffff", ar);
      sc.Form("%d",(Int_t)mark[3]);
      offset = 60 - gVirtualX->TextWidth(labelFS, sc.Data(), sc.Length());
      fImage->DrawText((Int_t)xc + offset, (Int_t)yc - 29, sc.Data(), 14, "#ffffff", ar);
      sc.Form("%d",(Int_t)mark[4]);
      offset = 52 - gVirtualX->TextWidth(labelFS, sc.Data(), sc.Length());
      fImage->DrawText((Int_t)xc + offset, (Int_t)yc + 30, sc.Data(), 14, "#ffffff", ar);
      // draw main label (two lines)
      fImage->DrawText((Int_t)xc + 13, (Int_t)yc - 17, fLabel1.Data(), 14, "#ffffff", ar);
      fImage->DrawText((Int_t)xc + 13, (Int_t)yc -  4, fLabel2.Data(), 12, "#ffffff", ar);
      if (fBase)
         gVirtualX->ShapeCombineMask(fId, 0, 0, fBase->GetMask());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TGSpeedo widget Destructor.

TGSpeedo::~TGSpeedo()
{
   if (fImage && fImage->IsValid())
      delete fImage;
   if (fImage2 && fImage2->IsValid())
      delete fImage2;
   if (fBase)
      fClient->FreePicture(fBase);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the widget.

TGDimension TGSpeedo::GetDefaultSize() const
{
   if (fBase)
      return TGDimension(fBase->GetWidth(), fBase->GetHeight());
   return TGDimension(100, 100);
}

////////////////////////////////////////////////////////////////////////////////
/// Make speedo glowing.

void TGSpeedo::Glow(EGlowColor col)
{
   static EGlowColor act_col = kNoglow;
   TImage *glowImage = 0;

   if (col == act_col)
      return;

   if (fImage && fImage->IsValid())
      delete fImage;

   switch (col) {
      case kNoglow:
         break;
      case kGreen:
         glowImage = TImage::Open("glow_green.png");
         if (!glowImage || !glowImage->IsValid()) {
            Error("TGSpeedo::Glow", "glow_green.png not found");
            glowImage = 0;
         }
         break;
      case kOrange:
         glowImage = TImage::Open("glow_orange.png");
         if (!glowImage || !glowImage->IsValid()) {
            Error("TGSpeedo::Glow", "glow_orange.png not found");
            glowImage = 0;
         }
         break;
      case kRed:
         glowImage = TImage::Open("glow_red.png");
         if (!glowImage || !glowImage->IsValid()) {
            Error("TGSpeedo::Glow", "glow_red.png not found");
            glowImage = 0;
         }
         break;
   }
   fImage = TImage::Open(fPicName);
   if (fImage && fImage->IsValid() && glowImage && glowImage->IsValid()) {
      fImage->Merge(glowImage);
      delete glowImage;
   }
   act_col = col;
   Build();
   DrawText();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event.

Bool_t TGSpeedo::HandleButton(Event_t *event)
{
   if (fBase) {
      int xc = (fBase->GetWidth() + 1) / 2;
      int yc = (fBase->GetHeight() + 1) / 2;
      if (event->fType == kButtonRelease && event->fCode == kButton1) {
         // check if in the selector area
         if ((event->fX > (xc - 26)) && (event->fX < (xc + 26)) &&
             (event->fY < (yc + 50)) && (event->fY > (yc + 28))) {
            OdoClicked();
         }
         // check if in the led area
         else if ((event->fX > (xc + 30)) && (event->fX < (xc + 40)) &&
                  (event->fY > (yc + 57)) && (event->fY < (yc + 67))) {
            LedClicked();
         }
      }
      if (event->fType == kButtonRelease && event->fCode == kButton3) {
         ResetPeakVal();
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set actual value of odo meter.

void TGSpeedo::SetOdoValue(Int_t val)
{
   // avoid useless redraw
   if (val == fCounter)
      return;
   fCounter = val;
   DrawText();
   DrawNeedle();
}

////////////////////////////////////////////////////////////////////////////////
/// Set small display text (two lines).

void TGSpeedo::SetDisplayText(const char *text1, const char *text2)
{
   if (!(fDisplay1.CompareTo(text1)) &&
       !(fDisplay2.CompareTo(text2)))
      return;
   fDisplay1 = text1;
   fDisplay2 = text2;
   DrawText();
   DrawNeedle();
}

////////////////////////////////////////////////////////////////////////////////
/// Set main label text (two lines).

void TGSpeedo::SetLabelText(const char *text1, const char *text2)
{
   if (fImage && fImage->IsValid())
      delete fImage;
   fLabel1 = text1;
   fLabel2 = text2;
   fImage = TImage::Open(fPicName);
   if (!fImage || !fImage->IsValid())
      Error("TGSpeedo::Build", "%s not found", fPicName.Data());
   Build();
   DrawText();
}

////////////////////////////////////////////////////////////////////////////////
/// Set min and max scale values.

void TGSpeedo::SetMinMaxScale(Float_t min, Float_t max)
{
   if (fImage && fImage->IsValid())
      delete fImage;
   fScaleMin = min;
   fScaleMax = max;
   fImage = TImage::Open(fPicName);
   if (!fImage || !fImage->IsValid())
      Error("TGSpeedo::Build", "%s not found", fPicName.Data());
   Build();
   DrawText();
}

////////////////////////////////////////////////////////////////////////////////
/// Set actual scale (needle position) value.

void TGSpeedo::SetScaleValue(Float_t val)
{
   // avoid useless redraw
   if (val == fValue)
      return;

   fValue = val;
   if (fValue > fScaleMax)
      fValue = fScaleMax;
   else if (fValue < fScaleMin)
      fValue = fScaleMin;

   if (fThresholdActive) {
      if (fValue < fThreshold[0])
         Glow(kNoglow);
      if (fValue >= fThreshold[0] && fValue < fThreshold[1])
         Glow(fThresholdColor[0]);
      if (fValue >= fThreshold[1] && fValue < fThreshold[2])
         Glow(fThresholdColor[1]);
      if (fValue >= fThreshold[2])
         Glow(fThresholdColor[2]);
   }
   if (fValue > fPeakVal)
      fPeakVal = fValue;

   fAngle = fAngleMin + (fValue / ((fScaleMax - fScaleMin) /
           (fAngleMax - fAngleMin)));

   if (fAngle > fAngleMax)
      fAngle = fAngleMax;
   else if (fAngle < fAngleMin)
      fAngle = fAngleMin;
   DrawNeedle();
}

////////////////////////////////////////////////////////////////////////////////
/// Set actual scale (needle position) value.

void TGSpeedo::SetScaleValue(Float_t val, Int_t damping)
{
   Float_t i;
   Float_t old_val = fValue;
   Float_t step, new_val = val;
   // avoid useless redraw
   if (val == fValue)
      return;

   if ((damping > 0) || (gVirtualX->InheritsFrom("TGX11")))
      step = 2.0;
   else
      step = 0.15;

   Float_t diff_angle = fAngleMax - fAngleMin;
   Float_t diff_scale = fScaleMax - fScaleMin;
   Float_t diff_ratio = diff_scale / diff_angle;
   Float_t old_angle  = fAngleMin + (old_val / diff_ratio);
   Float_t new_angle  = fAngleMin + (new_val / diff_ratio);

   if (new_angle > old_angle) {
      for (i=old_angle; i<new_angle; i+=step) {
         new_val = (i - fAngleMin) * diff_ratio;
         SetScaleValue(new_val);
         if (damping > 0)
            gSystem->Sleep(damping);
      }
   }
   if (new_angle < old_angle) {
      for (i=old_angle; i>new_angle; i-=step) {
         new_val = (i - fAngleMin) * diff_ratio;
         SetScaleValue(new_val);
         if (damping > 0)
            gSystem->Sleep(damping);
      }
   }
   // Last step
   SetScaleValue(val);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment/decrement scale (needle position) of "step" value.

void TGSpeedo::StepScale(Float_t step)
{
   SetScaleValue(fValue + step);
}

////////////////////////////////////////////////////////////////////////////////
/// Translate distance from center and angle to xy coordinates.

void TGSpeedo::Translate(Float_t val, Float_t angle, Int_t *x, Int_t *y)
{
   Float_t xc = (Float_t)(fBase ? (fBase->GetWidth() + 1) / 2 : 96.0);
   Float_t yc = (Float_t)(fBase ? (fBase->GetHeight() + 1) / 2 : 96.0);
   *x = (Int_t)(xc + val * sin(angle * TMath::Pi() / 180) + 0.5);
   *y = (Int_t)(yc - val * cos(angle * TMath::Pi() / 180) + 0.5);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw needle in speedo widget.

void TGSpeedo::DrawNeedle()
{
   Int_t xch0, xch1, ych0, ych1;
   Int_t xpk0, ypk0, xpk1, ypk1;
   Int_t xmn0, ymn0, xmn1, ymn1;
   fValue = (fAngle - fAngleMin) * ((fScaleMax - fScaleMin) /
            (fAngleMax - fAngleMin));

   // compute x/y position of the needle
   Translate(9.0, fAngle, &xch0, &ych0);
   Translate(73.0, fAngle, &xch1, &ych1);

   // compute x/y position of the peak mark
   Float_t angle = fAngleMin + (fPeakVal / ((fScaleMax - fScaleMin) /
                  (fAngleMax - fAngleMin)));
   Translate(80.0, angle, &xpk0, &ypk0);
   Translate(67.0, angle, &xpk1, &ypk1);

   // compute x/y position of the peak mark
   angle = fAngleMin + (fMeanVal / ((fScaleMax - fScaleMin) /
          (fAngleMax - fAngleMin)));
   Translate(80.0, angle, &xmn0, &ymn0);
   Translate(70.0, angle, &xmn1, &ymn1);

   if (fImage2 && fImage2->IsValid()) {
      // First clone original image.
      TImage *img = (TImage*)fImage2->Clone("img");
      if (!img || !img->IsValid()) return;
      if (fPeakMark) {
         img->DrawLine(xpk0, ypk0, xpk1, ypk1, "#00ff00", 3);
         img->DrawLine(xpk0, ypk0, xpk1, ypk1, "#ffffff", 1);
      }
      if (fMeanMark) {
         img->DrawLine(xmn0, ymn0, xmn1, ymn1, "#ffff00", 3);
         img->DrawLine(xmn0, ymn0, xmn1, ymn1, "#ff0000", 1);
      }
      // draw line (used to render the needle) directly on the image
      img->DrawLine(xch0, ych0, xch1, ych1, "#ff0000", 2);
      // finally paint image to the widget
      img->PaintImage(fId, 0, 0, 0, 0, 0, 0, "opaque");
      // and finally, to avoid memory leaks
      delete img;
   }
   gVirtualX->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text in speedo widget.

void TGSpeedo::DrawText()
{
   char sval[80];
   char dsval[80];
   Int_t strSize;

   // center of the image
   Float_t xc = fBase ? (fBase->GetWidth() + 1) / 2 : 96.0;
   Float_t yc = fBase ? (fBase->GetHeight() + 1) / 2 : 96.0;

   if (fImage && fImage->IsValid()) {
      // First clone original image.
      if (fImage2 && fImage2->IsValid())
         delete fImage2;
      fImage2 = (TImage*)fImage->Clone("fImage2");
      if (!fImage2 || !fImage2->IsValid()) return;
      TString ar = "arialbd.ttf";
      // format counter value
      Int_t nexe = 0;
      Int_t ww = fCounter;
      if (fCounter >= 10000) {
         while (1) {
            nexe++;
            ww /= 10;
            if (nexe%3 == 0 && ww < 10000) break;
         }
         fImage2->DrawText((Int_t)xc - 9, (Int_t)yc + 72, "x10", 10, "#ffffff", ar);
         snprintf(sval, 80, "%d", nexe);
         fImage2->DrawText((Int_t)xc + 9, (Int_t)yc + 69, sval, 8, "#ffffff", ar);
      }
      snprintf(sval, 80, "%04d", (int)ww);
      snprintf(dsval, 80, "%c %c %c %c", sval[0], sval[1], sval[2], sval[3]);
      // draw text in the counter
      if (gVirtualX->InheritsFrom("TGX11")) {
         // as there is a small difference between Windows and Linux...
         fImage2->DrawText((Int_t)xc - 18, (Int_t)yc + 55, dsval, 12, "#ffffff", ar);
      }
      else {
         fImage2->DrawText((Int_t)xc - 16, (Int_t)yc + 56, dsval, 12, "#ffffff", ar);
      }
      // compute the size of the string to draw in the small display box
      // first line
      strSize = gVirtualX->TextWidth(fTextFS, fDisplay1.Data(), fDisplay1.Length()) - 6;
      // draw text directly on the imaget_t)yc + 29, fDispla
      fImage2->DrawText((Int_t)xc - (strSize / 2), (Int_t)yc + 29, fDisplay1.Data(), 8, "#ffffff", ar);
      // second line
      strSize = gVirtualX->TextWidth(fTextFS, fDisplay2.Data(), fDisplay2.Length()) - 6;
      fImage2->DrawText((Int_t)xc - (strSize / 2), (Int_t)yc + 38, fDisplay2.Data(), 8, "#ffffff", ar);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw speedo widget.

void TGSpeedo::DoRedraw()
{
   char sval[80];
   char dsval[80];
   Int_t strSize;
   Int_t xch0, xch1, ych0, ych1;
   Int_t xpk0, ypk0, xpk1, ypk1;
   Int_t xmn0, ymn0, xmn1, ymn1;
   static Bool_t first = kTRUE;
   if (first) {
      TGFrame::DoRedraw();
      first = kFALSE;
   }
   fValue = (fAngle - fAngleMin) * ((fScaleMax - fScaleMin) /
            (fAngleMax - fAngleMin));

   // center of the image
   Float_t xc = fBase ? (fBase->GetWidth() + 1) / 2 : 96.0;
   Float_t yc = fBase ? (fBase->GetHeight() + 1) / 2 : 96.0;

   // compute x/y position of the needle
   Translate(9.0, fAngle, &xch0, &ych0);
   Translate(73.0, fAngle, &xch1, &ych1);

   // compute x/y position of the peak mark
   Float_t angle = fAngleMin + (fPeakVal / ((fScaleMax - fScaleMin) /
                  (fAngleMax - fAngleMin)));
   Translate(80.0, angle, &xpk0, &ypk0);
   Translate(67.0, angle, &xpk1, &ypk1);

   // compute x/y position of the peak mark
   angle = fAngleMin + (fMeanVal / ((fScaleMax - fScaleMin) /
          (fAngleMax - fAngleMin)));
   Translate(80.0, angle, &xmn0, &ymn0);
   Translate(70.0, angle, &xmn1, &ymn1);

   if (fImage && fImage->IsValid()) {
      // First clone original image.
      if (fImage2 && fImage2->IsValid())
         delete fImage2;
      fImage2 = (TImage*)fImage->Clone("fImage2");
      if (!fImage2 || !fImage2->IsValid()) return;
      TString ar = "arialbd.ttf";
      // format counter value
      Int_t nexe = 0;
      Int_t ww = fCounter;
      if (fCounter >= 10000) {
         while (1) {
            nexe++;
            ww /= 10;
            if (nexe%3 == 0 && ww < 10000) break;
         }
         fImage2->DrawText((Int_t)xc - 9, (Int_t)yc + 72, "x10", 10, "#ffffff", ar);
         snprintf(sval, 80, "%d", nexe);
         fImage2->DrawText((Int_t)xc + 9, (Int_t)yc + 69, sval, 8, "#ffffff", ar);
      }
      snprintf(sval, 80, "%04d", (int)ww);
      snprintf(dsval, 80, "%c %c %c %c", sval[0], sval[1], sval[2], sval[3]);
      // draw text in the counter
      if (gVirtualX->InheritsFrom("TGX11")) {
         // as there is a small difference between Windows and Linux...
         fImage2->DrawText((Int_t)xc - 18, (Int_t)yc + 55, dsval, 12, "#ffffff", ar);
      }
      else {
         fImage2->DrawText((Int_t)xc - 16, (Int_t)yc + 56, dsval, 12, "#ffffff", ar);
      }
      // compute the size of the string to draw in the small display box
      // first line
      strSize = gVirtualX->TextWidth(fTextFS, fDisplay1.Data(), fDisplay1.Length()) - 6;
      // draw text directly on the imaget_t)yc + 29, fDispla
      fImage2->DrawText((Int_t)xc - (strSize / 2), (Int_t)yc + 29, fDisplay1.Data(), 8, "#ffffff", ar);
      // second line
      strSize = gVirtualX->TextWidth(fTextFS, fDisplay2.Data(), fDisplay2.Length()) - 6;
      fImage2->DrawText((Int_t)xc - (strSize / 2), (Int_t)yc + 38, fDisplay2.Data(), 8, "#ffffff", ar);
      TImage *img = (TImage*)fImage2->Clone("img");
      if (!img || !img->IsValid()) return;
      if (fPeakMark) {
         img->DrawLine(xpk0, ypk0, xpk1, ypk1, "#00ff00", 3);
         img->DrawLine(xpk0, ypk0, xpk1, ypk1, "#ffffff", 1);
      }
      if (fMeanMark) {
         img->DrawLine(xmn0, ymn0, xmn1, ymn1, "#ffff00", 3);
         img->DrawLine(xmn0, ymn0, xmn1, ymn1, "#ff0000", 1);
      }
      // draw line (used to render the needle) directly on the image
      img->DrawLine(xch0, ych0, xch1, ych1, "#ff0000", 2);
      // finally paint image to the widget
      img->PaintImage(fId, 0, 0, 0, 0, 0, 0, "opaque");
      // and finally, to avoid memory leaks
      delete img;
   }
}
