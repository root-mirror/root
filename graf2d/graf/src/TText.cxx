// @(#)root/graf:$Id$
// Author: Nicolas Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TText.h"

#include "TROOT.h"
#include "TBuffer.h"
#include "TVirtualPad.h"
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H
#include "TTF.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TPoint.h"

#include <cwchar>
#include <cstdlib>
#include <iostream>

ClassImp(TText);


/** \class TText
\ingroup BasicGraphics

Base class for several text objects.

See TAttText for a list of text attributes or fonts,
and also for a discussion on text speed and font quality.

By default, the text is drawn in the pad coordinates system.
One can draw in NDC coordinates [0,1] if the function SetNDC
is called for a TText object.

Example:
Begin_Macro(source)
{
   TText *t = new TText(.5,.5,"Hello World !");
   t->SetTextAlign(22);
   t->SetTextColor(kRed+2);
   t->SetTextFont(43);
   t->SetTextSize(40);
   t->SetTextAngle(45);
   t->Draw();
}
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// Text normal constructor.

TText::TText(Double_t x, Double_t y, const char *text) : TNamed("",text), TAttText(), fWcsTitle(nullptr)
{
   fX = x;
   fY = y;
}

////////////////////////////////////////////////////////////////////////////////
/// Text normal constructor.

TText::TText(Double_t x, Double_t y, const wchar_t *text) : TNamed(), TAttText()
{
   fX = x;
   fY = y;
   fWcsTitle = new std::wstring(text);
   SetName("");
   SetMbTitle(text);
}

////////////////////////////////////////////////////////////////////////////////
/// Text default destructor.

TText::~TText()
{
   if (fWcsTitle) delete reinterpret_cast<std::wstring *>(fWcsTitle);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TText::TText(const TText &text) : TNamed(text), TAttText(text), TAttBBox2D(text), fWcsTitle(nullptr)
{
   text.TText::Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TText &TText::operator=(const TText &src)
{
   src.TText::Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this text to text.

void TText::Copy(TObject &obj) const
{
   ((TText&)obj).fX = fX;
   ((TText&)obj).fY = fY;
   TNamed::Copy(obj);
   TAttText::Copy(((TText&)obj));
   if (((TText&)obj).fWcsTitle != NULL) {
      if (fWcsTitle != NULL) {
         *reinterpret_cast<std::wstring*>(&((TText&)obj).fWcsTitle) = *reinterpret_cast<const std::wstring*>(&fWcsTitle);
      } else {
        delete reinterpret_cast<std::wstring*>(&((TText&)obj).fWcsTitle);
        ((TText&)obj).fWcsTitle = NULL;
      }
   } else {
      if (fWcsTitle != NULL) {
         ((TText&)(obj)).fWcsTitle = new std::wstring(*reinterpret_cast<const std::wstring*>(fWcsTitle));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the text as UNICODE.

const void *TText::GetWcsTitle(void) const
{
   if (fWcsTitle != NULL) {
      return reinterpret_cast<std::wstring *>(fWcsTitle)->c_str();
   } else {
      return NULL;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a string.
/// The rectangle surrounding this string is evaluated.
/// If the point (px,py) is in the rectangle, the distance is set to zero.

Int_t TText::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t ptx, pty;

   TAttText::Modify();  // change text attributes only if necessary

   if (TestBit(kTextNDC)) {
      ptx = gPad->UtoPixel(fX);
      pty = gPad->VtoPixel(fY);
   } else {
      ptx = gPad->XtoAbsPixel(gPad->XtoPad(fX));
      pty = gPad->YtoAbsPixel(gPad->YtoPad(fY));
   }

   // Get the text control box
   Int_t cBoxX[5], cBoxY[5];
   GetControlBox(ptx, pty, -fTextAngle, cBoxX, cBoxY);
   cBoxY[4] = cBoxY[0];
   cBoxX[4] = cBoxX[0];

   // Check if the point (px,py) is inside the text control box
   if (TMath::IsInside(px, py, 5, cBoxX, cBoxY)){
      return 0;
   } else {
      return 9999;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates.

TText *TText::DrawText(Double_t x, Double_t y, const char *text)
{
   TText *newtext = new TText(x, y, text);
   TAttText::Copy(*newtext);
   newtext->SetBit(kCanDelete);
   if (TestBit(kTextNDC)) newtext->SetNDC();
   newtext->AppendPad();
   return newtext;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates.

TText *TText::DrawText(Double_t x, Double_t y, const wchar_t *text)
{
   TText *newtext = new TText(x, y, text);
   TAttText::Copy(*newtext);
   newtext->SetBit(kCanDelete);
   if (TestBit(kTextNDC)) newtext->SetNDC();
   newtext->AppendPad();
   return newtext;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates in NDC.

TText *TText::DrawTextNDC(Double_t x, Double_t y, const char *text)
{
   TText *newtext = DrawText(x, y, text);
   newtext->SetNDC();
   return newtext;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates in NDC.

TText *TText::DrawTextNDC(Double_t x, Double_t y, const wchar_t *text)
{
   TText *newtext = DrawText(x, y, text);
   newtext->SetNDC();
   return newtext;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function must be implemented to realize the action
///  corresponding to the mouse click on the object in the window

void TText::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   static Int_t px1, py1, pxold, pyold, Size, height, width;
   static Bool_t resize,turn;
   Int_t dx, dy;
   const char *text = GetTitle();
   Int_t len = strlen(text);
   Double_t sizetowin = gPad->GetAbsHNDC()*Double_t(gPad->GetWh());
   Double_t fh = (fTextSize*sizetowin);
   Int_t h     = Int_t(fh/2);
   Int_t w     = h*len;
   Short_t halign = fTextAlign/10;
   Short_t valign = fTextAlign - 10*halign;
   Double_t co, si, dtheta, norm;
   static Bool_t right, ndcsav;
   static Double_t theta;
   Int_t ax, ay, bx, by, cx, cy;
   ax = ay = 0;
   Double_t lambda, x2,y2;
   Double_t dpx,dpy,xp1,yp1;
   Int_t cBoxX[4], cBoxY[4], part;
   Double_t div = 0;
   Bool_t opaque  = gPad->OpaqueMoving();

   if (!gPad->IsEditable()) return;
   switch (event) {

   case kArrowKeyPress:
   case kButton1Down:
      ndcsav = TestBit(kTextNDC);
      // No break !!!

   case kMouseMotion:
      if (TestBit(kTextNDC)) {
         px1 = gPad->UtoPixel(fX);
         py1 = gPad->VtoPixel(fY);
      } else {
         px1 = gPad->XtoAbsPixel(gPad->XtoPad(fX));
         py1 = gPad->YtoAbsPixel(gPad->YtoPad(fY));
      }
      theta  = fTextAngle;
      Size   = 0;
      pxold  = px;
      pyold  = py;
      co     = TMath::Cos(fTextAngle*0.017453293);
      si     = TMath::Sin(fTextAngle*0.017453293);
      resize = kFALSE;
      turn   = kFALSE;
      GetControlBox(px1, py1, -theta, cBoxX, cBoxY);
      div    = ((cBoxX[3]-cBoxX[0])*co-(cBoxY[3]-cBoxY[0])*si);
      if (TMath::Abs(div) > 1e-8) part = (Int_t)(3*((px-cBoxX[0])*co-(py-cBoxY[0])*si)/ div);
      else part = 0;
      switch (part) {
      case 0:
         if (halign == 3) {
            turn  = kTRUE;
            right = kTRUE;
            gPad->SetCursor(kRotate);
         } else {
            resize = kTRUE;
            height = valign;
            width  = halign;
            gPad->SetCursor(kArrowVer);
         }
         break;
      case 1:
         gPad->SetCursor(kMove);
         break;
      case 2:
         if (halign == 3) {
            resize = kTRUE;
            height = valign;
            width  = halign;
            gPad->SetCursor(kArrowVer);
         } else {
            turn  = kTRUE;
            right = kFALSE;
            gPad->SetCursor(kRotate);
         }
      }
      break;

   case kArrowKeyRelease:
   case kButton1Motion:
      if (!opaque) PaintControlBox(px1, py1, -theta);
      if (turn) {
         norm = TMath::Sqrt(Double_t((py-py1)*(py-py1)+(px-px1)*(px-px1)));
         if (norm>0) {
            theta = TMath::ACos((px-px1)/norm);
            dtheta= TMath::ASin((py1-py)/norm);
            if (dtheta<0) theta = -theta;
            theta = theta/TMath::ACos(-1)*180;
            if (theta<0) theta += 360;
            if (right) {theta = theta+180; if (theta>=360) theta -= 360;}
         }
      } else if (resize) {

         co = TMath::Cos(fTextAngle*0.017453293);
         si = TMath::Sin(fTextAngle*0.017453293);
         if (width == 1) {
            switch (valign) {
               case 1 : ax = px1; ay = py1; break;
               case 2 : ax = px1+Int_t(si*h/2); ay = py1+Int_t(co*h/2); break;
               case 3 : ax = px1+Int_t(si*h*3/2); ay = py1+Int_t(co*h*3/2); break;
            }
         }
         if (width == 2) {
            switch (valign) {
               case 1 : ax = px1-Int_t(co*w/2); ay = py1+Int_t(si*w/2); break;
               case 2 : ax = px1-Int_t(co*w/2+si*h/2); ay = py1+Int_t(si*w/2+co*h/2); break;
               case 3 : ax = px1-Int_t(co*w/2+si*h*3/2); ay = py1+Int_t(si*w/2+co*h*3/2); break;
            }
         }
         if (width == 3) {
            switch (valign) {
               case 1 : ax = px1-Int_t(co*w); ay = py1+Int_t(si*w); break;
               case 2 : ax = px1-Int_t(co*w+si*h/2); ay = py1+Int_t(si*w+co*h/2); break;
               case 3 : ax = px1-Int_t(co*w+si*h*3/2); ay = py1+Int_t(si*w+co*h*3/2); break;
            }
         }
         if (height == 3) {bx = ax-Int_t(si*h); by = ay-Int_t(co*h);}
         else {bx = ax; by = ay;}
         cx = bx+Int_t(co*w); cy = by-Int_t(si*w);
         lambda = Double_t(((px-bx)*(cx-bx)+(py-by)*(cy-by)))/Double_t(((cx-bx)*(cx-bx)+(cy-by)*(cy-by)));
         x2 = Double_t(px) - lambda*Double_t(cx-bx)-Double_t(bx);
         y2 = Double_t(py) - lambda*Double_t(cy-by)-Double_t(by);
         Size = Int_t(TMath::Sqrt(x2*x2+y2*y2)*2);
         if (Size<4) Size = 4;

         SetTextSize(Size/sizetowin);
         TAttText::Modify();
      } else {
         dx = px - pxold;  px1 += dx;
         dy = py - pyold;  py1 += dy;
      }
      if (opaque) {
         if (ndcsav) this->SetNDC(kFALSE);
         this->SetX(gPad->PadtoX(gPad->AbsPixeltoX(px1)));
         this->SetY(gPad->PadtoY(gPad->AbsPixeltoY(py1)));
         if (resize) gPad->ShowGuidelines(this, event, 't', false);
         if ((!resize)&&(!turn)) gPad->ShowGuidelines(this, event, 'i', true);
         gPad->ShowGuidelines(this, event, !resize&!turn);
         this->SetTextAngle(theta);
         gPad->Modified(kTRUE);
         gPad->Update();
      }
      if (!opaque) PaintControlBox(px1, py1, -theta);
      pxold = px;  pyold = py;
      break;

   case kButton1Up:
      if (opaque) {
         if (ndcsav && !this->TestBit(kTextNDC)) {
            this->SetX((fX - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1()));
            this->SetY((fY - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1()));
            this->SetNDC();
         }
         gPad->ShowGuidelines(this, event, !resize&!turn);
      } else {
         if (TestBit(kTextNDC)) {
            dpx  = gPad->GetX2() - gPad->GetX1();
            dpy  = gPad->GetY2() - gPad->GetY1();
            xp1  = gPad->GetX1();
            yp1  = gPad->GetY1();
            fX = (gPad->AbsPixeltoX(px1)-xp1)/dpx;
            fY = (gPad->AbsPixeltoY(py1)-yp1)/dpy;
         } else {
            fX = gPad->PadtoX(gPad->AbsPixeltoX(px1));
            fY = gPad->PadtoY(gPad->AbsPixeltoY(py1));
         }
         fTextAngle = theta;
      }
      gPad->Modified(kTRUE);
      break;

   case kButton1Locate:
      ExecuteEvent(kButton1Down, px, py);

      while (1) {
         px = py = 0;
         event = gVirtualX->RequestLocator(1, 1, px, py);

         ExecuteEvent(kButton1Motion, px, py);

         if (event != -1) {                     // button is released
            ExecuteEvent(kButton1Up, px, py);
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the text control box. The text position coordinates is (x,y) and
/// the text angle is theta. The control box coordinates are returned in cBoxX
/// and cBoxY.

void TText::GetControlBox(Int_t x, Int_t y, Double_t theta,
                          Int_t cBoxX[4], Int_t cBoxY[4])
{
   Short_t halign = fTextAlign/10;          // horizontal alignment
   Short_t valign = fTextAlign - 10*halign; // vertical alignment
   UInt_t cBoxW, cBoxH;                     // control box width and heigh
   UInt_t Dx = 0, Dy = 0;                   // delta along x and y to align the box

   GetBoundingBox(cBoxW, cBoxH);

   // compute the translations (Dx, Dy) required by the alignments
   switch (halign) {
      case 1 : Dx = 0      ; break;
      case 2 : Dx = cBoxW/2; break;
      case 3 : Dx = cBoxW  ; break;
   }
   switch (valign) {
      case 1 : Dy = 0      ; break;
      case 2 : Dy = cBoxH/2; break;
      case 3 : Dy = cBoxH  ; break;
   }

   // compute the control box coordinates before rotation
   cBoxX[0] = x-Dx;
   cBoxY[0] = y+Dy;
   cBoxX[1] = x-Dx;
   cBoxY[1] = y-cBoxH+Dy;
   cBoxX[2] = x+cBoxW-Dx;
   cBoxY[2] = y-cBoxH+Dy;
   cBoxX[3] = x+cBoxW-Dx;
   cBoxY[3] = y+Dy;

   // rotate the control box if needed
   if (theta) {
      Double_t cosTheta = TMath::Cos(theta*0.017453293);
      Double_t sinTheta = TMath::Sin(theta*0.017453293);
      for (int i=0; i<4 ; i++) {
         Int_t hcBoxX = cBoxX[i];
         Int_t hcBoxY = cBoxY[i];
         cBoxX[i] = (Int_t)((hcBoxX-x)*cosTheta-(hcBoxY-y)*sinTheta+x);
         cBoxY[i] = (Int_t)((hcBoxX-x)*sinTheta+(hcBoxY-y)*cosTheta+y);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return text size in pixels. By default the size returned does not take
/// into account the text angle (angle = kFALSE). If angle is set to kTRUE
/// w and h take the angle into account.

void TText::GetBoundingBox(UInt_t &w, UInt_t &h, Bool_t angle)
{
   const char *text = GetTitle();
   if (!text[0]) {
      w = h = 0;
      return;
   }

   if (angle) {
      Int_t cBoxX[4], cBoxY[4];
      Int_t ptx, pty;
      if (TestBit(kTextNDC)) {
         ptx = gPad->UtoPixel(fX);
         pty = gPad->VtoPixel(fY);
      } else {
         ptx = gPad->XtoAbsPixel(gPad->XtoPad(fX));
         pty = gPad->YtoAbsPixel(gPad->YtoPad(fY));
      }
      GetControlBox(ptx, pty, fTextAngle, cBoxX, cBoxY);
      Int_t x1 = cBoxX[0];
      Int_t x2 = cBoxX[0];
      Int_t y1 = cBoxY[0];
      Int_t y2 = cBoxY[0];
      for (Int_t i=1; i<4; i++) {
         if (cBoxX[i] < x1) x1 = cBoxX[i];
         if (cBoxX[i] > x2) x2 = cBoxX[i];
         if (cBoxY[i] < y1) y1 = cBoxY[i];
         if (cBoxY[i] > y2) y2 = cBoxY[i];
      }
      w = x2-x1;
      h = y2-y1;
   } else {
      if ((gVirtualX->HasTTFonts() && TTF::IsInitialized()) || gPad->IsBatch()) {
         TTF::GetTextExtent(w, h, (char*)GetTitle());
      } else {
         const Font_t oldFont = gVirtualX->GetTextFont();
         if (gVirtualX->InheritsFrom("TGCocoa"))
            gVirtualX->SetTextFont(fTextFont);
         gVirtualX->GetTextExtent(w, h, (char*)GetTitle());
         if (gVirtualX->InheritsFrom("TGCocoa"))
            gVirtualX->SetTextFont(oldFont);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return text ascent and descent for string text
///  - in a return total text ascent
///  - in d return text descent

void TText::GetTextAscentDescent(UInt_t &a, UInt_t &d, const char *text) const
{
   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Double_t tsize;
   if (wh < hh)  tsize = fTextSize*wh;
   else          tsize = fTextSize*hh;

   if (gVirtualX->HasTTFonts() || gPad->IsBatch()) {
      TTF::SetTextFont(fTextFont);
      TTF::SetTextSize(tsize);
      a = TTF::GetBox().yMax;
      d = TMath::Abs(TTF::GetBox().yMin);
   } else {
      const Font_t oldFont = gVirtualX->GetTextFont();
      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(fTextFont);
      gVirtualX->SetTextSize(tsize);
      a = gVirtualX->GetFontAscent(text);
      if (!a) {
         UInt_t w;
         gVirtualX->GetTextExtent(w, a, (char*)text);
      }
      d = gVirtualX->GetFontDescent(text);
      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(oldFont);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return text ascent and descent for string text
///  - in a return total text ascent
///  - in d return text descent

void TText::GetTextAscentDescent(UInt_t &a, UInt_t &d, const wchar_t *text) const
{
   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Double_t tsize;
   if (wh < hh)  tsize = fTextSize*wh;
   else          tsize = fTextSize*hh;

   if (gVirtualX->HasTTFonts() || gPad->IsBatch() || gVirtualX->InheritsFrom("TGCocoa")) {
      TTF::SetTextFont(fTextFont);
      TTF::SetTextSize(tsize);
      a = TTF::GetBox().yMax;
      d = TMath::Abs(TTF::GetBox().yMin);
   } else {
      gVirtualX->SetTextSize(tsize);
      a = gVirtualX->GetFontAscent();
      if (!a) {
         UInt_t w;
         gVirtualX->GetTextExtent(w, a, (wchar_t*)text);
      }
      d = gVirtualX->GetFontDescent();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return text extent for string text
///  - in w return total text width
///  - in h return text height

void TText::GetTextExtent(UInt_t &w, UInt_t &h, const char *text) const
{
   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Double_t tsize;
   if (wh < hh)  tsize = fTextSize*wh;
   else          tsize = fTextSize*hh;

   if (gVirtualX->HasTTFonts() || gPad->IsBatch()) {
      TTF::SetTextFont(fTextFont);
      TTF::SetTextSize(tsize);
      TTF::GetTextExtent(w, h, (char*)text);
   } else {
      const Font_t oldFont = gVirtualX->GetTextFont();
      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(fTextFont);
      gVirtualX->SetTextSize(tsize);
      gVirtualX->GetTextExtent(w, h, (char*)text);
      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(oldFont);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return text advance for string text
/// if kern is true (default) kerning is taken into account. If it is false
/// the kerning is not taken into account.

void TText::GetTextAdvance(UInt_t &a, const char *text, const Bool_t kern) const
{
   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Double_t tsize;
   if (wh < hh)  tsize = fTextSize*wh;
   else          tsize = fTextSize*hh;

   if (gVirtualX->HasTTFonts() || gPad->IsBatch()) {
      Bool_t kernsave = TTF::GetKerning();
      TTF::SetKerning(kern);
      TTF::SetTextFont(fTextFont);
      TTF::SetTextSize(tsize);
      TTF::GetTextAdvance(a, (char*)text);
      TTF::SetKerning(kernsave);
   } else {
      UInt_t h;
      const Font_t oldFont = gVirtualX->GetTextFont();
      //how do I know what to calculate without a font???
      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(fTextFont);

      gVirtualX->SetTextSize(tsize);
      gVirtualX->GetTextExtent(a, h, (char*)text);

      if (gVirtualX->InheritsFrom("TGCocoa"))
         gVirtualX->SetTextFont(oldFont);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return text extent for string text
///  - in w return total text width
///  - in h return text height

void TText::GetTextExtent(UInt_t &w, UInt_t &h, const wchar_t *text) const
{
   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Double_t tsize;
   if (wh < hh)  tsize = fTextSize*wh;
   else          tsize = fTextSize*hh;

   if (gVirtualX->HasTTFonts() || gPad->IsBatch() || gVirtualX->InheritsFrom("TGCocoa")) {
      TTF::SetTextFont(fTextFont);
      TTF::SetTextSize(tsize);
      TTF::GetTextExtent(w, h, (wchar_t*)text);
   } else {
      gVirtualX->SetTextSize(tsize);
      gVirtualX->GetTextExtent(w, h, (wchar_t*)text);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List this text with its attributes.

void TText::ls(Option_t *) const
{
   TROOT::IndentLevel();
   printf("Text  X=%f Y=%f Text=%s\n",fX,fY,GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this text with its current attributes.

void TText::Paint(Option_t *)
{
   TAttText::Modify();  //Change text attributes only if necessary
   if (TestBit(kTextNDC)) gPad->PaintTextNDC(fX,fY,GetTitle());
   else                   gPad->PaintText(gPad->XtoPad(fX),gPad->YtoPad(fY),GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the text control box. (x,y) are the coordinates where the control
/// box should be painted and theta is the angle of the box.

void TText::PaintControlBox(Int_t x, Int_t y, Double_t theta)
{
   Int_t cBoxX[4], cBoxY[4];
   Short_t halign = fTextAlign/10;               // horizontal alignment
   Short_t valign = fTextAlign - 10*halign;      // vertical alignment

   GetControlBox(x, y, theta, cBoxX, cBoxY);
   // Draw the text control box outline
   gVirtualX->SetLineStyle((Style_t)1);
   gVirtualX->SetLineWidth(1);
   gVirtualX->SetLineColor(1);
   gVirtualX->DrawLine(cBoxX[0], cBoxY[0], cBoxX[1], cBoxY[1]);
   gVirtualX->DrawLine(cBoxX[1], cBoxY[1], cBoxX[2], cBoxY[2]);
   gVirtualX->DrawLine(cBoxX[2], cBoxY[2], cBoxX[3], cBoxY[3]);
   gVirtualX->DrawLine(cBoxX[3], cBoxY[3], cBoxX[0], cBoxY[0]);

   // Draw a symbol at the text starting point
   TPoint p;
   Int_t ix = 0, iy = 0;
   switch (halign) {
      case 1 :
         switch (valign) {
            case 1 : ix = 0 ; iy = 0 ; break;
            case 2 : ix = 0 ; iy = 1 ; break;
            case 3 : ix = 1 ; iy = 1 ; break;
         }
      break;
      case 2 :
         switch (valign) {
            case 1 : ix = 0 ; iy = 3 ; break;
            case 2 : ix = 0 ; iy = 2 ; break;
            case 3 : ix = 1 ; iy = 2 ; break;
         }
      break;
      case 3 :
         switch (valign) {
            case 1 : ix = 3 ; iy = 3 ; break;
            case 2 : ix = 2 ; iy = 3 ; break;
            case 3 : ix = 2 ; iy = 2 ; break;
         }
      break;
   }
   p.fX = (cBoxX[ix]+cBoxX[iy])/2;
   p.fY = (cBoxY[ix]+cBoxY[iy])/2;
   gVirtualX->SetMarkerColor(1);
   gVirtualX->SetMarkerStyle(24);
   gVirtualX->SetMarkerSize(0.7);
   gVirtualX->DrawPolyMarker(1, &p);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates.

void TText::PaintText(Double_t x, Double_t y, const char *text)
{
   TAttText::Modify();  //Change text attributes only if necessary
   if (gPad) gPad->PaintText(x,y,text);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates.

void TText::PaintText(Double_t x, Double_t y, const wchar_t *text)
{
   TAttText::Modify();  //Change text attributes only if necessary
   if (gPad) gPad->PaintText(x,y,text);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates in NDC.

void TText::PaintTextNDC(Double_t u, Double_t v, const char *text)
{
   TAttText::Modify();  //Change text attributes only if necessary
   if (gPad) gPad->PaintTextNDC(u,v,text);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this text with new coordinates in NDC.

void TText::PaintTextNDC(Double_t u, Double_t v, const wchar_t *text)
{
   TAttText::Modify();  //Change text attributes only if necessary
   if (gPad) gPad->PaintTextNDC(u,v,text);
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this text with its attributes.

void TText::Print(Option_t *) const
{
   printf("Text  X=%f Y=%f Text=%s Font=%d Size=%f",fX,fY,GetTitle(),GetTextFont(),GetTextSize());
   if (GetTextColor() != 1 ) printf(" Color=%d",GetTextColor());
   if (GetTextAlign() != 10) printf(" Align=%d",GetTextAlign());
   if (GetTextAngle() != 0 ) printf(" Angle=%f",GetTextAngle());
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TText::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   if (gROOT->ClassSaved(TText::Class())) {
       out<<"   ";
   } else {
       out<<"   TText *";
   }
   TString s = GetTitle();
   s.ReplaceAll("\"","\\\"");
   out<<"text = new TText("<<fX<<","<<fY<<","<<quote<<s.Data()<<quote<<");"<<std::endl;
   if (TestBit(kTextNDC)) out<<"   text->SetNDC();"<<std::endl;

   SaveTextAttributes(out,"text",11,0,1,62,0.05);

   out<<"   text->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NDC mode on if isNDC = kTRUE, off otherwise

void TText::SetNDC(Bool_t isNDC)
{
   ResetBit(kTextNDC);
   if (isNDC) SetBit(kTextNDC);
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the title of the TNamed.

void TText::SetMbTitle(const wchar_t *title)
{
   char *mb_title = new char[MB_CUR_MAX * wcslen(title) + 1]();
   char *p = mb_title;
   size_t length = wcslen(title);
   for (size_t i = 0; i < length; i++) {
      const int n = wctomb(p, title[i]);
      if (n >= 0) p += n;
   }
   fTitle = mb_title;
   delete [] mb_title;
   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TText.

void TText::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TText::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      TAttText::Streamer(R__b);
      Float_t x,y;
      R__b >> x; fX = x;
      R__b >> y; fY = y;
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TText::Class(),this);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Return the "bounding Box" of the Box

Rectangle_t TText::GetBBox()
{
   UInt_t w, h;
   Int_t Dx, Dy;
   Dx = Dy = 0;
   GetBoundingBox(w, h, false);

   Short_t halign = fTextAlign/10;
   Short_t valign = fTextAlign - 10*halign;

   switch (halign) {
      case 1 : Dx = 0      ; break;
      case 2 : Dx = w/2   ; break;
      case 3 : Dx = w     ; break;
   }
   switch (valign) {
      case 1 : Dy = h     ; break;
      case 2 : Dy = h/2   ; break;
      case 3 : Dy = 0      ; break;
   }

   Rectangle_t BBox;
   BBox.fX = gPad->XtoPixel(fX)-Dx;
   BBox.fY = gPad->YtoPixel(fY)-Dy;
   BBox.fWidth  = w;
   BBox.fHeight = h;
   return (BBox);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the point given by Alignment as 'center'

TPoint TText::GetBBoxCenter()
{
   TPoint p;
   p.SetX(gPad->XtoPixel(fX));
   p.SetY(gPad->YtoPixel(fY));
   return(p);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the point given by Alignment as 'center'

void TText::SetBBoxCenter(const TPoint &p)
{
   this->SetX(gPad->PixeltoX(p.GetX()));
   this->SetY(gPad->PixeltoY(p.GetY()-gPad->VtoPixel(0)));
}

////////////////////////////////////////////////////////////////////////////////
/// Set X coordinate of the point given by Alignment as 'center'

void TText::SetBBoxCenterX(const Int_t x)
{
   this->SetX(gPad->PixeltoX(x));
}

////////////////////////////////////////////////////////////////////////////////
/// Set Y coordinate of the point given by Alignment as 'center'

void TText::SetBBoxCenterY(const Int_t y)
{
   this->SetY(gPad->PixeltoY(y - gPad->VtoPixel(0)));
}

////////////////////////////////////////////////////////////////////////////////
/// Set left hand side of BoundingBox to a value
/// (resize in x direction on left)

void TText::SetBBoxX1(const Int_t /*x*/)
{
   //NOT IMPLEMENTED
}

////////////////////////////////////////////////////////////////////////////////
/// Set right hand side of BoundingBox to a value
/// (resize in x direction on right)

void TText::SetBBoxX2(const Int_t /*x*/)
{
   //NOT IMPLEMENTED
}

////////////////////////////////////////////////////////////////////////////////
/// Set top of BoundingBox to a value (resize in y direction on top)

void TText::SetBBoxY1(const Int_t /*y*/)
{
   //NOT IMPLEMENTED
}

////////////////////////////////////////////////////////////////////////////////
/// Set bottom of BoundingBox to a value
/// (resize in y direction on bottom)

void TText::SetBBoxY2(const Int_t /*y*/)
{
   //NOT IMPLEMENTED
}
