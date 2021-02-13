// @(#)root/gui:$Id$
// Author: Fons Rademakers   10/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGProgressBar, TGHProgressBar and TGVProgressBar                     //
//                                                                      //
// The classes in this file implement progress bars. Progress bars can  //
// be used to show progress of tasks taking more then a few seconds.    //
// TGProgressBar is an abstract base class, use either TGHProgressBar   //
// or TGVProgressBar. TGHProgressBar can in addition show the position  //
// as text in the bar.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGProgressBar.h"
#include "TGResourcePool.h"
#include "TColor.h"
#include "TVirtualX.h"

#include <iostream>

const TGFont *TGProgressBar::fgDefaultFont = nullptr;
TGGC         *TGProgressBar::fgDefaultGC = nullptr;


ClassImp(TGProgressBar);
ClassImp(TGHProgressBar);
ClassImp(TGVProgressBar);

////////////////////////////////////////////////////////////////////////////////
/// Create progress bar.

TGProgressBar::TGProgressBar(const TGWindow *p, UInt_t w, UInt_t h,
                             ULong_t back, ULong_t barcolor, GContext_t norm,
                             FontStruct_t font, UInt_t options) :
   TGFrame(p, w, h, options | kOwnBackground, back)
{
   fMin        = 0;
   fMax        = 100;
   fPos        = 0;
   fPosPix     = 0;
   fFillType   = kSolidFill;
   fBarType    = kStandard;
   fShowPos    = kFALSE;
   fPercent    = kTRUE;
   fNormGC     = norm;
   fFontStruct = font;
   fBarColorGC.SetFillStyle(kFillSolid);
   fBarColorGC.SetForeground(barcolor);
   fBarWidth   = kProgressBarStandardWidth;
   fDrawBar    = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set min and max of progress bar.

void TGProgressBar::SetRange(Float_t min, Float_t max)
{
   if (min >= max) {
      Error("SetRange", "max must be > min");
      return;
   }

   Bool_t draw = kFALSE;
   if (fPos > fMin) {
      // already in progress... rescale
      if (fPos < min) fPos = min;
      if (fPos > max) fPos = max;
      draw = kTRUE;
   } else
      fPos = min;

   fMin = min;
   fMax = max;

   if (draw)
      DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Set progress position between [min,max].

void TGProgressBar::SetPosition(Float_t pos)
{
   if (pos < fMin) pos = fMin;
   if (pos > fMax) pos = fMax;

   if (fPos == pos)
      return;

   fPos = pos;

   //fClient->NeedRedraw(this);
   fDrawBar = kTRUE;
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Increment progress position.

void TGProgressBar::Increment(Float_t inc)
{
   if (fPos == fMax)
      return;

   fPos += inc;
   if (fPos > fMax) fPos = fMax;

   //fClient->NeedRedraw(this);
   fDrawBar = kTRUE;
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset progress bar (i.e. set pos to 0).

void TGProgressBar::Reset()
{
   fPos = 0;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill type.

void TGProgressBar::SetFillType(EFillType type)
{
   fFillType = type;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set bar type.

void TGProgressBar::SetBarType(EBarType type)
{
   fBarType = type;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set progress bar color.

void TGProgressBar::SetBarColor(ULong_t color)
{
   fBarColorGC.SetForeground(color);

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set progress bar color.

void TGProgressBar::SetBarColor(const char *color)
{
   ULong_t ic;
   fClient->GetColorByName(color, ic);
   fBarColorGC.SetForeground(ic);
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set format for displaying a value.

void TGProgressBar::Format(const char *format)
{
   fFormat = format;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure in use.

FontStruct_t TGProgressBar::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context in use.

const TGGC &TGProgressBar::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = new TGGC(*gClient->GetResourcePool()->GetFrameGC());
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Change text color drawing.

void TGProgressBar::SetForegroundColor(Pixel_t pixel)
{
   TGGC *gc = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);

   if (!gc) {
      return;
   }
   gc->SetForeground(pixel);
   fNormGC = gc->GetGC();

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Horizontal progress bar constructor.

TGHProgressBar::TGHProgressBar(const TGWindow *p, UInt_t w, UInt_t h,
                              Pixel_t back, Pixel_t barcolor,
                              GContext_t norm, FontStruct_t font, UInt_t options) :
      TGProgressBar(p, w, h, back, barcolor, norm, font, options)
{
   fBarWidth = h;
   fEditDisabled = kEditDisableHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Simple constructor allow you to create either a standard progress
/// bar, or a more fancy progress bar (fancy means: double sized border,
/// white background and a bit wider to allow for text to be printed
/// in the bar.

TGHProgressBar::TGHProgressBar(const TGWindow *p, EBarType type, UInt_t w)
   : TGProgressBar(p, w, type == kStandard ? kProgressBarStandardWidth :
                   kProgressBarTextWidth, type == kStandard ? GetDefaultFrameBackground() :
                   fgWhitePixel, fgDefaultSelectedBackground, GetDefaultGC()(),
                   GetDefaultFontStruct(),
                   type == kStandard ? kSunkenFrame : kDoubleBorder | kSunkenFrame)
{
   fBarType  = type;
   fBarWidth = (type == kStandard) ? kProgressBarStandardWidth : kProgressBarTextWidth;
   fEditDisabled = kEditDisableHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Show postion text, either in percent or formatted according format.

void TGHProgressBar::ShowPosition(Bool_t set, Bool_t percent, const char *format)
{
   fShowPos = set;
   fPercent = percent;
   fFormat  = format;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw horizontal progress bar.

void TGHProgressBar::DoRedraw()
{
   if (!fDrawBar) {
      // calls TGProgressBar::DrawBorder()
      TGFrame::DoRedraw();
   }

   fPosPix = (Int_t)((fWidth - 2.0*fBorderWidth ) *(fPos - fMin) / (fMax - fMin)) + fBorderWidth;

   Int_t pospix = fPosPix;

   if (fFillType == kSolidFill)
      gVirtualX->FillRectangle(fId, fBarColorGC(), fBorderWidth,
                               fBorderWidth, fPosPix - fBorderWidth, fBarWidth -
                               (fBorderWidth << 1));
   else {
      Int_t blocksize = kBlockSize;
      Int_t delta     = kBlockSpace;
      Int_t pos       = fBorderWidth;
      while (pos < fPosPix) {
         if (pos + blocksize > Int_t(fWidth)-fBorderWidth)
            blocksize = fWidth-fBorderWidth-pos;
         gVirtualX->FillRectangle(fId, fBarColorGC(), pos,
                                  fBorderWidth, blocksize, fBarWidth -
                                  (fBorderWidth << 1));
         if (fDrawBar && fShowPos)
            gVirtualX->ClearArea(fId, pos+blocksize, fBorderWidth,
                                 delta, fBarWidth - (fBorderWidth << 1));

         pos += blocksize + delta;
      }
      pospix = pos - delta;
   }

   if (fShowPos) {
      TString buf;
      if (fPercent)
         buf = TString::Format("%d%%", Int_t((fPos-fMin)/(fMax-fMin)*100.));
      else
         buf = TString::Format(fFormat.Data(), fPos);

      Int_t x, y, max_ascent, max_descent;
      UInt_t twidth  = gVirtualX->TextWidth(fFontStruct, buf.Data(), buf.Length());
      gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
      UInt_t theight = max_ascent + max_descent;

      x = (Int_t)((fWidth - twidth)*0.5);
      y = (Int_t)((fHeight - theight)*0.5);

      if (fDrawBar && fPosPix < Int_t(x+twidth))
         gVirtualX->ClearArea(fId, pospix, fBorderWidth,
                              fWidth - pospix - fBorderWidth,
                              fBarWidth - (fBorderWidth << 1));

      gVirtualX->DrawString(fId, fNormGC, x, y + max_ascent, buf.Data(), buf.Length());
   }

   fDrawBar = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// cconstructor

TGVProgressBar::TGVProgressBar(const TGWindow *p, UInt_t w, UInt_t h,
                              Pixel_t back, Pixel_t barcolor, GContext_t norm,
                              FontStruct_t font,UInt_t options) :
      TGProgressBar(p, w, h, back, barcolor, norm, font, options)
{
   fBarWidth = w;
   fEditDisabled = kEditDisableWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Simple constructor allow you to create either a standard progress
/// bar, or a more fancy progress bar (fancy means: double sized border,
/// white background and a bit wider to allow for text to be printed
/// in the bar.

TGVProgressBar::TGVProgressBar(const TGWindow *p, EBarType type, UInt_t h)
   : TGProgressBar(p, type == kStandard ? kProgressBarStandardWidth :
                   kProgressBarTextWidth, h, type == kStandard ? GetDefaultFrameBackground() :
                   fgWhitePixel, fgDefaultSelectedBackground, GetDefaultGC()(),
                   GetDefaultFontStruct(),
                   type == kStandard ? kSunkenFrame : kDoubleBorder | kSunkenFrame)
{
   fBarType  = type;
   fBarWidth = (type == kStandard) ? kProgressBarStandardWidth : kProgressBarTextWidth;
   fDrawBar  = kFALSE;
   fEditDisabled = kEditDisableWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw vertical progress bar.

void TGVProgressBar::DoRedraw()
{
   if (!fDrawBar) {
      // calls TGProgressBar::DrawBorder()
      TGFrame::DoRedraw();
   }

    fPosPix = (Int_t)((fHeight - 2.0f*fBorderWidth) *(fPos - fMin) / (fMax - fMin)) + fBorderWidth;

   if (fFillType == kSolidFill)
      gVirtualX->FillRectangle(fId, fBarColorGC(), fBorderWidth,
                               fHeight - fPosPix, fBarWidth - (fBorderWidth << 1),
                               fPosPix - fBorderWidth);
   else {
      Int_t blocksize = kBlockSize;
      Int_t delta     = kBlockSpace;
      Int_t pos       = fBorderWidth;
      while (pos < fPosPix) {
         if (pos + blocksize > Int_t(fHeight)-fBorderWidth)
            blocksize = fHeight-fBorderWidth-pos;
         gVirtualX->FillRectangle(fId, fBarColorGC(), fBorderWidth,
                                  fHeight - pos - blocksize, fBarWidth - (fBorderWidth << 1),
                                  blocksize);
         pos += blocksize + delta;
      }
   }

   if (fShowPos) {
      // not text shown for vertical progress bars
   }

   fDrawBar = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save progress bar parameters as a C++ statement(s) on output stream out.

void TGProgressBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   const char *barcolor;
   char quote = '"';
   switch (fBarType) {
      case kFancy:
         if (GetOptions() != (kSunkenFrame | kDoubleBorder | kOwnBackground))
            out << "   " << GetName() << "->ChangeOptions(" << GetOptionString()
                << ");" << std::endl;
         if (GetBackground() != GetWhitePixel()) {
            SaveUserColor(out, option);
            out << "   " << GetName() << "->SetBackgroundColor(ucolor);" << std::endl;
         }
         break;

      case kStandard:
         if (GetOptions() != (kSunkenFrame | kOwnBackground))
            out << "   " << GetName() << "->ChangeOptions(" << GetOptionString()
                << ");" << std::endl;
         if (GetBackground() != GetDefaultFrameBackground()) {
            SaveUserColor(out, option);
            out << "   " << GetName() << "->SetBackgroundColor(ucolor);" << std::endl;
         }
         break;
   }

   if (fBarColorGC.GetForeground() != GetDefaultSelectedBackground()) {
      barcolor = TColor::PixelAsHexString(fBarColorGC.GetForeground());
      out << "   " << GetName() <<"->SetBarColor(" << quote << barcolor << quote
          << ");"  << std::endl;
   }

   if (fMin != 0 && fMax != 100)
      out << "   " << GetName() << "->SetRange(" << fMin << "," << fMax << ");" << std::endl;

   out <<"   "<< GetName() <<"->SetPosition("<< fPos <<");"<< std::endl;

}

////////////////////////////////////////////////////////////////////////////////
/// Save a vertical progress bar as a C++ statement(s) on output stream out.

void TGVProgressBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{

   out << "   TGVProgressBar *";
   out << GetName() << " = new TGVProgressBar(" << fParent->GetName();

   if ((fBarType == kFancy) && (fBarWidth == kProgressBarTextWidth)) {
      out << ",TGProgressBar::kFancy";
   } else if ((fBarType == kStandard) && (fBarWidth == kProgressBarStandardWidth)){
      out << ",TGProgressBar::kStandard";
   }

   out << "," << GetHeight() <<");" << std::endl;

   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (GetFillType() == kBlockFill)
      out << "   " << GetName() <<"->SetFillType(TGProgressBar::kBlockFill);"<< std::endl;

   TGProgressBar::SavePrimitive(out, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a horizontal progress bar as a C++ statement(s) on output stream out

void TGHProgressBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   out <<"   TGHProgressBar *";
   out << GetName() <<" = new TGHProgressBar("<< fParent->GetName();

   if ((fBarType == kFancy) && (fBarWidth == kProgressBarTextWidth)) {
      out << ",TGProgressBar::kFancy";
   } else if ((fBarType == kStandard) && (fBarWidth == kProgressBarStandardWidth)){
      out << ",TGProgressBar::kStandard";
   }

   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   out << "," << GetWidth() << ");" << std::endl;

   if (GetFillType() == kBlockFill)
      out << "   " << GetName() <<"->SetFillType(TGProgressBar::kBlockFill);"<< std::endl;

   if (GetShowPos()) {
      out << "   " << GetName() <<"->ShowPosition(kTRUE,";
      if (UsePercent()) {
         out << "kTRUE,";
      } else {
         out << "kFALSE,";
      }
      out << quote << GetFormat() << quote << ");"<< std::endl;

   } else if (UsePercent() && !GetFillType()) {
      out << "   " << GetName() <<"->ShowPosition();" << std::endl;
   }
   TGProgressBar::SavePrimitive(out, option);
}
