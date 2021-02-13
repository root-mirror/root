// Author:  Sergey Linev, GSI  10/04/2017

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TWebPadPainter
#define ROOT_TWebPadPainter

#include "TVirtualPadPainter.h"

#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttText.h"
#include "TAttMarker.h"
#include <string>

#include "TWebPainting.h"

class TWebCanvas;

class TWebPadPainter : public TVirtualPadPainter, public TAttLine, public TAttFill, public TAttText, public TAttMarker {

friend class TWebCanvas;

protected:

   TWebPainting *fPainting{nullptr};      ///!< object to store all painting, owned by TWebPS object

   enum { attrLine = 0x1, attrFill = 0x2, attrMarker = 0x4, attrText = 0x8, attrAll = 0xf };

   Float_t *StoreOperation(const std::string &oper, unsigned attrkind, int opersize = 0);

public:

   TWebPadPainter() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300

   void SetPainting(TWebPainting *p) { fPainting = p; }

   //Final overrides for TVirtualPadPainter pure virtual functions.
   //1. Part, which simply catch attributes.
   //Line attributes.
   Color_t  GetLineColor() const override { return TAttLine::GetLineColor(); }
   Style_t  GetLineStyle() const override { return TAttLine::GetLineStyle(); }
   Width_t  GetLineWidth() const override { return TAttLine::GetLineWidth(); }

   void     SetLineColor(Color_t lcolor) override { TAttLine::SetLineColor(lcolor); }
   void     SetLineStyle(Style_t lstyle) override { TAttLine::SetLineStyle(lstyle); }
   void     SetLineWidth(Width_t lwidth) override { TAttLine::SetLineWidth(lwidth); }

   //Fill attributes.
   Color_t  GetFillColor() const override { return TAttFill::GetFillColor(); }
   Style_t  GetFillStyle() const override { return TAttFill::GetFillStyle(); }
   Bool_t   IsTransparent() const override { return TAttFill::IsTransparent(); }

   void     SetFillColor(Color_t fcolor) override { TAttFill::SetFillColor(fcolor); }
   void     SetFillStyle(Style_t fstyle) override { TAttFill::SetFillStyle(fstyle); }
   void     SetOpacity(Int_t percent) override { TAttFill::SetFillStyle(4000 + percent); }

   //Text attributes.
   Short_t  GetTextAlign() const override { return TAttText::GetTextAlign(); }
   Float_t  GetTextAngle() const override { return TAttText::GetTextAngle(); }
   Color_t  GetTextColor() const override { return TAttText::GetTextColor(); }
   Font_t   GetTextFont()  const override { return TAttText::GetTextFont(); }
   Float_t  GetTextSize()  const override { return TAttText::GetTextSize(); }
   Float_t  GetTextMagnitude() const override { return  0; }

   void     SetTextAlign(Short_t align) override { TAttText::SetTextAlign(align); }
   void     SetTextAngle(Float_t tangle) override { TAttText::SetTextAngle(tangle); }
   void     SetTextColor(Color_t tcolor) override { TAttText::SetTextColor(tcolor); }
   void     SetTextFont(Font_t tfont) override { TAttText::SetTextFont(tfont); }
   void     SetTextSize(Float_t tsize) override { TAttText::SetTextSize(tsize); }
   void     SetTextSizePixels(Int_t npixels) override { TAttText::SetTextSizePixels(npixels); }

   //2. "Off-screen management" part.
   Int_t    CreateDrawable(UInt_t, UInt_t) override { return -1; }
   void     ClearDrawable() override {}
   void     CopyDrawable(Int_t, Int_t, Int_t) override {}
   void     DestroyDrawable(Int_t) override {}
   void     SelectDrawable(Int_t) override {}

   //jpg, png, bmp, gif output.
   void     SaveImage(TVirtualPad *, const char *, Int_t) const override;

   //TASImage support (noop for a non-gl pad).
   void     DrawPixels(const unsigned char *pixelData, UInt_t width, UInt_t height,
                       Int_t dstX, Int_t dstY, Bool_t enableAlphaBlending) override;

   void     DrawLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2) override;
   void     DrawLineNDC(Double_t u1, Double_t v1, Double_t u2, Double_t v2) override;

   void     DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, EBoxMode mode) override;
   //TPad needs double and float versions.
   void     DrawFillArea(Int_t n, const Double_t *x, const Double_t *y) override;
   void     DrawFillArea(Int_t n, const Float_t *x, const Float_t *y) override;

   //TPad needs both double and float versions of DrawPolyLine.
   void     DrawPolyLine(Int_t n, const Double_t *x, const Double_t *y) override;
   void     DrawPolyLine(Int_t n, const Float_t *x, const Float_t *y) override;
   void     DrawPolyLineNDC(Int_t n, const Double_t *u, const Double_t *v) override;

   //TPad needs both versions.
   void     DrawPolyMarker(Int_t n, const Double_t *x, const Double_t *y) override;
   void     DrawPolyMarker(Int_t n, const Float_t *x, const Float_t *y) override;

   void     DrawText(Double_t x, Double_t y, const char *text, ETextMode mode) override;
   void     DrawText(Double_t x, Double_t y, const wchar_t *text, ETextMode mode) override;
   void     DrawTextNDC(Double_t u, Double_t v, const char *text, ETextMode mode) override;
   void     DrawTextNDC(Double_t u, Double_t v, const wchar_t *text, ETextMode mode) override;

private:
   //Let's make this clear:
   TWebPadPainter(const TWebPadPainter &rhs) = delete;
   TWebPadPainter(TWebPadPainter && rhs) = delete;
   TWebPadPainter & operator = (const TWebPadPainter &rhs) = delete;
   TWebPadPainter & operator = (TWebPadPainter && rhs) = delete;

   ClassDefOverride(TWebPadPainter, 0) //Abstract interface for painting in TPad
};

#endif
