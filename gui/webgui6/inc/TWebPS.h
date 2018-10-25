// Author:  Sergey Linev, GSI  23/10/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TWebPS
#define ROOT_TWebPS

#include "TVirtualPS.h"

#include "TWebPainting.h"

#include "TWebPadPainter.h"

class TWebPS : public TVirtualPS {

   TWebPainting *fPainting{nullptr};      ///!< object to store all painting

   enum EAttrKinds { attrLine = 0x1, attrFill = 0x2, attrMarker = 0x4, attrText = 0x8 };

   Float_t *StoreOperation(const std::string &oper, unsigned attrkind, int opersize = 0);

public:
   TWebPS() {}
   virtual ~TWebPS();

   TWebPainting *TakePainting();
   void ResetPainting();

   //Redirect calls to WebPainter
   //Line attributes.
/*   Color_t  GetLineColor() const override { return fPainter.GetLineColor(); }
   Style_t  GetLineStyle() const override { return fPainter.GetLineStyle(); }
   Width_t  GetLineWidth() const override { return fPainter.GetLineWidth(); }

   void     SetLineColor(Color_t lcolor) override { fPainter.SetLineColor(lcolor); }
   void     SetLineStyle(Style_t lstyle) override { fPainter.SetLineStyle(lstyle); }
   void     SetLineWidth(Width_t lwidth) override { fPainter.SetLineWidth(lwidth); }

   //Fill attributes.
   Color_t  GetFillColor() const override { return fPainter.GetFillColor(); }
   Style_t  GetFillStyle() const override { return fPainter.GetFillStyle(); }
   Bool_t   IsTransparent() const override { return fPainter.IsTransparent(); }

   void     SetFillColor(Color_t fcolor)  override { fPainter.SetFillColor(fcolor); }
   void     SetFillStyle(Style_t fstyle)  override { fPainter.SetFillStyle(fstyle); }
   void     SetOpacity(Int_t percent) { fPainter.SetOpacity(percent); }

   //Text attributes.
   Short_t  GetTextAlign() const override { return fPainter.GetTextAlign(); }
   Float_t  GetTextAngle() const override { return fPainter.GetTextAngle(); }
   Color_t  GetTextColor() const override { return fPainter.GetTextColor(); }
   Font_t   GetTextFont()  const override { return fPainter.GetTextFont(); }
   Float_t  GetTextSize()  const override { return fPainter.GetTextSize(); }
   Float_t  GetTextMagnitude() const { return fPainter.GetTextMagnitude(); }

   void     SetTextAlign(Short_t align) override { fPainter.SetTextAlign(align); }
   void     SetTextAngle(Float_t tangle) override { fPainter.SetTextAngle(tangle); }
   void     SetTextColor(Color_t tcolor) override { fPainter.SetTextColor(tcolor); }
   void     SetTextFont(Font_t tfont) override { fPainter.SetTextFont(tfont); }
   void     SetTextSize(Float_t tsize) override { fPainter.SetTextSize(tsize); }
   void     SetTextSizePixels(Int_t npixels) override { fPainter.SetTextSizePixels(npixels); }

   //MISSING in base class - Marker attributes

   Color_t   GetMarkerColor() const override  { return fPainter.GetMarkerColor(); }
   Size_t    GetMarkerSize() const override  { return fPainter.GetMarkerSize(); }
   Style_t   GetMarkerStyle() const override  { return fPainter.GetMarkerStyle(); }

   void      SetMarkerColor(Color_t cindex) override  { fPainter.SetMarkerColor(cindex); }
   void      SetMarkerSize(Float_t markersize) override  { fPainter.SetMarkerSize(markersize); }
   void      SetMarkerStyle(Style_t markerstyle) override  { fPainter.SetMarkerStyle(markerstyle); }

*/

   /// not yet implemented

   void CellArrayBegin(Int_t, Int_t, Double_t, Double_t, Double_t, Double_t) override  {}
   void CellArrayFill(Int_t, Int_t, Int_t)  override {}
   void CellArrayEnd()  override  {}
   void Close(Option_t * = "")  override  {}
   void DrawFrame(Double_t, Double_t, Double_t, Double_t, Int_t, Int_t, Int_t, Int_t) override {}
   void NewPage() override {}
   void Open(const char *, Int_t = -111) override {}
   void SetColor(Float_t, Float_t, Float_t) override {}


   // overwritten methods
   void DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2) override;
   void DrawPolyMarker(Int_t n, Float_t *x, Float_t *y) override;
   void DrawPolyMarker(Int_t n, Double_t *x, Double_t *y) override;
   void DrawPS(Int_t n, Float_t *xw, Float_t *yw) override;
   void DrawPS(Int_t n, Double_t *xw, Double_t *yw) override;
   void Text(Double_t x, Double_t y, const char *str) override;
   void Text(Double_t x, Double_t y, const wchar_t *str) override;

   ClassDefOverride(TWebPS, 0) // Redirection of VirtualPS to Web painter
};

#endif
