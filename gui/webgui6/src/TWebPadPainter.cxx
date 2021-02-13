// Author:  Sergey Linev, GSI  10/04/2017

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebPadPainter.h"
#include "TError.h"
#include "TImage.h"
#include "TVirtualPad.h"
#include "TWebCanvas.h"
#include "TBufferJSON.h"

#include <ROOT/RMakeUnique.hxx>


/** \class TWebPadPainter
\ingroup webgui6
\brief Implement TVirtualPadPainter which abstracts painting operations.

TWebPadPainter tries to support old Paint methods of the ROOT classes.
Main classes (like histograms or graphs) should be painted on JavaScript side

*/


//////////////////////////////////////////////////////////////////////////
/// Store operation identifier with appropriate attributes

Float_t *TWebPadPainter::StoreOperation(const std::string &oper, unsigned attrkind, int opersize)
{
   if (!fPainting) return nullptr;

   if (attrkind & attrLine)
      fPainting->AddLineAttr(*this);

   if (attrkind & attrFill)
      fPainting->AddFillAttr(*this);

   if (attrkind & attrMarker)
      fPainting->AddMarkerAttr(*this);

   if (attrkind & attrText)
      fPainting->AddTextAttr(*this);

   fPainting->AddOper(oper);

   return fPainting->Reserve(opersize);
}

////////////////////////////////////////////////////////////////////////////////
///Noop, for non-gl pad TASImage calls gVirtualX->CopyArea.

void TWebPadPainter::DrawPixels(const unsigned char * /*pixelData*/, UInt_t /*width*/, UInt_t /*height*/,
                             Int_t /*dstX*/, Int_t /*dstY*/, Bool_t /*enableAlphaBlending*/)
{

}


////////////////////////////////////////////////////////////////////////////////
/// Paint a simple line.

void TWebPadPainter::DrawLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
   if (GetLineWidth() <= 0)
      return;

   auto buf = StoreOperation("l2", attrLine, 4);
   if (buf) {
      buf[0] = x1;
      buf[1] = y1;
      buf[2] = x2;
      buf[3] = y2;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Paint a simple line in normalized coordinates.

void TWebPadPainter::DrawLineNDC(Double_t u1, Double_t v1, Double_t u2, Double_t v2)
{
   if (GetLineWidth()<=0) return;

   ::Error("DrawLineNDC", "Not supported correctly");

   auto buf = StoreOperation("l2", attrLine, 4);
   if (buf) {
      buf[0] = u1;
      buf[1] = v1;
      buf[2] = u2;
      buf[3] = v2;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Paint a simple box.

void TWebPadPainter::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, EBoxMode mode)
{
   if (GetLineWidth()<=0 && mode == TVirtualPadPainter::kHollow) return;

   Float_t *buf = nullptr;

   if (mode == TVirtualPadPainter::kHollow)
      buf = StoreOperation("r", attrLine, 4); // only border
   else
      buf = StoreOperation("b", attrFill, 4); // only fill

   if (buf) {
      buf[0] = x1;
      buf[1] = y1;
      buf[2] = x2;
      buf[3] = y2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint filled area.

void TWebPadPainter::DrawFillArea(Int_t nPoints, const Double_t *xs, const Double_t *ys)
{
   if ((GetFillStyle() <= 0) || (nPoints < 3))
      return;

   auto buf = StoreOperation("f" + std::to_string(nPoints), attrFill, nPoints * 2);
   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = xs[n];
         buf[n * 2 + 1] = ys[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint filled area.

void TWebPadPainter::DrawFillArea(Int_t nPoints, const Float_t *xs, const Float_t *ys)
{
   if ((GetFillStyle() <= 0) || (nPoints < 3))
      return;

   auto buf = StoreOperation("f" + std::to_string(nPoints), attrFill, nPoints * 2);
   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = xs[n];
         buf[n * 2 + 1] = ys[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint Polyline.

void TWebPadPainter::DrawPolyLine(Int_t nPoints, const Double_t *xs, const Double_t *ys)
{
   if ((GetLineWidth() <= 0) || (nPoints < 2))
      return;

   auto buf = StoreOperation("l" + std::to_string(nPoints), attrLine, nPoints * 2);
   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = xs[n];
         buf[n * 2 + 1] = ys[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint polyline.

void TWebPadPainter::DrawPolyLine(Int_t nPoints, const Float_t *xs, const Float_t *ys)
{
   if ((GetLineWidth() <= 0) || (nPoints < 2))
      return;

   auto buf = StoreOperation("l" + std::to_string(nPoints), attrLine, nPoints * 2);
   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = xs[n];
         buf[n * 2 + 1] = ys[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint polyline in normalized coordinates.

void TWebPadPainter::DrawPolyLineNDC(Int_t nPoints, const Double_t *u, const Double_t *v)
{
   if ((GetLineWidth() <= 0) || (nPoints < 2))
      return;

   ::Error("DrawPolyLineNDC", "Not supported correctly");

   auto buf = StoreOperation("l" + std::to_string(nPoints), attrLine, nPoints * 2);
   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = u[n];
         buf[n * 2 + 1] = v[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint polymarker.

void TWebPadPainter::DrawPolyMarker(Int_t nPoints, const Double_t *x, const Double_t *y)
{
   if (nPoints < 1)
      return;

   auto buf = StoreOperation(std::string("m") + std::to_string(nPoints), attrLine | attrMarker, nPoints * 2);

   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = x[n];
         buf[n * 2 + 1] = y[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint polymarker.

void TWebPadPainter::DrawPolyMarker(Int_t nPoints, const Float_t *x, const Float_t *y)
{
   if (nPoints < 1)
      return;

   auto buf = StoreOperation(std::string("m") + std::to_string(nPoints), attrLine | attrMarker, nPoints * 2);

   if (buf)
      for (Int_t n = 0; n < nPoints; ++n) {
         buf[n * 2] = x[n];
         buf[n * 2 + 1] = y[n];
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint text.

void TWebPadPainter::DrawText(Double_t x, Double_t y, const char *text, ETextMode /*mode*/)
{
   auto buf = StoreOperation(TWebPainting::MakeTextOper(text), attrText, 2);
   if (buf) {
      buf[0] = x;
      buf[1] = y;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Special version working with wchar_t and required by TMathText.

void TWebPadPainter::DrawText(Double_t x, Double_t y, const wchar_t * /*text*/, ETextMode /*mode*/)
{
   auto buf = StoreOperation(TWebPainting::MakeTextOper("wchar_t"), attrText, 2);
   if (buf) {
      buf[0] = x;
      buf[1] = y;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint text in normalized coordinates.

void TWebPadPainter::DrawTextNDC(Double_t u, Double_t v, const char *text, ETextMode /*mode*/)
{
   ::Error("DrawTextNDC", "Not supported correctly");

   auto buf = StoreOperation(TWebPainting::MakeTextOper(text), attrText, 2);

   if (buf) {
      buf[0] = u;
      buf[1] = v;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint text in normalized coordinates.

void TWebPadPainter::DrawTextNDC(Double_t  u , Double_t v, const wchar_t * /*text*/, ETextMode /*mode*/)
{
   ::Error("DrawTextNDC", "Not supported correctly");

   auto buf = StoreOperation(TWebPainting::MakeTextOper("wchar_t"), attrText, 2);

   if (buf) {
      buf[0] = u;
      buf[1] = v;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Produce image from WebPadPainter

void TWebPadPainter::SaveImage(TVirtualPad *pad, const char *fileName, Int_t gtype) const
{
   if ((gtype == TImage::kPng) || (gtype == TImage::kJpeg))
      TWebCanvas::ProduceImage(pad->GetCanvas(), fileName);
}

