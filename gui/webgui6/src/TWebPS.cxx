// Author:  Sergey Linev, GSI  23/10/2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebPS.h"

#include <ROOT/RMakeUnique.hxx>

/** \class TWebPS
\ingroup webgui6

*/


TWebPS::TWebPS()
{
   CreatePainting();
}

void TWebPS::CreatePainting()
{
   fPainting = std::make_unique<TWebPainting>();
}


Float_t *TWebPS::StoreOperation(const std::string &oper, unsigned attrkind, int opersize)
{
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

void TWebPS::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
   Float_t *buf = (GetFillStyle() > 0) ? StoreOperation("b", attrFill, 4) : StoreOperation("r", attrLine, 4);

   buf[0] = x1;
   buf[1] = y1;
   buf[2] = x2;
   buf[3] = y2;
}

void TWebPS::DrawPolyMarker(Int_t nPoints, Float_t *x, Float_t *y)
{
   if (nPoints < 1) return;

   Float_t *buf = StoreOperation(std::string("m") + std::to_string(nPoints), attrLine | attrMarker, nPoints*2);

   for (Int_t n=0;n<nPoints;++n) {
      buf[n*2] = x[n];
      buf[n*2+1] = y[n];
   }
}

void TWebPS::DrawPolyMarker(Int_t nPoints, Double_t *x, Double_t *y)
{
   if (nPoints < 1) return;

   Float_t *buf = StoreOperation(std::string("m") + std::to_string(nPoints), attrLine | attrMarker, nPoints*2);

   for (Int_t n=0;n<nPoints;++n) {
      buf[n*2] = x[n];
      buf[n*2+1] = y[n];
   }
}

void TWebPS::DrawPS(Int_t nPoints, Float_t *xw, Float_t *yw)
{
   Float_t *buf = nullptr;
   if (nPoints < 0) {
      nPoints = -nPoints;
      if ((GetFillStyle() <= 0) || (nPoints < 3))  return;
      buf = StoreOperation(std::string("f") + std::to_string(nPoints), attrFill, nPoints*2);
   } else {
      if ((GetLineWidth() <= 0) || (nPoints < 2))  return;
      buf = StoreOperation(std::string("l") + std::to_string(nPoints), attrLine, nPoints*2);
   }
   for (Int_t n=0;n<nPoints;++n) {
      buf[n*2] = xw[n];
      buf[n*2+1] = yw[n];
   }
}

void TWebPS::DrawPS(Int_t nPoints, Double_t *xw, Double_t *yw)
{
   Float_t *buf = nullptr;
   if (nPoints < 0) {
      nPoints = -nPoints;
      if ((GetFillStyle() <= 0) || (nPoints < 3))  return;
      buf = StoreOperation(std::string("f") + std::to_string(nPoints), attrFill, nPoints*2);
   } else {
      if ((GetLineWidth() <= 0) || (nPoints < 2))  return;
      buf = StoreOperation(std::string("l") + std::to_string(nPoints), attrLine, nPoints*2);
   }
   for (Int_t n=0;n<nPoints;++n) {
      buf[n*2] = xw[n];
      buf[n*2+1] = yw[n];
   }
}

void TWebPS::Text(Double_t x, Double_t y, const char *str)
{
   Float_t *buf = StoreOperation(TWebPainting::MakeTextOper(str), attrText, 2);
   buf[0] = x;
   buf[1] = y;
}


void TWebPS::Text(Double_t x, Double_t y, const wchar_t *)
{
   Float_t *buf = StoreOperation(TWebPainting::MakeTextOper("wchar_t"), attrText, 2);
   buf[0] = x;
   buf[1] = y;
}
