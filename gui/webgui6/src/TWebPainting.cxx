// Author:  Sergey Linev, GSI  10/04/2017

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebPainting.h"

#include "TColor.h"

/** \class TWebPainting
\ingroup webgui6

Object used to store paint operations and deliver them to JSROOT

*/


///////////////////////////////////////////////////////////////////////////////////////
/// Constructor

TWebPainting::TWebPainting()
{
   fLastFill.SetFillStyle(9999);
   fLastLine.SetLineWidth(-123);
   fLastMarker.SetMarkerStyle(9999);
}

///////////////////////////////////////////////////////////////////////////////////////
/// Add next custom operator to painting
/// Operations are separated by semicolons
/// Following operations are supported:
///   t - text
///   h - text coded into simple hex
///   r - rectangle
///   b - rectangular fill region
///   l - polyline
///   f - poly fill region
///   m - poly marker
///   z - line attributes
///   y - fill attributes
///   x - marker attributes
///   o - text attributes
///  After operation code optional arguments can be append like length of operation or coded text
///  Each operation may use data from binary float buffer

void TWebPainting::AddOper(const std::string &oper)
{
   if (!fOper.empty())
      fOper.append(";");
   fOper.append(oper);
}

///////////////////////////////////////////////////////////////////////////////////////
/// Create text operation
/// If text include special symbols - use simple hex coding

std::string TWebPainting::MakeTextOper(const char *str)
{
   bool isany_special = false;
   if (!str)
      str = "";
   for (auto p = str; *p; ++p) {
      if ((*p < 32) || (*p > 126) || (*p == ';') || (*p == '\'') || (*p == '\"') || (*p == '%')) {
         isany_special = true;
         break;
      }
   }

   if (!isany_special)
      return std::string("t") + str;

   // use simple hex coding while special symbols are hard to handle
   std::string oper("h");
   static const char *digits = "0123456789abcdef";
   for (auto p = str; *p; ++p) {
      unsigned code = (unsigned)*p;
      oper.append(1, digits[(code >> 4) & 0xF]);
      oper.append(1, digits[code & 0xF]);
   }
   return oper;
}

///////////////////////////////////////////////////////////////////////////////////////
/// Reserve place in the float buffer
/// Returns pointer on first element in reserved area

Float_t *TWebPainting::Reserve(Int_t sz)
{
   if (sz <= 0)
      return nullptr;

   if (fSize + sz > fBuf.GetSize()) {
      Int_t nextsz = fBuf.GetSize() + TMath::Max(1024, (sz/128 + 1) * 128);
      fBuf.Set(nextsz);
   }

   Float_t *res = fBuf.GetArray() + fSize;
   fSize += sz;
   return res; // return size where drawing can start
}

///////////////////////////////////////////////////////////////////////////////////////
/// Store line attributes
/// If attributes were not changed - ignore operation

void TWebPainting::AddLineAttr(const TAttLine &attr)
{
   if ((attr.GetLineColor() == fLastLine.GetLineColor()) &&
       (attr.GetLineStyle() == fLastLine.GetLineStyle()) &&
       (attr.GetLineWidth() == fLastLine.GetLineWidth())) return;

   fLastLine = attr;

   AddOper(std::string("z") +
           std::to_string((int) attr.GetLineColor()) + ":" +
           std::to_string((int) attr.GetLineStyle()) + ":" +
           std::to_string((int) attr.GetLineWidth()));
}

///////////////////////////////////////////////////////////////////////////////////////
/// Store fill attributes
/// If attributes were not changed - ignore operation

void TWebPainting::AddFillAttr(const TAttFill &attr)
{
   if ((fLastFill.GetFillColor() == attr.GetFillColor()) &&
       (fLastFill.GetFillStyle() == attr.GetFillStyle())) return;

   fLastFill = attr;

   AddOper(std::string("y") +
           std::to_string((int) attr.GetFillColor()) + ":" +
           std::to_string((int) attr.GetFillStyle()));
}

///////////////////////////////////////////////////////////////////////////////////////
/// Store text attributes
/// If attributes were not changed - ignore operation

void TWebPainting::AddTextAttr(const TAttText &attr)
{
   AddOper(std::string("o") +
           std::to_string((int) attr.GetTextColor()) + ":" +
           std::to_string((int) attr.GetTextFont()) + ":" +
           std::to_string((int) (attr.GetTextSize()>=1 ? attr.GetTextSize() : -1000*attr.GetTextSize())) + ":" +
           std::to_string((int) attr.GetTextAlign()) + ":" +
           std::to_string((int) attr.GetTextAngle()));
}

///////////////////////////////////////////////////////////////////////////////////////
/// Store marker attributes
/// If attributes were not changed - ignore operation

void TWebPainting::AddMarkerAttr(const TAttMarker &attr)
{
   if ((attr.GetMarkerColor() == fLastMarker.GetMarkerColor()) &&
       (attr.GetMarkerStyle() == fLastMarker.GetMarkerStyle()) &&
       (attr.GetMarkerSize() == fLastMarker.GetMarkerSize())) return;

   fLastMarker = attr;

   AddOper(std::string("x") +
           std::to_string((int) attr.GetMarkerColor()) + ":" +
           std::to_string((int) attr.GetMarkerStyle()) + ":" +
           std::to_string((int) attr.GetMarkerSize()));
}

///////////////////////////////////////////////////////////////////////////////////////
/// Add custom color to operations

void TWebPainting::AddColor(Int_t indx, TColor *col)
{
   if (!col) return;

   TString code;

   if (col->GetAlpha() == 1)
      code.Form("%d:%d,%d,%d", indx, (int) (255*col->GetRed()), (int) (255*col->GetGreen()), (int) (255*col->GetBlue()));
   else
      code.Form("%d=%d,%d,%d,%5.3f", indx, (int) (255*col->GetRed()), (int) (255*col->GetGreen()), (int) (255*col->GetBlue()), col->GetAlpha());

   AddOper(code.Data());
}

