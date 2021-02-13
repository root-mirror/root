// @(#)root/gui:$Id$
// Author: Fons Rademakers   20/9/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGGC and TGGCPool                                                    //
//                                                                      //
// Encapsulate a graphics context used in the low level graphics.       //
// TGGCPool provides a pool of graphics contexts.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGClient.h"
#include "TGGC.h"
#include "TVirtualX.h"
#include "THashTable.h"
#include "TColor.h"
#include "TROOT.h"

#include <iostream>
#include <cstring>


ClassImp(TGGC);

////////////////////////////////////////////////////////////////////////////////
/// Create a graphics context (only called via TGGCPool::GetGC()).

TGGC::TGGC(GCValues_t *values, Bool_t)
{
   fContext = 0;
   if (values) {
      fValues = *values;
      fContext = gVirtualX->CreateGC(gVirtualX->GetDefaultRootWindow(), values);
      if (values->fMask & kGCDashList) {
         if (values->fDashLen > (Int_t)sizeof(fValues.fDashes))
            Warning("TGGC", "dash list can have only up to %ld elements",
                    (Long_t)sizeof(fValues.fDashes));
         fValues.fDashLen = TMath::Min(values->fDashLen, (Int_t)sizeof(fValues.fDashes));
         gVirtualX->SetDashes(fContext, fValues.fDashOffset, fValues.fDashes,
                              fValues.fDashLen);
      }
   } else {
      fValues = {};
      fContext = 0;
   }
   SetRefCount(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a graphics context, registers GC in GCPool.

TGGC::TGGC(GCValues_t *values)
{
   fContext = 0;
   // case of default ctor at program startup before gClient exists
   if (!values) {
      fValues = {};
      fContext = 0;
      SetRefCount(1);
      return;
   }

   if (gClient)
      gClient->GetGC(values, kTRUE);
   else {
      fContext = 0;
      Error("TGGC", "TGClient not yet initialized, should never happen");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a graphics context.

TGGC::TGGC(const TGGC &g) : TObject(g), TRefCnt()
{
   fValues = g.fValues;
   if (g.fContext) {
      fContext = gVirtualX->CreateGC(gVirtualX->GetDefaultRootWindow(), &fValues);
      if (fValues.fMask & kGCDashList)
         gVirtualX->SetDashes(fContext, fValues.fDashOffset, fValues.fDashes,
                              fValues.fDashLen);
   } else
      fContext = 0;
   SetRefCount(1);

   if (gClient)
      gClient->GetGCPool()->fList->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete graphics context.

TGGC::~TGGC()
{
   if (gClient)
      gClient->GetGCPool()->ForceFreeGC(this);

   if (fContext)
      gVirtualX->DeleteGC(fContext);
}

////////////////////////////////////////////////////////////////////////////////
/// Graphics context assignment operator.

TGGC &TGGC::operator=(const TGGC &rhs)
{
   if (this != &rhs) {
      if (!fContext && gClient) {
         TGGC *gc = gClient->GetGCPool()->FindGC(this);
         if (!gc)
            gClient->GetGCPool()->fList->Add(this);
      }
      if (fContext)
         gVirtualX->DeleteGC(fContext);
      TObject::operator=(rhs);
      fValues  = rhs.fValues;
      fContext = gVirtualX->CreateGC(gVirtualX->GetDefaultRootWindow(), &fValues);
      if (fValues.fMask & kGCDashList)
         gVirtualX->SetDashes(fContext, fValues.fDashOffset, fValues.fDashes,
                              fValues.fDashLen);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Not inline due to a bug in g++ 2.96 20000731 (Red Hat Linux 7.0).

GContext_t TGGC::operator()() const
{
   return fContext;
}

////////////////////////////////////////////////////////////////////////////////
/// Update values + mask.

void TGGC::UpdateValues(GCValues_t *values)
{
   fValues.fMask |= values->fMask;

   for (Mask_t bit = 1; bit <= fValues.fMask; bit <<= 1) {
      switch (bit & values->fMask) {
         default:
         case 0:
            continue;
         case kGCFunction:
            fValues.fFunction = values->fFunction;
            break;
         case kGCPlaneMask:
            fValues.fPlaneMask = values->fPlaneMask;
            break;
         case kGCForeground:
            fValues.fForeground = values->fForeground;
            break;
         case kGCBackground:
            fValues.fBackground = values->fBackground;
            break;
         case kGCLineWidth:
            fValues.fLineWidth = values->fLineWidth;
            break;
         case kGCLineStyle:
            fValues.fLineStyle = values->fLineStyle;
            break;
         case kGCCapStyle:
            fValues.fCapStyle = values->fCapStyle;
            break;
         case kGCJoinStyle:
            fValues.fJoinStyle = values->fJoinStyle;
            break;
         case kGCFillStyle:
            fValues.fFillStyle = values->fFillStyle;
            break;
         case kGCFillRule:
            fValues.fFillRule = values->fFillRule;
            break;
         case kGCTile:
            fValues.fTile = values->fTile;
            break;
         case kGCStipple:
            fValues.fStipple = values->fStipple;
            break;
         case kGCTileStipXOrigin:
            fValues.fTsXOrigin = values->fTsXOrigin;
            break;
         case kGCTileStipYOrigin:
            fValues.fTsYOrigin = values->fTsYOrigin;
            break;
         case kGCFont:
            fValues.fFont = values->fFont;
            break;
         case kGCSubwindowMode:
            fValues.fSubwindowMode = values->fSubwindowMode;
            break;
         case kGCGraphicsExposures:
            fValues.fGraphicsExposures = values->fGraphicsExposures;
            break;
         case kGCClipXOrigin:
            fValues.fClipXOrigin = values->fClipXOrigin;
            break;
         case kGCClipYOrigin:
            fValues.fClipYOrigin = values->fClipYOrigin;
            break;
         case kGCClipMask:
            fValues.fClipMask = values->fClipMask;
            break;
         case kGCDashOffset:
            fValues.fDashOffset = values->fDashOffset;
            break;
         case kGCDashList:
            if (values->fDashLen > (Int_t)sizeof(fValues.fDashes))
               Warning("UpdateValues", "dash list can have only up to %ld elements",
                       (Long_t)sizeof(fValues.fDashes));
            fValues.fDashLen = TMath::Min(values->fDashLen, (Int_t)sizeof(fValues.fDashes));
            memcpy(fValues.fDashes, values->fDashes, fValues.fDashLen);
            break;
         case kGCArcMode:
            fValues.fArcMode = values->fArcMode;
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set attributes as specified in the values structure.

void TGGC::SetAttributes(GCValues_t *values)
{
   if (!fContext && gClient) {
      TGGC *gc = gClient->GetGCPool()->FindGC(this);
      if (!gc)
         gClient->GetGCPool()->fList->Add(this);
   }

   if (fContext)
      gVirtualX->ChangeGC(fContext, values);
   else
      fContext = gVirtualX->CreateGC(gVirtualX->GetDefaultRootWindow(), values);
   UpdateValues(values);
   if (values->fMask & kGCDashList)
      gVirtualX->SetDashes(fContext, fValues.fDashOffset, fValues.fDashes,
                           fValues.fDashLen);
}

////////////////////////////////////////////////////////////////////////////////
/// Set graphics context drawing function.

void TGGC::SetFunction(EGraphicsFunction v)
{
   GCValues_t values;
   values.fFunction = v;
   values.fMask     = kGCFunction;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set plane mask.

void TGGC::SetPlaneMask(ULong_t v)
{
   GCValues_t values;
   values.fPlaneMask = v;
   values.fMask      = kGCPlaneMask;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set foreground color.

void TGGC::SetForeground(ULong_t v)
{
   GCValues_t values;
   values.fForeground = v;
   values.fMask       = kGCForeground;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set background color.

void TGGC::SetBackground(ULong_t v)
{
   GCValues_t values;
   values.fBackground = v;
   values.fMask       = kGCBackground;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line width.

void TGGC::SetLineWidth(Int_t v)
{
   GCValues_t values;
   values.fLineWidth = v;
   values.fMask      = kGCLineWidth;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line style (kLineSolid, kLineOnOffDash, kLineDoubleDash).

void TGGC::SetLineStyle(Int_t v)
{
   GCValues_t values;
   values.fLineStyle = v;
   values.fMask      = kGCLineStyle;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set cap style (kCapNotLast, kCapButt, kCapRound, kCapProjecting).

void TGGC::SetCapStyle(Int_t v)
{
   GCValues_t values;
   values.fCapStyle = v;
   values.fMask     = kGCCapStyle;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line join style (kJoinMiter, kJoinRound, kJoinBevel).

void TGGC::SetJoinStyle(Int_t v)
{
   GCValues_t values;
   values.fJoinStyle = v;
   values.fMask      = kGCJoinStyle;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill style (kFillSolid, kFillTiled, kFillStippled,
/// kFillOpaeueStippled).

void TGGC::SetFillStyle(Int_t v)
{
   GCValues_t values;
   values.fFillStyle = v;
   values.fMask      = kGCFillStyle;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill rule (kEvenOddRule, kWindingRule).

void TGGC::SetFillRule(Int_t v)
{
   GCValues_t values;
   values.fFillRule = v;
   values.fMask     = kGCFillRule;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tile pixmap for tiling operations.

void TGGC::SetTile(Pixmap_t v)
{
   GCValues_t values;
   values.fTile = v;
   values.fMask = kGCTile;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set 1 plane pixmap for stippling.

void TGGC::SetStipple(Pixmap_t v)
{
   GCValues_t values;
   values.fStipple = v;
   values.fMask    = kGCStipple;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// X offset for tile or stipple operations.

void TGGC::SetTileStipXOrigin(Int_t v)
{
   GCValues_t values;
   values.fTsXOrigin = v;
   values.fMask      = kGCTileStipXOrigin;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Y offset for tile or stipple operations.

void TGGC::SetTileStipYOrigin(Int_t v)
{
   GCValues_t values;
   values.fTsYOrigin = v;
   values.fMask      = kGCTileStipYOrigin;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set font.

void TGGC::SetFont(FontH_t v)
{
   GCValues_t values;
   values.fFont = v;
   values.fMask = kGCFont;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set sub window mode (kClipByChildren, kIncludeInferiors).

void TGGC::SetSubwindowMode(Int_t v)
{
   GCValues_t values;
   values.fSubwindowMode = v;
   values.fMask          = kGCSubwindowMode;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// True if graphics exposure should be generated.

void TGGC::SetGraphicsExposures(Bool_t v)
{
   GCValues_t values;
   values.fGraphicsExposures = v;
   values.fMask              = kGCGraphicsExposures;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// X origin for clipping.

void TGGC::SetClipXOrigin(Int_t v)
{
   GCValues_t values;
   values.fClipXOrigin = v;
   values.fMask        = kGCClipXOrigin;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Y origin for clipping.

void TGGC::SetClipYOrigin(Int_t v)
{
   GCValues_t values;
   values.fClipYOrigin = v;
   values.fMask        = kGCClipYOrigin;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Bitmap for clipping.

void TGGC::SetClipMask(Pixmap_t v)
{
   GCValues_t values;
   values.fClipMask = v;
   values.fMask     = kGCClipMask;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Patterned/dashed line offset.

void TGGC::SetDashOffset(Int_t v)
{
   GCValues_t values;
   values.fDashOffset = v;
   values.fMask       = kGCDashOffset;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dash pattern. First use SetDashOffset() if not 0.

void TGGC::SetDashList(const char v[], Int_t len)
{
   GCValues_t values;
   if (len > (Int_t)sizeof(values.fDashes))
      Warning("SetDashList", "dash list can have only up to %ld elements",
              (Long_t)sizeof(values.fDashes));
   values.fDashLen = TMath::Min(len, (Int_t)sizeof(values.fDashes));
   memcpy(values.fDashes, v, values.fDashLen);
   values.fMask    = kGCDashList;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Set arc mode (kArcChord, kArcPieSlice).

void TGGC::SetArcMode(Int_t v)
{
   GCValues_t values;
   values.fArcMode = v;
   values.fMask    = kGCArcMode;
   SetAttributes(&values);
}

////////////////////////////////////////////////////////////////////////////////
/// Print graphics contexts info.

void TGGC::Print(Option_t *) const
{
   Printf("TGGC: mask = %x, handle = %lx, ref cnt = %u", fValues.fMask,
          fContext, References());
}

////////////////////////////////////////////////////////////////////////////////
/// Returns GC mask as a string - used in SavePrimitive().

TString TGGC::GetMaskString() const
{
   TString mask;

   Mask_t fmask = GetMask();

   if (fmask & kGCFunction) {
      if (mask.Length() == 0) mask  = "kGCFunction";
      else                    mask += " | kGCFunction";
   }
   if (fmask & kGCPlaneMask) {
      if (mask.Length() == 0) mask  = "kGCPlaneMask";
      else                    mask += " | kGCPlaneMask";
   }
   if (fmask & kGCForeground) {
      if (mask.Length() == 0) mask  = "kGCForeground";
      else                    mask += " | kGCForeground";
   }
   if (fmask & kGCBackground) {
      if (mask.Length() == 0) mask  = "kGCBackground";
      else                    mask += " | kGCBackground";
   }
   if (fmask & kGCLineWidth) {
      if (mask.Length() == 0) mask  = "kGCLineWidth";
      else                    mask += " | kGCLineWidth";
   }
   if (fmask & kGCLineStyle) {
      if (mask.Length() == 0) mask  = "kGCLineStyle";
      else                    mask += " | kGCLineStyle";
   }
   if (fmask & kGCCapStyle) {
      if (mask.Length() == 0) mask  = "kGCCapStyle";
      else                    mask += " | kGCCapStyle";
   }
   if (fmask & kGCJoinStyle) {
      if (mask.Length() == 0) mask  = "kGCJoinStyle";
      else                    mask += " | kGCJoinStyle";
   }
   if (fmask & kGCFillStyle) {
      if (mask.Length() == 0) mask  = "kGCFillStyle";
      else                    mask += " | kGCFillStyle";
   }
   if (fmask & kGCFillRule) {
      if (mask.Length() == 0) mask  = "kGCFillRule";
      else                    mask += " | kGCFillRule";
   }
   if (fmask & kGCTile) {
      if (mask.Length() == 0) mask  = "kGCTile";
      else                    mask += " | kGCTile";
   }
   if (fmask & kGCStipple) {
      if (mask.Length() == 0) mask  = "kGCStipple";
      else                    mask += " | kGCStipple";
   }
   if (fmask & kGCTileStipXOrigin) {
      if (mask.Length() == 0) mask  = "kGCTileStipXOrigin";
      else                    mask += " | kGCTileStipXOrigin";
   }
   if (fmask & kGCTileStipYOrigin) {
      if (mask.Length() == 0) mask  = "kGCTileStipYOrigin";
      else                    mask += " | kGCTileStipYOrigin";
   }
   if (fmask & kGCFont) {
      if (mask.Length() == 0) mask  = "kGCFont";
      else                    mask += " | kGCFont";
   }
   if (fmask & kGCSubwindowMode) {
      if (mask.Length() == 0) mask  = "kGCSubwindowMode";
      else                    mask += " | kGCSubwindowMode";
   }
   if (fmask & kGCGraphicsExposures) {
      if (mask.Length() == 0) mask  = "kGCGraphicsExposures";
      else                    mask += " | kGCGraphicsExposures";
   }
   if (fmask & kGCClipXOrigin) {
      if (mask.Length() == 0) mask  = "kGCClipXOrigin";
      else                    mask += " | kGCClipXOrigin";
   }
   if (fmask & kGCClipYOrigin) {
      if (mask.Length() == 0) mask  = "kGCClipYOrigin";
      else                    mask += " | kGCClipYOrigin";
   }
   if (fmask & kGCClipMask) {
      if (mask.Length() == 0) mask  = "kGCClipMask";
      else                    mask += " | kGCClipMask";
   }
   if (fmask & kGCDashOffset) {
      if (mask.Length() == 0) mask  = "kGCDashOffset";
      else                    mask += " | kGCDashOffset";
   }
   if (fmask & kGCDashList) {
      if (mask.Length() == 0) mask  = "kGCDashList";
      else                    mask += " | kGCDashList";
   }
   if (fmask & kGCArcMode) {
      if (mask.Length() == 0) mask  = "kGCArcMode";
      else                    mask += " | kGCArcMode";
   }
   return mask;
}

////////////////////////////////////////////////////////////////////////////////
/// Save graphics context info as a C++ statement(s) on output stream out

void TGGC::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (gROOT->ClassSaved(TGGC::Class())) {
      out << std::endl;
   } else {
      //  declare graphics context object to reflect required user changes
      out << std::endl;
      out << "   TGGC   *uGC;           // will reflect user GC changes" << std::endl;
   }

   Mask_t fmask = GetMask();

   const char *colorname;
   TString valname;
   char quote ='"';
   ULong_t color;

   valname = TString::Format("val%s", option);

   out << "   // graphics context changes" << std::endl;
   //out << "   TGGC *uGC" << option << ";" << std::endl;
   out << "   GCValues_t " << valname.Data() << ";" << std::endl;
   out << "   " << valname.Data() << ".fMask = " << GetMaskString() << ";" << std::endl;

   for (Mask_t bit = 1; bit <= fmask; bit <<= 1) {
      switch (bit & fmask) {
         default:
         case 0:
            continue;
         case kGCFunction:
            out << "   " << valname.Data() << ".fFunction = ";
            switch (GetFunction()) {
               case kGXclear:
                  out << "kGXclear";
                  break;
               case kGXand:
                  out << "kGXand";
                  break;
               case kGXandReverse:
                  out << "kGXandReverse";
                  break;
               case kGXcopy:
                  out << "kGXcopy";
                  break;
               case kGXandInverted:
                  out << "kGXandInverted";
                  break;
               case kGXnoop:
                  out << "kGXnoop";
                  break;
               case kGXxor:
                  out << "kGXxor";
                  break;
               case kGXor:
                  out << "kGXor";
                  break;
               case kGXnor:
                  out << "kGXnor";
                  break;
               case kGXequiv:
                  out << "kGXequiv";
                  break;
               case kGXinvert:
                  out << "kGXinvert";
                  break;
               case kGXorReverse:
                  out << "kGXorReverse";
                  break;
               case kGXcopyInverted:
                  out << "kGXcopyInverted";
                  break;
               case kGXorInverted:
                  out << "kGXorInverted";
                  break;
               case kGXnand:
                  out << "kGXnand";
                  break;
               case kGXset:
                  out << "kGXset";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCPlaneMask:
            out << "   " << valname.Data() << ".fPlaneMask = " << GetPlaneMask() << ";" << std::endl;
            break;
         case kGCForeground:
            color = GetForeground();
            colorname = TColor::PixelAsHexString(color);
            out << "   gClient->GetColorByName(" << quote << colorname << quote
                << "," << valname.Data() << ".fForeground);" << std::endl;
            break;
         case kGCBackground:
            color = GetBackground();
            colorname = TColor::PixelAsHexString(color);
            out << "   gClient->GetColorByName(" << quote << colorname << quote
                << "," << valname.Data() << ".fBackground);" << std::endl;
            break;
         case kGCLineWidth:
            out << "   " << valname.Data() << ".fLineWidth = " << GetLineWidth() << ";" << std::endl;
            break;
         case kGCLineStyle:
            out << "   " << valname.Data() << ".fLineStyle = ";
            switch (GetLineStyle()) {
               case kLineSolid:
                  out << "kLineSolid";
                  break;
               case kLineOnOffDash:
                  out << "kLineOnOffDash";
                  break;
               case kLineDoubleDash:
                  out << "kLineDoubleDash";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCCapStyle:
            out << "   " << valname.Data() << ".fCapStyle = ";
            switch (GetCapStyle()) {
               case kCapNotLast:
                  out << "kCapNotLast";
                  break;
               case kCapButt:
                  out << "kCapButt";
                  break;
               case kCapRound:
                  out << "kCapRound";
                  break;
               case kCapProjecting:
                  out << "kCapProjecting";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCJoinStyle:
            out << "   " << valname.Data() << ".fJoinStyle = ";
            switch (GetJoinStyle()) {
               case kJoinMiter:
                  out << "kJoinMiter";
                  break;
               case kJoinRound:
                  out << "kJoinRound";
                  break;
               case kJoinBevel:
                  out << "kJoinBevel";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCFillStyle:
            out << "   " << valname.Data() << ".fFillStyle = ";
            switch (GetFillStyle()) {
               case kFillSolid:
                  out << "kFillSolid";
                  break;
               case kFillTiled:
                  out << "kFillTiled";
                  break;
               case kFillStippled:
                  out << "kFillStippled";
                  break;
               case kFillOpaqueStippled:
                  out << "kFillOpaqueStippled";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCFillRule:
            out << "   " << valname.Data() << ".fFillRule = ";
            switch (GetFillRule()) {
               case kEvenOddRule:
                  out << "kEvenOddRule";
                  break;
               case kWindingRule:
                  out << "kWindingRule";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCTile:
            out << "   " << valname.Data() << ".fTile = " << GetTile() << ";" << std::endl;
            break;
         case kGCStipple:
            out << "   " << valname.Data() << ".fStipple = " << GetStipple() << ";" << std::endl;
            break;
         case kGCTileStipXOrigin:
            out << "   " << valname.Data() << ".fTsXOrigin = " << GetTileStipXOrigin() << ";" << std::endl;
            break;
         case kGCTileStipYOrigin:
            out << "   " << valname.Data() << ".fTsYOrigin = " << GetTileStipYOrigin() << ";" << std::endl;
            break;
         case kGCFont:
            out << "   " << valname.Data() << ".fFont = ufont->GetFontHandle();" << std::endl;
            break;
         case kGCSubwindowMode:
            out << "   " << valname.Data() << ".fSubwindowMode = ";
            switch (GetSubwindowMode()) {
               case kClipByChildren:
                  out << "kClipByChildren";
                  break;
               case kIncludeInferiors:
                  out << "kIncludeInferiors";
                  break;
            }
            out << ";" << std::endl;
            break;
         case kGCGraphicsExposures:
            out << "   " << valname.Data() << ".fGraphicsExposures = ";
            if (GetGraphicsExposures())
               out << "kTRUE";
            else
               out << "kFALSE";
            out << ";" << std::endl;
            break;
         case kGCClipXOrigin:
            out << "   " << valname.Data() << ".fClipXOrigin = " << GetClipXOrigin() << ";" << std::endl;
            break;
         case kGCClipYOrigin:
            out << "   " << valname.Data() << ".fClipYOrigin = " << GetClipYOrigin() << ";" << std::endl;
            break;
         case kGCClipMask:
            out << "   " << valname.Data() << ".fClipMask = " << GetClipMask() << ";" << std::endl;
            break;
         case kGCDashOffset:
            out << "   " << valname.Data() << ".fDashOffset = " << GetDashOffset() << ";" << std::endl;
            break;
         case kGCDashList:
            if (GetDashLen() > (Int_t)sizeof(GetDashes()))
               Warning("TGGC::SavePrimitive", "dash list can have only up to %ld elements",
                       (Long_t)sizeof(GetDashes()));
            out << "   " << valname.Data() << ".fDashLen = "
                << TMath::Min(GetDashLen(),(Int_t)sizeof(GetDashes())) << ";" << std::endl;
            out << "   memcpy(GetDashes()," << valname.Data() << ".fDashes,"
                                            << valname.Data() << ".fDashLen);" << std::endl;
            break;
         case kGCArcMode:
            out << "   " << valname.Data() << ".fArcMode = ";
            switch (GetArcMode()) {
               case kArcChord:
                  out << "kArcChord";
                  break;
               case kArcPieSlice:
                  out << "kArcPieSlice";
                  break;
            }
            out << ";" << std::endl;
            break;
      }
   }
   out << "   uGC = gClient->GetGC(&" << valname.Data() << ", kTRUE);" << std::endl;
}


ClassImp(TGGCPool);

////////////////////////////////////////////////////////////////////////////////
/// Create graphics context pool.

TGGCPool::TGGCPool(TGClient *client)
{
   fClient = client;
   fList   = new THashTable;
   fList->SetOwner();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete graphics context pool.

TGGCPool::~TGGCPool()
{
   delete fList;
}

////////////////////////////////////////////////////////////////////////////////
/// Force remove graphics context from list. Is only called via ~TGGC().

void TGGCPool::ForceFreeGC(const TGGC *gct)
{
   TGGC *gc = (TGGC *) fList->FindObject(gct);

   if (gc) {
      if (gc->References() > 1)
         Error("ForceFreeGC", "removed a shared graphics context\n"
               "best to use graphics contexts via the TGGCPool()");
      fList->Remove(gc);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete graphics context if it is not used anymore.

void TGGCPool::FreeGC(const TGGC *gct)
{
   TGGC *gc = (TGGC *) fList->FindObject(gct);

   if (gc) {
      if (gc->RemoveReference() == 0) {
         fList->Remove(gc);
         delete gc;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete graphics context if it is not used anymore.

void TGGCPool::FreeGC(GContext_t gct)
{
   TIter next(fList);

   while (TGGC *gc = (TGGC *) next()) {
      if (gc->fContext == gct) {
         if (gc->RemoveReference() == 0) {
            fList->Remove(gc);
            delete gc;
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find graphics context. Returns 0 in case gc is not found.

TGGC *TGGCPool::FindGC(const TGGC *gct)
{
   return (TGGC*) fList->FindObject(gct);
}

////////////////////////////////////////////////////////////////////////////////
/// Find graphics context based on its GContext_t handle. Returns 0
/// in case gc is not found.

TGGC *TGGCPool::FindGC(GContext_t gct)
{
   TIter next(fList);

   while (TGGC *gc = (TGGC *) next()) {
      if (gc->fContext == gct)
         return gc;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// returns graphics context based on its GContext_t handle.

TGGC *TGGCPool::GetGC(GContext_t gct)
{
   GCValues_t gval;
   gVirtualX->GetGCValues(gct, gval);
   return GetGC(&gval, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the best matching graphics context depending on values.
/// If rw is false only a readonly, not modifiable graphics context
/// is returned. If rw is true a new modifiable graphics context is
/// returned.

TGGC *TGGCPool::GetGC(GCValues_t *values, Bool_t rw)
{
   TGGC *gc, *best_match = 0;
   Int_t matching_bits, best_matching_bits = -1;
   Bool_t exact = kFALSE;

   if (!values)
      rw = kTRUE;

   if (!rw) {

      // First, try to find an exact matching GC.
      // If no one found, then use the closest one.

      TIter next(fList);

      while ((gc = (TGGC *) next())) {
         matching_bits = MatchGC(gc, values);
         if (matching_bits > best_matching_bits) {
            best_matching_bits = matching_bits;
            best_match = gc;
            if ((gc->fValues.fMask & values->fMask) == values->fMask) {
               exact = kTRUE;
               break;
            }
         }
      }

      if (best_match) {
         if (gDebug > 0)
            Printf("<TGGCPool::GetGC>: %smatching GC found\n", exact ? "exact " : "");
         best_match->AddReference();
         if (!exact) {
            // add missing values to the best_match'ing GC...
            UpdateGC(best_match, values);
         }
         return best_match;
      }
   }

   gc = new TGGC(values, kTRUE);

   fList->Add(gc);

   return gc;
}

////////////////////////////////////////////////////////////////////////////////
/// Try to find matching graphics context. On success returns the amount
/// of matching bits (which may be zero if masks have no common bits),
/// -1 on failure (when there are common bits but not a single match).

Int_t TGGCPool::MatchGC(const TGGC *gc, GCValues_t *values)
{
   Mask_t bit, common_bits;
   Int_t  matching_bits = -1;
   Bool_t match = kFALSE;
   const GCValues_t *gcv = &gc->fValues;

   common_bits = values->fMask & gcv->fMask;

   if (common_bits == 0) return 0;  // no common bits, a possible
                                    // candidate anyway.

   // Careful, check first the tile and stipple mask bits, as these
   // influence nearly all other GC functions... (do the same for
   // some other such bits as GCFunction, etc...). Perhaps we should
   // allow only exact GC matches.

   if (gcv->fMask & kGCTile)
      if ((gcv->fTile != kNone) && !(values->fMask & kGCTile)) return -1;
   if (values->fMask & kGCTile)
      if ((values->fTile != kNone) && !(gcv->fMask & kGCTile)) return -1;
   if (gcv->fMask & kGCStipple)
      if ((gcv->fStipple != kNone) && !(values->fMask & kGCStipple)) return -1;
   if (values->fMask & kGCStipple)
      if ((values->fStipple != kNone) && !(gcv->fMask & kGCStipple)) return -1;

   for (bit = 1; bit <= common_bits; bit <<= 1) {
      switch (bit & common_bits) {
         default:
         case 0:
            continue;
         case kGCFunction:
            match = (values->fFunction == gcv->fFunction);
            break;
         case kGCPlaneMask:
            match = (values->fPlaneMask == gcv->fPlaneMask);
            break;
         case kGCForeground:
            match = (values->fForeground == gcv->fForeground);
            break;
         case kGCBackground:
            match = (values->fBackground == gcv->fBackground);
            break;
         case kGCLineWidth:
            match = (values->fLineWidth == gcv->fLineWidth);
            break;
         case kGCLineStyle:
            match = (values->fLineStyle == gcv->fLineStyle);
            break;
         case kGCCapStyle:
            match = (values->fCapStyle == gcv->fCapStyle);
            break;
         case kGCJoinStyle:
            match = (values->fJoinStyle == gcv->fJoinStyle);
            break;
         case kGCFillStyle:
            match = (values->fFillStyle == gcv->fFillStyle);
            break;
         case kGCFillRule:
            match = (values->fFillRule == gcv->fFillRule);
            break;
         case kGCTile:
            match = (values->fTile == gcv->fTile);
            break;
         case kGCStipple:
            match = (values->fStipple == gcv->fStipple);
            break;
         case kGCTileStipXOrigin:
            match = (values->fTsXOrigin == gcv->fTsXOrigin);
            break;
         case kGCTileStipYOrigin:
            match = (values->fTsYOrigin == gcv->fTsYOrigin);
            break;
         case kGCFont:
            match = (values->fFont == gcv->fFont);
            break;
         case kGCSubwindowMode:
            match = (values->fSubwindowMode == gcv->fSubwindowMode);
            break;
         case kGCGraphicsExposures:
            match = (values->fGraphicsExposures == gcv->fGraphicsExposures);
            break;
         case kGCClipXOrigin:
            match = (values->fClipXOrigin == gcv->fClipXOrigin);
            break;
         case kGCClipYOrigin:
            match = (values->fClipYOrigin == gcv->fClipYOrigin);
            break;
         case kGCClipMask:
            match = (values->fClipMask == gcv->fClipMask);
            break;
         case kGCDashOffset:
            match = (values->fDashOffset == gcv->fDashOffset);
            break;
         case kGCDashList:
            if (values->fDashLen == gcv->fDashLen)
               match = (strncmp(values->fDashes, gcv->fDashes, gcv->fDashLen) == 0);
            break;
         case kGCArcMode:
            match = (values->fArcMode == gcv->fArcMode);
            break;
      }
      if (!match)
         return -1;
      matching_bits++;
      match = kFALSE;
   }

   return matching_bits;
}

////////////////////////////////////////////////////////////////////////////////
/// Update graphics context with the values spcified in values->fMask.

void TGGCPool::UpdateGC(TGGC *gc, GCValues_t *values)
{
   gc->SetAttributes(values);
}

////////////////////////////////////////////////////////////////////////////////
/// List all graphics contexts in the pool.

void TGGCPool::Print(Option_t *) const
{
   fList->Print();
}
