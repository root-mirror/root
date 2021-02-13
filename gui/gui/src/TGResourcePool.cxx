// @(#)root/gui:$Id$
// Author: Fons Rademakers   19/5/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGResourcePool                                                       //
//                                                                      //
// This class implements a pool for the default GUI resource set,       //
// like GC's, colors, fonts, etc..                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TGResourcePool.h"
#include "TGWindow.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TVirtualX.h"
#include "TGGC.h"
#include "TGPicture.h"
#include "TGFont.h"
#include "TGMimeTypes.h"


const int gray_width  = 8;
const int gray_height = 8;
static unsigned char gray_bits[] = {
   0x55, 0xaa, 0x55, 0xaa,
   0x55, 0xaa, 0x55, 0xaa
};


ClassImp(TGResourcePool);

////////////////////////////////////////////////////////////////////////////////
/// Create the global GUI resource pool manager.

TGResourcePool::TGResourcePool(TGClient *client)
{
   fClient = client;
   fDefaultColormap = gVirtualX->GetColormap();

   // Get GUI defaults from [system].rootrc
   TString default_font   = gEnv->GetValue("Gui.DefaultFont",  "-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
   TString menu_font      = gEnv->GetValue("Gui.MenuFont",     "-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
   TString menu_hi_font   = gEnv->GetValue("Gui.MenuHiFont",   "-*-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
   TString doc_fixed_font = gEnv->GetValue("Gui.DocFixedFont", "-*-courier-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
   TString doc_prop_font  = gEnv->GetValue("Gui.DocPropFont",  "-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
   TString icon_font      = gEnv->GetValue("Gui.IconFont",     "-*-helvetica-medium-r-*-*-10-*-*-*-*-*-iso8859-1");
   TString status_font    = gEnv->GetValue("Gui.StatusFont",   "-*-helvetica-medium-r-*-*-10-*-*-*-*-*-iso8859-1");
   TString backcolor      = gEnv->GetValue("Gui.BackgroundColor",         "#e0e0e0");
   TString forecolor      = gEnv->GetValue("Gui.ForegroundColor",         "black");
   TString highlightcolor = gEnv->GetValue("Gui.HighLightColor",          "#d0d0d0");
   TString selbackcolor   = gEnv->GetValue("Gui.SelectBackgroundColor",   "#86abd9");
   TString selforecolor   = gEnv->GetValue("Gui.SelectForegroundColor",   "white");
   TString docbackcolor   = gEnv->GetValue("Gui.DocumentBackgroundColor", "white");
   TString docforecolor   = gEnv->GetValue("Gui.DocumentForegroundColor", "black");
   TString tipbackcolor   = gEnv->GetValue("Gui.TooltipBackgroundColor",  "LightYellow");
   TString tipforecolor   = gEnv->GetValue("Gui.TooltipForegroundColor",  "black");
   TString framebgpixmap  = gEnv->GetValue("Gui.FrameBackgroundPixmap", "");
   TString docbgpixmap    = gEnv->GetValue("Gui.DocumentBackgroundPixmap", "");

   TString icon_path = gEnv->GetValue("Gui.IconPath", "");
   if (icon_path.IsNull()) {
      icon_path = "icons";
      gSystem->PrependPathName(gSystem->HomeDirectory(), icon_path);
#ifndef R__WIN32
      icon_path = ".:" + icon_path + ":" + TROOT::GetIconPath() + ":" + EXTRAICONPATH;
#else
      icon_path = ".;" + icon_path + ";" + TROOT::GetIconPath() + ";" + EXTRAICONPATH;
#endif
   }

   TString mime_file = ".root.mimes";
   gSystem->PrependPathName(gSystem->HomeDirectory(), mime_file);
   mime_file = gEnv->GetValue("Gui.MimeTypeFile", mime_file.Data());
   gSystem->ExpandPathName(mime_file);

   if (gSystem->AccessPathName(mime_file, kReadPermission)) {
      mime_file = "root.mimes";
      gSystem->PrependPathName(TROOT::GetEtcDir(), mime_file);
   }

   // Setup colors...
   fClient->GetColorByName("white", fWhite);  // white and black always exist
   fClient->GetColorByName("black", fBlack);
   if (!fClient->GetColorByName(backcolor, fBackColor))
      fBackColor = fWhite;
   if (!fClient->GetColorByName(forecolor, fForeColor))
      fForeColor = fBlack;
   if (!fClient->GetColorByName(highlightcolor, fHighLightColor))
      fHighLightColor = fWhite;
   fHilite = fClient->GetHilite(fBackColor);
   fShadow = fClient->GetShadow(fBackColor);
   if (!fClient->GetColorByName(selbackcolor, fSelBackColor))
      fSelBackColor = fBlack;
   if (!fClient->GetColorByName(selforecolor, fSelForeColor))
      fSelForeColor = fWhite;
   if (!fClient->GetColorByName(docbackcolor, fDocBackColor))
      fDocBackColor = fWhite;
   if (!fClient->GetColorByName(docforecolor, fDocForeColor))
      fDocForeColor = fBlack;
   if (!fClient->GetColorByName(tipbackcolor, fTipBackColor))
      fTipBackColor = fWhite;
   if (!fClient->GetColorByName(tipforecolor, fTipForeColor))
      fTipForeColor = fBlack;

   // Setup checkered pix/bit-maps...
   fCheckered = gVirtualX->CreatePixmap(fClient->GetDefaultRoot()->GetId(),
                                        (const char *)gray_bits,
                                        gray_width, gray_height,
                                        fBackColor, fWhite /*fHilite*/,
                                        gVirtualX->GetDepth());
   fCheckeredBitmap = gVirtualX->CreatePixmap(fClient->GetDefaultRoot()->GetId(),
                                              (const char *)gray_bits,
                                              gray_width, gray_height,
                                              1, 0, 1);

   // Create picture pool, GC pool, font pool, mime type list, etc.

   // Create picture pool and pre-load some pictures...
   fPicturePool = new TGPicturePool(fClient, icon_path);

   fDefaultBackPicture    = 0;
   fDefaultDocBackPicture = 0;
   if (!framebgpixmap.IsNull())
      fDefaultBackPicture = fPicturePool->GetPicture(framebgpixmap);
   if (!docbgpixmap.IsNull())
      fDefaultDocBackPicture = fPicturePool->GetPicture(docbgpixmap);

   // Create font pool and pre-load some fonts...
   fFontPool = new TGFontPool(fClient);
   fDefaultFont  = fFontPool->GetFont(default_font);
   fMenuFont     = fFontPool->GetFont(menu_font);
   fMenuHiFont   = fFontPool->GetFont(menu_hi_font);
   fDocFixedFont = fFontPool->GetFont(doc_fixed_font);
   fDocPropFont  = fFontPool->GetFont(doc_prop_font);
   fIconFont     = fFontPool->GetFont(icon_font);
   fStatusFont   = fFontPool->GetFont(status_font);

   // Create GC pool and pre-load some GCs...
   fGCPool = new TGGCPool(fClient);

   GCValues_t gval;
   gval.fMask = kGCForeground | kGCBackground | kGCFont |
                kGCFillStyle  | kGCGraphicsExposures;
   gval.fGraphicsExposures = kFALSE;
   gval.fFillStyle  = kFillSolid;
   gval.fFont       = fDefaultFont ? fDefaultFont->GetFontHandle() : FontH_t(0);
   gval.fBackground = fBackColor;
   gval.fForeground = fBlack;
   fBlackGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fWhite;
   fWhiteGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fHilite;
   fHiliteGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fShadow;
   fShadowGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fBackColor;
   fBckgndGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fForeColor;
   fFrameGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fSelBackColor;
   fSelbgndGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fSelForeColor;
   gval.fBackground = fSelBackColor;
   fSelGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fFont       = fDocPropFont ? fDocPropFont->GetFontHandle() : FontH_t(0);
   gval.fForeground = fDocForeColor;
   gval.fBackground = fDocBackColor;
   fDocGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fForeground = fDocBackColor;
   fDocbgndGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fFont       = fStatusFont ? fStatusFont->GetFontHandle() : FontH_t(0);
   gval.fForeground = fTipForeColor;
   gval.fBackground = fTipBackColor;
   fTipGC = fGCPool->GetGC(&gval, kTRUE);

   gval.fMask = kGCForeground | kGCBackground | kGCFillStyle |
                kGCLineWidth  | kGCSubwindowMode |
                kGCGraphicsExposures | kGCStipple;
   gval.fForeground = fBlack;
   gval.fBackground = fSelBackColor;
   gval.fLineWidth  = 0;
   gval.fFillStyle  = kFillStippled;
   gval.fSubwindowMode = kIncludeInferiors;
   gval.fGraphicsExposures = kFALSE;
   gval.fStipple = fCheckeredBitmap;
   fFocusGC = fGCPool->GetGC(&gval, kTRUE);

   // Setup cursors...
   fDefaultCursor = kNone;
   fGrabCursor    = gVirtualX->CreateCursor(kArrowRight);
   fTextCursor    = gVirtualX->CreateCursor(kCaret);
   fWaitCursor    = gVirtualX->CreateCursor(kWatch);

   // Read in mime type...
   fMimeTypeList  = new TGMimeTypes(fClient, mime_file);

   // Clipboard handle...
#ifndef R__WIN32
   fClipboardAtom = gVirtualX->InternAtom("_ROOT_CLIPBOARD", kFALSE);
#else
   fClipboardAtom = gVirtualX->InternAtom("CLIPBOARD", kFALSE);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup the resource pool...

TGResourcePool::~TGResourcePool()
{
   delete fPicturePool;
   delete fFontPool;
   delete fGCPool;
   delete fMimeTypeList;

   gVirtualX->DeletePixmap(fCheckered);
   gVirtualX->DeletePixmap(fCheckeredBitmap);
}
