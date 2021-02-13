// @(#)root/gui:$Id$
// Author: Fons Rademakers   05/01/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
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
// TGIcon                                                               //
//                                                                      //
// This class handles GUI icons.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGIcon.h"
#include "TGPicture.h"
#include "TSystem.h"
#include "TImage.h"
#include "TVirtualDragManager.h"
#include "TVirtualX.h"

#include <iostream>


ClassImp(TGIcon);

////////////////////////////////////////////////////////////////////////////////
/// Create icon.

TGIcon::TGIcon(const TGWindow *p, const char *image) : TGFrame(p, 1, 1)
{
   fPic = nullptr;

   if (!image)
      image = "bld_rgb.xpm";

   char *path = StrDup(image);

   fPath = gSystem->GetDirName(path);

   fImage = TImage::Open(path);
   if (fImage) {
      fPic = fClient->GetPicturePool()->GetPicture(gSystem->BaseName(path),
                                   fImage->GetPixmap(), fImage->GetMask());
      TGFrame::Resize(fImage->GetWidth(), fImage->GetHeight());
   }
   SetWindowName();
   delete [] path;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete icon and free picture.

TGIcon::~TGIcon()
{
   if (fPic) fClient->FreePicture(fPic);
}

////////////////////////////////////////////////////////////////////////////////
/// Set icon picture.

void TGIcon::SetPicture(const TGPicture *pic)
{
   fPic = pic;
   gVirtualX->ClearWindow(fId);
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set icon image.

void TGIcon::SetImage(const char *img)
{
   //delete fImage;
   TImage *i = TImage::Open(img);
   fPath = gSystem->GetDirName(img);

   SetImage(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Change icon image.

void TGIcon::SetImage(TImage *img)
{
   if (!img) {
      return;
   }

   delete fImage; //  !! mem.leak!!
   fImage = img;

   Resize(fImage->GetWidth(), fImage->GetHeight());
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return size of icon.

TGDimension TGIcon::GetDefaultSize() const
{
   return TGDimension((fPic) ? fPic->GetWidth()  : fWidth,
                      (fPic) ? fPic->GetHeight() : fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw picture.

void TGIcon::DoRedraw()
{
   Bool_t border = (GetOptions() & kRaisedFrame) ||
                   (GetOptions() & kSunkenFrame) ||
                   (GetOptions() & kDoubleBorder);

   if (fPic) fPic->Draw(fId, GetBckgndGC()(), border, border);
   if (border)  DrawBorder();
}

////////////////////////////////////////////////////////////////////////////////
/// Resize.

void TGIcon::Resize(UInt_t w, UInt_t h)
{
   TGFrame::Resize(w, h);

   // allow scaled resize for icons with TImage
   if (!fImage) {
      return;
   }

   gVirtualX->ClearWindow(fId);

   if (fPic) {
      fClient->FreePicture(fPic);
   }
   Bool_t border = (GetOptions() & kRaisedFrame) ||
                   (GetOptions() & kSunkenFrame) ||
                   (GetOptions() & kDoubleBorder);

   fImage->Scale(w - 2*border, h - 2*border);
   fPic = fClient->GetPicturePool()->GetPicture(fImage->GetName(),
                                                fImage->GetPixmap(), fImage->GetMask());
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Move icon to (x,y) and resize it to (w,h).

void TGIcon::MoveResize(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   Move(x, y);
   Resize(w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset icon to original image. It can be used only via context menu.

void TGIcon::Reset()
{
   if (!fImage || !fClient->IsEditable()) return;

   TString name = fImage->GetName();
   name.Chop();
   char *path = gSystem->ConcatFileName(fPath.Data(), name.Data());
   SetImage(path);

   delete [] path;
}

////////////////////////////////////////////////////////////////////////////////
/// Set directory where image is located

void TGIcon::SetImagePath(const char *path)
{
   if (!path) {
      return;
   }

   fPath = gSystem->UnixPathName(path);
   gSystem->ExpandPathName(fPath);
}

////////////////////////////////////////////////////////////////////////////////
/// Save an icon widget as a C++ statement(s) on output stream out.

void TGIcon::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   if (!fPic) {
      Error("SavePrimitive()", "icon pixmap not found ");
      return;
   }

   TString picname = gSystem->UnixPathName(fPic->GetName());
   gSystem->ExpandPathName(picname);

   out <<"   TGIcon *";
   if (!fImage) {
      out << GetName() << " = new TGIcon(" << fParent->GetName()
         << ",gClient->GetPicture(" << quote
         << picname   // if no path
         << quote << ")" << "," << GetWidth() << "," << GetHeight();
      if (fBackground == GetDefaultFrameBackground()) {
         if (!GetOptions()) {
            out <<");" << std::endl;
         } else {
            out << "," << GetOptionString() <<");" << std::endl;
         }
      } else {
         out << "," << GetOptionString() << ",ucolor);" << std::endl;
      }
   } else {
      TString name = fPath;
      name += "/";
      name += fImage->GetName();
      name.Chop();
      out << GetName() << " = new TGIcon(" << fParent->GetName()  << ","
          << quote << name.Data() << quote << ");" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
}
