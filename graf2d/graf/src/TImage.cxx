// @(#)root/graf:$Id$
// Author: Fons Rademakers   15/10/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TImage.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSystem.h"

ClassImp(TImage);

/** \class TImage
\ingroup BasicGraphics

An abstract interface to image processing library.

It allows for the reading and writing of images in different formats, several
image manipulations (scaling, tiling, merging, etc.) and displaying in pads.

The concrete implementation of this class is done by the TASImage class. The
methods are documented in that class.
*/

////////////////////////////////////////////////////////////////////////////////
/// Create an image.
/// Use ReadImage() or SetImage() to initialize the image.

TImage *TImage::Create()
{
   static TPluginHandler *h = 0;

   if (!h) {
      h = gROOT->GetPluginManager()->FindHandler("TImage");
      if (!h) return 0;
      if (h->LoadPlugin() == -1) {
         h = 0;   // try to reload plugin next time
         return 0;
      }
   }
   TImage *img = (TImage *) h->ExecPlugin(0);
   if (img) img->SetName("dummy_image");

   return img;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the image type for the extension specified in filename.
/// Case of the extension is ignored. E.g. for a filename "myimg.GIF",
/// kGif is returned.
/// kAnimGif is returned if the file extension is ".anim.gif".

TImage::EImageFileTypes TImage::GetImageFileTypeFromFilename(const char* filename)
{
   if (!filename) return kUnknown;

   TString sFilename(filename);
   if (sFilename.EndsWith(".xpm.gz", TString::kIgnoreCase))
      return kGZCompressedXpm;
   else if (sFilename.EndsWith(".xpm.z", TString::kIgnoreCase))
      return kZCompressedXpm;
   else if (sFilename.EndsWith(".png", TString::kIgnoreCase))
      return kPng;
   else if (sFilename.EndsWith(".jpeg", TString::kIgnoreCase))
      return kJpeg;
   else if (sFilename.EndsWith(".jpg", TString::kIgnoreCase))
      return kJpeg;
   else if (sFilename.EndsWith(".xcf", TString::kIgnoreCase))
      return kXcf;
   else if (sFilename.EndsWith(".ppm", TString::kIgnoreCase))
      return kPpm;
   else if (sFilename.EndsWith(".pnm", TString::kIgnoreCase))
      return kPnm;
   else if (sFilename.EndsWith(".bmp", TString::kIgnoreCase))
      return kBmp;
   else if (sFilename.EndsWith(".ico", TString::kIgnoreCase))
      return kIco;
   else if (sFilename.EndsWith(".cur", TString::kIgnoreCase))
      return kCur;
   else if (sFilename.EndsWith(".gif", TString::kIgnoreCase))
      return kGif;
   else if (sFilename.EndsWith(".tiff", TString::kIgnoreCase))
      return kTiff;
   else if (sFilename.EndsWith(".tif", TString::kIgnoreCase))
      return kTiff;
   else if (sFilename.EndsWith(".xbm", TString::kIgnoreCase))
      return kXbm;
   else if (sFilename.EndsWith(".fits", TString::kIgnoreCase))
      return kFits;
   else if (sFilename.EndsWith(".tga", TString::kIgnoreCase))
      return kTga;
   else if (sFilename.EndsWith(".xml", TString::kIgnoreCase))
      return kXml;
   else if (sFilename.EndsWith(".anim.gif", TString::kIgnoreCase))
      return kAnimGif;

   return kUnknown;
}

////////////////////////////////////////////////////////////////////////////////
/// List this image with its attributes.

void TImage::ls(Option_t *) const
{
   TROOT::IndentLevel();
   printf("TImage: \"%s\"\n", GetName() );
}

////////////////////////////////////////////////////////////////////////////////
/// Open a specified image file.

TImage *TImage::Open(const char *file, EImageFileTypes type)
{
   TImage *img = Create();
   char *fullname = gSystem->ExpandPathName(file);

   if (img)
      img->ReadImage(fullname, type);

   delete [] fullname;

   return img;
}

////////////////////////////////////////////////////////////////////////////////
/// Open an image with the specified data in a Double_t array.

TImage *TImage::Open(const char *name, const Double_t *imageData, UInt_t width,
                     UInt_t height, TImagePalette *palette)
{
   TImage *img = Create();

   if (img) {
      img->SetImage(imageData, width, height, palette);
      img->SetName(name);
   }
   return img;
}

////////////////////////////////////////////////////////////////////////////////
/// Open an image with the specified data in a TArrayD.

TImage *TImage::Open(const char *name, const TArrayD &imageData, UInt_t width,
                     TImagePalette *palette)
{
   TImage *img = Create();

   if (img) {
      img->SetImage(imageData, width, palette);
      img->SetName(name);
   }
   return img;
}

////////////////////////////////////////////////////////////////////////////////
/// Open an image with the specified data in a TVectorD.

TImage *TImage::Open(const char *name, const TVectorD &imageData, UInt_t width,
                     TImagePalette *palette)
{
   TImage *img = Create();

   if (img) {
      img->SetImage(imageData, width, palette);
      img->SetName(name);
   }
   return img;
}

////////////////////////////////////////////////////////////////////////////////
/// Create image from XPM data array.

TImage *TImage::Open(char **data)
{
   TImage *img = Create();

   if (img) {
      img->SetImageBuffer(data, TImage::kXpm);
      img->SetName("XPM_image");
   }
   return img;
}


TImage operator+(const TImage &i1, const TImage &i2) { TImage ret(i1); ret.Append(&i2, "+"); return ret; }
TImage operator/(const TImage &i1, const TImage &i2) { TImage ret(i1); ret.Append(&i2, "/"); return ret; }
