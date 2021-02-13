// @(#)root/asimage:$Id: TASImage.cxx,v 1.54 2006/03/13 15:18:56 rdm E
// Author: Fons Rademakers, Reiner Rohlfs, Valeriy Onuchin   28/11/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************
 * Some parts of this source are based on libAfterImage 2.00.00
 *          (http://www.afterstep.org/)
 *
 * Copyright (c) 2002 Sasha Vasko <sasha@aftercode.net>
 * Copyright (c) 1998, 1999 Ethan Fischer <allanon@crystaltokyo.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 **************************************************************************/

/** \class TASImage
\ingroup asimage

Image class.

TASImage is the concrete interface to the image processing library
libAfterImage.

It allows reading and writing of images in different formats, several image
manipulations (scaling, tiling, merging, etc.) and displaying in pads. The size
of the image on the screen does not depend on the original size of the image but
on the size of the pad. Therefore it is very easy to resize the image on the
screen by resizing the pad.

Besides reading an image from a file an image can be defined by a two
dimensional array of values. A palette defines the color of each value.

The image can be zoomed by defining a rectangle with the mouse. The color
palette can be modified with a GUI, just select StartPaletteEditor() from the
context menu.

Several examples showing how to use this class are available in the
ROOT tutorials: `$ROOTSYS/tutorials/image/`
*/

#  include <ft2build.h>
#  include FT_FREETYPE_H
#  include FT_GLYPH_H
#include "TASImage.h"
#include "TASImagePlugin.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVirtualX.h"
#include "TVirtualPad.h"
#include "TArrayD.h"
#include "TVectorD.h"
#include "TVirtualPS.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TArrayL.h"
#include "TPoint.h"
#include "TFrame.h"
#include "TTF.h"
#include "TRandom.h"
#include <iostream>
#include "THashTable.h"
#include "TPluginManager.h"
#include "TEnv.h"
#include "TStyle.h"
#include "TText.h"
#include "RConfigure.h"
#include "TVirtualPadPainter.h"
#include "snprintf.h"

#ifndef WIN32
#ifndef R__HAS_COCOA
#   include <X11/Xlib.h>
#endif
#else
#   include "Windows4root.h"
#endif
#ifndef WIN32
#ifdef R__HAS_COCOA
#   define X_DISPLAY_MISSING 1
#endif
#   include <afterbase.h>
#else
#   include <win32/config.h>
#   include <win32/afterbase.h>
#   define X_DISPLAY_MISSING 1
#endif
#   include <afterimage.h>
#   include <bmp.h>
extern "C" {
#   include <draw.h>
}

// auxiliary functions for general polygon filling
#include "TASPolyUtils.c"


ASVisual *TASImage::fgVisual = 0;
Bool_t TASImage::fgInit = kFALSE;

static ASFontManager *gFontManager = 0;
static unsigned long kAllPlanes = ~0;
THashTable *TASImage::fgPlugList = new THashTable(50);

// default icon paths
static char *gIconPaths[7] = {0, 0, 0, 0, 0, 0, 0};

// To scale fonts to the same size as the old TT version
const Float_t kScale = 0.985;

///////////////////////////// alpha-blending macros ///////////////////////////////

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3)) && !__INTEL_COMPILER
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#ifdef R__BYTESWAP
typedef struct {
   unsigned char b;
   unsigned char g;
   unsigned char r;
   unsigned char a;
} __argb32__;
#else
typedef struct {
   unsigned char a;
   unsigned char r;
   unsigned char g;
   unsigned char b;
} __argb32__;
#endif


//______________________________________________________________________________
#define _alphaBlend(bot, top) {\
   __argb32__ *T = (__argb32__*)(top);\
   __argb32__ *B = (__argb32__*)(bot);\
   int aa = 255-T->a;\
   if (!aa) {\
      *bot = *top;\
   } else { \
      B->a = ((B->a*aa)>>8) + T->a;\
      B->r = (B->r*aa + T->r*T->a)>>8;\
      B->g = (B->g*aa + T->g*T->a)>>8;\
      B->b = (B->b*aa + T->b*T->a)>>8;\
   }\
}\


ClassImp(TASImage);
ClassImp(TASImagePlugin);

////////////////////////////////////////////////////////////////////////////////
/// Destroy image.

void TASImage::DestroyImage()
{
   if (fImage) {
      destroy_asimage(&fImage);
   }

   if (fIsGray && fGrayImage) {
      destroy_asimage(&fGrayImage);
   }

   fIsGray     = kFALSE;
   fGrayImage  = 0;
   fImage      = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default parameters.

void TASImage::SetDefaults()
{
   fImage         = 0;
   fScaledImage   = 0;
   fMaxValue      = 1;
   fMinValue      = 0;
   fEditable      = kFALSE;
   fPaintMode     = 1;
   fZoomOffX      = 0;
   fZoomOffY      = 0;
   fZoomWidth     = 0;
   fZoomHeight    = 0;
   fZoomUpdate    = kZoomOps;

   fGrayImage     = 0;
   fIsGray        = kFALSE;
   fPaletteEnabled = kFALSE;

   if (!fgInit) {
      set_application_name((char*)(gProgName ? gProgName : "ROOT"));
      fgInit = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default image constructor.

TASImage::TASImage()
{
   SetDefaults();
}

////////////////////////////////////////////////////////////////////////////////
/// Create an empty image.

TASImage::TASImage(UInt_t w, UInt_t h) : TImage(w, h)
{
   SetDefaults();
   fImage = create_asimage(w ? w : 20, h ? h : 20, 0);
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image object and read from specified file.
/// For more information see description of function ReadImage()
/// which is called by this constructor.

TASImage::TASImage(const char *file, EImageFileTypes) : TImage(file)
{
   SetDefaults();
   TString fname = file;
   gSystem->ExpandPathName(fname);
   ReadImage(fname.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image depending on the values of imageData.
/// For more information see function SetImage() which is called
/// by this constructor.

TASImage::TASImage(const char *name, const Double_t *imageData, UInt_t width,
                   UInt_t height, TImagePalette *palette) : TImage(name)
{
   SetDefaults();
   SetImage(imageData, width, height, palette);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image depending on the values of imageData.
/// The size of the image is width X (imageData.fN / width).
/// For more information see function SetImage() which is called by
/// this constructor.

TASImage::TASImage(const char *name, const TArrayD &imageData, UInt_t width,
                   TImagePalette *palette) : TImage(name)
{
   SetDefaults();
   SetImage(imageData, width, palette);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image depending on the values of imageData.
/// The size of the image is width X (imageData.fN / width).
/// For more information see function SetImage() which is called by
/// this constructor.

TASImage::TASImage(const char *name, const TVectorD &imageData, UInt_t width,
                   TImagePalette *palette) : TImage(name)
{
   SetDefaults();
   SetImage(imageData, width, palette);
}

////////////////////////////////////////////////////////////////////////////////
/// Image copy constructor.

TASImage::TASImage(const TASImage &img) : TImage(img)
{
   SetDefaults();

   if (img.IsValid()) {
      fImage = clone_asimage(img.fImage, SCL_DO_ALL);
      fScaledImage   = fScaledImage ? (TASImage*)img.fScaledImage->Clone("") : 0;
      fGrayImage     = fGrayImage ? clone_asimage(img.fGrayImage, SCL_DO_ALL) : 0;

      if (img.fImage->alt.vector) {
         Int_t size = img.fImage->width * img.fImage->height * sizeof(double);
         fImage->alt.vector = (double*)malloc(size);
         memcpy(fImage->alt.vector, img.fImage->alt.vector, size);
      }

      fZoomUpdate = kNoZoom;
      fZoomOffX   = img.fZoomOffX;
      fZoomOffY   = img.fZoomOffY;
      fZoomWidth  = img.fZoomWidth;
      fZoomHeight = img.fZoomHeight;
      fEditable   = img.fEditable;
      fIsGray     = img.fIsGray;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Image assignment operator.

TASImage &TASImage::operator=(const TASImage &img)
{
   if (this != &img && img.IsValid()) {
      TImage::operator=(img);

      DestroyImage();
      delete fScaledImage;
      fImage = clone_asimage(img.fImage, SCL_DO_ALL);
      fScaledImage = fScaledImage ? (TASImage*)img.fScaledImage->Clone("") : 0;
      fGrayImage = fGrayImage ? clone_asimage(img.fGrayImage, SCL_DO_ALL) : 0;

      if (img.fImage->alt.vector) {
         Int_t size = img.fImage->width * img.fImage->height * sizeof(double);
         fImage->alt.vector = (double*)malloc(size);
         memcpy(fImage->alt.vector, img.fImage->alt.vector, size);
      }

      fScaledImage = img.fScaledImage ? (TASImage*)img.fScaledImage->Clone("") : 0;
      fZoomUpdate = kNoZoom;
      fZoomOffX   = img.fZoomOffX;
      fZoomOffY   = img.fZoomOffY;
      fZoomWidth  = img.fZoomWidth;
      fZoomHeight = img.fZoomHeight;
      fEditable   = img.fEditable;
      fIsGray     = img.fIsGray;
      fPaintMode  = 1;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Image destructor, clean up image and visual.

TASImage::~TASImage()
{
   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set icons paths.

static void init_icon_paths()
{
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

   Int_t cnt = 0;
   Ssiz_t from = 0;
   TString token;
#ifndef R__WIN32
   const char *delim = ":";
#else
   const char *delim = ";";
#endif
   while (icon_path.Tokenize(token, from, delim) && cnt < 6) {
      char *path = gSystem->ExpandPathName(token.Data());
      if (path) {
         gIconPaths[cnt] = path;
         cnt++;
      }
   }
   gIconPaths[cnt] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Guess the file type from the first byte of file.

const char *TASImage::TypeFromMagicNumber(const char *file)
{
   UChar_t magic;
   FILE *fp = fopen(file, "rb");
   const char *ret = "";

   if (!fp) return 0;

   if (!fread(&magic, 1, 1, fp)) {
      fclose(fp);
      return 0;
   }

   switch (magic) {
      case 0x00:
      {
         if (!fread(&magic, 1, 1, fp)) {
            fclose(fp);
            return 0;
         }
         if (!fread(&magic, 1, 1, fp)) {
            fclose(fp);
            return 0;
         }

         ret = (magic == 1) ? "ico" : "cur";
         break;
      }
      case 0x25:
      {
         if (!fread(&magic, 1, 1, fp)) {
            fclose(fp);
            return 0;
         }

         if (magic == 0x21) ret = "ps";
         else if (magic == 0x50) ret = "pdf";
         break;
      }
      case 0x42:
         ret = "bmp";
         break;
      case 0x47:
         ret = "gif";
         break;
      case 0x54:
         ret = "tga";
         break;
      case 0x49:
         ret = "tiff";
         break;
      case 0x89:
         ret = "png";
         break;
      case 0xff:
         ret = "jpg";
         break;
      default:
         ret = "";
   }

   fclose(fp);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Read specified image file.
/// The file type is determined by the file extension (the type argument is
/// ignored). It will attempt to append .gz and then .Z to the filename and
/// find such a file. If the filename ends with extension consisting of digits
/// only, it will attempt to find the file with this extension stripped
/// off. On success this extension will be used to load subimage from
/// the file with that number. Subimage is supported for GIF files
/// (ICO, BMP, CUR, TIFF, XCF to be supported in future).
///  For example,
/// ~~~ {.cpp}
///    i1 = TImage::Open("anim.gif.0"); // read the first subimage
///    i4 = TImage::Open("anim.gif.3"); // read the forth subimage
/// ~~~
/// It is also possible to put XPM raw string (see also SetImageBuffer) as
/// the first input parameter ("filename"), such string  is returned by
/// GetImageBuffer method.

void TASImage::ReadImage(const char *filename, EImageFileTypes /*type*/)
{
   if (!InitVisual()) {
      Warning("Scale", "Visual not initiated");
      return;
   }

   Bool_t xpm = filename && (filename[0] == '/' &&
                filename[1] == '*') && filename[2] == ' ';

   if (xpm) {  // XPM strings in-memory array
      SetImageBuffer((char**)&filename, TImage::kXpm);
      fName = "XPM_image";
      return;
   }

   if (!gIconPaths[0]) {
      init_icon_paths();
   }
   // suppress the "root : looking for image ..." messages
   set_output_threshold(0);

   static ASImageImportParams iparams;
   iparams.flags = 0;
   iparams.width = 0;
   iparams.height = 0;
   iparams.filter = SCL_DO_ALL;
   iparams.gamma = SCREEN_GAMMA;
   iparams.gamma_table = NULL;
   iparams.compression = GetImageCompression();
   iparams.format = ASA_ASImage;
   iparams.search_path = gIconPaths;
   iparams.subimage = 0;
   iparams.return_animation_delay = -1;

   TString ext;
   const char *dot;
   if (filename) dot = strrchr(filename, '.');
   else          dot = 0;
   ASImage *image = 0;
   TString fname = filename;

   if (!dot) {
      if (filename) ext = TypeFromMagicNumber(filename);
      else ext = dot + 1;
   } else {
      ext = dot + 1;
   }

   if (!ext.IsNull() && ext.IsDigit()) { // read subimage
      iparams.subimage = ext.Atoi();
      fname = fname(0, fname.Length() - ext.Length() - 1);
      ext = strrchr(fname.Data(), '.') + 1;
   }

   image = file2ASImage_extra(fname.Data(), &iparams);

   if (image) { // it's OK
      goto end;
   } else {  // try to read it via plugin
      if (ext.IsNull()) {
         return;
      }
      ext.ToLower();
      ext.Strip();
      UInt_t w = 0;
      UInt_t h = 0;
      unsigned char *bitmap = 0;

      TImagePlugin *plug = (TImagePlugin*)fgPlugList->FindObject(ext.Data());

      if (!plug) {
         TPluginHandler *handler = gROOT->GetPluginManager()->FindHandler("TImagePlugin", ext);
         if (!handler || ((handler->LoadPlugin() == -1))) {
            return;
         }
         plug = (TImagePlugin*)handler->ExecPlugin(1, ext.Data());

         if (!plug) {
            return;
         }

         fgPlugList->Add(plug);
      }

      if (plug) {
         if (plug->InheritsFrom(TASImagePlugin::Class())) {
            image = ((TASImagePlugin*)plug)->File2ASImage(fname.Data());
            if (image) goto end;
         }
         bitmap = plug->ReadFile(fname.Data(), w, h);
         if (bitmap) {
            image = bitmap2asimage(bitmap, w, h, 0, 0);
         }
         if (!image) {
            return;
         }
      }
   }

end:
   fName.Form("%s.", gSystem->BaseName(fname.Data()));

   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   fImage      = image;
   fZoomUpdate = kNoZoom;
   fEditable   = kFALSE;
   fZoomOffX   = 0;
   fZoomOffY   = 0;
   fZoomWidth  = fImage->width;
   fZoomHeight = fImage->height;
   fPaintMode  = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Write image to specified file.
///
/// If there is no file extension or if the file extension is unknown, the
/// type argument will be used to determine the file type. The quality and
/// compression is derived from the TAttImage values.
///
/// It's possible to write image into an animated GIF file by specifying file
/// name as "myfile.gif+" or "myfile.gif+NN", where NN is the delay of displaying
/// subimages during animation in 10ms seconds units. NN is not restricted
/// to two digits. If NN is omitted the delay between subimages is zero.
/// For an animation that stops after last subimage is reached, one has to
/// write the last image as .gif+ (zero delay of last image) or .gif+NN
/// (NN*10ms delay of last image).
///
/// For repeated animation (looping), the last subimage must be specified as:
///  - "myfile.gif++NN++" if you want an infinite looping gif with NN*10ms
///    delay of the last image.
///  - "myfile.gif++" for an infinite loop with zero delay of last image.
///  - "myfile.gif+NN++RR" if you want a finite looping gif with NN*10ms
///    delay of the last image and the animation to be stopped after RR
///    repeats. RR is not restricted to two digits.
///
/// A deprecated version for saving the last subimage of a looping gif animation is:
///  - "myfile.gif++NN" for a finite loop where NN is number of repetitions
///    and NN*10ms the delay of last image. (No separate control of repeats and delay).
///    Note: If the file "myfile.gif" already exists, the new frames are appended at
///    the end of the file. To avoid this, delete it first with gSystem->Unlink(myfile.gif);
///
/// The following macro creates animated gif from jpeg images with names
///  - imageNN.jpg, where 1<= NN <= 10
///  - The delays are set to 10*10ms.
/// ~~~ {.cpp}
/// {
///    TImage *img = 0;
///    gSystem->Unlink("anim.gif");  // delete existing file
///
///    for (int i = 1; i <= 10; i++) {
///       delete img; // delete previous image
///
///       // Read image data. Image can be in any format, e.g. png, gif, etc.
///       img = TImage::Open(Form("image%d.jpg", i));
///
///       if (i < 10) {
///          img->WriteImage("anim.gif+10"); // 10 centiseconds delay
///       } else { // the last image written.  "++" stands for infinit animation.
///          img->WriteImage("anim.gif++10++"); // 10 centiseconds delay of last image
///       }
///    }
/// }
/// ~~~

void TASImage::WriteImage(const char *file, EImageFileTypes type)
{
   if (!IsValid()) {
      Error("WriteImage", "no image loaded");
      return;
   }

   if (!file || !*file) {
      Error("WriteImage", "no file name specified");
      return;
   }

   const char *s;
   if ((s = strrchr(file, '.'))) {
      s++;
      EImageFileTypes t = GetFileType(s);
      if (t == kUnknown) {
         Error("WriteImage", "cannot determine a valid file type");
         return;
      }
      if (t != kUnknown)
         type = t;
   }

   if (type == kUnknown) {
      Error("WriteImage", "not a valid file type was specified");
      return;
   }

   UInt_t mytype;
   MapFileTypes(type, mytype);
   ASImageFileTypes atype = (ASImageFileTypes)mytype;

   UInt_t aquality;
   EImageQuality quality = GetImageQuality();
   MapQuality(quality, aquality);

   static TString fname;
   fname = file;
   static ASImageExportParams parms;
   ASImage *im = fScaledImage ? fScaledImage->fImage : fImage;

   switch (type) {
   case kXpm:
      parms.xpm.type = atype;
      parms.xpm.flags = EXPORT_ALPHA;
      parms.xpm.dither = 4;
      parms.xpm.opaque_threshold = 127;
      parms.xpm.max_colors = 512;
      break;
   case kBmp:
      ASImage2bmp(im, fname.Data(), 0);
      return;
   case kXcf:
      ASImage2xcf(im, fname.Data(), 0);
      return;
   case kPng:
      parms.png.type = atype;
      parms.png.flags = EXPORT_ALPHA;
      parms.png.compression = !GetImageCompression() ? -1 : int(GetImageCompression());
      break;
   case kJpeg:
      parms.jpeg.type = atype;
      parms.jpeg.flags = 0;
      parms.jpeg.quality = aquality;
      break;
   case kGif:
      parms.gif.type = atype;
      parms.gif.flags = EXPORT_ALPHA;
      parms.gif.dither = 0;
      parms.gif.opaque_threshold = 0;
      break;
   case kAnimGif:
   {
      parms.gif.type = atype;
      parms.gif.flags = EXPORT_ALPHA | EXPORT_APPEND;
      parms.gif.dither = 0;
      parms.gif.opaque_threshold = 0;
      parms.gif.animate_repeats = 0;

      s += 4; // skip "gif+"
      int delay = 0;

      const TString sufix = s; // we denote as suffix as everything that is after .gif+
      const UInt_t sLength = sufix.Length();

      if (sufix=="+") {
         // .gif++ implies that this is the last image of the animation
         // and that the gif will loop forever (infinite animation)
         // and that the delay of last image is 0ms (backward compatibility reasons)
         delay = 0;
         parms.gif.flags |= EXPORT_ANIMATION_REPEATS;// activate repetition
         parms.gif.animate_repeats = 0;// 0 is code for looping forever (if EXPORT_ANIMATION_REPEATS is also set)
      } else if(sufix=="") {
         // .gif+ implies that this is a subimage of the animation with zero delay
         // or the last image of an animation that will not repeat.
         // Last image delay is zero because atoi("")=0.
         delay = atoi(s);
         //Nothing else needed here
      } else if(!sufix.Contains("+")) {
         // .gif+NN implies that this is a subimage of the animation
         // with NN*10ms delay (latency) until the next one.
         // You can also use this option on the last image if you do not want the gif to replay
         delay = atoi(s);
         //Nothing else needed here
      } else if(sLength>1 && sufix.BeginsWith("+") && sufix.CountChar('+')==1) {
         // .gif++NN implies that this is the last image of the animation
         // and that it will loop NN number of times (finite animation)
         // and that the delay of last image is NN*10ms (backward compatibility reasons).
         delay = atoi(s);// atoi is smart enough to ignore the "+" sign before.
         parms.gif.flags |= EXPORT_ANIMATION_REPEATS;// activate repetition
         parms.gif.animate_repeats = atoi(s);// loops only NN times, then it stops. atoi discards + sign.
      } else if(sLength>3 && sufix.BeginsWith("+") && sufix.EndsWith("++") && !TString(sufix(1,sLength-3)).Contains("+")) {
         // .gif++NN++ implies that this is the last image of the animation
         // and that the gif will loop forever (infinite animation)
         // and that the delay of last image is NN*10ms.
         // In contrast, .gif++ is an infinite loop but with 0 delay, whereas the option
         // .gif++NN is a loop repeated NN times (not infinite) with NN*10ms delay
         // between last and first loop images.
         delay = atoi(s);// atoi discards the three plus signs
         parms.gif.flags |= EXPORT_ANIMATION_REPEATS;// activate repetition
         parms.gif.animate_repeats = 0;// 0 is code for looping forever (if EXPORT_ANIMATION_REPEATS is also set)
      } else if(sLength>3 && sufix.CountChar('+')==2 && TString(sufix(1,sLength-2)).Contains("++")) {
         // .gif+NN++RR implies that this is the last image animation
         // and that the gif will loop RR number of times (finite animation)
         // and that the delay of last image is NN*10ms.
         const TString sDelay = sufix(0,sufix.First('+'));
         const TString sRepeats = sufix(sufix.First('+')+2,sLength-(sufix.First('+')+2));
         delay = atoi(sDelay);
         parms.gif.flags |= EXPORT_ANIMATION_REPEATS;// activate repetition
         parms.gif.animate_repeats = atoi(sRepeats);// loops NN times.
      } else {
         Error("WriteImage", "gif suffix %s not yet supported", s);
         return;
      }

      parms.gif.animate_delay = delay;

      int i1 = fname.Index("gif+");
      if (i1 != kNPOS) {
         fname = fname(0, i1 + 3);
      }
      else {
         Error("WriteImage", "unexpected gif extension structure %s", fname.Data());
         return;
      }
      break;
   }
   case kTiff:
      parms.tiff.type = atype;
      parms.tiff.flags = EXPORT_ALPHA;
      parms.tiff.rows_per_strip = 0;
      parms.tiff.compression_type = aquality <= 50 ? TIFF_COMPRESSION_JPEG :
                                                     TIFF_COMPRESSION_NONE;
      parms.tiff.jpeg_quality = 100;
      parms.tiff.opaque_threshold = 0;
      break;
   default:
      Error("WriteImage", "file type %s not yet supported", s);
      return;
   }

   if (!ASImage2file(im, 0, fname.Data(), atype, &parms)) {
      Error("WriteImage", "error writing file %s", file);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return file type depending on specified extension.
/// Protected method.

TImage::EImageFileTypes TASImage::GetFileType(const char *ext)
{
   TString s(ext);
   s.Strip();
   s.ToLower();

   if (s == "xpm")
      return kXpm;
   if (s == "png")
      return kPng;
   if (s == "jpg" || s == "jpeg")
      return kJpeg;
   if (s == "xcf")
      return kXcf;
   if (s == "ppm")
      return kPpm;
   if (s == "pnm")
      return kPnm;
   if (s == "bmp")
      return kBmp;
   if (s == "ico")
      return kIco;
   if (s == "cur")
      return kCur;
   if (s == "gif")
      return kGif;
   if (s.Contains("gif+"))
      return kAnimGif;
   if (s == "tiff")
      return kTiff;
   if (s == "xbm")
      return kXbm;
   if (s == "tga")
      return kTga;
   if (s == "xml")
      return kXml;

   return kUnknown;
}

////////////////////////////////////////////////////////////////////////////////
/// Map file type to/from AfterImage types.
/// Protected method.

void TASImage::MapFileTypes(EImageFileTypes &type, UInt_t &astype, Bool_t toas)
{
   if (toas) {
      switch (type) {
         case kXpm:
            astype = ASIT_Xpm; break;
         case kZCompressedXpm:
            astype = ASIT_ZCompressedXpm; break;
         case kGZCompressedXpm:
            astype = ASIT_GZCompressedXpm; break;
         case kPng:
            astype = ASIT_Png; break;
         case kJpeg:
            astype = ASIT_Jpeg; break;
         case kXcf:
            astype = ASIT_Xcf; break;
         case kPpm:
            astype = ASIT_Ppm; break;
         case kPnm:
            astype = ASIT_Pnm; break;
         case kBmp:
            astype = ASIT_Bmp; break;
         case kIco:
            astype = ASIT_Ico; break;
         case kCur:
            astype = ASIT_Cur; break;
         case kGif:
            astype = ASIT_Gif; break;
         case kAnimGif:
            astype = ASIT_Gif; break;
         case kTiff:
            astype = ASIT_Tiff; break;
         case kXbm:
            astype = ASIT_Xbm; break;
         case kTga:
            astype = ASIT_Targa; break;
         case kXml:
            astype = ASIT_XMLScript; break;
         default:
            astype = ASIT_Unknown;
      }
   } else {
      switch (astype) {
         case ASIT_Xpm:
            type = kXpm; break;
         case ASIT_ZCompressedXpm:
            type = kZCompressedXpm; break;
         case ASIT_GZCompressedXpm:
            type = kGZCompressedXpm; break;
         case ASIT_Png:
            type = kPng; break;
         case ASIT_Jpeg:
            type = kJpeg; break;
         case ASIT_Xcf:
            type = kXcf; break;
         case ASIT_Ppm:
            type = kPpm; break;
         case ASIT_Pnm:
            type = kPnm; break;
         case ASIT_Bmp:
            type = kBmp; break;
         case ASIT_Ico:
            type = kIco; break;
         case ASIT_Cur:
            type = kCur; break;
         case ASIT_Gif:
            type = kGif; break;
         case ASIT_Tiff:
            type = kTiff; break;
         case ASIT_Xbm:
            type = kXbm; break;
         case ASIT_XMLScript:
            type = kXml; break;
         case ASIT_Targa:
            type = kTga; break;
         default:
            type = kUnknown;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map quality to/from AfterImage quality.
/// Protected method.

void TASImage::MapQuality(EImageQuality &quality, UInt_t &asquality, Bool_t toas)
{
   if (toas) {
      switch (quality) {
         case kImgPoor:
            asquality = 25; break;
         case kImgFast:
            asquality = 75; break;
         case kImgGood:
            asquality = 50; break;
         case kImgBest:
            asquality = 100; break;
         default:
            asquality = 0;
      }
   } else {
      quality = kImgDefault;
      if (asquality > 0  && asquality <= 25)
         quality = kImgPoor;
      if (asquality > 26 && asquality <= 50)
         quality = kImgFast;
      if (asquality > 51 && asquality <= 75)
         quality = kImgGood;
      if (asquality > 76 && asquality <= 100)
         quality = kImgBest;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the old image and creates a new image depending on the values
/// of imageData. The size of the image is width X height.
///
/// The color of each pixel depends on the imageData of the corresponding
/// pixel. The palette is used to convert an image value into its color.
/// If palette is not defined (palette = 0) a default palette is used.
/// Any previously defined zooming is reset.

void TASImage::SetImage(const Double_t *imageData, UInt_t width, UInt_t height,
                        TImagePalette *palette)
{
   TAttImage::SetPalette(palette);

   if (!InitVisual()) {
      Warning("SetImage", "Visual not initiated");
      return;
   }

   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   // get min and max value of image
   fMinValue = fMaxValue = *imageData;
   for (Int_t pixel = 1; pixel < Int_t(width * height); pixel++) {
      if (fMinValue > *(imageData + pixel)) fMinValue = *(imageData + pixel);
      if (fMaxValue < *(imageData + pixel)) fMaxValue = *(imageData + pixel);
   }

   // copy ROOT palette to asImage palette
   const TImagePalette &pal = GetPalette();

   ASVectorPalette asPalette;

   asPalette.npoints = pal.fNumPoints;
   Int_t col;
   for (col = 0; col < 4; col++)
      asPalette.channels[col] = new UShort_t[asPalette.npoints];

   memcpy(asPalette.channels[0], pal.fColorBlue,  pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[1], pal.fColorGreen, pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[2], pal.fColorRed,   pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[3], pal.fColorAlpha, pal.fNumPoints * sizeof(UShort_t));

   asPalette.points = new Double_t[asPalette.npoints];
   for (Int_t point = 0; point < Int_t(asPalette.npoints); point++)
      asPalette.points[point] = fMinValue + (fMaxValue - fMinValue) * pal.fPoints[point];

   fImage = create_asimage_from_vector(fgVisual, (Double_t*)imageData, width,
                                       height, &asPalette, ASA_ASImage,
                                       GetImageCompression(), GetImageQuality());

   delete [] asPalette.points;
   for (col = 0; col < 4; col++)
      delete [] asPalette.channels[col];

   fZoomUpdate = 0;
   fZoomOffX   = 0;
   fZoomOffY   = 0;
   fZoomWidth  = width;
   fZoomHeight = height;
   fPaletteEnabled = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the old image and creates a new image depending on the values
/// of imageData. The size of the image is width X (imageData.fN / width).
/// The color of each pixel depends on the imageData of the corresponding
/// pixel. The palette is used to convert an image value into its color.
/// If palette is not defined (palette = 0) a default palette is used.
/// Any previously defined zooming is reset.

void TASImage::SetImage(const TArrayD &imageData, UInt_t width, TImagePalette *palette)
{
   SetImage(imageData.GetArray(), width, imageData.GetSize() / width, palette);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the old image and creates a new image depending on the values
/// of imageData. The size of the image is width X (imageData.fN / width).
/// The color of each pixel depends on the imageData of the corresponding
/// pixel. The palette is used to convert an image value into its color.
/// If palette is not defined (palette = 0) a default palette is used.
/// Any previously defined zooming is reset.

void TASImage::SetImage(const TVectorD &imageData, UInt_t width, TImagePalette *palette)
{
   SetImage(imageData.GetMatrixArray(), width,
            imageData.GetNoElements() / width, palette);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image from the given pad, afterwards this image can be
/// saved in any of the supported image formats.

void TASImage::FromPad(TVirtualPad *pad, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   if (!pad) {
      Error("FromPad", "pad cannot be 0");
      return;
   }

   if (!InitVisual()) {
      Warning("FromPad", "Visual not initiated");
      return;
   }

   SetName(pad->GetName());

   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   if (gROOT->IsBatch()) { // in batch mode
      TVirtualPS *psave = gVirtualPS;
      gVirtualPS = (TVirtualPS*)gROOT->ProcessLineFast("new TImageDump()");
      gVirtualPS->Open(pad->GetName(), 114); // in memory
      gVirtualPS->SetBit(BIT(11)); //kPrintingPS

      TASImage *itmp = (TASImage*)gVirtualPS->GetStream();

      if (itmp && itmp->fImage) {
         itmp->BeginPaint();
      }

      TVirtualPad *sav = gPad;
      gPad = pad;
      pad->Paint();
      gPad = sav;

      if (itmp && itmp->fImage && (itmp != this)) {
         fImage = clone_asimage(itmp->fImage, SCL_DO_ALL);
         if (itmp->fImage->alt.argb32) {
            UInt_t sz = itmp->fImage->width*itmp->fImage->height;
            fImage->alt.argb32 = (ARGB32*)safemalloc(sz*sizeof(ARGB32));
            memcpy(fImage->alt.argb32, itmp->fImage->alt.argb32, sz*4);
         }
      }
      delete gVirtualPS;
      gVirtualPS = psave;
      return;
   }

   // X11 Synchronization
   gVirtualX->Update(1);
   if (!gThreadXAR) {
      gSystem->Sleep(100);
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
      gSystem->ProcessEvents();
   }

   TVirtualPad *canvas = (TVirtualPad*)pad->GetCanvas();
   Int_t wid = (pad == canvas) ? pad->GetCanvasID() : pad->GetPixmapID();
   gVirtualX->SelectWindow(wid);

   Window_t wd = (Window_t)gVirtualX->GetCurrentWindow();
   if (!wd) return;

   if (w == 0) w = TMath::Abs(pad->UtoPixel(1.));
   if (h == 0) h = pad->VtoPixel(0.);

   static int x11 = -1;
   if (x11 < 0) x11 = gVirtualX->InheritsFrom("TGX11");

   if (x11) { //use built-in optimized version
      fImage = pixmap2asimage(fgVisual, wd, x, y, w, h, kAllPlanes, 0, 0);
   } else {
      unsigned char *bits = gVirtualX->GetColorBits(wd, 0, 0, w, h);

      if (!bits) { // error
         return;
      }
      fImage = bitmap2asimage(bits, w, h, 0, 0);
      delete [] bits;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw image.
/// Support the following drawing options:
///  - "T[x,y[,tint]]" : tile image (use specified offset and tint),
///                      e.g. "T100,100,#556655"
///                      with this option the zooming is not possible
///                      and disabled
///  - "N"             : display in new canvas (of original image size)
///  - "X"             : image is drawn expanded to pad size
///  - "Z"             : image is vectorized and image palette is drawn
///
/// The default is to display the image in the current gPad.

void TASImage::Draw(Option_t *option)
{
   if (!fImage) {
      Error("Draw", "no image set");
      return;
   }

   TString opt = option;
   opt.ToLower();
   if (opt.Contains("n") || !gPad || !gPad->IsEditable()) {
      Int_t w = -64;
      Int_t h = 64;
      w = (fImage->width > 64) ? (Int_t)fImage->width : w;
      h = (fImage->height > 64) ? (Int_t)fImage->height : h;

      Float_t cx = 1./gStyle->GetScreenFactor();
      w = Int_t(w*cx) + 4;
      h = Int_t(h*cx) + 28;
      TString rname = GetName();
      rname.ReplaceAll(".", "");
      rname += Form("\", \"%s (%d x %d)", rname.Data(), fImage->width, fImage->height);
      rname = "new TCanvas(\"" + rname + Form("\", %d, %d);", w, h);
      gROOT->ProcessLineFast(rname.Data());
   }

   if (!opt.Contains("x")) {
      Double_t left = gPad->GetLeftMargin();
      Double_t right = gPad->GetRightMargin();
      Double_t top = gPad->GetTopMargin();
      Double_t bottom = gPad->GetBottomMargin();

      gPad->Range(-left / (1.0 - left - right),
                  -bottom / (1.0 - top - bottom),
                  1 + right / (1.0 - left - right),
                  1 + top / ( 1.0 - top - bottom));
      gPad->RangeAxis(0, 0, 1, 1);
   }

   TFrame *frame = gPad->GetFrame();
   if (frame) {
      frame->SetBorderMode(0);
      frame->SetFillColor(gPad->GetFillColor());
      frame->SetLineColor(gPad->GetFillColor());
      frame->Draw();
   }

   TObject::Draw(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw asimage on drawable.

void TASImage::Image2Drawable(ASImage *im, Drawable_t wid, Int_t x, Int_t y,
                              Int_t xsrc, Int_t ysrc, UInt_t wsrc, UInt_t hsrc,
                              Option_t *opt)
{
   if (!im) return;

   wsrc = wsrc ? wsrc : im->width;
   hsrc = hsrc ? hsrc : im->height;

   static int x11 = -1;
   if (x11 < 0) x11 = gVirtualX->InheritsFrom("TGX11");

   Pixmap_t mask = kNone;

   if (x11) {
      UInt_t hh = hsrc;
      UInt_t ow = wsrc%8;
      UInt_t ww = wsrc - ow + (ow ? 8 : 0);

      UInt_t bit = 0;
      int i = 0;
      UInt_t yy = 0;
      UInt_t xx = 0;

      char *bits = new char[ww*hh]; //an array of bits

      ASImageDecoder *imdec = start_image_decoding(fgVisual, im, SCL_DO_ALPHA,
                                                   xsrc, ysrc, ww, 0, 0);
      if (imdec) {
         for (yy = 0; yy < hh; yy++) {
            imdec->decode_image_scanline(imdec);
            CARD32 *a = imdec->buffer.alpha;

            for (xx = 0; xx < ww; xx++) {
               if (a[xx]) {
                  SETBIT(bits[i], bit);
               } else {
                  CLRBIT(bits[i], bit);
               }
               bit++;
               if (bit == 8) {
                  bit = 0;
                  i++;
               }
            }
         }
      }

      stop_image_decoding(&imdec);

      mask = gVirtualX->CreateBitmap(gVirtualX->GetDefaultRootWindow(),
                                          (const char *)bits, ww, hh);
      delete [] bits;
   }

   GCValues_t gv;
   static GContext_t gc = 0;

   gv.fMask = kGCClipMask | kGCClipXOrigin | kGCClipYOrigin;
   gv.fClipMask = mask;
   gv.fClipXOrigin = x;
   gv.fClipYOrigin = y;

   if (!gc) {
      gc = gVirtualX->CreateGC(gVirtualX->GetDefaultRootWindow(), &gv);
   } else {
      gVirtualX->ChangeGC(gc, &gv);
   }

   if (x11 && (!gPad || gPad->GetGLDevice() == -1)) { //use built-in optimized version
      asimage2drawable(fgVisual, wid, im, (GC)gc, xsrc, ysrc, x, y, wsrc, hsrc, 1);
   } else {
      ASImage *img = 0;
      unsigned char *bits = (unsigned char *)im->alt.argb32;
      if (!bits) {
         img = tile_asimage(fgVisual, im, xsrc, ysrc, wsrc, hsrc,
                            0, ASA_ARGB32, 0, ASIMAGE_QUALITY_DEFAULT);
         if (img)
            bits = (unsigned char *)img->alt.argb32;
      }

      if (bits) {
         TString option(opt);
         option.ToLower();

         if (gPad && gPad->GetGLDevice() != -1) {
            if (TVirtualPadPainter *painter = gPad->GetPainter())
               painter->DrawPixels(bits, wsrc, hsrc, x, y, !option.Contains("opaque"));
         } else {
            Pixmap_t pic = gVirtualX->CreatePixmapFromData(bits, wsrc, hsrc);
            if (pic) {
               if (!option.Contains("opaque")) {
                  SETBIT(wsrc,31);
                  SETBIT(hsrc,31);
               }
               gVirtualX->CopyArea(pic, wid, gc, 0, 0, wsrc, hsrc, x, y);
               gVirtualX->DeletePixmap(pic);
            }
         }
      }

      if (img) {
         destroy_asimage(&img);
      }
   }

   // free mask pixmap
   if (gv.fClipMask != kNone) gVirtualX->DeletePixmap(gv.fClipMask);

   gv.fMask = kGCClipMask;
   gv.fClipMask = kNone;
   if (gc) gVirtualX->ChangeGC(gc, &gv);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw image on the drawable wid (pixmap, window) at x,y position.
///
/// \param[in] wid        : Drawable (pixmap or window) on which image is drawn.
/// \param[in] x,y        : Window coordinates where image is drawn.
/// \param[in] xsrc, ysrc : X and Y coordinates of an image area to be drawn.
/// \param[in] wsrc, hsrc : Width and height image area to be drawn.

void TASImage::PaintImage(Drawable_t wid, Int_t x, Int_t y, Int_t xsrc, Int_t ysrc,
                          UInt_t wsrc, UInt_t hsrc, Option_t *opt)
{
   Image2Drawable(fScaledImage ? fScaledImage->fImage : fImage, wid, x, y,
                  xsrc, ysrc, wsrc, hsrc, opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint image.
/// Support the following drawing options:
///  - "T[x,y[,tint]]" : tile image (use specified offset and tint),
///                      e.g. "T100,100,#556655"
///                      with this option the zooming is not possible
///                      and disabled
///  - "N"             : display in new canvas (of original image size)
///  - "X"             : image is drawn expanded to pad size
///  - "Z"             : image is vectorized and image palette is drawn
///
/// The default is to display the image in the current gPad.

void TASImage::Paint(Option_t *option)
{
   if (!fImage) {
      Error("Paint", "no image set");
      return;
   }

   if (!InitVisual()) {
      Warning("Paint", "Visual not initiated");
      return;
   }

   Int_t   tile_x = 0, tile_y = 0;
   CARD32  tile_tint = 0;
   Bool_t  tile = kFALSE;
   Bool_t  expand = kFALSE;

   TString opt = option;
   opt.ToLower();

   if (opt.Contains("t")) {
      char stint[64];
      if (sscanf(opt.Data() + opt.Index("t"), "t%d,%d,%s", &tile_x, &tile_y,
                 stint) <= 3) {
         tile = kTRUE;
         if (parse_argb_color(stint, (CARD32*)&tile_tint) == stint)
            tile_tint = 0;
      } else {
         Error("Paint", "tile option error");
      }
   } else if (opt.Contains("x")) {
      expand = kTRUE;
      fConstRatio = kFALSE;
   } else if (opt.Contains("z")) {
      fPaletteEnabled = kTRUE;

      if (!fImage->alt.vector) {
         Vectorize(256);
      }
   }

   ASImage *image = fImage;

   // Get geometry of pad
   Int_t to_w = gPad->UtoPixel(1.);
   Int_t to_h = gPad->VtoPixel(0.);

   // remove the size by the margin of the pad
   if (!expand) {
      to_h  = (Int_t)(to_h * (1.0 - gPad->GetBottomMargin() - gPad->GetTopMargin() ) + 0.5);
      to_w  = (Int_t)(to_w * (1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin() ) + 0.5);
   }

   if ((to_w < 25 || to_h < 25) && !expand) {
      Error("Paint", "pad too small to display an image");
      return;
   }

   if (GetConstRatio()) {
      if ((Double_t)to_w / (Double_t)fZoomWidth <
          (Double_t)to_h / (Double_t)fZoomHeight)
         to_h = Int_t(Double_t(fZoomHeight) * to_w / fZoomWidth);
      else
         to_w = Int_t(Double_t(fZoomWidth) * to_h / fZoomHeight);
   }
   // upper left corner and size of the palette in pixels
   Int_t pal_Ax = to_w + gPad->UtoAbsPixel(gPad->GetLeftMargin()) +
                 (gPad->UtoAbsPixel(gPad->GetRightMargin()) / 10);
   Int_t pal_Ay = gPad->YtoAbsPixel(1.0);
   Int_t pal_x = to_w + gPad->UtoPixel(gPad->GetLeftMargin()) +
                 (gPad->UtoPixel(gPad->GetRightMargin()) / 10);
   Int_t pal_y = gPad->YtoPixel(1.0);
   Int_t pal_w = gPad->UtoPixel(gPad->GetRightMargin()) / 3;
   Int_t pal_h = to_h;

   ASImage  *grad_im = 0;

   if (fImage->alt.vector && fPaletteEnabled) {
      // draw the palette
      ASGradient grad;
      const TImagePalette &pal = GetPalette();

      grad.npoints = pal.fNumPoints;
      grad.type    = GRADIENT_Top2Bottom;
      grad.color   = new ARGB32[grad.npoints];
      grad.offset  = new double[grad.npoints];

      for (Int_t pt = 0; pt < grad.npoints; pt++) {
         Int_t oldPt = grad.npoints - pt -1;
         grad.offset[pt] = 1 - pal.fPoints[oldPt];
         grad.color[pt] = (((ARGB32)(pal.fColorBlue[oldPt]  & 0xff00)) >>  8) |
                          (((ARGB32)(pal.fColorGreen[oldPt] & 0xff00))      ) |
                          (((ARGB32)(pal.fColorRed[oldPt]   & 0xff00)) <<  8) |
                          (((ARGB32)(pal.fColorAlpha[oldPt] & 0xff00)) << 16);
      }

      grad_im = make_gradient(fgVisual, &grad , UInt_t(pal_w),
                              pal_h, SCL_DO_COLOR,
                              ASA_ARGB32, GetImageCompression(), GetImageQuality());

      delete [] grad.color;
      delete [] grad.offset;
   }

   if (tile) {
      delete fScaledImage;
      fScaledImage = (TASImage*)TImage::Create();
      if (!fScaledImage) return;
      fScaledImage->fImage = tile_asimage(fgVisual, fImage, tile_x, tile_y,
                                          to_w, to_h, tile_tint, ASA_ASImage,
                                          GetImageCompression(), GetImageQuality());
      image = fScaledImage->fImage;

   } else if (fZoomUpdate == kZoomOps) {
      image = fImage;

   } else {
      // Scale and zoom image if needed
      if (Int_t(fImage->width) != to_w || Int_t(fImage->height) != to_h ||
          fImage->width != fZoomWidth || fImage->height != fZoomHeight) {

         if (fScaledImage && (Int_t(fScaledImage->GetWidth()) != to_w ||
                Int_t(fScaledImage->GetHeight()) != to_h ||
                fZoomUpdate)) {

            delete fScaledImage;
            fScaledImage = 0;
         }

         if (!fScaledImage) {
            fScaledImage = (TASImage*)TImage::Create();
            if (!fScaledImage) return;

            if (fZoomWidth && fZoomHeight &&
                ((fImage->width != fZoomWidth) || (fImage->height != fZoomHeight))) {
               // zoom and scale image
               ASImage *tmpImage = 0;

               tmpImage = tile_asimage(fgVisual, fImage, fZoomOffX,
                                          fImage->height - fZoomHeight - fZoomOffY,
                                          fZoomWidth, fZoomHeight, 0, ASA_ASImage,
                                          GetImageCompression(), GetImageQuality());

               if (tmpImage) {
                  fScaledImage->fImage = scale_asimage(fgVisual, tmpImage, to_w, to_h,
                                                       ASA_ASImage, GetImageCompression(),
                                                      GetImageQuality());
                  destroy_asimage(&tmpImage);
               }
            } else {
               // scale image, no zooming
               fScaledImage->fImage = scale_asimage(fgVisual, fImage, to_w, to_h,
                                                    ASA_ASImage, GetImageCompression(),
                                                    GetImageQuality());
            }
         }
         image = fScaledImage->fImage;
      }
   }
   fZoomUpdate = 0;

   if (!image) {
      Error("Paint", "image could not be rendered to display");
      return;
   }

   int tox = expand  ? 0 : int(gPad->UtoPixel(1.) * gPad->GetLeftMargin());
   int toy = expand  ? 0 : int(gPad->VtoPixel(0.) * gPad->GetTopMargin());

   if (!gROOT->IsBatch()) {
      Window_t wid = (Window_t)gVirtualX->GetWindowID(gPad->GetPixmapID());
      Image2Drawable(fScaledImage ? fScaledImage->fImage : fImage, wid, tox, toy);

      if (grad_im && fPaletteEnabled) {
         // draw color bar
         Image2Drawable(grad_im, wid, pal_x, pal_y);

         // values of palette
         TGaxis axis;
         Int_t ndiv = 510;
         double min = fMinValue;
         double max = fMaxValue;
         axis.SetLineColor(0);       // draw white ticks
         Double_t pal_Xpos = gPad->AbsPixeltoX(pal_Ax + pal_w);
         axis.PaintAxis(pal_Xpos, gPad->PixeltoY(pal_Ay + pal_h - 1),
                        pal_Xpos, gPad->PixeltoY(pal_Ay),
                        min, max, ndiv, "+LU");
         min = fMinValue;
         max = fMaxValue;
         axis.SetLineColor(1);       // draw black ticks
         axis.PaintAxis(pal_Xpos, gPad->AbsPixeltoY(pal_Ay + pal_h),
                        pal_Xpos, gPad->AbsPixeltoY(pal_Ay + 1),
                        min, max, ndiv, "+L");
      }
   }

   // loop over pixmap and draw image to PostScript
   if (gVirtualPS) {
      if (gVirtualPS->InheritsFrom("TImageDump")) { // PostScript is asimage
         TImage *dump = (TImage *)gVirtualPS->GetStream();
         if (!dump) return;
         dump->Merge(fScaledImage ? fScaledImage : this, "alphablend",
                     gPad->XtoAbsPixel(0), gPad->YtoAbsPixel(1));

         if (grad_im) {
            TASImage tgrad;
            tgrad.fImage = grad_im;
            dump->Merge(&tgrad, "alphablend", pal_Ax, pal_Ay);

            // values of palette
            TGaxis axis;
            Int_t ndiv = 510;
            double min = fMinValue;
            double max = fMaxValue;
            axis.SetLineColor(1);       // draw black ticks
            Double_t pal_Xpos = gPad->AbsPixeltoX(pal_Ax + pal_w);
            axis.PaintAxis(pal_Xpos, gPad->AbsPixeltoY(pal_Ay + pal_h),
                           pal_Xpos, gPad->AbsPixeltoY(pal_Ay + 1),
                           min, max, ndiv, "+L");
         }
         return;
      } else if (gVirtualPS->InheritsFrom("TPDF")) {
         Warning("Paint", "PDF not implemented yet");
         return;
      } else if (gVirtualPS->InheritsFrom("TSVG")) {
         Warning("Paint", "SVG not implemented yet");
         return;
      }

      // get special color cell to be reused during image printing
      TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();
      TColor *color = 0;
      // Look for color by name
      if ((color = (TColor*)colors->FindObject("Image_PS")) == 0)
         color = new TColor(colors->GetEntries(), 1., 1., 1., "Image_PS");

      gVirtualPS->SetFillColor(color->GetNumber());
      gVirtualPS->SetFillStyle(1001);

      Double_t dx = gPad->GetX2()-gPad->GetX1();
      Double_t dy = gPad->GetY2()-gPad->GetY1();
      Double_t x1,x2,y1,y2;

      if (expand) {
         x1 = gPad->GetX1();
         x2 = x1+dx/image->width;
         y1 = gPad->GetY2();
         y2 = y1+dy/image->height;
      } else {
         x1 = gPad->GetX1()+dx*gPad->GetLeftMargin();
         x2 = x1+(dx*(1-gPad->GetRightMargin()-gPad->GetLeftMargin()))/image->width;
         y1 = gPad->GetY2()-dy*gPad->GetTopMargin();
         y2 = y1+(dy*(1-gPad->GetTopMargin()-gPad->GetBottomMargin()))/image->height;
      }

      gVirtualPS->CellArrayBegin(image->width, image->height, x1, x2, y1, y2);

      ASImageDecoder *imdec = start_image_decoding(fgVisual, image, SCL_DO_ALL,
                                                   0, 0, image->width, image->height, 0);
      if (!imdec) return;
      for (Int_t yt = 0; yt < (Int_t)image->height; yt++) {
         imdec->decode_image_scanline(imdec);
         for (Int_t xt = 0; xt < (Int_t)image->width; xt++)
            gVirtualPS->CellArrayFill(imdec->buffer.red[xt],
                                      imdec->buffer.green[xt],
                                      imdec->buffer.blue[xt]);
      }
      stop_image_decoding(&imdec);
      gVirtualPS->CellArrayEnd();

      // print the color bar
      if (grad_im) {
         Double_t xconv = (gPad->AbsPixeltoX(pal_Ax + pal_w) - gPad->AbsPixeltoX(pal_Ax)) / grad_im->width;
         Double_t yconv = (gPad->AbsPixeltoY(pal_Ay - pal_h) - gPad->AbsPixeltoY(pal_Ay)) / grad_im->height;
         x1 = gPad->AbsPixeltoX(pal_Ax);
         x2 = x1 + xconv;
         y2 = gPad->AbsPixeltoY(pal_Ay);
         y1 = y2 - yconv;
         gVirtualPS->CellArrayBegin(grad_im->width, grad_im->height,
                                    x1, x2, y1, y2);

         imdec = start_image_decoding(fgVisual, grad_im, SCL_DO_ALL,
                                      0, 0, grad_im->width, grad_im->height, 0);
         if (imdec) {
            for (Int_t yt = 0; yt < (Int_t)grad_im->height; yt++) {
               imdec->decode_image_scanline(imdec);
               for (Int_t xt = 0; xt < (Int_t)grad_im->width; xt++)
                  gVirtualPS->CellArrayFill(imdec->buffer.red[xt],
                                            imdec->buffer.green[xt],
                                            imdec->buffer.blue[xt]);
            }
         }
         stop_image_decoding(&imdec);
         gVirtualPS->CellArrayEnd();

         // values of palette
         TGaxis axis;
         Int_t ndiv = 510;
         double min = fMinValue;
         double max = fMaxValue;
         axis.SetLineColor(1);       // draw black ticks
         Double_t pal_Xpos = gPad->AbsPixeltoX(pal_Ax + pal_w);
         axis.PaintAxis(pal_Xpos, gPad->AbsPixeltoY(pal_Ay + pal_h),
                        pal_Xpos, gPad->AbsPixeltoY(pal_Ay + 1),
                        min, max, ndiv, "+L");

      }
   }

   if (grad_im) {
      destroy_asimage(&grad_im);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Is the mouse in the image ?

Int_t TASImage::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t pxl, pyl, pxt, pyt;

   Int_t px1 = gPad->XtoAbsPixel(0);
   Int_t py1 = gPad->YtoAbsPixel(0);
   Int_t px2 = gPad->XtoAbsPixel(1);
   Int_t py2 = gPad->YtoAbsPixel(1);

   if (px1 < px2) {pxl = px1; pxt = px2;}
   else           {pxl = px2; pxt = px1;}
   if (py1 < py2) {pyl = py1; pyt = py2;}
   else           {pyl = py2; pyt = py1;}

   if ((px > pxl && px < pxt) && (py > pyl && py < pyt))
      return 0;

   return 999999;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse events.

void TASImage::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   static TBox *ZoomBox;

   if (!gPad) return;

   if (IsEditable()) {
      gPad->ExecuteEvent(event, px, py);
      return;
   }

   gPad->SetCursor(kCross);

   static Int_t px1old, py1old, px2old, py2old;
   static Int_t px1, py1, px2, py2, pxl, pyl, pxt, pyt;

   if (!IsValid()) return;

   if (event == kButton1Motion || event == kButton1Down  ||
       event == kButton1Up) {

      // convert to image pixel on screen
      Int_t imgX = px - gPad->XtoAbsPixel(0);
      Int_t imgY = py - gPad->YtoAbsPixel(1);

      if (imgX < 0)  px = px - imgX;
      if (imgY < 0)  py = py - imgY;

      ASImage *image = fImage;
      if (fScaledImage) image = fScaledImage->fImage;

      if (imgX >= (int)image->width)  px = px - imgX + image->width - 1;
      if (imgY >= (int)image->height) py = py - imgY + image->height - 1;

      switch (event) {

         case kButton1Down:
            px1 = gPad->XtoAbsPixel(gPad->GetX1());
            py1 = gPad->YtoAbsPixel(gPad->GetY1());
            px2 = gPad->XtoAbsPixel(gPad->GetX2());
            py2 = gPad->YtoAbsPixel(gPad->GetY2());
            px1old = px; py1old = py;
            break;

         case kButton1Motion:
            px2old = px;
            px2old = TMath::Max(px2old, px1);
            px2old = TMath::Min(px2old, px2);
            py2old = py;
            py2old = TMath::Max(py2old, py2);
            py2old = TMath::Min(py2old, py1);
            pxl = TMath::Min(px1old, px2old);
            pxt = TMath::Max(px1old, px2old);
            pyl = TMath::Max(py1old, py2old);
            pyt = TMath::Min(py1old, py2old);

            if (ZoomBox) {
               ZoomBox->SetX1(gPad->AbsPixeltoX(pxl));
               ZoomBox->SetY1(gPad->AbsPixeltoY(pyl));
               ZoomBox->SetX2(gPad->AbsPixeltoX(pxt));
               ZoomBox->SetY2(gPad->AbsPixeltoY(pyt));
            }
            else {
               ZoomBox = new TBox(pxl, pyl, pxt, pyt);
               ZoomBox->SetFillStyle(0);
               ZoomBox->Draw("l*");
            }
            gPad->Modified(kTRUE);
            gPad->Update();
            break;

         case kButton1Up:
            // do nothing if zoom area is too small
            if ( TMath::Abs(pxl - pxt) < 5 || TMath::Abs(pyl - pyt) < 5)
               return;

            pxl = 0;
            pxt = 0;
            pyl = 0;
            pyt = 0;

            Double_t xfact = (fScaledImage) ? (Double_t)fScaledImage->fImage->width  / fZoomWidth  : 1;
            Double_t yfact = (fScaledImage) ? (Double_t)fScaledImage->fImage->height / fZoomHeight : 1;

            Int_t imgX1 = px1old - gPad->XtoAbsPixel(0);
            Int_t imgY1 = py1old - gPad->YtoAbsPixel(1);
            Int_t imgX2 = px  - gPad->XtoAbsPixel(0);
            Int_t imgY2 = py  - gPad->YtoAbsPixel(1);

            imgY1 = image->height - 1 - imgY1;
            imgY2 = image->height - 1 - imgY2;
            imgX1 = (Int_t)(imgX1 / xfact) + fZoomOffX;
            imgY1 = (Int_t)(imgY1 / yfact) + fZoomOffY;
            imgX2 = (Int_t)(imgX2 / xfact) + fZoomOffX;
            imgY2 = (Int_t)(imgY2 / yfact) + fZoomOffY;

            Zoom((imgX1 < imgX2) ? imgX1 : imgX2, (imgY1 < imgY2) ? imgY1 : imgY2,
                 TMath::Abs(imgX1 - imgX2) + 1, TMath::Abs(imgY1 - imgY2) + 1);

            if (ZoomBox) {
               ZoomBox->Delete();
               ZoomBox = 0;
            }
            gPad->Modified(kTRUE);
            gPad->Update();
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get image pixel coordinates and the pixel value at the mouse pointer.

char *TASImage::GetObjectInfo(Int_t px, Int_t py) const
{
   static char info[64];
   info[0] = 0;

   if (!IsValid()) return info;

   // convert to image pixel on screen
   px -= gPad->XtoAbsPixel(0);
   py -= gPad->YtoAbsPixel(1);

   // no info if mouse is outside of image
   if (px < 0 || py < 0)  return info;

   ASImage *image = fImage;
   if (fScaledImage) image = fScaledImage->fImage;
   if (px >= (int)image->width || py >= (int)image->height)
      return info;

   py = image->height - 1 - py;
   // convert to original image size and take zooming into account
   if (fScaledImage) {
      px = (Int_t)(px / (Double_t)fScaledImage->fImage->width  * fZoomWidth ) + fZoomOffX;
      py = (Int_t)(py / (Double_t)fScaledImage->fImage->height * fZoomHeight) + fZoomOffY;
   }

   if (fImage->alt.vector) {
      snprintf(info,64,"x: %d  y: %d   %.5g",
              px, py, fImage->alt.vector[px + py * fImage->width]);
   } else {
      snprintf(info,64,"x: %d  y: %d", px, py);
   }

   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a new palette to an image.
/// Only images that were created with the SetImage() functions can be
/// modified with this function. The previously used palette is destroyed.

void TASImage::SetPalette(const TImagePalette *palette)
{
   TAttImage::SetPalette(palette);

   if (!InitVisual()) {
      Warning("SetPalette", "Visual not initiated");
      return;
   }

   if (!IsValid()) {
      Warning("SetPalette", "Image not valid");
      return;
   }

   if (fImage->alt.vector == 0)
      return;

   // copy ROOT palette to asImage palette
   const TImagePalette &pal = GetPalette();

   ASVectorPalette asPalette;
   asPalette.npoints = pal.fNumPoints;
   asPalette.channels[0] = new CARD16 [asPalette.npoints];
   asPalette.channels[1] = new CARD16 [asPalette.npoints];
   asPalette.channels[2] = new CARD16 [asPalette.npoints];
   asPalette.channels[3] = new CARD16 [asPalette.npoints];
   memcpy(asPalette.channels[0], pal.fColorBlue,  pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[1], pal.fColorGreen, pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[2], pal.fColorRed,   pal.fNumPoints * sizeof(UShort_t));
   memcpy(asPalette.channels[3], pal.fColorAlpha, pal.fNumPoints * sizeof(UShort_t));

   asPalette.points = new double[asPalette.npoints];
   for (Int_t point = 0; point < Int_t(asPalette.npoints); point++)
      asPalette.points[point] = fMinValue + (fMaxValue - fMinValue) * pal.fPoints[point];

   // use the new palette in this image
   colorize_asimage_vector(fgVisual, fImage, &asPalette, ASA_ASImage, GetImageQuality());

   delete [] asPalette.points;
   for (Int_t col = 0; col < 4; col++)
      delete [] asPalette.channels[col];


   delete fScaledImage;
   fScaledImage = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Scale the original image.
/// The size of the image on the screen does not change because it is defined
/// by the size of the pad.
/// This function can be used to change the size of an image before writing
/// it into a file. The colors of the new pixels are interpolated.
/// An image created with the SetImage() functions cannot be modified with
/// the function SetPalette() any more after a call of this function!

void TASImage::Scale(UInt_t toWidth, UInt_t toHeight)
{
   if (!IsValid()) {
      Warning("Scale", "Image not initiated");
      return;
   }

   if (!InitVisual()) {
      Warning("Scale", "Visual not initiated");
      return;
   }

   if (toWidth < 1)
       toWidth = 1;
   if (toHeight < 1 )
      toHeight = 1;
   if (toWidth > 30000)
      toWidth = 30000;
   if (toHeight > 30000)
      toHeight = 30000;

   ASImage *img = scale_asimage(fgVisual, fImage, toWidth, toHeight,
                                ASA_ASImage, GetImageCompression(),
                                GetImageQuality());
   DestroyImage();
   fImage = img;
   UnZoom();
   fZoomUpdate = kZoomOps;
}

////////////////////////////////////////////////////////////////////////////////
/// Another method of enlarging images where corners remain unchanged,
/// but middle part gets tiled.

void TASImage::Slice(UInt_t xStart, UInt_t xEnd, UInt_t yStart,  UInt_t yEnd,
                     UInt_t toWidth, UInt_t toHeight)
{
   if (!IsValid()) {
      Warning("Scale", "Image not initiated");
      return;
   }

   if (!InitVisual()) {
      Warning("Scale", "Visual not initiated");
      return;
   }

   if (toWidth < 1)
       toWidth = 1;
   if (toHeight < 1 )
      toHeight = 1;
   if (toWidth > 30000)
      toWidth = 30000;
   if (toHeight > 30000)
      toHeight = 30000;

   ASImage *img = slice_asimage(fgVisual, fImage, xStart, xEnd,
                                yStart, yEnd, toWidth, toHeight,
                                ASA_ASImage, GetImageCompression(),
                                GetImageQuality());

   DestroyImage();
   fImage = img;
   UnZoom();
   fZoomUpdate = kZoomOps;
}

////////////////////////////////////////////////////////////////////////////////
/// Tile the original image.

void TASImage::Tile(UInt_t toWidth, UInt_t toHeight)
{
   if (!IsValid()) {
      Warning("Tile", "Image not initiated");
      return;
   }

   if (!InitVisual()) {
      Warning("Tile", "Visual not initiated");
      return;
   }

   if (toWidth < 1)
       toWidth = 1;
   if (toHeight < 1 )
      toHeight = 1;
   if (toWidth > 30000)
      toWidth = 30000;
   if (toHeight > 30000)
      toHeight = 30000;

   ASImage *img = tile_asimage(fgVisual, fImage, 0, 0, toWidth, toHeight, 0,
                                ASA_ASImage, GetImageCompression(), GetImageQuality());
   DestroyImage();
   fImage = img;
   UnZoom();
   fZoomUpdate = kZoomOps;
}

////////////////////////////////////////////////////////////////////////////////
/// The area of an image displayed in a pad is defined by this function.
/// Note: the size on the screen is defined by the size of the pad.
/// The original image is not modified by this function.
/// If width or height is larger than the original image they are reduced to
/// the width and height of the image.
/// If the off values are too large (off + width > image width) than the off
/// values are decreased. For example: offX = image width - width
/// Note: the parameters are always relative to the original image not to the
/// size of an already zoomed image.

void TASImage::Zoom(UInt_t offX, UInt_t offY, UInt_t width, UInt_t height)
{
   if (!IsValid()) {
      Warning("Zoom", "Image not valid");
      return;
   }
   fZoomUpdate = kZoom;

   fZoomWidth  = (width == 0) ? 1 : ((width > fImage->width) ? fImage->width : width);
   fZoomHeight = (height == 0) ? 1 : ((height > fImage->height) ? fImage->height : height);
   fZoomOffX   = offX;
   if (fZoomOffX + fZoomWidth > fImage->width)
      fZoomOffX = fImage->width - fZoomWidth;
   fZoomOffY   = offY;
   if (fZoomOffY + fZoomHeight > fImage->height)
      fZoomOffY = fImage->height - fZoomHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Un-zoom the image to original size.
/// UnZoom() - performs undo for Zoom,Crop,Scale actions

void TASImage::UnZoom()
{
   if (!IsValid()) {
      Warning("UnZoom", "Image not valid");
      return;
   }
   fZoomUpdate = kZoom;
   fZoomOffX   = 0;
   fZoomOffY   = 0;
   fZoomWidth  = fImage->width;
   fZoomHeight = fImage->height;

   delete fScaledImage;
   fScaledImage = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Flip image in place.
///
/// Flip is either 90, 180, 270, 180 is default.
/// This function manipulates the original image and destroys the
/// scaled and zoomed image which will be recreated at the next call of
/// the Draw function. If the image is zoomed the zoom - coordinates are
/// now relative to the new image.
/// This function cannot be used for images which were created with the
/// SetImage() functions, because the original pixel values would be
/// destroyed.

void TASImage::Flip(Int_t flip)
{
   if (!IsValid()) {
      Warning("Flip", "Image not valid");
      return;
   }
   if (!InitVisual()) {
      Warning("Flip", "Visual not initiated");
      return;
   }

   if (fImage->alt.vector) {
      Warning("Flip", "flip does not work for data images");
      return;
   }

   Int_t rflip = flip/90;

   UInt_t w = fImage->width;
   UInt_t h = fImage->height;

   if (rflip & 1) {
      w = fImage->height;
      h = fImage->width;
   }

   ASImage *img = flip_asimage(fgVisual, fImage, 0, 0, w, h, rflip,
                               ASA_ASImage, GetImageCompression(),
                               GetImageQuality());
   DestroyImage();
   fImage = img;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Mirror image in place.
///
/// If vert is true mirror in vertical axis, horizontal otherwise.
/// Vertical is default.
/// This function manipulates the original image and destroys the
/// scaled and zoomed image which will be recreated at the next call of
/// the Draw function. If the image is zoomed the zoom - coordinates are
/// now relative to the new image.
/// This function cannot be used for images which were created with the
/// SetImage() functions, because the original pixel values would be
/// destroyed.

void TASImage::Mirror(Bool_t vert)
{
   if (!IsValid()) {
      Warning("Mirror", "Image not valid");
      return;
   }

   if (!InitVisual()) {
      Warning("Mirror", "Visual not initiated");
      return;
   }

   if (fImage->alt.vector) {
      Warning("Mirror", "mirror does not work for data images");
      return;
   }

   ASImage *img = mirror_asimage(fgVisual, fImage, 0, 0,
                                 fImage->width, fImage->height, vert,
                                 ASA_ASImage, GetImageCompression(),
                                 GetImageQuality());
   DestroyImage();
   fImage = img;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Return width of original image not of the displayed image.
/// (Number of image pixels)

UInt_t TASImage::GetWidth() const
{
   return fImage ? fImage->width : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return height of original image not of the displayed image.
/// (Number of image pixels)

UInt_t TASImage::GetHeight() const
{
   return fImage ? fImage->height : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return width of the displayed image not of the original image.
/// (Number of screen pixels)

UInt_t TASImage::GetScaledWidth() const
{
   return fScaledImage ? fScaledImage->fImage->width : GetWidth();
}

////////////////////////////////////////////////////////////////////////////////
/// Return height of the displayed image not of the original image.
/// (Number of screen pixels)

UInt_t TASImage::GetScaledHeight() const
{
   return fScaledImage ? fScaledImage->fImage->height : GetHeight();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the zoom parameters.
/// This is useful when the zoom has been done interactively using the mouse.

void TASImage::GetZoomPosition(UInt_t &x, UInt_t &y, UInt_t &w, UInt_t &h) const
{
   x = fZoomOffX;
   y = fZoomOffY;
   w = fZoomWidth;
   h = fZoomHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to initialize the ASVisual.

Bool_t TASImage::InitVisual()
{
   Display *disp;

   Bool_t inbatch = fgVisual && (fgVisual->dpy == (void*)1); // was in batch
   Bool_t noX = gROOT->IsBatch() || gVirtualX->InheritsFrom("TGWin32");

   // was in batch, but switched to gui
   if (inbatch && !noX) {
      destroy_asvisual(fgVisual, kFALSE);
      fgVisual = 0;
   }

   if (fgVisual && fgVisual->dpy) { // already initialized
      return kTRUE;
   }

   // batch or win32 mode
   if (!fgVisual && noX) {
      disp = 0;
      fgVisual = create_asvisual(0, 0, 0, 0);
      fgVisual->dpy = (Display*)1; //fake (not used)
      return kTRUE;
   }

#ifndef WIN32
#ifdef R__HAS_COCOA
   fgVisual = create_asvisual(0, 0, 0, 0);
   fgVisual->dpy = (Display*)1; //fake (not used)
#else
   disp = (Display*) gVirtualX->GetDisplay();
   Int_t screen  = gVirtualX->GetScreen();
   Int_t depth   = gVirtualX->GetDepth();
   Visual *vis   = (Visual*) gVirtualX->GetVisual();
   Colormap cmap = (Colormap) gVirtualX->GetColormap();

   if (vis == 0 || cmap == 0) {
      fgVisual = create_asvisual(0, 0, 0, 0);
   } else {
      fgVisual = create_asvisual_for_id(disp, screen, depth,
                                        XVisualIDFromVisual(vis), cmap, 0);
   }
#endif
#else
   fgVisual = create_asvisual(0, 0, 0, 0);
   fgVisual->dpy = (Display*)1; //fake (not used)
#endif

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Start palette editor.

void TASImage::StartPaletteEditor()
{
   if (!IsValid()) {
      Warning("StartPaletteEditor", "Image not valid");
      return;
   }
   if (fImage->alt.vector == 0) {
      Warning("StartPaletteEditor", "palette can be modified only for data images");
      return;
   }

   // Opens a GUI to edit the color palette
   TAttImage::StartPaletteEditor();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns image pixmap.
/// The pixmap must deleted by user.

Pixmap_t TASImage::GetPixmap()
{
   if (!InitVisual()) {
      Warning("GetPixmap", "Visual not initiated");
      return 0;
   }

   Pixmap_t ret;

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;

   static int x11 = -1;
   if (x11 < 0) x11 = gVirtualX->InheritsFrom("TGX11");

   if (x11) {   // use builtin version
      ret = (Pixmap_t)asimage2pixmap(fgVisual, gVirtualX->GetDefaultRootWindow(),
                                       img, 0, kTRUE);
   } else {
      if (!fImage->alt.argb32) {
         BeginPaint();
      }
      ret = gVirtualX->CreatePixmapFromData((unsigned char*)fImage->alt.argb32,
                                             fImage->width, fImage->height);
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns image mask pixmap (alpha channel).
/// The pixmap must deleted by user.

Pixmap_t TASImage::GetMask()
{
   Pixmap_t pxmap = 0;

   if (!InitVisual()) {
      Warning("GetMask", "Visual not initiated");
      return pxmap;
   }

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;

   if (!img) {
      Warning("GetMask", "No image");
      return pxmap;
   }

   UInt_t hh = img->height;
   UInt_t ow = img->width%8;
   UInt_t ww = img->width - ow + (ow ? 8 : 0);

   UInt_t bit = 0;
   int i = 0;
   UInt_t y = 0;
   UInt_t x = 0;

   char *bits = new char[ww*hh]; //an array of bits

   ASImageDecoder *imdec = start_image_decoding(fgVisual, img, SCL_DO_ALPHA,
                                                0, 0, ww, 0, 0);
   if (!imdec) {
      delete [] bits;
      return 0;
   }

   for (y = 0; y < hh; y++) {
      imdec->decode_image_scanline(imdec);
      CARD32 *a = imdec->buffer.alpha;

      for (x = 0; x < ww; x++) {
         if (a[x]) {
            SETBIT(bits[i], bit);
         } else {
            CLRBIT(bits[i], bit);
         }
         bit++;
         if (bit == 8) {
            bit = 0;
            i++;
         }
      }
   }

   stop_image_decoding(&imdec);
   pxmap = gVirtualX->CreateBitmap(gVirtualX->GetDefaultRootWindow(), (const char *)bits,
                                   ww, hh);
   delete [] bits;
   return pxmap;
}

////////////////////////////////////////////////////////////////////////////////
/// Create image from pixmap.

void TASImage::SetImage(Pixmap_t pxm, Pixmap_t mask)
{
   if (!InitVisual()) {
      Warning("SetImage", "Visual not initiated");
      return;
   }

   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   Int_t xy;
   UInt_t w, h;
   gVirtualX->GetWindowSize(pxm, xy, xy, w, h);

   if (fName.IsNull()) fName.Form("img_%dx%d",w, h);

   static int x11 = -1;
   if (x11 < 0) x11 = gVirtualX->InheritsFrom("TGX11");

   if (x11) { //use built-in optimized version
      fImage = picture2asimage(fgVisual, pxm, mask, 0, 0, w, h, kAllPlanes, 1, 0);
   } else {
      unsigned char *bits = gVirtualX->GetColorBits(pxm, 0, 0, w, h);
      if (!bits) {   // error
         return;
      }

      // no mask
      if (!mask) {
         fImage = bitmap2asimage(bits, w, h, 0, 0);
         delete [] bits;
         return;
      }
      unsigned char *mask_bits = gVirtualX->GetColorBits(mask, 0, 0, w, h);
      fImage = bitmap2asimage(bits, w, h, 0, mask_bits);
      delete [] mask_bits;
      delete [] bits;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return 2D array of machine dependent pixel values.

TArrayL *TASImage::GetPixels(Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   if (!fImage) {
      Warning("GetPixels", "Wrong Image");
      return 0;
   }

   ASImage *img =  fScaledImage ? fScaledImage->fImage : fImage;
   ASImageDecoder *imdec;

   width = !width  ? img->width : width;
   height = !height ? img->height : height;

   if (x < 0) {
      width -= x;
      x = 0 ;
   }
   if (y < 0) {
      height -= y;
      y = 0;
   }

   if ((x >= (int)img->width) || (y >= (int)img->height)) {
      return 0;
   }

   if ((int)(x + width) > (int)img->width) {
      width = img->width - x;
   }

   if ((int)(y + height) > (int)img->height) {
      height = img->height - y;
   }

   if ((imdec = start_image_decoding(0, fImage, SCL_DO_ALL, 0, y,
                                     img->width, height, 0)) == 0) {
      Warning("GetPixels", "Failed to create image decoder");
      return 0;
   }

   TArrayL *ret = new TArrayL(width * height);
   Int_t r = 0;
   Int_t g = 0;
   Int_t b = 0;
   Long_t p = 0;

   for (UInt_t k = 0; k < height; k++) {
      imdec->decode_image_scanline(imdec);

      for (UInt_t i = 0; i < width; ++i)  {
         if ((r == (Int_t)imdec->buffer.red[i]) &&
             (g == (Int_t)imdec->buffer.green[i]) &&
             (b == (Int_t)imdec->buffer.blue[i])) {
         } else {
            r = (Int_t)imdec->buffer.red[i];
            g = (Int_t)imdec->buffer.green[i];
            b = (Int_t)imdec->buffer.blue[i];
            p = (Long_t)TColor::RGB2Pixel(r, g, b);
         }
         ret->AddAt(p, k*width + i);
      }
   }

   stop_image_decoding(&imdec);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to internal array[width x height] of double values [0,1].
/// This array is directly accessible. That allows to manipulate/change the
/// image.

Double_t *TASImage::GetVecArray()
{
   if (!fImage) {
      Warning("GetVecArray", "Bad Image");
      return 0;
   }
   if (fImage->alt.vector) {
      return fImage->alt.vector;
   }
   // vectorize
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// In case of vectorized image return an associated array of doubles
/// otherwise this method creates and returns a 2D array of doubles corresponding to palette.
/// If palette is ZERO a color converted to double value [0, 1] according to formula
/// ~~~ {.cpp}
///   Double_t((r << 16) + (g << 8) + b)/0xFFFFFF
/// ~~~
/// The returned array must be deleted after usage.

TArrayD *TASImage::GetArray(UInt_t w, UInt_t h, TImagePalette *palette)
{
   if (!fImage) {
      Warning("GetArray", "Bad Image");
      return 0;
   }

   TArrayD *ret;

   if (fImage->alt.vector) {
      ret = new TArrayD(fImage->width*fImage->height, fImage->alt.vector);
      return ret;
   }

   ASImageDecoder *imdec;

   w = w ? w : fImage->width;
   h = h ? h : fImage->height;

   if ((fImage->width != w) || (fImage->height != h)) {
      Scale(w, h);
   }

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;

   if ((imdec = start_image_decoding(0, img, SCL_DO_ALL, 0, 0,
                                     img->width, 0, 0)) == 0) {
      Warning("GetArray", "Failed to create image decoder");
      return 0;
   }

   ret = new TArrayD(w * h);
   CARD32 r = 0;
   CARD32 g = 0;
   CARD32 b = 0;
   Int_t p = 0;
   Double_t v = 0;

   for (UInt_t k = 0; k < h; k++) {
      imdec->decode_image_scanline(imdec);

      for (UInt_t i = 0; i < w; ++i)  {
         if ((r == imdec->buffer.red[i]) &&
             (g == imdec->buffer.green[i]) &&
             (b == imdec->buffer.blue[i])) {
         } else {
            r = imdec->buffer.red[i];
            g = imdec->buffer.green[i];
            b = imdec->buffer.blue[i];
            if (palette) p = palette->FindColor(r, g, b);
         }
         v = palette ? palette->fPoints[p] : Double_t((r << 16) + (g << 8) + b)/0xFFFFFF;
         ret->AddAt(v, (h-k-1)*w + i);
      }
   }

   stop_image_decoding(&imdec);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text of size (in pixels for TrueType fonts)
/// at position (x, y) with color  specified by hex string.
///
///  - font_name: TrueType font's filename or X font spec or alias.
///               3D style of text is one of the following:
///    * 0  plain 2D text,
///    * 1  embossed,
///    * 2  sunken,
///    * 3  shade above,
///    * 4  shade below,
///    * 5  embossed thick,
///    * 6  sunken thick.
///    * 7  outline above,
///    * 8  ouline below,
///    * 9  full ouline.
///  - fore_file specifies foreground texture of text.

void TASImage::DrawText(Int_t x, Int_t y, const char *text, Int_t size,
                        const char *color, const char *font_name,
                        EText3DType type, const char *fore_file, Float_t angle)
{
   UInt_t width=0, height=0;
   ARGB32 text_color = ARGB32_Black;
   ASImage *fore_im = 0;
   ASImage *text_im = 0;
   Bool_t ttfont = kFALSE;

   if (!InitVisual()) {
      Warning("DrawText", "Visual not initiated");
      return;
   }

   TString fn = font_name;
   fn.Strip();

   // This is for backward compatibility...
   if (fn.Last('/') == 0) fn = fn(1, fn.Length() - 1);

   const char *ttpath = gEnv->GetValue("Root.TTFontPath",
                                       TROOT::GetTTFFontDir());
   char *tmpstr = gSystem->Which(ttpath, fn, kReadPermission);
   fn = tmpstr;
   delete [] tmpstr;

   if (fn.EndsWith(".pfa") || fn.EndsWith(".PFA") || fn.EndsWith(".pfb") || fn.EndsWith(".PFB") || fn.EndsWith(".ttf") || fn.EndsWith(".TTF") || fn.EndsWith(".otf") || fn.EndsWith(".OTF")) {
      ttfont = kTRUE;
   }

   if (color) {
      parse_argb_color(color, &text_color);
   }

   if (fImage && fImage->alt.argb32 && ttfont) {
      DrawTextTTF(x, y, text, size, text_color, fn.Data(), angle);
      return;
   }

   if (!gFontManager) {
      gFontManager = create_font_manager(fgVisual->dpy, 0, 0);
   }

   if (!gFontManager) {
      Warning("DrawText", "cannot create Font Manager");
      return;
   }

   ASFont *font = get_asfont(gFontManager, fn.Data(), 0, size, ASF_GuessWho);

   if (!font) {
      font = get_asfont(gFontManager, "fixed", 0, size, ASF_GuessWho);
      if (!font) {
         Warning("DrawText", "cannot find a font %s", font_name);
         return;
      }
   }

   get_text_size(text, font, (ASText3DType)type, &width, &height);

   if (!fImage) {
      fImage = create_asimage(width, height, 0);
      fill_asimage(fgVisual, fImage, 0, 0, width, height, 0xFFFFFFFF);
   }

   text_im = draw_text(text, font, (ASText3DType)type, 0);

   ASImage *rimg = fImage;

   if (fore_file) {
      ASImage *tmp = file2ASImage(fore_file, 0xFFFFFFFF, SCREEN_GAMMA, 0, 0);
      if (tmp) {
         if ((tmp->width != width) || (tmp->height != height)) {
            fore_im = tile_asimage(fgVisual, tmp, 0, 0, width, height, 0,
                                   ASA_ASImage, GetImageCompression(), GetImageQuality());
         }
         destroy_asimage(&tmp);
      } else {
         fore_im = tmp;
      }
   }

   if (fore_im) {
      move_asimage_channel(fore_im, IC_ALPHA, text_im, IC_ALPHA);
      destroy_asimage(&text_im);
   } else {
      fore_im = text_im ;
   }

   release_font(font);

   if (fore_im) {
      ASImage *rendered_im;
      ASImageLayer layers[2];

      init_image_layers(&(layers[0]), 2);
      fore_im->back_color = text_color;
      layers[0].im = rimg;
      layers[0].dst_x = 0;
      layers[0].dst_y = 0;
      layers[0].clip_width = rimg->width;
      layers[0].clip_height = rimg->height;
      layers[0].bevel = 0;
      layers[1].im = fore_im;
      layers[1].dst_x = x;
      layers[1].dst_y = y;
      layers[1].clip_width = fore_im->width;
      layers[1].clip_height = fore_im->height;

      rendered_im = merge_layers(fgVisual, &(layers[0]), 2, rimg->width, rimg->height,
                                 ASA_ASImage, GetImageCompression(), GetImageQuality());

      destroy_asimage(&fore_im);
      DestroyImage();
      fImage = rendered_im;
      UnZoom();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Merge two images.
///
/// op is string which specifies overlay operation. Supported operations are:
///
///  -  add            - color addition with saturation
///  -  alphablend     - alpha-blending
///  -  allanon        - color values averaging
///  -  colorize       - hue and saturate bottom image same as top image
///  -  darken         - use lowest color value from both images
///  -  diff           - use absolute value of the color difference between two images
///  -  dissipate      - randomly alpha-blend images
///  -  hue            - hue bottom image same as top image
///  -  lighten        - use highest color value from both images
///  -  overlay        - some weird image overlaying(see GIMP)
///  -  saturate       - saturate bottom image same as top image
///  -  screen         - another weird image overlaying(see GIMP)
///  -  sub            - color substraction with saturation
///  -  tint           - tinting image with image
///  -  value          - value bottom image same as top image

void TASImage::Merge(const TImage *im, const char *op, Int_t x, Int_t y)
{
   if (!im) return;

   if (!InitVisual()) {
      Warning("Merge", "Visual not initiated");
      return;
   }

   ASImage *rendered_im;
   ASImageLayer layers[2];

   init_image_layers(&(layers[0]), 2);
   layers[0].im = fImage;
   layers[0].dst_x = 0;
   layers[0].dst_y = 0;
   layers[0].clip_width = fImage->width;
   layers[0].clip_height = fImage->height;
   layers[0].bevel = 0;
   layers[1].im = ((TASImage*)im)->fImage;
   layers[1].dst_x = x;
   layers[1].dst_y = y;
   layers[1].clip_width  = im->GetWidth();
   layers[1].clip_height = im->GetHeight();
   layers[1].merge_scanlines = blend_scanlines_name2func(op ? op : "add");

   rendered_im = merge_layers(fgVisual, &(layers[0]), 2, fImage->width, fImage->height,
                              ASA_ASImage, GetImageCompression(), GetImageQuality());

   DestroyImage();
   fImage = rendered_im;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Perform Gaussian blur of the image (useful for drop shadows).
///  -  hr         - horizontal radius of the blur
///  -  vr         - vertical radius of the blur

void TASImage::Blur(Double_t hr, Double_t vr)
{
   if (!InitVisual()) {
      Warning("Blur", "Visual not initiated");
      return;
   }

   if (!fImage) {
      fImage = create_asimage(100, 100, 0);

      if (!fImage) {
         Warning("Blur", "Failed to create image");
         return;
      }

      fill_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height, ARGB32_White);
   }

   ASImage *rendered_im = blur_asimage_gauss(fgVisual, fImage, hr > 0 ? hr : 3,
                                             vr > 0 ? vr : 3, SCL_DO_ALL,
                                             ASA_ASImage, GetImageCompression(), GetImageQuality());
   DestroyImage();
   fImage = rendered_im;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Clone image.

TObject *TASImage::Clone(const char *newname) const
{
   if (!InitVisual() || !fImage) {
      Warning("Clone", "Image not initiated");
      return 0;
   }

   TASImage *im = (TASImage*)TImage::Create();

   if (!im) {
      Warning("Clone", "Failed to create image");
      return 0;
   }

   im->SetName(newname);

   im->fImage = clone_asimage(fImage, SCL_DO_ALL);
   im->fMaxValue = fMaxValue;
   im->fMinValue = fMinValue;
   im->fZoomOffX = fZoomOffX;
   im->fZoomOffY = fZoomOffY;
   im->fZoomWidth = fZoomWidth;
   im->fZoomHeight = fZoomHeight;
   im->fZoomUpdate = fZoomUpdate;
   im->fScaledImage = fScaledImage ? (TASImage*)fScaledImage->Clone("") : 0;

   if (fImage->alt.argb32) {
      UInt_t sz = fImage->width * fImage->height;
      im->fImage->alt.argb32 = (ARGB32*)safemalloc(sz*sizeof(ARGB32));
      memcpy(im->fImage->alt.argb32, fImage->alt.argb32, sz * sizeof(ARGB32));
   }

   return im;
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce color-depth of an image and fills vector of "scientific data"
/// [0...1]
///
/// Colors are reduced by allocating color cells to most used colors first,
/// and then approximating other colors with those allocated.
///
/// \param[in] max_colors       - maximum size of the colormap.
/// \param[in] dither           - number of bits to strip off the color data ( 0...7 )
/// \param[in] opaque_threshold - alpha channel threshold at which pixel should be treated as opaque

Double_t *TASImage::Vectorize(UInt_t max_colors, UInt_t dither, Int_t opaque_threshold)
{
   if (!InitVisual()) {
      Warning("Vectorize", "Visual not initiated");
      return 0;
   }

   if (!fImage) {
      fImage = create_asimage(100, 100, 0);

      if (!fImage) {
         Warning("Vectorize", "Failed to create image");
         return 0;
      }

      fill_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height, ARGB32_White);
   }

   ASColormap cmap;
   int *res;
   UInt_t r=0, g=0, b=0;

   dither = dither > 7 ? 7 : dither;

   res = colormap_asimage(fImage, &cmap, max_colors, dither, opaque_threshold);

   Double_t *vec = new Double_t[fImage->height*fImage->width];
   UInt_t v;
   Double_t tmp;
   fMinValue = 2;
   fMaxValue = -1;

   for (UInt_t y = 0; y < fImage->height; y++) {
      for (UInt_t x = 0; x < fImage->width; x++) {
         int i = y*fImage->width + x;
         if (res) {
            g = INDEX_SHIFT_GREEN(cmap.entries[res[i]].green);
            b = INDEX_SHIFT_BLUE(cmap.entries[res[i]].blue);
            r = INDEX_SHIFT_RED(cmap.entries[res[i]].red);
         }
         v = MAKE_INDEXED_COLOR24(r,g,b);
         v = (v>>12)&0x0FFF;
         tmp = Double_t(v)/0x0FFF;
         vec[(fImage->height - y - 1)*fImage->width + x] = tmp;
         if (fMinValue > tmp) fMinValue = tmp;
         if (fMaxValue < tmp) fMaxValue = tmp;
      }
   }
   TImagePalette *pal = new TImagePalette(cmap.count);

   for (UInt_t j = 0; j < cmap.count; j++) {
      g = INDEX_SHIFT_GREEN(cmap.entries[j].green);
      b = INDEX_SHIFT_BLUE(cmap.entries[j].blue);
      r = INDEX_SHIFT_RED(cmap.entries[j].red);
      v = MAKE_INDEXED_COLOR24(r,g,b);

      v = (v>>12) & 0x0FFF;
      pal->fPoints[j] = Double_t(v)/0x0FFF;

      pal->fColorRed[j] = cmap.entries[j].red << 8;
      pal->fColorGreen[j] = cmap.entries[j].green << 8;
      pal->fColorBlue[j] = cmap.entries[j].blue << 8;
      pal->fColorAlpha[j] = 0xFF00;
   }

   destroy_colormap(&cmap, kTRUE);

   fPalette = *pal;
   fImage->alt.vector = vec;
   UnZoom();
   // ROOT-7647: res is allocated with `safemalloc` by colormap_asimage
   if (res) safefree(res);
   return (Double_t*)fImage->alt.vector;
}

////////////////////////////////////////////////////////////////////////////////
/// This function will tile original image to specified size with offsets
/// requested, and then it will go though it and adjust hue, saturation and
/// value of those pixels that have specific hue, set by affected_hue/
/// affected_radius parameters. When affected_radius is greater then 180
/// entire image will be adjusted. Note that since grayscale colors have
/// no hue - the will not get adjusted. Only saturation and value will be
/// adjusted in gray pixels.
///
/// Hue is measured as an angle on a 360 degree circle, The following is
/// relationship of hue values to regular color names :
/// - red      - 0
/// - yellow   - 60
/// - green    - 120
/// - cyan     - 180
/// - blue     - 240
/// - magenta  - 300
/// - red      - 360
///
/// All the hue values in parameters will be adjusted to fall within 0-360 range.
///
/// \param[in] hue   hue in degrees in range 0-360. This allows to limit
///                  impact of color adjustment to affect only limited range of hues.
///
/// \param[in]  radius   value in degrees to be used in order to
///                      calculate the range of affected hues. Range is determined by
///                      substracting and adding this value from/to affected_hue.
///
/// \param[in] H      value by which to change hues in affected range.
/// \param[in] S     value by which to change saturation of the pixels in affected hue range.
/// \param[in] V     value by which to change Value(brightness) of pixels in affected hue range.
///
/// \param[in] x,y   position on infinite surface tiled with original image, of the
///                  left-top corner of the area to be used for new image.
///
/// \param[in] width, height   size of the area of the original image to be used for new image.
///                            Default is current width, height of the image.

void TASImage::HSV(UInt_t hue, UInt_t radius, Int_t H, Int_t S, Int_t V,
                   Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   if (!InitVisual()) {
      Warning("HSV", "Visual not initiated");
      return;
   }

   if (!fImage) {
      fImage = create_asimage(width ? width : 20, height ? height : 20, 0);

      if (!fImage) {
         Warning("HSV", "Failed to create image");
         return;
      }

      x = 0;
      y = 0;
      fill_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height, ARGB32_White);
   }

   width = !width ? fImage->width : width;
   height = !height ? fImage->height : height;

   ASImage *rendered_im = 0;

   if (H || S || V) {
      rendered_im = adjust_asimage_hsv(fgVisual, fImage, x, y, width, height,
                                       hue, radius, H, S, V, ASA_ASImage, 100,
                                       ASIMAGE_QUALITY_TOP);
   }
   if (!rendered_im) {
      Warning("HSV", "Failed to create rendered image");
      return;
   }

   DestroyImage();
   fImage = rendered_im;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Render multipoint gradient inside rectangle of size (width, height)
/// at position (x,y) within the existing image.
///
/// \param[in] angle    Given in degrees.  Default is 0.  This is the
///          direction of the gradient.  Currently the only supported
///          values are 0, 45, 90, 135, 180, 225, 270, 315.  0 means left
///          to right, 90 means top to bottom, etc.
///
/// \param[in] colors   Whitespace-separated list of colors.  At least two
///          colors are required.  Each color in this list will be visited
///          in turn, at the intervals given by the offsets attribute.
///
/// \param[in] offsets  Whitespace-separated list of floating point values
///          ranging from 0.0 to 1.0.  The colors from the colors attribute
///          are given these offsets, and the final gradient is rendered
///          from the combination of the two.  If both colors and offsets
///          are given but the number of colors and offsets do not match,
///          the minimum of the two will be used, and the other will be
///          truncated to match.  If offsets are not given, a smooth
///          stepping from 0.0 to 1.0 will be used.

void TASImage::Gradient(UInt_t angle, const char *colors, const char *offsets,
                        Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   if (!InitVisual()) {
      Warning("Gradient", "Visual not initiated");
      return;
   }

   ASImage *rendered_im = 0;
   ASGradient gradient;

   int reverse = 0, npoints1 = 0, npoints2 = 0;
   char *p;
   char *pb;
   char ch;
   TString str = colors;
   TString col;

   if ((angle > 2 * 180 * 15 / 16) || (angle < 2 * 180 * 1 / 16)) {
      gradient.type = GRADIENT_Left2Right;
   } else if (angle < 2 * 180 * 3 / 16) {
      gradient.type = GRADIENT_TopLeft2BottomRight;
   } else if (angle < 2 * 180 * 5 / 16) {
      gradient.type = GRADIENT_Top2Bottom;
   } else if (angle < 2 * 180 * 7 / 16) {
      gradient.type = GRADIENT_BottomLeft2TopRight; reverse = 1;
   } else if (angle < 2 * 180 * 9 / 16) {
      gradient.type = GRADIENT_Left2Right; reverse = 1;
   } else if (angle < 2 * 180 * 11 / 16) {
      gradient.type = GRADIENT_TopLeft2BottomRight; reverse = 1;
   } else if (angle < 2 * 180 * 13 / 16) {
      gradient.type = GRADIENT_Top2Bottom; reverse = 1;
   } else {
      gradient.type = GRADIENT_BottomLeft2TopRight;
   }

   for (p = (char*)colors; isspace((int)*p); p++) { }

   for (npoints1 = 0; *p; npoints1++) {
      if (*p) {
         for ( ; *p && !isspace((int)*p); p++) { }
      }
      for ( ; isspace((int)*p); p++) { }
   }
   if (offsets) {
      for (p = (char*)offsets; isspace((int)*p); p++) { }

      for (npoints2 = 0; *p; npoints2++) {
         if (*p) {
            for ( ; *p && !isspace((int)*p); p++) { }
         }
         for ( ; isspace((int)*p); p++) { }
      }
   }
   if (npoints1 > 1) {
      int i;
      if (offsets && (npoints1 > npoints2)) npoints1 = npoints2;

      if (!width) {
         width = fImage ? fImage->width : 20;
      }
      if (!height) {
         height = fImage ? fImage->height : 20;
      }

      gradient.color = new ARGB32[npoints1];
      gradient.offset = new double[npoints1];

      for (p = (char*)colors; isspace((int)*p); p++) { }

      for (npoints1 = 0; *p; ) {
         pb = p;

         if (*p) {
            for ( ; *p && !isspace((int)*p); p++) { }
         }
         for ( ; isspace((int)*p); p++) { }

         col = str(pb - colors, p - pb);

         if (parse_argb_color(col.Data(), gradient.color + npoints1) != col) {
            npoints1++;
         } else {
            Warning("Gradient", "Failed to parse color [%s] - defaulting to black", pb);
         }
      }

      if (offsets) {
         for (p = (char*)offsets; isspace((int)*p); p++) { }

         for (npoints2 = 0; *p; ) {
            pb = p;

            if (*p) {
               for ( ; *p && !isspace((int)*p); p++) { }
            }
            ch = *p; *p = '\0';
            gradient.offset[npoints2] = strtod(pb, &pb);

            if (pb == p) npoints2++;
            *p = ch;
            for ( ; isspace((int)*p); p++) { }
         }
      } else {
         for (npoints2 = 0; npoints2 < npoints1; npoints2++) {
            gradient.offset[npoints2] = (double)npoints2 / (npoints1 - 1);
         }
      }
      gradient.npoints = npoints1;

      if (npoints2 && (gradient.npoints > npoints2)) {
         gradient.npoints = npoints2;
      }
      if (reverse) {
         for (i = 0; i < gradient.npoints/2; i++) {
            int i2 = gradient.npoints - 1 - i;
            ARGB32 c = gradient.color[i];
            double o = gradient.offset[i];
            gradient.color[i] = gradient.color[i2];
            gradient.color[i2] = c;
            gradient.offset[i] = gradient.offset[i2];
            gradient.offset[i2] = o;
         }
         for (i = 0; i < gradient.npoints; i++) {
            gradient.offset[i] = 1.0 - gradient.offset[i];
         }
      }
      rendered_im = make_gradient(fgVisual, &gradient, width, height, SCL_DO_ALL,
                                  ASA_ASImage, GetImageCompression(), GetImageQuality());

      delete [] gradient.color;
      delete [] gradient.offset;
   }

   if (!rendered_im) {  // error
      Warning("Gradient", "Failed to create gradient image");
      return;
   }

   if (!fImage) {
      fImage = rendered_im;
      return;
   }

   ASImageLayer layers[2];

   init_image_layers(&(layers[0]), 2);
   layers[0].im = fImage;
   layers[0].dst_x = 0;
   layers[0].dst_y = 0;
   layers[0].clip_width = fImage->width;
   layers[0].clip_height = fImage->height;
   layers[0].bevel = 0;
   layers[1].im = rendered_im;
   layers[1].dst_x = x;
   layers[1].dst_y = y;
   layers[1].clip_width = width;
   layers[1].clip_height = height;
   layers[1].merge_scanlines = alphablend_scanlines;

   ASImage *merge_im = merge_layers(fgVisual, &(layers[0]), 2, fImage->width, fImage->height,
                                    ASA_ASImage, GetImageCompression(), GetImageQuality());
   if (!merge_im) {
      Warning("Gradient", "Failed to create merged image");
      return;
   }

   destroy_asimage(&rendered_im);
   DestroyImage();
   fImage = merge_im;
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Make component hilite.
/// (used internally)

static CARD8 MakeComponentHilite(int cmp)
{
   if (cmp < 51) {
      cmp = 51;
   }
   cmp = (cmp * 12) / 10;

   return (cmp > 255) ? 255 : cmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate highlite color.
/// (used internally)

static ARGB32 GetHilite(ARGB32 background)
{
   return ((MakeComponentHilite((background>>24) & 0x000000FF) << 24) & 0xFF000000) |
           ((MakeComponentHilite((background & 0x00FF0000) >> 16) << 16) & 0x00FF0000) |
           ((MakeComponentHilite((background & 0x0000FF00) >> 8) << 8) & 0x0000FF00) |
           ((MakeComponentHilite((background & 0x000000FF))) & 0x000000FF);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate shadow color.
/// (used internally)

static ARGB32 GetShadow(ARGB32 background)
{
   return (background >> 1) & 0x7F7F7F7F;
}

////////////////////////////////////////////////////////////////////////////////
/// Get average.
/// (used internally)

static ARGB32 GetAverage(ARGB32 foreground, ARGB32 background)
{
   CARD16   a, r, g, b;

   a = ARGB32_ALPHA8(foreground) + ARGB32_ALPHA8(background);
   a = (a<<3)/10;
   r = ARGB32_RED8(foreground) + ARGB32_RED8(background);
   r = (r<<3)/10;
   g = ARGB32_GREEN8(foreground) + ARGB32_GREEN8(background);
   g = (g<<3)/10;
   b = ARGB32_BLUE8(foreground) + ARGB32_BLUE8(background);
   b = (b<<3)/10;

   return MAKE_ARGB32(a, r, g, b);
}


////////////////////////////////////////////////////////////////////////////////
/// Bevel is used to create 3D effect while drawing buttons, or any other
/// image that needs to be framed. Bevel is drawn using 2 primary colors:
/// one for top and left sides - hi color, and another for bottom and
/// right sides - low color. Bevel can be drawn over existing image or
/// as newly created,  as it is shown in code below:
/// ~~~ {.cpp}
///  TImage *img = TImage::Create();
///  img->Bevel(0, 0, 400, 300, "#dddddd", "#000000", 3);
/// ~~~

void TASImage::Bevel(Int_t x, Int_t y, UInt_t width, UInt_t height,
                     const char *hi_color, const char *lo_color, UShort_t thick,
                     Bool_t reverse)
{
   if (!InitVisual()) {
      Warning("Bevel", "Visual not initiated");
      return;
   }

   ASImageBevel bevel;
   bevel.type = 0;

   ARGB32 hi=ARGB32_White, lo=ARGB32_White;
   parse_argb_color(hi_color, &hi);
   parse_argb_color(lo_color, &lo);

   if (reverse) {
      bevel.lo_color = hi;
      bevel.lolo_color = GetHilite(hi);
      bevel.hi_color = lo;
      bevel.hihi_color = GetShadow(lo);
   } else {
      bevel.hi_color = hi;
      bevel.hihi_color = GetHilite(hi);
      bevel.lo_color = lo;
      bevel.lolo_color = GetShadow(lo);
   }
   bevel.hilo_color = GetAverage(hi, lo);

   int extra_hilite = 2;
   bevel.left_outline = bevel.top_outline = bevel.right_outline = bevel.bottom_outline = thick;
   bevel.left_inline = bevel.top_inline = bevel.right_inline = bevel.bottom_inline = extra_hilite + 1;

   if (bevel.top_outline > 1) {
      bevel.top_inline += bevel.top_outline - 1;
   }

   if (bevel.left_outline > 1) {
      bevel.left_inline += bevel.left_outline - 1;
   }

   if (bevel.right_outline > 1) {
      bevel.right_inline += bevel.right_outline - 1;
   }

   if (bevel.bottom_outline > 1) {
      bevel.bottom_inline += bevel.bottom_outline - 1;
   }

   ASImage *merge_im;
   ARGB32 fill = ((hi>>24) != 0xff) || ((lo>>24) != 0xff) ? bevel.hilo_color : (bevel.hilo_color | 0xff000000);

   if (!fImage) {
      fImage = create_asimage(width ? width : 20, height ? height : 20, 0);

      if (!fImage) {
         Warning("Bevel", "Failed to create image");
         return;
      }

      x = 0;
      y = 0;
      fill_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height, fill);
   }

   width = !width ? fImage->width : width;
   height = !height ? fImage->height : height;

   ASImageLayer layers[2];
   init_image_layers(&(layers[0]), 2);

   layers[0].im = fImage;
   layers[0].dst_x = 0;
   layers[0].dst_y = 0;
   layers[0].clip_width = fImage->width;
   layers[0].clip_height = fImage->height;
   layers[0].bevel = 0;

   UInt_t w = width - (bevel.left_outline + bevel.right_outline);
   UInt_t h = height - (bevel.top_outline + bevel.bottom_outline);
   ASImage *bevel_im = create_asimage(w, h, 0);

   if (!bevel_im) {
      Warning("Bevel", "Failed to create bevel image");
      return;
   }

   layers[1].im = bevel_im;
   fill_asimage(fgVisual, bevel_im, 0, 0, w, h, fill);

   layers[1].dst_x = x;
   layers[1].dst_y = y;
   layers[1].clip_width = width;
   layers[1].clip_height = height;
   layers[1].bevel = &bevel;
   layers[1].merge_scanlines = alphablend_scanlines;

   merge_im = merge_layers(fgVisual, &(layers[0]), 2, fImage->width, fImage->height,
                           ASA_ASImage, GetImageCompression(), GetImageQuality());
   destroy_asimage(&bevel_im);

   if (!merge_im) {
      Warning("Bevel", "Failed to image");
      return;
   }

   DestroyImage();
   fImage = merge_im;
   UnZoom();
}


////////////////////////////////////////////////////////////////////////////////
/// Enlarge image, padding it with specified color on each side in
/// accordance with requested geometry.

void TASImage::Pad(const char *col, UInt_t l, UInt_t r, UInt_t t, UInt_t b)
{
   Int_t x, y;
   UInt_t w, h;

   if (!InitVisual()) {
      Warning("Pad", "Visual not initiated");
      return;
   }

   if (!fImage) {
      fImage = create_asimage(100, 100, 0);

      if (!fImage) {
         Warning("Pad", "Failed to create image");
         return;
      }

      x = 0;
      y = 0;
      fill_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height, ARGB32_White);
   }

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   x = l;
   y = t;
   w = l + fImage->width + r;
   h = t + fImage->height + b;

   ASImage *img = pad_asimage(fgVisual, fImage, x, y, w, h, color,
                              ASA_ASImage, GetImageCompression(), GetImageQuality());

   if (!img) {
      Warning("Pad", "Failed to create output image");
      return;
   }

   DestroyImage();
   fImage = img;
   UnZoom();
   fZoomUpdate = kZoomOps;
}


////////////////////////////////////////////////////////////////////////////////
/// Crop an image.

void TASImage::Crop(Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   if (!InitVisual()) {
      Warning("Crop", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("Crop", "No image");
      return;
   }

   x = x < 0 ? 0 : x;
   y = y < 0 ? 0 : y;

   width = x + width > fImage->width ? fImage->width - x : width;
   height = y + height > fImage->height ? fImage->height - y : height;

   if ((width == fImage->width) && (height == fImage->height)) {
      Warning("Crop", "input size larger than image");
      return;
   }
   ASImageDecoder *imdec = start_image_decoding(fgVisual, fImage, SCL_DO_ALL,
                                                x, y, width, height, 0);

   if (!imdec) {
      Warning("Crop", "Failed to start image decoding");
      return;
   }

   ASImage *img = create_asimage(width, height, 0);

   if (!img) {
      delete [] imdec;
      Warning("Crop", "Failed to create image");
      return;
   }

   ASImageOutput *imout = start_image_output(fgVisual, img, ASA_ASImage,
                                             GetImageCompression(), GetImageQuality());

   if (!imout) {
      Warning("Crop", "Failed to start image output");
      destroy_asimage(&img);
      if (imdec) delete [] imdec;
      return;
   }

#ifdef HAVE_MMX
   mmx_init();
#endif

   for (UInt_t i = 0; i < height; i++) {
      imdec->decode_image_scanline(imdec);
      imout->output_image_scanline(imout, &(imdec->buffer), 1);
   }

   stop_image_decoding(&imdec);
   stop_image_output(&imout);

#ifdef HAVE_MMX
   mmx_off();
#endif

   DestroyImage();
   fImage = img;
   UnZoom();
   fZoomUpdate = kZoomOps;
}

////////////////////////////////////////////////////////////////////////////////
/// Append image.
///
/// option:
///  - "+" - appends to the right side
///  - "/" - appends to the bottom

void TASImage::Append(const TImage *im, const char *option, const char *color )
{
   if (!im) return;

   if (!InitVisual()) {
      Warning("Append", "Visual not initiated");
      return;
   }

   if (!fImage) {
      fImage = ((TASImage*)im)->fImage;
      return;
   }

   TString opt = option;
   opt.Strip();

   UInt_t width = fImage->width;
   UInt_t height = fImage->height;

   if (opt == "+") {
      Pad(color, 0, im->GetWidth(), 0, 0);
      Merge(im, "alphablend", width, 0);
   } else if (opt == "/") {
      Pad(color, 0, 0, 0, im->GetHeight());
      Merge(im, "alphablend", 0, height);
   } else {
      return;
   }

   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// BeginPaint initializes internal array[width x height] of ARGB32 pixel
/// values.
///
/// That provides quick access to image during paint operations.
/// To RLE compress image one needs to call EndPaint method when painting
/// is over.

void TASImage::BeginPaint(Bool_t mode)
{
   if (!InitVisual()) {
      Warning("BeginPaint", "Visual not initiated");
      return;
   }

   if (!fImage) {
      return;
   }

   fPaintMode = mode;

   if (!fPaintMode || fImage->alt.argb32) {
      return;
   }

   ASImage *img = tile_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height,
                               0, ASA_ARGB32, 0, ASIMAGE_QUALITY_DEFAULT);

   if (!img) {
      Warning("BeginPaint", "Failed to create image");
      return;
   }

   DestroyImage();
   fImage = img;
}

////////////////////////////////////////////////////////////////////////////////
/// EndPaint does internal RLE compression of image data.

void TASImage::EndPaint()
{
   if (!fImage) {
      Warning("EndPaint", "no image");
      return;
   }

   if (!fImage->alt.argb32) return;

   ASImage *img = tile_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height,
                               0, ASA_ASImage, 0, ASIMAGE_QUALITY_DEFAULT);

   if (!img) {
      Warning("EndPaint", "Failed to create image");
      return;
   }

   fPaintMode = kFALSE;
   DestroyImage();
   fImage = img;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to internal array[width x height] of ARGB32 values
/// This array is directly accessible. That allows to manipulate/change the
/// image.

UInt_t *TASImage::GetArgbArray()
{
   if (!fImage) {
      Warning("GetArgbArray", "no image");
      return 0;
   }

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;
   if (!img) return 0;

   if (!img->alt.argb32) {
      if (fScaledImage) {
         fScaledImage->BeginPaint();
         img = fScaledImage->fImage;
      } else {
         BeginPaint();
         img = fImage;
      }
   }

   return (UInt_t *)img->alt.argb32;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to an array[width x height] of RGBA32 values.
/// This array is created from internal ARGB32 array,
/// must be deleted after usage.

UInt_t *TASImage::GetRgbaArray()
{
   if (!fImage) {
      Warning("GetRgbaArray", "no image");
      return 0;
   }

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;
   if (!img) return 0;

   if (!img->alt.argb32) {
      if (fScaledImage) {
         fScaledImage->BeginPaint();
         img = fScaledImage->fImage;
      } else {
         BeginPaint();
         img = fImage;
      }
   }

   UInt_t i, j;
   Int_t y = 0;
   Int_t idx = 0;
   UInt_t a, rgb, rgba, argb;
   y = 0;

   UInt_t *ret = new UInt_t[img->width*img->height];

   for (i = 0; i < img->height; i++) {
      for (j = 0; j < img->width; j++) {
         idx = Idx(y + j);
         argb = img->alt.argb32[idx];
         a = argb >> 24;
         rgb =  argb & 0x00ffffff;
         rgba = (rgb <<  8) + a;
         ret[idx] = rgba;
      }
      y += img->width;
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to scan-line.

UInt_t *TASImage::GetScanline(UInt_t y)
{
   if (!fImage) {
      Warning("GetScanline", "no image");
      return 0;
   }

   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;
   CARD32 *ret = new CARD32[img->width];

   ASImageDecoder *imdec = start_image_decoding(fgVisual, img, SCL_DO_ALL,
                                                0, y, img->width, 1, 0);

   if (!imdec) {
      delete [] ret;
      Warning("GetScanline", "Failed to start image decoding");
      return 0;
   }

#ifdef HAVE_MMX
   mmx_init();
#endif

   imdec->decode_image_scanline(imdec);
   memcpy(imdec->buffer.buffer, ret, img->width*sizeof(CARD32));
   stop_image_decoding(&imdec);

#ifdef HAVE_MMX
   mmx_off();
#endif

   return (UInt_t*)ret;
}


//______________________________________________________________________________
//
// Vector graphics
// a couple of macros which can be "assembler accelerated"

#if defined(R__GNU) && defined(__i386__) && !defined(__sun)
#define _MEMSET_(dst, lng, val)   __asm__("movl  %0,%%eax \n"\
                                      "movl  %1,%%edi \n"              \
                                      "movl  %2,%%ecx \n"              \
                                      "cld \n"                         \
                                      "rep \n"                         \
                                      "stosl \n"                       \
                                      : /* no output registers */      \
                                      :"g" (val), "g" (dst), "g" (lng) \
                                      :"eax","edi","ecx"               \
                                     )

#else
 #define _MEMSET_(dst, lng, val) do {\
 for( UInt_t j=0; j < lng; j++) *((dst)+j) = val; } while (0)

#endif

#define FillSpansInternal(npt, ppt, widths, color) do {\
   UInt_t yy = ppt[0].fY*fImage->width;\
   for (UInt_t i = 0; i < npt; i++) {\
      _MEMSET_(&fImage->alt.argb32[Idx(yy + ppt[i].fX)], widths[i], color);\
      yy += ((i+1 < npt) && (ppt[i].fY != ppt[i+1].fY) ? fImage->width : 0);\
   }\
} while (0)

////////////////////////////////////////////////////////////////////////////////
/// Fill rectangle of size (width, height) at position (x,y)
/// within the existing image with specified color.

void TASImage::FillRectangleInternal(UInt_t col, Int_t x, Int_t y, UInt_t width, UInt_t height)
{

   if (!InitVisual()) {
      Warning("FillRectangle", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("FillRectangle", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("FillRectangle", "Failed to get pixel array");
      return;
   }

   ARGB32 color = (ARGB32)col;

   if (width  == 0) width = 1;
   if (height == 0) height = 1;

   if (x < 0) {
      width += x;
      x = 0;
   }
   if (y < 0) {
      height += y;
      y = 0;
   }

   Bool_t has_alpha = (color & 0xff000000) != 0xff000000;

   x = x > (int)fImage->width ? (Int_t)fImage->width : x;
   y = y > (int)fImage->height ? (Int_t)fImage->height : y;

   width = x + width > fImage->width ? fImage->width - x : width;
   height = y + height > fImage->height ? fImage->height - y : height;

   if (!fImage->alt.argb32) {
      fill_asimage(fgVisual, fImage, x, y, width, height, color);
   } else {
      int yyy = y*fImage->width;
      if (!has_alpha) { // use faster memset
         ARGB32 *p0 = fImage->alt.argb32 + yyy + x;
         ARGB32 *p = p0;
         for (UInt_t i = 0; i < height; i++) {
            _MEMSET_(p, width, color);
            p += fImage->width;
         }
      } else {
         for (UInt_t i = y; i < y + height; i++) {
            int j = x + width;
            while (j > x) {
               j--;
               _alphaBlend(&fImage->alt.argb32[Idx(yyy + j)], &color);
            }
            yyy += fImage->width;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill rectangle of size (width, height) at position (x,y)
/// within the existing image with specified color.
///
/// To create new image with Fill method the following code can be used:
/// ~~~ {.cpp}
///  TImage *img = TImage::Create();
///  img->Fill("#FF00FF", 0, 0, 400, 300);
/// ~~~

void TASImage::FillRectangle(const char *col, Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   if (!InitVisual()) {
      Warning("Fill", "Visual not initiated");
      return;
   }

   ARGB32 color = ARGB32_White;

   if (col) {
      parse_argb_color(col, &color);
   }

   if (!fImage) {
      fImage = create_asimage(width ? width : 20, height ? height : 20, 0);
      x = 0;
      y = 0;
   }

   FillRectangleInternal((UInt_t)color, x, y, width, height);
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a vertical line.

void TASImage::DrawVLine(UInt_t x, UInt_t y1, UInt_t y2, UInt_t col, UInt_t thick)
{
   ARGB32 color = (ARGB32)col;
   UInt_t half = 0;

   if (!thick)  thick = 1;

   if (thick > 1) {
      half = thick >> 1;
      if (x > half) {
         x =  x - half;
      } else {
         x = 0;
         thick += (x - half);
      }
   }

   y2 = y2 >= fImage->height ? fImage->height - 1 : y2;
   y1 = y1 >= fImage->height ? fImage->height - 1 : y1;
   x = x + thick >= fImage->width ? fImage->width - thick - 1 : x;

   int yy = y1*fImage->width;
   for (UInt_t y = y1; y <= y2; y++) {
      for (UInt_t w = 0; w < thick; w++) {
         if (x + w < fImage->width) {
            _alphaBlend(&fImage->alt.argb32[Idx(yy + (x + w))], &color);
         }
      }
      yy += fImage->width;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw an horizontal line.

void TASImage::DrawHLine(UInt_t y, UInt_t x1, UInt_t x2, UInt_t col, UInt_t thick)
{
   ARGB32 color = (ARGB32)col;
   UInt_t half = 0;

   if (!thick)  thick = 1;

   if (thick > 1) {
      half = thick >> 1;
      if (y > half) {
         y =  y - half;
      } else {
         y = 0;
         thick += (y - half);
      }
   }

   y = y + thick >= fImage->height ? fImage->height - thick - 1 : y;
   x2 = x2 >= fImage->width ? fImage->width - 1 : x2;
   x1 = x1 >= fImage->width ? fImage->width - 1 : x1;

   int yy = y*fImage->width;
   for (UInt_t w = 0; w < thick; w++) {
      for (UInt_t x = x1; x <= x2; x++) {
         if (y + w < fImage->height) {
            _alphaBlend(&fImage->alt.argb32[Idx(yy + x)], &color);
         }
      }
      yy += fImage->width;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line.

void TASImage::DrawLine(UInt_t x1,  UInt_t y1, UInt_t x2, UInt_t y2,
                        const char *col, UInt_t thick)
{
   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);
   DrawLineInternal(x1, y1, x2, y2, (UInt_t)color, thick);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal line drawing.

void TASImage::DrawLineInternal(UInt_t x1, UInt_t y1, UInt_t x2, UInt_t y2,
                                UInt_t col, UInt_t thick)
{
   int dx, dy, d;
   int i1, i2;
   int x, y, xend, yend;
   int xdir, ydir;
   int q;
   int idx;
   int yy;

   if (!InitVisual()) {
      Warning("DrawLine", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("DrawLine", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("DrawLine", "Failed to get pixel array");
      return;
   }

   ARGB32 color = (ARGB32)col;

   dx = TMath::Abs(Int_t(x2) - Int_t(x1));
   dy = TMath::Abs(Int_t(y2) - Int_t(y1));

   if (!dx && !dy) return; // invisible line

   if (!dx) {
      DrawVLine(x1, y2 > y1 ? y1 : y2,
                    y2 > y1 ? y2 : y1, color, thick);
      return;
   }

   if (!dy) {
      DrawHLine(y1, x2 > x1 ? x1 : x2,
                    x2 > x1 ? x2 : x1, color, thick);
      return;
   }

   if (thick > 1) {
      DrawWideLine(x1, y1, x2, y2, color, thick);
      return;
   }

   if (dy <= dx) {
      UInt_t ddy = dy << 1;
      i1 = ddy;
      i2 = i1 - (dx << 1);
      d = i1 - dx;

      if (x1 > x2) {
         x = x2;
         y = y2;
         ydir = -1;
         xend = x1;
      } else {
         x = x1;
         y = y1;
         ydir = 1;
         xend = x2;
      }

      yy = y*fImage->width;
      _alphaBlend(&fImage->alt.argb32[Idx(yy + x)], &color);
      q = (y2 - y1) * ydir;

      if (q > 0) {
         while (x < xend) {

            idx = Idx(yy + x);
            _alphaBlend(&fImage->alt.argb32[idx], &color);
            x++;

            if (d >= 0) {
               yy += fImage->width;
               d += i2;
            } else {
               d += i1;
            }
         }
      } else {
         while (x < xend) {
            idx = Idx(yy + x);
            _alphaBlend(&fImage->alt.argb32[idx], &color);
            x++;

            if (d >= 0) {
               yy -= fImage->width;
               d += i2;
            } else {
               d += i1;
            }
         }
      }
   } else {
      UInt_t ddx = dx << 1;
      i1 = ddx;
      i2 = i1 - (dy << 1);
      d = i1 - dy;

      if (y1 > y2) {
         y = y2;
         x = x2;
         yend = y1;
         xdir = -1;
      } else {
         y = y1;
         x = x1;
         yend = y2;
         xdir = 1;
      }

      yy = y*fImage->width;
      _alphaBlend(&fImage->alt.argb32[Idx(yy + x)], &color);
      q = (x2 - x1) * xdir;

      if (q > 0) {
         while (y < yend) {
            idx = Idx(yy + x);
            _alphaBlend(&fImage->alt.argb32[idx], &color);
            y++;
            yy += fImage->width;

            if (d >= 0) {
               x++;
               d += i2;
            } else {
               d += i1;
            }
         }
      } else {
         while (y < yend) {
            idx = Idx(yy + x);
            _alphaBlend(&fImage->alt.argb32[idx], &color);
            y++;
            yy += fImage->width;

            if (d >= 0) {
               x--;
               d += i2;
            } else {
               d += i1;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a rectangle.

void TASImage::DrawRectangle(UInt_t x, UInt_t y, UInt_t w, UInt_t h,
                             const char *col, UInt_t thick)
{
   if (!InitVisual()) {
      Warning("DrawRectangle", "Visual not initiated");
      return;
   }

   if (!fImage) {
      w = w ? w : 20;
      h = h ? h : 20;
      x = 0;
      y = 0;
      fImage = create_asimage(w, h, 0);
      FillRectangle(col, 0, 0, w, h);
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("DrawRectangle", "Failed to get pixel array");
      return;
   }

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   DrawHLine(y, x, x + w, (UInt_t)color, thick);
   DrawVLine(x + w, y, y + h, (UInt_t)color, thick);
   DrawHLine(y + h, x, x + w, (UInt_t)color, thick);
   DrawVLine(x, y, y + h, (UInt_t)color, thick);
   UnZoom();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a box.

void TASImage::DrawBox(Int_t x1, Int_t y1, Int_t x2, Int_t y2, const char *col,
                       UInt_t thick, Int_t mode)
{
   Int_t x = TMath::Min(x1, x2);
   Int_t y = TMath::Min(y1, y2);
   Int_t w = TMath::Abs(x2 - x1);
   Int_t h = TMath::Abs(y2 - y1);

   ARGB32 color = ARGB32_White;

   if (!fImage) {
      w = w ? x+w : x+20;
      h = h ? y+h : y+20;
      fImage = create_asimage(w, h, 0);
      FillRectangle(col, 0, 0, w, h);
      return;
   }

   if (x1 == x2) {
      parse_argb_color(col, &color);
      DrawVLine(x1, y1, y2, color, 1);
      return;
   }

   if (y1 == y2) {
      parse_argb_color(col, &color);
      DrawHLine(y1, x1, x2, color, 1);
      return;
   }


   switch (mode) {
      case TVirtualX::kHollow:
         DrawRectangle(x, y, w, h, col, thick);
         break;

      case TVirtualX::kFilled:
         FillRectangle(col, x, y, w, h);
         break;

      default:
         FillRectangle(col, x, y, w, h);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a dashed horizontal line.

void TASImage::DrawDashHLine(UInt_t y, UInt_t x1, UInt_t x2, UInt_t nDash,
                             const char *pDash, UInt_t col, UInt_t thick)
{
   UInt_t iDash = 0;    // index of current dash
   int i = 0;

   ARGB32 color = (ARGB32)col;

   UInt_t half = 0;

   if (thick > 1) {
      half = thick >> 1;
      if (y > half) {
         y =  y - half;
      } else {
         y = 0;
         thick += (y - half);
      }
   }
   thick = thick <= 0 ? 1 : thick;

   y = y + thick >= fImage->height ? fImage->height - thick - 1 : y;
   x2 = x2 >= fImage->width ? fImage->width - 1 : x2;
   x1 = x1 >= fImage->width ? fImage->width - 1 : x1;

   // switch x1, x2
   UInt_t tmp = x1;
   x1 = x2 < x1 ? x2 : x1;
   x2 = x2 < tmp ? tmp : x2;

   for (UInt_t x = x1; x <= x2; x++) {
      for (UInt_t w = 0; w < thick; w++) {
         if (y + w < fImage->height) {
            if ((iDash%2)==0) {
               _alphaBlend(&fImage->alt.argb32[Idx((y + w)*fImage->width + x)], &color);
            }
         }
      }
      i++;

      if (i >= pDash[iDash]) {
         iDash++;
         i = 0;
      }
      if (iDash >= nDash) {
         iDash = 0;
         i = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a dashed vertical line.

void TASImage::DrawDashVLine(UInt_t x, UInt_t y1, UInt_t y2, UInt_t nDash,
                             const char *pDash, UInt_t col, UInt_t thick)
{
   UInt_t iDash = 0;    // index of current dash
   int i = 0;

   ARGB32 color = (ARGB32)col;

   UInt_t half = 0;

   if (thick > 1) {
      half = thick >> 1;
      if (x > half) {
         x =  x - half;
      } else {
         x = 0;
         thick += (x - half);
      }
   }
   thick = thick <= 0 ? 1 : thick;

   y2 = y2 >= fImage->height ? fImage->height - 1 : y2;
   y1 = y1 >= fImage->height ? fImage->height - 1 : y1;

   // switch x1, x2
   UInt_t tmp = y1;
   y1 = y2 < y1 ? y2 : y1;
   y2 = y2 < tmp ? tmp : y2;

   x = x + thick >= fImage->width ? fImage->width - thick - 1 : x;

   int yy = y1*fImage->width;
   for (UInt_t y = y1; y <= y2; y++) {
      for (UInt_t w = 0; w < thick; w++) {
         if (x + w < fImage->width) {
            if ((iDash%2)==0) {
               _alphaBlend(&fImage->alt.argb32[Idx(yy + (x + w))], &color);
            }
         }
      }
      i++;

      if (i >= pDash[iDash]) {
         iDash++;
         i = 0;
      }
      if (iDash >= nDash) {
         iDash = 0;
         i = 0;
      }
      yy += fImage->width;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a dashed line with one pixel width.

void TASImage::DrawDashZLine(UInt_t x1, UInt_t y1, UInt_t x2, UInt_t y2,
                             UInt_t nDash, const char *tDash, UInt_t color)
{
   int dx, dy, d;
   int i, i1, i2;
   int x, y, xend, yend;
   int xdir, ydir;
   int q;
   UInt_t iDash = 0;    // index of current dash
   int yy;
   int idx;

   dx = TMath::Abs(Int_t(x2) - Int_t(x1));
   dy = TMath::Abs(Int_t(y2) - Int_t(y1));

   char *pDash = new char[nDash];

   if (dy <= dx) {
      double ac = TMath::Cos(TMath::ATan2(dy, dx));

      for (i = 0; i < (int)nDash; i++) {
         pDash[i] = TMath::Nint(tDash[i] * ac);
      }

      UInt_t ddy = dy << 1;
      i1 = ddy;
      i2 = i1 - (dx << 1);
      d = i1 - dx;
      i = 0;

      if (x1 > x2) {
         x = x2;
         y = y2;
         ydir = -1;
         xend = x1;
      } else {
         x = x1;
         y = y1;
         ydir = 1;
         xend = x2;
      }

      yy = y*fImage->width;
      _alphaBlend(&fImage->alt.argb32[Idx(y*fImage->width + x)], &color);
      q = (y2 - y1) * ydir;

      if (q > 0) {
         while (x < xend) {
            idx = Idx(yy + x);
            if ((iDash%2) == 0) {
               _alphaBlend(&fImage->alt.argb32[idx], &color);
            }
            x++;
            if (d >= 0) {
               yy += fImage->width;
               d += i2;
            } else {
               d += i1;
            }

            i++;
            if (i >= pDash[iDash]) {
               iDash++;
               i = 0;
            }
            if (iDash >= nDash) {
               iDash = 0;
               i = 0;
            }
         }
      } else {
         while (x < xend) {
            idx = Idx(yy + x);
            if ((iDash%2) == 0) {
               _alphaBlend(&fImage->alt.argb32[idx], &color);
            }
            x++;
            if (d >= 0) {
               yy -= fImage->width;
               d += i2;
            } else {
               d += i1;
            }

            i++;
            if (i >= pDash[iDash]) {
               iDash++;
               i = 0;
            }
            if (iDash >= nDash) {
               iDash = 0;
               i = 0;
            }
         }
      }
   } else {
      double as = TMath::Sin(TMath::ATan2(dy, dx));

      for (i = 0; i < (int)nDash; i++) {
         pDash[i] = TMath::Nint(tDash[i] * as);
      }

      UInt_t ddx = dx << 1;
      i1 = ddx;
      i2 = i1 - (dy << 1);
      d = i1 - dy;
      i = 0;

      if (y1 > y2) {
         y = y2;
         x = x2;
         yend = y1;
         xdir = -1;
      } else {
         y = y1;
         x = x1;
         yend = y2;
         xdir = 1;
      }

      yy = y*fImage->width;
      _alphaBlend(&fImage->alt.argb32[Idx(y*fImage->width + x)], &color);
      q = (x2 - x1) * xdir;

      if (q > 0) {
         while (y < yend) {
            idx = Idx(yy + x);
            if ((iDash%2) == 0) {
               _alphaBlend(&fImage->alt.argb32[idx], &color);
            }
            y++;
            yy += fImage->width;

            if (d >= 0) {
               x++;
               d += i2;
            } else {
               d += i1;
            }

            i++;
            if (i >= pDash[iDash]) {
               iDash++;
               i = 0;
            }
            if (iDash >= nDash) {
               iDash = 0;
               i = 0;
            }
         }
      } else {
         while (y < yend) {
            idx = Idx(yy + x);
            if ((iDash%2) == 0) {
               _alphaBlend(&fImage->alt.argb32[idx], &color);
            }
            y++;
            yy += fImage->width;

            if (d >= 0) {
               x--;
               d += i2;
            } else {
               d += i1;
            }

            i++;
            if (i >= pDash[iDash]) {
               iDash++;
               i = 0;
            }
            if (iDash >= nDash) {
               iDash = 0;
               i = 0;
            }
         }
      }
   }
   delete [] pDash;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a dashed line with thick pixel width.

void TASImage::DrawDashZTLine(UInt_t x1, UInt_t y1, UInt_t x2, UInt_t y2,
                             UInt_t nDash, const char *tDash, UInt_t color, UInt_t thick)
{
   int dx, dy;
   int i;
   double x, y, xend=0, yend=0, x0, y0;
   int xdir, ydir;
   int q;
   UInt_t iDash = 0;    // index of current dash

   dx = TMath::Abs(Int_t(x2) - Int_t(x1));
   dy = TMath::Abs(Int_t(y2) - Int_t(y1));

   double *xDash = new double[nDash];
   double *yDash = new double[nDash];
   double a = TMath::ATan2(dy, dx);
   double ac = TMath::Cos(a);
   double as = TMath::Sin(a);

   for (i = 0; i < (int)nDash; i++) {
      xDash[i] = tDash[i] * ac;
      yDash[i] = tDash[i] * as;

      // dirty trick (must be fixed)
      if ((i%2) == 0) {
         xDash[i] = xDash[i]/2;
         yDash[i] = yDash[i]/2;
      } else {
         xDash[i] = xDash[i]*2;
         yDash[i] = yDash[i]*2;
      }
   }

   if (dy <= dx) {
      if (x1 > x2) {
         x = x2;
         y = y2;
         ydir = -1;
         xend = x1;
      } else {
         x = x1;
         y = y1;
         ydir = 1;
         xend = x2;
      }

      q = (y2 - y1) * ydir;
      x0 = x;
      y0 = y;
      iDash = 0;
      yend = y + q;

      if (q > 0) {
         while ((x < xend) && (y < yend)) {
            x += xDash[iDash];
            y += yDash[iDash];

            if ((iDash%2) == 0) {
              DrawWideLine(TMath::Nint(x0), TMath::Nint(y0),
                           TMath::Nint(x), TMath::Nint(y), color, thick);
            } else {
               x0 = x;
               y0 = y;
            }

            iDash++;

            if (iDash >= nDash) {
               iDash = 0;
            }
        }
      } else {
         while ((x < xend) && (y > yend)) {
            x += xDash[iDash];
            y -= yDash[iDash];

            if ((iDash%2) == 0) {
               DrawWideLine(TMath::Nint(x0), TMath::Nint(y0),
                            TMath::Nint(x), TMath::Nint(y), color, thick);
            } else {
               x0 = x;
               y0 = y;
            }

            iDash++;

            if (iDash >= nDash) {
               iDash = 0;
            }
         }
      }
   } else {

      if (y1 > y2) {
         y = y2;
         x = x2;
         yend = y1;
         xdir = -1;
      } else {
         y = y1;
         x = x1;
         yend = y2;
         xdir = 1;
      }

      q = (x2 - x1) * xdir;
      x0 = x;
      y0 = y;
      iDash = 0;
      xend = x + q;

      if (q > 0) {
         while ((x < xend) && (y < yend)) {
            x += xDash[iDash];
            y += yDash[iDash];

            if ((iDash%2) == 0) {
               DrawWideLine(TMath::Nint(x0), TMath::Nint(y0),
                            TMath::Nint(x), TMath::Nint(y), color, thick);
            } else {
               x0 = x;
               y0 = y;
            }

            iDash++;

            if (iDash >= nDash) {
               iDash = 0;
            }
         }
      } else {
         while ((x > xend) && (y < yend)) {
            x -= xDash[iDash];
            y += yDash[iDash];

            if ((iDash%2) == 0) {
               DrawWideLine(TMath::Nint(x0), TMath::Nint(y0),
                            TMath::Nint(x), TMath::Nint(y), color, thick);
            } else {
               x0 = x;
               y0 = y;
            }

            iDash++;

            if (iDash >= nDash) {
               iDash = 0;
            }
         }
      }
   }
   delete [] xDash;
   delete [] yDash;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a dashed line.

void TASImage::DrawDashLine(UInt_t x1,  UInt_t y1, UInt_t x2, UInt_t y2, UInt_t nDash,
                            const char *pDash, const char *col, UInt_t thick)

{
   if (!InitVisual()) {
      Warning("DrawDashLine", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("DrawDashLine", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("DrawDashLine", "Failed to get pixel array");
      return;
   }

   if ((nDash < 2) || !pDash || (nDash%2)) {
      Warning("DrawDashLine", "Wrong input parameters n=%d %ld", nDash, (Long_t)sizeof(pDash)-1);
      return;
   }

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (x1 == x2) {
      DrawDashVLine(x1, y1, y2, nDash, pDash, (UInt_t)color, thick);
   } else if (y1 == y2) {
      DrawDashHLine(y1, x1, x2, nDash, pDash, (UInt_t)color, thick);
   } else {
      if (thick < 2) DrawDashZLine(x1, y1, x2, y2, nDash, pDash, (UInt_t)color);
      else DrawDashZTLine(x1, y1, x2, y2, nDash, pDash, (UInt_t)color, thick);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a polyline.

void TASImage::DrawPolyLine(UInt_t nn, TPoint *xy, const char *col, UInt_t thick,
                            TImage::ECoordMode mode)
{
   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   Int_t x0 = xy[0].GetX();
   Int_t y0 = xy[0].GetY();
   Int_t x = 0;
   Int_t y = 0;

   for (UInt_t i = 1; i < nn; i++) {
      x = (mode == kCoordModePrevious) ? x + xy[i].GetX() : xy[i].GetX();
      y = (mode == kCoordModePrevious) ? y + xy[i].GetY() : xy[i].GetY();

      DrawLineInternal(x0, y0, x, y, (UInt_t)color, thick);

      x0 = x;
      y0 = y;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a point at the specified position.

void TASImage::PutPixel(Int_t x, Int_t y, const char *col)
{
   if (!InitVisual()) {
      Warning("PutPixel", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("PutPixel", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("PutPixel", "Failed to get pixel array");
      return;
   }

   ARGB32 color;
   parse_argb_color(col, &color);

   if ((x < 0) || (y < 0) || (x >= (int)fImage->width) || (y >= (int)fImage->height)) {
      Warning("PutPixel", "Out of range width=%d x=%d, height=%d y=%d",
               fImage->width, x, fImage->height, y);
      return;
   }
   _alphaBlend(&fImage->alt.argb32[Idx(y*fImage->width + x)], &color);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a poly point.

void TASImage::PolyPoint(UInt_t npt, TPoint *ppt, const char *col, TImage::ECoordMode mode)
{
   if (!InitVisual()) {
      Warning("PolyPoint", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("PolyPoint", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("PolyPoint", "Failed to get pixel array");
      return;
   }

   if (!npt || !ppt) {
      Warning("PolyPoint", "No points specified");
      return;
   }

   TPoint *ipt = 0;
   UInt_t i = 0;
   ARGB32 color;
   parse_argb_color(col, &color);

   //make pointlist origin relative
   if (mode == kCoordModePrevious) {
      ipt = new TPoint[npt];

      for (i = 0; i < npt; i++) {
         ipt[i].fX += ppt[i].fX;
         ipt[i].fY += ppt[i].fY;
      }
   }
   int x, y;

   for (i = 0; i < npt; i++) {
      x = ipt ? ipt[i].fX : ppt[i].fX;
      y = ipt ? ipt[i].fY : ppt[i].fY;

      if ((x < 0) || (y < 0) || (x >= (int)fImage->width) || (y >= (int)fImage->height)) {
         continue;
      }
      _alphaBlend(&fImage->alt.argb32[Idx(y*fImage->width + x)], &color);
   }

   if (ipt) {
      delete [] ipt;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw segments.

void TASImage::DrawSegments(UInt_t nseg, Segment_t *seg, const char *col, UInt_t thick)
{
   if (!nseg || !seg) {
      Warning("DrawSegments", "Invalid data nseg=%d seg=0x%lx", nseg, (Long_t)seg);
      return;
   }

   TPoint pt[2];

   for (UInt_t i = 0; i < nseg; i++) {
      pt[0].fX = seg->fX1;
      pt[1].fX = seg->fX2;
      pt[0].fY = seg->fY1;
      pt[1].fY = seg->fY2;

      DrawPolyLine(2, pt, col, thick, kCoordModeOrigin);
      seg++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill spans with specified color or/and stipple.

void TASImage::FillSpans(UInt_t npt, TPoint *ppt, UInt_t *widths, const char *col,
                         const char *stipple, UInt_t w, UInt_t h)
{
   if (!InitVisual()) {
      Warning("FillSpans", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("FillSpans", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("FillSpans", "Failed to get pixel array");
      return;
   }

   if (!npt || !ppt || !widths || (stipple && (!w || !h))) {
      Warning("FillSpans", "Invalid input data npt=%d ppt=0x%lx col=%s widths=0x%lx stipple=0x%lx w=%d h=%d",
              npt, (Long_t)ppt, col, (Long_t)widths, (Long_t)stipple, w, h);
      return;
   }

   ARGB32 color;
   parse_argb_color(col, &color);
   Int_t idx = 0;
   UInt_t x = 0;
   UInt_t yy;

   for (UInt_t i = 0; i < npt; i++) {
      yy = ppt[i].fY*fImage->width;
      for (UInt_t j = 0; j < widths[i]; j++) {
         if ((ppt[i].fX >= (Int_t)fImage->width) || (ppt[i].fX < 0) ||
             (ppt[i].fY >= (Int_t)fImage->height) || (ppt[i].fY < 0)) continue;

         x = ppt[i].fX + j;
         idx = Idx(yy + x);

         if (!stipple) {
            _alphaBlend(&fImage->alt.argb32[idx], &color);
         } else {
            Int_t ii = (ppt[i].fY%h)*w + x%w;

            if (stipple[ii >> 3] & (1 << (ii%8))) {
               _alphaBlend(&fImage->alt.argb32[idx], &color);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill spans with tile image.

void TASImage::FillSpans(UInt_t npt, TPoint *ppt, UInt_t *widths, TImage *tile)
{
   if (!InitVisual()) {
      Warning("FillSpans", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("FillSpans", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("FillSpans", "Failed to get pixel array");
      return;
   }

   if (!npt || !ppt || !widths || !tile) {
      Warning("FillSpans", "Invalid input data npt=%d ppt=0x%lx widths=0x%lx tile=0x%lx",
              npt, (Long_t)ppt, (Long_t)widths, (Long_t)tile);
      return;
   }

   Int_t idx = 0;
   Int_t ii = 0;
   UInt_t x = 0;
   UInt_t *arr = tile->GetArgbArray();
   if (!arr) return;
   UInt_t xx = 0;
   UInt_t yy = 0;
   UInt_t yyy = 0;

   for (UInt_t i = 0; i < npt; i++) {
      yyy = ppt[i].fY*fImage->width;

      for (UInt_t j = 0; j < widths[i]; j++) {
         if ((ppt[i].fX >= (Int_t)fImage->width) || (ppt[i].fX < 0) ||
             (ppt[i].fY >= (Int_t)fImage->height) || (ppt[i].fY < 0)) continue;
         x = ppt[i].fX + j;
         idx = Idx(yyy + x);
         xx = x%tile->GetWidth();
         yy = ppt[i].fY%tile->GetHeight();
         ii = yy*tile->GetWidth() + xx;
         _alphaBlend(&fImage->alt.argb32[idx], &arr[ii]);
      }
      yyy += fImage->width;;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Crop spans.

void TASImage::CropSpans(UInt_t npt, TPoint *ppt, UInt_t *widths)
{
   if (!InitVisual()) {
      Warning("CropSpans", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("CropSpans", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("CropSpans", "Failed to get pixel array");
      return;
   }

   if (!npt || !ppt || !widths) {
      Warning("CropSpans", "No points specified npt=%d ppt=0x%lx widths=0x%lx", npt, (Long_t)ppt, (Long_t)widths);
      return;
   }

   int y0 = ppt[0].fY;
   int y1 = ppt[npt-1].fY;
   UInt_t y = 0;
   UInt_t x = 0;
   UInt_t i = 0;
   UInt_t idx = 0;
   UInt_t sz = fImage->width*fImage->height;
   UInt_t yy = y*fImage->width;

   for (y = 0; (int)y < y0; y++) {
      for (x = 0; x < fImage->width; x++) {
         idx = Idx(yy + x);
         if (idx < sz) fImage->alt.argb32[idx] = 0;
      }
      yy += fImage->width;
   }

   for (i = 0; i < npt; i++) {
      for (x = 0; (int)x < ppt[i].fX; x++) {
         idx = Idx(ppt[i].fY*fImage->width + x);
         if (idx < sz) fImage->alt.argb32[idx] = 0;
      }
      for (x = ppt[i].fX + widths[i] + 1; x < fImage->width; x++) {
         idx = Idx(ppt[i].fY*fImage->width + x);
         if (idx < sz) fImage->alt.argb32[idx] = 0;
      }
   }

   yy = y1*fImage->width;
   for (y = y1; y < fImage->height; y++) {
      for (x = 0; x < fImage->width; x++) {
         idx = Idx(yy + x);
         if (idx < sz) fImage->alt.argb32[idx] = 0;
      }
      yy += fImage->width;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy source region to the destination image. Copy is done according
/// to specified function:
/// ~~~ {.cpp}
/// enum EGraphicsFunction {
///    kGXclear = 0,               // 0
///    kGXand,                     // src AND dst
///    kGXandReverse,              // src AND NOT dst
///    kGXcopy,                    // src (default)
///    kGXandInverted,             // NOT src AND dst
///    kGXnoop,                    // dst
///    kGXxor,                     // src XOR dst
///    kGXor,                      // src OR dst
///    kGXnor,                     // NOT src AND NOT dst
///    kGXequiv,                   // NOT src XOR dst
///    kGXinvert,                  // NOT dst
///    kGXorReverse,               // src OR NOT dst
///    kGXcopyInverted,            // NOT src
///    kGXorInverted,              // NOT src OR dst
///    kGXnand,                    // NOT src OR NOT dst
///    kGXset                      // 1
/// };
/// ~~~

void TASImage::CopyArea(TImage *dst, Int_t xsrc, Int_t ysrc, UInt_t w,  UInt_t h,
                        Int_t xdst, Int_t ydst, Int_t gfunc, EColorChan)
{
   if (!InitVisual()) {
      Warning("CopyArea", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("CopyArea", "no image");
      return;
   }
   if (!dst) return;

   ASImage *out = ((TASImage*)dst)->GetImage();

   int x = 0;
   int y = 0;
   int idx = 0;
   int idx2 = 0;
   xsrc = xsrc < 0 ? 0 : xsrc;
   ysrc = ysrc < 0 ? 0 : ysrc;

   if ((xsrc >= (int)fImage->width) || (ysrc >= (int)fImage->height)) return;

   w = xsrc + w > fImage->width ? fImage->width - xsrc : w;
   h = ysrc + h > fImage->height ? fImage->height - ysrc : h;
   UInt_t yy = (ysrc + y)*fImage->width;

   if (!fImage->alt.argb32) {
      BeginPaint();
   }
   if (!out->alt.argb32) {
      dst->BeginPaint();
      out = ((TASImage*)dst)->GetImage();
   }

   if (fImage->alt.argb32 && out->alt.argb32) {
      for (y = 0; y < (int)h; y++) {
         for (x = 0; x < (int)w; x++) {
            idx = Idx(yy + x + xsrc);
            if ((x + xdst < 0) || (ydst + y < 0) ||
                (x + xdst >= (int)out->width) || (y + ydst >= (int)out->height) ) continue;

            idx2 = Idx((ydst + y)*out->width + x + xdst);

            switch ((EGraphicsFunction)gfunc) {
               case kGXclear:
                  out->alt.argb32[idx2] = 0;
                  break;
               case kGXand:
                  out->alt.argb32[idx2] &= fImage->alt.argb32[idx];
                  break;
               case kGXandReverse:
                  out->alt.argb32[idx2] = fImage->alt.argb32[idx] & (~out->alt.argb32[idx2]);
                  break;
               case kGXandInverted:
                  out->alt.argb32[idx2] &= ~fImage->alt.argb32[idx];
                  break;
               case kGXnoop:
                  break;
               case kGXxor:
                  out->alt.argb32[idx2] ^= fImage->alt.argb32[idx];
                  break;
               case kGXor:
                  out->alt.argb32[idx2] |= fImage->alt.argb32[idx];
                  break;
               case kGXnor:
                  out->alt.argb32[idx2] = (~fImage->alt.argb32[idx]) & (~out->alt.argb32[idx2]);
                  break;
               case kGXequiv:
                  out->alt.argb32[idx2] ^= ~fImage->alt.argb32[idx];
                  break;
               case kGXinvert:
                  out->alt.argb32[idx2] = ~out->alt.argb32[idx2];
                  break;
               case kGXorReverse:
                  out->alt.argb32[idx2] = fImage->alt.argb32[idx] | (~out->alt.argb32[idx2]);
                  break;
               case kGXcopyInverted:
                  out->alt.argb32[idx2] = ~fImage->alt.argb32[idx];
                  break;
               case kGXorInverted:
                  out->alt.argb32[idx2] |= ~fImage->alt.argb32[idx];
                  break;
               case kGXnand:
                  out->alt.argb32[idx2] = (~fImage->alt.argb32[idx]) | (~out->alt.argb32[idx2]);
                  break;
               case kGXset:
                  out->alt.argb32[idx2] = 0xFFFFFFFF;
                  break;
               case kGXcopy:
               default:
                  out->alt.argb32[idx2] = fImage->alt.argb32[idx];
                  break;
            }
         }
         yy += fImage->width;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a cell array.
///
/// \param[in] x1,y1        : left down corner
/// \param[in] x2,y2        : right up corner
/// \param[in] nx,ny        : array size
/// \param[in] ic           : array of ARGB32 colors
///
/// Draw a cell array. The drawing is done with the pixel precision
/// if (X2-X1)/NX (or Y) is not a exact pixel number the position of
/// the top right corner may be wrong.

void TASImage::DrawCellArray(Int_t x1, Int_t y1, Int_t x2, Int_t y2, Int_t nx,
                             Int_t ny, UInt_t *ic)
{
   int i, j, ix, iy, w, h;

   ARGB32 color = 0xFFFFFFFF;
   ARGB32 icol;

   w  = TMath::Max((x2-x1)/(nx),1);
   h  = TMath::Max((y1-y2)/(ny),1);
   ix = x1;

   for (i = 0; i < nx; i++) {
      iy = y1 - h;
      for (j = 0; j < ny; j++) {
         icol = (ARGB32)ic[i + (nx*j)];
         if (icol != color) {
            color = icol;
         }
         FillRectangleInternal((UInt_t)color, ix, iy, w, h);
         iy = iy - h;
      }
      ix = ix + w;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return alpha-blended value computed from bottom and top pixel values.

UInt_t TASImage::AlphaBlend(UInt_t bot, UInt_t top)
{
   UInt_t ret = bot;

   _alphaBlend(&ret, &top);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Return visual.

const ASVisual *TASImage::GetVisual()
{
   return fgVisual;
}

////////////////////////////////////////////////////////////////////////////////
/// Get poly bounds along Y.

static int GetPolyYBounds(TPoint *pts, int n, int *by, int *ty)
{
   TPoint *ptMin;
   int ymin, ymax;
   TPoint *ptsStart = pts;

   ptMin = pts;
   ymin = ymax = (pts++)->fY;

   while (--n > 0) {
      if (pts->fY < ymin) {
         ptMin = pts;
         ymin = pts->fY;
      }
      if (pts->fY > ymax) {
         ymax = pts->fY;
      }
      pts++;
    }

    *by = ymin;
    *ty = ymax;
    return (ptMin - ptsStart);
}

////////////////////////////////////////////////////////////////////////////////
/// The code is based on Xserver/mi/mipolycon.c
///    "Copyright 1987, 1998  The Open Group"

Bool_t TASImage::GetPolygonSpans(UInt_t npt, TPoint *ppt, UInt_t *nspans,
                                 TPoint **outPoint, UInt_t **outWidth)
{
   int xl = 0;                   // x vals of leftedges
   int xr = 0;                   // x vals of right edges
   int dl = 0;                   // decision variables
   int dr = 0;                   // decision variables
   int ml = 0;                   // left edge slope
   int m1l = 0;                  // left edge slope+1
   int mr = 0, m1r = 0;          // right edge slope and slope+1
   int incr1l = 0, incr2l = 0;   // left edge error increments
   int incr1r = 0, incr2r = 0;   // right edge error increments
   int dy;                       // delta y
   int y;                        // current scanline
   int left, right;              // indices to first endpoints
   int i;                        // loop counter
   int nextleft, nextright;      // indices to second endpoints
   TPoint *ptsOut;               // output buffer
   UInt_t *width;                // output buffer
   TPoint *firstPoint=0;
   UInt_t *firstWidth=0;
   int imin;                     // index of smallest vertex (in y)
   int ymin;                     // y-extents of polygon
   int ymax;
   Bool_t  ret = kTRUE;

   *nspans = 0;

   if (!InitVisual()) {
      Warning("GetPolygonSpans", "Visual not initiated");
      return kFALSE;
   }

   if (!fImage) {
      Warning("GetPolygonSpans", "no image");
      return kFALSE;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("GetPolygonSpans", "Failed to get pixel array");
      return kFALSE;
   }

   if ((npt < 3) || !ppt) {
      Warning("GetPolygonSpans", "No points specified npt=%d ppt=0x%lx", npt, (Long_t)ppt);
      return kFALSE;
   }

   //  find leftx, bottomy, rightx, topy, and the index
   //  of bottomy. Also translate the points.
   imin = GetPolyYBounds(ppt, npt, &ymin, &ymax);

   dy = ymax - ymin + 1;
   if ((npt < 3) || (dy < 0)) return kFALSE;

   ptsOut = firstPoint = new TPoint[dy];
   width = firstWidth = new UInt_t[dy];
   ret = kTRUE;

   nextleft = nextright = imin;
   y = ppt[nextleft].fY;

   //  loop through all edges of the polygon
   do {
      // add a left edge if we need to
      if (ppt[nextleft].fY == y) {
         left = nextleft;

         //  find the next edge, considering the end
         //  conditions of the array.
         nextleft++;
         if (nextleft >= (int)npt) {
            nextleft = 0;
         }

         // now compute all of the random information
         // needed to run the iterative algorithm.
         BRESINITPGON(ppt[nextleft].fY - ppt[left].fY,
                      ppt[left].fX, ppt[nextleft].fX,
                      xl, dl, ml, m1l, incr1l, incr2l);
      }

      // add a right edge if we need to
      if (ppt[nextright].fY == y) {
         right = nextright;

         // find the next edge, considering the end
         // conditions of the array.
         nextright--;
         if (nextright < 0) {
            nextright = npt-1;
         }

         //  now compute all of the random information
         //  needed to run the iterative algorithm.
         BRESINITPGON(ppt[nextright].fY - ppt[right].fY,
                      ppt[right].fX, ppt[nextright].fX,
                      xr, dr, mr, m1r, incr1r, incr2r);
      }

      // generate scans to fill while we still have
      //  a right edge as well as a left edge.
      i = TMath::Min(ppt[nextleft].fY, ppt[nextright].fY) - y;

      // in case of non-convex polygon
      if (i < 0) {
         delete [] firstWidth;
         delete [] firstPoint;
         return kTRUE;
      }

      while (i-- > 0)  {
         ptsOut->fY = y;

         // reverse the edges if necessary
         if (xl < xr) {
            *(width++) = xr - xl;
            (ptsOut++)->fX = xl;
         } else {
            *(width++) = xl - xr;
            (ptsOut++)->fX = xr;
         }
         y++;

         // increment down the edges
         BRESINCRPGON(dl, xl, ml, m1l, incr1l, incr2l);
         BRESINCRPGON(dr, xr, mr, m1r, incr1r, incr2r);
      }
   }  while (y != ymax);

   *nspans = UInt_t(ptsOut - firstPoint);
   *outPoint = firstPoint;
   *outWidth = firstWidth;

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a convex polygon with background color or bitmap.
/// For non convex polygon one must use DrawFillArea method

void TASImage::FillPolygon(UInt_t npt, TPoint *ppt, const char *col,
                           const char *stipple, UInt_t w, UInt_t h)
{
   UInt_t  nspans = 0;
   TPoint *firstPoint = 0;   // output buffer
   UInt_t *firstWidth = 0;   // output buffer

   Bool_t del = GetPolygonSpans(npt, ppt, &nspans, &firstPoint, &firstWidth);
   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (nspans) {
      if (!stipple && ((color & 0xff000000)==0xff000000)) { //no stipple no alpha
         FillSpansInternal(nspans, firstPoint, firstWidth, color);
      } else {
         FillSpans(nspans, firstPoint, firstWidth, col, stipple, w, h);
      }

      if (del) {
         delete [] firstWidth;
         delete [] firstPoint;
      }
   } else {
      if (firstWidth) delete [] firstWidth;
      if (firstPoint) delete [] firstPoint;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a convex polygon with background image.
/// For non convex polygon one must use DrawFillArea method

void TASImage::FillPolygon(UInt_t npt, TPoint *ppt, TImage *tile)
{
   UInt_t  nspans = 0;
   TPoint *firstPoint = 0;   // output buffer
   UInt_t *firstWidth = 0;   // output buffer

   Bool_t del = GetPolygonSpans(npt, ppt, &nspans, &firstPoint, &firstWidth);

   if (nspans) {
      FillSpans(nspans, firstPoint, firstWidth, tile);

      if (del) {
         delete [] firstWidth;
         delete [] firstPoint;
      }
   } else {
      if (firstWidth) delete [] firstWidth;
      if (firstPoint) delete [] firstPoint;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Crop a convex polygon.

void TASImage::CropPolygon(UInt_t npt, TPoint *ppt)
{
   UInt_t  nspans = 0;
   TPoint *firstPoint = 0;
   UInt_t *firstWidth = 0;

   Bool_t del = GetPolygonSpans(npt, ppt, &nspans, &firstPoint, &firstWidth);

   if (nspans) {
      CropSpans(nspans, firstPoint, firstWidth);

      if (del) {
         delete [] firstWidth;
         delete [] firstPoint;
      }
   } else {
         if (firstWidth) delete [] firstWidth;
         if (firstPoint) delete [] firstPoint;
   }
}

static const UInt_t NUMPTSTOBUFFER = 512;

////////////////////////////////////////////////////////////////////////////////
/// Fill a polygon (any type convex, non-convex).

void TASImage::DrawFillArea(UInt_t count, TPoint *ptsIn, const char *col,
                           const char *stipple, UInt_t w, UInt_t h)
{
   if (!InitVisual()) {
      Warning("DrawFillArea", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("DrawFillArea", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("DrawFillArea", "Failed to get pixel array");
      return;
   }

   if ((count < 3) || !ptsIn) {
      Warning("DrawFillArea", "No points specified npt=%d ppt=0x%lx", count, (Long_t)ptsIn);
      return;
   }

   if (count < 5) {
      FillPolygon(count, ptsIn, col, stipple, w, h);
      return;
   }

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   EdgeTableEntry *pAET;  // the Active Edge Table
   int y;                 // the current scanline
   UInt_t nPts = 0;          // number of pts in buffer

   ScanLineList *pSLL;   // Current ScanLineList
   TPoint *ptsOut;       // ptr to output buffers
   UInt_t *width;
   TPoint firstPoint[NUMPTSTOBUFFER];  // the output buffers
   UInt_t firstWidth[NUMPTSTOBUFFER];
   EdgeTableEntry *pPrevAET;       // previous AET entry
   EdgeTable ET;                   // Edge Table header node
   EdgeTableEntry AET;             // Active ET header node
   EdgeTableEntry *pETEs;          // Edge Table Entries buff
   ScanLineListBlock SLLBlock;     // header for ScanLineList
   Bool_t del = kTRUE;

   static const UInt_t gEdgeTableEntryCacheSize = 200;
   static EdgeTableEntry gEdgeTableEntryCache[gEdgeTableEntryCacheSize];

   if (count < gEdgeTableEntryCacheSize) {
      pETEs = (EdgeTableEntry*)&gEdgeTableEntryCache;
      del = kFALSE;
   } else {
      pETEs = new EdgeTableEntry[count];
      del = kTRUE;
   }

   ptsOut = firstPoint;
   width = firstWidth;
   CreateETandAET(count, ptsIn, &ET, &AET, pETEs, &SLLBlock);
   pSLL = ET.scanlines.next;

   for (y = ET.ymin; y < ET.ymax; y++) {
      if (pSLL && y == pSLL->scanline) {
         loadAET(&AET, pSLL->edgelist);
         pSLL = pSLL->next;
      }
      pPrevAET = &AET;
      pAET = AET.next;

      while (pAET) {
         ptsOut->fX = pAET->bres.minor_axis;
         ptsOut->fY = y;
         ptsOut++;
         nPts++;

         *width++ = pAET->next->bres.minor_axis - pAET->bres.minor_axis;

         if (nPts == NUMPTSTOBUFFER) {
            if (!stipple && ((color & 0xff000000)==0xff000000)) { //no stipple, no alpha
               FillSpansInternal(nPts, firstPoint, firstWidth, color);
            } else {
               FillSpans(nPts, firstPoint, firstWidth, col, stipple, w, h);
            }
            ptsOut = firstPoint;
            width = firstWidth;
            nPts = 0;
         }
         EVALUATEEDGEEVENODD(pAET, pPrevAET, y)
         EVALUATEEDGEEVENODD(pAET, pPrevAET, y)
      }
      InsertionSort(&AET);
   }

   if (nPts) {
      if (!stipple && ((color & 0xff000000)==0xff000000)) {  //no stipple, no alpha
         FillSpansInternal(nPts, firstPoint, firstWidth, color);
      } else {
         FillSpans(nPts, firstPoint, firstWidth, col, stipple, w, h);
      }
   }

   if (del) delete [] pETEs;
   FreeStorage(SLLBlock.next);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a polygon (any type convex, non-convex).

void TASImage::DrawFillArea(UInt_t count, TPoint *ptsIn, TImage *tile)
{
   if (!InitVisual()) {
      Warning("DrawFillArea", "Visual not initiated");
      return;
   }

   if (!fImage) {
      Warning("DrawFillArea", "no image");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!fImage->alt.argb32) {
      Warning("DrawFillArea", "Failed to get pixel array");
      return;
   }

   if ((count < 3) || !ptsIn) {
      Warning("DrawFillArea", "No points specified npt=%d ppt=0x%lx", count, (Long_t)ptsIn);
      return;
   }

   if (count < 5) {
      FillPolygon(count, ptsIn, tile);
      return;
   }

   EdgeTableEntry *pAET;   // the Active Edge Table
   int y;                  // the current scanline
   UInt_t nPts = 0;       // number of pts in buffer

   ScanLineList *pSLL;    // Current ScanLineList
   TPoint *ptsOut;        // ptr to output buffers
   UInt_t *width;
   TPoint firstPoint[NUMPTSTOBUFFER]; // the output buffers
   UInt_t firstWidth[NUMPTSTOBUFFER];
   EdgeTableEntry *pPrevAET;       // previous AET entry
   EdgeTable ET;                   // Edge Table header node
   EdgeTableEntry AET;             // Active ET header node
   EdgeTableEntry *pETEs;          // Edge Table Entries buff
   ScanLineListBlock SLLBlock;     // header for ScanLineList

   pETEs = new EdgeTableEntry[count];

   ptsOut = firstPoint;
   width = firstWidth;
   CreateETandAET(count, ptsIn, &ET, &AET, pETEs, &SLLBlock);
   pSLL = ET.scanlines.next;

   for (y = ET.ymin; y < ET.ymax; y++) {
      if (pSLL && y == pSLL->scanline) {
         loadAET(&AET, pSLL->edgelist);
         pSLL = pSLL->next;
      }
      pPrevAET = &AET;
      pAET = AET.next;

      while (pAET) {
         ptsOut->fX = pAET->bres.minor_axis;
         ptsOut->fY = y;
         ptsOut++;
         nPts++;

         *width++ = pAET->next->bres.minor_axis - pAET->bres.minor_axis;

         if (nPts == NUMPTSTOBUFFER) {
            FillSpans(nPts, firstPoint, firstWidth, tile);
            ptsOut = firstPoint;
            width = firstWidth;
            nPts = 0;
         }
         EVALUATEEDGEEVENODD(pAET, pPrevAET, y)
         EVALUATEEDGEEVENODD(pAET, pPrevAET, y)
      }
      InsertionSort(&AET);
   }
   FillSpans(nPts, firstPoint, firstWidth, tile);

   delete [] pETEs;
   FreeStorage(SLLBlock.next);
}

////////////////////////////////////////////////////////////////////////////////
/// Create draw context.

static ASDrawContext *create_draw_context_argb32(ASImage *im, ASDrawTool *brush)
{
   ASDrawContext *ctx = new ASDrawContext;

   ctx->canvas_width = im->width;
   ctx->canvas_height = im->height;
   ctx->canvas = im->alt.argb32;
   ctx->scratch_canvas = 0;

   ctx->flags = ASDrawCTX_CanvasIsARGB;
   asim_set_custom_brush_colored( ctx, brush);
   return ctx;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy asdraw context32.

static void destroy_asdraw_context32( ASDrawContext *ctx )
{
   if (ctx) {
      if (ctx->scratch_canvas) free(ctx->scratch_canvas);
      delete ctx;
   }
}

static const UInt_t kBrushCacheSize = 20;
static CARD32 gBrushCache[kBrushCacheSize*kBrushCacheSize];

////////////////////////////////////////////////////////////////////////////////
/// Draw wide line.

void TASImage::DrawWideLine(UInt_t x1, UInt_t y1, UInt_t x2, UInt_t y2,
                            UInt_t color, UInt_t thick)
{
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = thick < kBrushCacheSize;

   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   };

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.width = thick;
   brush.height = thick;
   brush.center_y = brush.center_x = thick/2;

   // When the first or last point of a wide line is exactly on the
   // window limit the line is drawn vertically or horizontally.
   // see https://sft.its.cern.ch/jira/browse/ROOT-8021
   UInt_t xx1 = x1;
   UInt_t yy1 = y1;
   UInt_t xx2 = x2;
   UInt_t yy2 = y2;
   if (xx1 == fImage->width)  --xx1;
   if (yy1 == fImage->height) --yy1;
   if (xx2 == fImage->width)  --xx2;
   if (yy2 == fImage->height) --yy2;
   ASDrawContext *ctx = create_draw_context_argb32(fImage, &brush);
   asim_move_to(ctx, xx1, yy1);
   asim_line_to(ctx, xx2, yy2);

   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw glyph bitmap.

void TASImage::DrawGlyph(void *bitmap, UInt_t color, Int_t bx, Int_t by)
{
   static UInt_t col[5];
   Int_t x, y, yy, y0, xx;
   Bool_t has_alpha = (color & 0xff000000) != 0xff000000;

   ULong_t r, g, b;
   int idx = 0;
   FT_Bitmap *source = (FT_Bitmap*)bitmap;
   UChar_t d = 0, *s = source->buffer;

   Int_t dots = Int_t(source->width * source->rows);
   r = g = b = 0;
   Int_t bxx, byy;

   yy = y0 = by > 0 ? by * fImage->width : 0;
   for (y = 0; y < (int) source->rows; y++) {
      byy = by + y;
      if ((byy >= (int)fImage->height) || (byy <0)) continue;

      for (x = 0; x < (int) source->width; x++) {
         bxx = bx + x;
         if ((bxx >= (int)fImage->width) || (bxx < 0)) continue;

         idx = Idx(bxx + yy);
         r += ((fImage->alt.argb32[idx] & 0xff0000) >> 16);
         g += ((fImage->alt.argb32[idx] & 0x00ff00) >> 8);
         b += (fImage->alt.argb32[idx] & 0x0000ff);
      }
      yy += fImage->width;
   }
   if (dots != 0) {
      r /= dots;
      g /= dots;
      b /= dots;
   }

   col[0] = (r << 16) + (g << 8) + b;
   col[4] = color;
   Int_t col4r = (col[4] & 0xff0000) >> 16;
   Int_t col4g = (col[4] & 0x00ff00) >> 8;
   Int_t col4b = (col[4] & 0x0000ff);

   // interpolate between fore and background colors
   for (x = 3; x > 0; x--) {
      xx = 4-x;
      Int_t colxr = (col4r*x + r*xx) >> 2;
      Int_t colxg = (col4g*x + g*xx) >> 2;
      Int_t colxb = (col4b*x + b*xx) >> 2;
      col[x] = (colxr << 16) + (colxg << 8) + colxb;
   }

   yy = y0;
   ARGB32 acolor;

   Int_t clipx1=0, clipx2=0, clipy1=0, clipy2=0;
   Bool_t noClip = kTRUE;

   if (gPad) {
      Float_t is = gStyle->GetImageScaling();
      clipx1 = gPad->XtoAbsPixel(gPad->GetX1())*is;
      clipx2 = gPad->XtoAbsPixel(gPad->GetX2())*is;
      clipy1 = gPad->YtoAbsPixel(gPad->GetY1())*is;
      clipy2 = gPad->YtoAbsPixel(gPad->GetY2())*is;
      noClip = kFALSE;
   }

   for (y = 0; y < (int) source->rows; y++) {
      byy = by + y;

      for (x = 0; x < (int) source->width; x++) {
         bxx = bx + x;

         d = *s++ & 0xff;
         d = ((d + 10) * 5) >> 8;
         if (d > 4) d = 4;

         if (d) {
            if ( noClip || ((x < (int) source->width) &&
                 (bxx <  (int)clipx2) && (bxx >= (int)clipx1) &&
                 (byy >= (int)clipy2) && (byy <  (int)clipy1) )) {
               idx    = Idx(bxx + yy);
               acolor = (ARGB32)col[d];
               if (has_alpha) {
                  _alphaBlend(&fImage->alt.argb32[idx], &acolor);
               } else {
                  fImage->alt.argb32[idx] = acolor;
               }
            }
         }
      }
      yy += fImage->width;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text at the pixel position (x,y).

void TASImage::DrawText(TText *text, Int_t x, Int_t y)
{
   if (!text)   return;
   if (!fImage) return;
   if (!gPad)   return;

   if (!InitVisual()) {
      Warning("DrawText", "Visual not initiated");
      return;
   }

   if (!fImage->alt.argb32) {
      BeginPaint();
   }

   if (!TTF::IsInitialized()) TTF::Init();

   // set text font
   TTF::SetTextFont(text->GetTextFont());

   Int_t wh = gPad->XtoPixel(gPad->GetX2());
   Int_t hh = gPad->YtoPixel(gPad->GetY1());

   // set text size
   Float_t ttfsize;
   if (wh < hh) {
      ttfsize = text->GetTextSize()*wh;
   } else {
      ttfsize = text->GetTextSize()*hh;
   }
   TTF::SetTextSize(ttfsize*kScale);

   // set text angle
   TTF::SetRotationMatrix(text->GetTextAngle());

   // set text
   const wchar_t *wcsTitle = reinterpret_cast<const wchar_t *>(text->GetWcsTitle());
   if (wcsTitle != NULL) {
      TTF::PrepareString(wcsTitle);
   } else {
      TTF::PrepareString(text->GetTitle());
   }
   TTF::LayoutGlyphs();

   // color
   TColor *col = gROOT->GetColor(text->GetTextColor());
   if (!col) { // no color, make it black
      col = gROOT->GetColor(1);
      if (!col) return;
   }
   ARGB32 color = ARGB32_White;
   parse_argb_color(col->AsHexString(), &color);

   // Align()
   Int_t align = 0;
   Int_t txalh = text->GetTextAlign()/10;
   Int_t txalv = text->GetTextAlign()%10;

   switch (txalh) {
      case 0 :
      case 1 :
         switch (txalv) {  //left
            case 1 :
               align = 7;   //bottom
               break;
            case 2 :
               align = 4;   //center
               break;
            case 3 :
               align = 1;   //top
               break;
         }
         break;
      case 2 :
         switch (txalv) { //center
            case 1 :
               align = 8;   //bottom
               break;
            case 2 :
               align = 5;   //center
               break;
            case 3 :
               align = 2;   //top
               break;
         }
         break;
      case 3 :
         switch (txalv) {  //right
            case 1 :
               align = 9;   //bottom
               break;
            case 2 :
               align = 6;   //center
               break;
            case 3 :
               align = 3;   //top
               break;
         }
         break;
   }

   FT_Vector ftal;

   // vertical alignment
   if (align == 1 || align == 2 || align == 3) {
      ftal.y = TTF::GetAscent();
   } else if (align == 4 || align == 5 || align == 6) {
      ftal.y = TTF::GetAscent()/2;
   } else {
      ftal.y = 0;
   }

   // horizontal alignment
   if (align == 3 || align == 6 || align == 9) {
      ftal.x = TTF::GetWidth();
   } else if (align == 2 || align == 5 || align == 8) {
      ftal.x = TTF::GetWidth()/2;
   } else {
      ftal.x = 0;
   }

   FT_Vector_Transform(&ftal, TTF::GetRotMatrix());
   ftal.x = (ftal.x >> 6);
   ftal.y = (ftal.y >> 6);

   TTF::TTGlyph *glyph = TTF::GetGlyphs();

   for (int n = 0; n < TTF::GetNumGlyphs(); n++, glyph++) {
      if (FT_Glyph_To_Bitmap(&glyph->fImage, ft_render_mode_normal, 0, 1 )) continue;

      FT_BitmapGlyph bitmap = (FT_BitmapGlyph)glyph->fImage;
      FT_Bitmap *source = &bitmap->bitmap;

      Int_t bx = x - ftal.x + bitmap->left;
      Int_t by = y + ftal.y - bitmap->top;

      DrawGlyph(source, color, bx, by);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text using TrueType fonts.

void TASImage::DrawTextTTF(Int_t x, Int_t y, const char *text, Int_t size,
                           UInt_t color, const char *font_name, Float_t angle)
{
   if (!TTF::IsInitialized()) TTF::Init();

   TTF::SetTextFont(font_name);
   TTF::SetTextSize(size);
   TTF::SetRotationMatrix(angle);
   TTF::PrepareString(text);
   TTF::LayoutGlyphs();

   TTF::TTGlyph *glyph = TTF::GetGlyphs();

   // compute the size and position  that will contain the text
   // Int_t Xoff = 0; if (TTF::GetBox().xMin < 0) Xoff = -TTF::GetBox().xMin;
   Int_t Yoff = 0; if (TTF::GetBox().yMin < 0) Yoff = -TTF::GetBox().yMin;
   Int_t h    = TTF::GetBox().yMax + Yoff;

   for (int n = 0; n < TTF::GetNumGlyphs(); n++, glyph++) {
      if (FT_Glyph_To_Bitmap(&glyph->fImage, ft_render_mode_normal, 0, 1 )) continue;

      FT_BitmapGlyph bitmap = (FT_BitmapGlyph)glyph->fImage;
      FT_Bitmap *source = &bitmap->bitmap;

      Int_t bx = x + bitmap->left;
      Int_t by = y + h - bitmap->top;
      DrawGlyph(source, color, bx, by);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return in-memory buffer compressed according image type.
/// Buffer must be deallocated after usage.
/// This method can be used for sending images over network.

void TASImage::GetImageBuffer(char **buffer, int *size, EImageFileTypes type)
{
   static ASImageExportParams params;
   Bool_t ret = kFALSE;
   int   isize = 0;
   char *ibuff = 0;
   ASImage *img = fScaledImage ? fScaledImage->fImage : fImage;

   if (!img) return;

   switch (type) {
      case TImage::kXpm:
         ret = ASImage2xpmRawBuff(img, (CARD8 **)buffer, size, 0);
         break;
      default:
         ret = ASImage2PNGBuff(img, (CARD8 **)buffer, size, &params);
   }

   if (!ret) {
      *size = isize;
      *buffer = ibuff;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create image from compressed buffer.
/// Supported formats:
///
///  -  PNG - by default
///  -  XPM - two options exist:
///      1.  xpm as a single string (raw buffer). Such string
///          is returned by GetImageBuffer method.
///          For example:
/// ~~~ {.cpp}
///       char *buf;
///       int sz;
///       im1->GetImageBuffer(&buf, &int, TImage::kXpm); /*raw buffer*/
///       TImage *im2 = TImage::Create();
///       im2->SetImageBuffer(&buf, TImage::kXpm);
/// ~~~
///      2.  xpm as an array of strings (pre-parsed)
/// ~~~ {.cpp}
///    For example:
///       char *xpm[] = {
///          "64 28 58 1",
///          "  c #0A030C",
///          ". c #1C171B"
///             ...
///    TImage *im = TImage::Create();
///    im->SetImageBuffer(xpm, TImage::kXpm);
/// ~~~

Bool_t TASImage::SetImageBuffer(char **buffer, EImageFileTypes type)
{
   DestroyImage();

   static ASImageImportParams params;
   params.flags = 0;
   params.width = 0;
   params.height = 0 ;
   params.filter = SCL_DO_ALL;
   params.gamma = 0;
   params.gamma_table = 0;
   params.compression = 0;
   params.format = ASA_ASImage;
   params.search_path = 0;
   params.subimage = 0;

   switch (type) {
      case TImage::kXpm:
      {
         char *ptr = buffer[0];
         while (isspace((int)*ptr)) ++ptr;
         if (atoi(ptr)) {  // pre-parsed and preloaded data
            fImage = xpm_data2ASImage((const char**)buffer, &params);
         } else {
            fImage = xpmRawBuff2ASImage((const char*)*buffer, &params);
         }
         break;
      }
      default:
         fImage = PNGBuff2ASimage((CARD8 *)*buffer, &params);
         break;
   }

   if (!fImage) {
      return kFALSE;
   }

   if (fName.IsNull()) {
      fName.Form("img_%dx%d.%d", fImage->width, fImage->height, gRandom->Integer(1000));
   }
   UnZoom();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create image thumbnail.

void TASImage::CreateThumbnail()
{
   int size;
   const int sz = 64;

   if (!fImage) {
      return;
   }

   if (!InitVisual()) {
      return;
   }

   static char *buf = 0;
   int w, h;
   ASImage *img = 0;

   if (fImage->width > fImage->height) {
      w = sz;
      h = (fImage->height*sz)/fImage->width;
   } else {
      h = sz;
      w = (fImage->width*sz)/fImage->height;
   }

   w = w < 8 ? 8 : w;
   h = h < 8 ? 8 : h;

   img = scale_asimage(fgVisual, fImage, w, h, ASA_ASImage,
                       GetImageCompression(), GetImageQuality());
   if (!img) {
      return;
   }

   // contrasting
   ASImage *rendered_im;
   ASImageLayer layers[2];
   init_image_layers(&(layers[0]), 2);
   layers[0].im = img;
   layers[0].dst_x = 0;
   layers[0].dst_y = 0;
   layers[0].clip_width = img->width;
   layers[0].clip_height = img->height;
   layers[0].bevel = 0;
   layers[1].im = img;
   layers[1].dst_x = 0;
   layers[1].dst_y = 0;
   layers[1].clip_width = img->width;
   layers[1].clip_height = img->height;
   layers[1].merge_scanlines = blend_scanlines_name2func("tint");
   rendered_im = merge_layers(fgVisual, &(layers[0]), 2, img->width, img->height,
                              ASA_ASImage, GetImageCompression(), GetImageQuality());
   destroy_asimage(&img);
   img = rendered_im;

   // pad image
   ASImage *padimg = 0;
   int d = 0;

   if (w == sz) {
      d = (sz - h) >> 1;
      padimg = pad_asimage(fgVisual, img, 0, d, sz, sz, 0x00ffffff,
                           ASA_ASImage, GetImageCompression(), GetImageQuality());
   } else {
      d = (sz - w) >> 1;
      padimg = pad_asimage(fgVisual, img, d, 0, sz, sz, 0x00ffffff,
                           ASA_ASImage, GetImageCompression(), GetImageQuality());
   }

   if (!padimg) {
      destroy_asimage(&img);
      return;
   }

   void *ptr = &buf;
   ASImage2xpmRawBuff(padimg, (CARD8 **)ptr, &size, 0);
   fTitle = buf;

   destroy_asimage(&padimg);
}

////////////////////////////////////////////////////////////////////////////////
/// Streamer for ROOT I/O.

void TASImage::Streamer(TBuffer &b)
{
   Bool_t image_type = 0;
   char *buffer = 0;
   int size = 0;
   int w, h;
   UInt_t R__s, R__c;

   if (b.IsReading()) {
      Version_t version = b.ReadVersion(&R__s, &R__c);
      if (version == 0) { //dumb prototype for schema evolution
         return;
      }

      if ( version == 1 ) {
         Int_t fileVersion = b.GetVersionOwner();
         if (fileVersion > 0 && fileVersion < 50000 ) {
            TImage::Streamer(b);
            b >> fMaxValue;
            b >> fMinValue;
            b >> fZoomOffX;
            b >> fZoomOffY;
            b >> fZoomWidth;
            b >> fZoomHeight;
            if ( fileVersion < 40200 ) {
               Bool_t zoomUpdate;
               b >> zoomUpdate;
               fZoomUpdate = zoomUpdate;
            } else {
               b >> fZoomUpdate;
               b >> fEditable;
               Bool_t paintMode;
               b >> paintMode;
               fPaintMode = paintMode;
            }
            b.CheckByteCount(R__s, R__c, TASImage::IsA());
            return;
         }
      }

      TNamed::Streamer(b);
      b >> image_type;

      if (image_type != 0) {     // read PNG compressed image
         b >> size;
         buffer = new char[size];
         b.ReadFastArray(buffer, size);
         SetImageBuffer(&buffer, TImage::kPng);
         delete [] buffer;
      } else {                   // read vector with palette
         TAttImage::Streamer(b);
         b >> w;
         b >> h;
         size = w*h;
         Double_t *vec = new Double_t[size];
         b.ReadFastArray(vec, size);
         SetImage(vec, w, h, &fPalette);
         delete [] vec;
      }
      b.CheckByteCount(R__s, R__c, TASImage::IsA());
   } else {
      if (!fImage) {
         return;
      }
      R__c = b.WriteVersion(TASImage::IsA(), kTRUE);

      if (fName.IsNull()) {
         fName.Form("img_%dx%d.%d", fImage->width, fImage->height, gRandom->Integer(1000));
      }
      TNamed::Streamer(b);

      image_type = fImage->alt.vector ? 0 : 1;
      b << image_type;

      if (image_type != 0) {     // write PNG compressed image
         GetImageBuffer(&buffer, &size, TImage::kPng);
         b << size;
         b.WriteFastArray(buffer, size);
         delete buffer;
      } else {                   // write vector  with palette
         TAttImage::Streamer(b);
         b << fImage->width;
         b << fImage->height;
         b.WriteFastArray(fImage->alt.vector, fImage->width*fImage->height);
      }
      b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Browse image.

void TASImage::Browse(TBrowser *)
{
   if (fImage->alt.vector) {
      Draw("n");
   } else {
      Draw("nxxx");
   }
   CreateThumbnail();
}

////////////////////////////////////////////////////////////////////////////////
/// Title is used to keep 32x32 xpm image's thumbnail.

const char *TASImage::GetTitle() const
{
   if (!gDirectory || !gDirectory->IsWritable()) {
      return 0;
   }

   TASImage *mutble = (TASImage *)this;

   if (fTitle.IsNull()) {
      mutble->SetTitle(fName.Data());
   }

   return fTitle.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Set a title for an image.

void TASImage::SetTitle(const char *title)
{
   if (fTitle.IsNull()) {
      CreateThumbnail();
   }

   if (fTitle.IsNull()) {
      return;
   }

   int start = fTitle.Index("/*") + 3;
   int stop = fTitle.Index("*/") - 1;

   if ((start > 0) && (stop - start > 0)) {
      fTitle.Replace(start, stop - start, title);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a cubic bezier line.

void TASImage::DrawCubeBezier(Int_t x1, Int_t y1, Int_t x2, Int_t y2,
                             Int_t x3, Int_t y3, const char *col, UInt_t thick)
{
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = thick < kBrushCacheSize;

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   };

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.width = thick;
   brush.height = thick;
   brush.center_y = brush.center_x = thick/2;

   ASDrawContext *ctx = 0;

   ctx = create_draw_context_argb32(fImage, &brush);
   asim_cube_bezier(ctx, x1, y1, x2, y2, x3, y3);

   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a straight ellipse.
/// If thick < 0 - draw filled ellipse.

void TASImage::DrawStraightEllips(Int_t x, Int_t y, Int_t rx, Int_t ry,
                                  const char *col, Int_t thick)
{
   thick = !thick ? 1 : thick;
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = (thick > 0) && ((UInt_t)thick < kBrushCacheSize);

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   };

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.width = thick > 0 ? thick : 1;
   brush.height = thick > 0 ? thick : 1;
   brush.center_y = brush.center_x = thick > 0 ? thick/2 : 0;

   ASDrawContext *ctx = create_draw_context_argb32(fImage, &brush);
   asim_straight_ellips(ctx, x, y, rx, ry, thick < 0);

   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a circle.
/// If thick < 0 - draw filled circle

void TASImage::DrawCircle(Int_t x, Int_t y, Int_t r, const char *col, Int_t thick)
{
   thick = !thick ? 1 : thick;
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = (thick > 0) && ((UInt_t)thick < kBrushCacheSize);

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

///matrix = new CARD32[sz];
   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   }

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.height = brush.width = thick > 0 ? thick : 1;
   brush.center_y = brush.center_x = thick > 0 ? thick/2 : 0;

   ASDrawContext *ctx = create_draw_context_argb32(fImage, &brush);
   asim_circle(ctx, x,  y, r, thick < 0);

///free (matrix);
   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw an ellipse.
/// If thick < 0 - draw filled ellips

void TASImage::DrawEllips(Int_t x, Int_t y, Int_t rx, Int_t ry, Int_t angle,
                           const char *col, Int_t thick)
{
   thick = !thick ? 1 : thick;
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = (thick > 0) && ((UInt_t)thick < kBrushCacheSize);

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   };

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.width = thick > 0 ? thick : 1;
   brush.height = thick > 0 ? thick : 1;
   brush.center_y = brush.center_x = thick > 0 ? thick/2 : 0;

   ASDrawContext *ctx = create_draw_context_argb32(fImage, &brush);
   asim_ellips(ctx, x, y, rx, ry, angle, thick < 0);

   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw an ellipse.
/// If thick < 0 - draw filled ellipse.

void TASImage::DrawEllips2(Int_t x, Int_t y, Int_t rx, Int_t ry, Int_t angle,
                           const char *col, Int_t thick)
{
   thick = !thick ? 1 : thick;
   Int_t sz = thick*thick;
   CARD32 *matrix;
   Bool_t use_cache = (thick > 0) && ((UInt_t)thick < kBrushCacheSize);

   ARGB32 color = ARGB32_White;
   parse_argb_color(col, &color);

   if (use_cache) {
      matrix = gBrushCache;
   } else {
      matrix = new CARD32[sz];
   }

   for (int i = 0; i < sz; i++) {
      matrix[i] = (CARD32)color;
   };

   ASDrawTool brush;
   brush.matrix = matrix;
   brush.width = thick > 0 ? thick : 1;
   brush.height = thick > 0 ? thick : 1;
   brush.center_y = brush.center_x = thick > 0 ? thick/2 : 0;

   ASDrawContext *ctx = create_draw_context_argb32(fImage, &brush);
   asim_ellips2(ctx, x, y, rx, ry, angle, thick < 0);

   if (!use_cache) {
      delete [] matrix;
   }
   destroy_asdraw_context32(ctx);
}

////////////////////////////////////////////////////////////////////////////////
/// Flood fill.

void TASImage::FloodFill(Int_t /*x*/, Int_t /*y*/, const char * /*col*/,
                         const char * /*minc*/, const char * /*maxc*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Convert RGB image to Gray image and vice versa.

void TASImage::Gray(Bool_t on)
{
   if (fIsGray == on) {
      return;
   }

   if (!IsValid()) {
      Warning("Gray", "Image not initiated");
      return;
   }

   if (!InitVisual()) {
      Warning("Gray", "Visual not initiated");
      return;
   }

   if (!fGrayImage && !on) {
      return;
   }
   ASImage *sav = 0;
   delete fScaledImage;
   fScaledImage = 0;

   if (fGrayImage)  {
      sav = fImage;
      fImage = fGrayImage;
      fGrayImage = sav;
      fIsGray = on;
      return;
   }

   if (!on) return;

   UInt_t l, r, g, b, idx;
   int y = 0;
   UInt_t i, j;

   if (fImage->alt.argb32) {
      fGrayImage = tile_asimage(fgVisual, fImage, 0, 0, fImage->width, fImage->height,
                                0, ASA_ARGB32, 0, ASIMAGE_QUALITY_DEFAULT);

      for (i = 0; i < fImage->height; i++) {
         for (j = 0; j < fImage->width; j++) {
            idx = Idx(y + j);

            r = ((fImage->alt.argb32[idx] & 0xff0000) >> 16);
            g = ((fImage->alt.argb32[idx] & 0x00ff00) >> 8);
            b = (fImage->alt.argb32[idx] & 0x0000ff);
            l = (57*r + 181*g + 18*b)/256;
            fGrayImage->alt.argb32[idx] = (l << 16) + (l << 8) + l;
         }
         y += fImage->width;
      }
   } else {
      fGrayImage = create_asimage(fImage->width, fImage->height, 0);

      ASImageDecoder *imdec = start_image_decoding(fgVisual, fImage, SCL_DO_ALL,
                                                   0, 0, fImage->width, fImage->height, 0);

      if (!imdec) {
         return;
      }
#ifdef HAVE_MMX
   mmx_init();
#endif
      ASImageOutput *imout = start_image_output(fgVisual, fGrayImage, ASA_ASImage,
                                                GetImageCompression(), GetImageQuality());
      if (!imout) {
         Warning("ToGray", "Failed to start image output");
         delete fScaledImage;
         fScaledImage = 0;
         delete [] imdec;
         return;
      }

      CARD32 *aa = imdec->buffer.alpha;
      CARD32 *rr = imdec->buffer.red;
      CARD32 *gg = imdec->buffer.green;
      CARD32 *bb = imdec->buffer.blue;

      ASScanline result;
      prepare_scanline(fImage->width, 0, &result, fgVisual->BGR_mode);

      for (i = 0; i < fImage->height; i++) {
         imdec->decode_image_scanline(imdec);
         result.flags = imdec->buffer.flags;
         result.back_color = imdec->buffer.back_color;

         for (j = 0; j < fImage->width; j++) {
            l = (57*rr[j] + 181*gg[j]+ 18*bb[j])/256;
            result.alpha[j] = aa[j];
            result.red[j] = l;
            result.green[j] = l;
            result.blue[j] = l;
         }
         imout->output_image_scanline(imout, &result, 1);
      }

      stop_image_decoding(&imdec);
      stop_image_output(&imout);
#ifdef HAVE_MMX
   mmx_off();
#endif
   }

   sav = fImage;
   fImage = fGrayImage;
   fGrayImage = sav;
   fIsGray = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image (screenshot) from specified window.

void TASImage::FromWindow(Drawable_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   Int_t xy;

   x = x < 0 ? 0 : x;
   y = y < 0 ? 0 : y;

   // X11 Synchronization
   gVirtualX->Update(1);
   if (!gThreadXAR) {
      gSystem->Sleep(10);
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
      gSystem->ProcessEvents();
   }

   if (!w || !h) {
      gVirtualX->GetWindowSize(wid, xy, xy, w, h);
   }

   if ((x >= (Int_t)w) || (y >= (Int_t)h)) {
      return;
   }

   if (!InitVisual()) {
      Warning("FromWindow", "Visual not initiated");
      return;
   }

   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   static int x11 = -1;
   if (x11 < 0) x11 = gVirtualX->InheritsFrom("TGX11");

   if (x11) { //use built-in optimized version
      fImage = pixmap2asimage(fgVisual, wid, x, y, w, h, kAllPlanes, 0, 0);
   } else {
      unsigned char *bits = gVirtualX->GetColorBits(wid, 0, 0, w, h);

      if (!bits) { // error
         return;
      }
      fImage = bitmap2asimage(bits, w, h, 0, 0);
      delete [] bits;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates an image (screenshot) from a RGBA buffer.

void TASImage::FromGLBuffer(UChar_t* buf, UInt_t w, UInt_t h)
{
   DestroyImage();
   delete fScaledImage;
   fScaledImage = 0;

   UChar_t* xx = new UChar_t[4*w];
   for (UInt_t i = 0; i < h/2; ++i) {
      memcpy(xx, buf + 4*w*i, 4*w);
      memcpy(buf + 4*w*i, buf + 4*w*(h-i-1), 4*w);
      memcpy(buf + 4*w*(h-i-1), xx, 4*w);
   }
   delete [] xx;

   fImage = bitmap2asimage(buf, w, h, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch on/off the image palette.
/// That also invokes calling vectorization of image.

void TASImage::SetPaletteEnabled(Bool_t on)
{
   if (!fImage) {
      return;
   }

   if (!fImage->alt.vector && on) {
      Vectorize();
   }
   fPaletteEnabled = on;

   if (on) {
      Double_t left = gPad->GetLeftMargin();
      Double_t right = gPad->GetRightMargin();
      Double_t top = gPad->GetTopMargin();
      Double_t bottom = gPad->GetBottomMargin();

      gPad->Range(-left / (1.0 - left - right),
                  -bottom / (1.0 - top - bottom),
                  1 + right / (1.0 - left - right),
                  1 + top / ( 1.0 - top - bottom));
      gPad->RangeAxis(0, 0, 1, 1);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TASImage::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char *buf = 0;
   int sz;

   UInt_t w = GetWidth();
   UInt_t h = GetHeight();

   if (w > 500) { // workaround CINT limitations
      w = 500;
      Double_t scale = 500./GetWidth();
      h = TMath::Nint(GetHeight()*scale);
      Scale(w, h);
   }

   GetImageBuffer(&buf, &sz, TImage::kXpm);

   TString name = GetName();
   name.ReplaceAll(".", "_");
   TString str = buf;
   static int ii = 0;
   ii++;

   str.ReplaceAll("static", "const");
   TString xpm = "xpm_";
   xpm += name;
   xpm += ii;
   str.ReplaceAll("asxpm", xpm.Data());
   out << std::endl << str << std::endl << std::endl;

   out << "   TImage *";
   out << xpm << "_img = TImage::Create();" << std::endl;
   out << "   " << xpm << "_img->SetImageBuffer( (char **)" << xpm << ", TImage::kXpm);" << std::endl;
   out << "   " << xpm << "_img->Draw();" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set an image printing resolution in Dots Per Inch units.
///
///  \param[in] name - the name of jpeg file.
///  \param[in] set - dpi resolution.
///
/// Returns kFALSE in case of error.

Bool_t TASImage::SetJpegDpi(const char *name, UInt_t set)
{
   static char buf[32];
   FILE *fp = fopen(name, "rb+");

   if (!fp) {
      printf("file %s : failed to open\n", name);
      return kFALSE;
   }

   if (!fread(buf, 1, 20, fp)) {
      fclose(fp);
      return kFALSE;
   }

   char dpi1 = (set & 0xffff) >> 8;
   char dpi2 = set & 0xff;

   int i = 0;

   int dpi = 0; // start of dpi data
   for (i = 0; i < 20; i++) {
      if ((buf[i] == 0x4a) && (buf[i+1] == 0x46) &&  (buf[i+2] == 0x49) &&
          (buf[i+3] == 0x46) && (buf[i+4] == 0x00) ) {
         dpi = i + 7;
         break;
      }
   }

   if (i == 20 || dpi+4 >= 20) { // jpeg maker was not found
      fclose(fp);
      printf("file %s : wrong JPEG format\n", name);
      return kFALSE;
   }

   buf[dpi] = 1;   // format specified in  dots per inch

   // set x density in dpi units
   buf[dpi + 1] = dpi1;
   buf[dpi + 2] = dpi2;

   // set y density in dpi units
   buf[dpi + 3] = dpi1;
   buf[dpi + 4] = dpi2;

   rewind(fp);
   fwrite(buf, 1, 20, fp);
   fclose(fp);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a valid index in fImage tables to avoid seg-fault by accessing out of
/// indices out of array's ranges.

Int_t TASImage::Idx(Int_t idx)
{
   // The size of arrays like fImage->alt.argb32 is fImage->width*fImage->height
   return TMath::Min(idx,(Int_t)(fImage->width*fImage->height));
}

