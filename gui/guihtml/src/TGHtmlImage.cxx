// $Id: TGHtmlImage.cxx,v 1.2 2007/05/07 15:28:48 brun Exp $
// Author:  Valeriy Onuchin   03/05/2007

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************

    HTML widget for xclass. Based on tkhtml 1.28
    Copyright (C) 1997-2000 D. Richard Hipp <drh@acm.org>
    Copyright (C) 2002-2003 Hector Peraza.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

**************************************************************************/

// Routines used for processing <IMG> markup

#include <cstring>
#include <cstdlib>

#include "TGHtml.h"
#include "snprintf.h"
#include "TImage.h"
#include "TUrl.h"
#include "TSocket.h"
#include "TSystem.h"
#include "TError.h"
#ifdef R__SSL
#include "TSSLSocket.h"
#endif

////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGHtmlImage::TGHtmlImage(TGHtml *htm, const char *url, const char *width,
                       const char *height)
{
   fHtml = htm;
   fZUrl = StrDup(url);
   fZWidth = StrDup(width);
   fZHeight = StrDup(height);
   fImage = NULL;
   fPNext = NULL;
   fPList = NULL;
   fW = 0;
   fH = 0;
   fTimer = NULL;
}

////////////////////////////////////////////////////////////////////////////////
///  dtor.

TGHtmlImage::~TGHtmlImage()
{
   delete [] fZUrl;
   delete [] fZWidth;
   delete [] fZHeight;

   if (fImage) delete fImage;
   if (fTimer) delete fTimer;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the alignment for an image

int TGHtml::GetImageAlignment(TGHtmlElement *p)
{
   const char *z;
   int i;
   int result;

   static struct {
      const char *zName;
      int iValue;
   } aligns[] = {
      { "bottom",    IMAGE_ALIGN_Bottom    },
      { "baseline",  IMAGE_ALIGN_Bottom    },
      { "middle",    IMAGE_ALIGN_Middle    },
      { "top",       IMAGE_ALIGN_Top       },
      { "absbottom", IMAGE_ALIGN_AbsBottom },
      { "absmiddle", IMAGE_ALIGN_AbsMiddle },
      { "texttop",   IMAGE_ALIGN_TextTop   },
      { "left",      IMAGE_ALIGN_Left      },
      { "right",     IMAGE_ALIGN_Right     },
   };

   z = p->MarkupArg("align", 0);
   result = IMAGE_ALIGN_Bottom;
   if (z) {
      for (i = 0; i < int(sizeof(aligns) / sizeof(aligns[0])); i++) {
         if (strcasecmp(aligns[i].zName, z) == 0) {
            result = aligns[i].iValue;
            break;
         }
      }
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine is called when an image changes. If the size of the
/// images changes, then we need to completely redo the layout. If
/// only the appearance changes, then this works like an expose event.
///
/// pImage    - Pointer to an TGHtmlImage object
/// newWidth  - New width of the image
/// newHeight - New height of the image

void TGHtml::ImageChanged(TGHtmlImage *pImage, int newWidth, int newHeight)
{
   TGHtmlImageMarkup *pElem;

   if (pImage->fW != newWidth || pImage->fH != newHeight) {
      // We have to completely redo the layout after adjusting the size
      // of the images
      for (pElem = pImage->fPList; pElem; pElem = pElem->fINext) {
         pElem->fW = newWidth;
         pElem->fH = newHeight;
      }
      fFlags |= RELAYOUT;
      pImage->fW = newWidth;
      pImage->fH = newHeight;
      RedrawEverything();
   } else {
#if 0
      for (pElem = pImage->fPList; pElem; pElem = pElem->fINext) {
         pElem->fRedrawNeeded = 1;
      }
      fFlags |= REDRAW_IMAGES;
      ScheduleRedraw();
#else
      for (pElem = pImage->fPList; pElem; pElem = pElem->fINext) {
         pElem->fRedrawNeeded = 1;
         DrawRegion(pElem->fX, pElem->fY - pElem->fAscent, pElem->fW, pElem->fH);
      }
#endif
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Given an <IMG> markup, find or create an appropriate TGHtmlImage
/// object and return a pointer to that object. NULL might be returned.

TGHtmlImage *TGHtml::GetImage(TGHtmlImageMarkup *p)
{
   const char *zWidth;
   const char *zHeight;
   const char *zSrc;
   TGHtmlImage *pImage;

   if (p->fType != Html_IMG) { CANT_HAPPEN; return 0; }

   zSrc = p->MarkupArg("src", 0);
   if (zSrc == 0) return 0;

   zSrc = ResolveUri(zSrc);
   if (zSrc == 0) return 0;

   zWidth = p->MarkupArg("width", "");
   zHeight = p->MarkupArg("height", "");

   //p->w = atoi(fZWidth);
   //p->h = atoi(zHeight);

   for (pImage = fImageList; pImage; pImage = pImage->fPNext) {
      if (strcmp(pImage->fZUrl, zSrc) == 0
          &&  strcmp(pImage->fZWidth, zWidth) == 0
          &&  strcmp(pImage->fZHeight, zHeight) == 0) {
         delete [] zSrc;
         return pImage;
      }
   }

   TImage *img = LoadImage(zSrc, atoi(zWidth), atoi(zHeight));

   if (img) {
      pImage = new TGHtmlImage(this, zSrc, zWidth, zHeight);
      pImage->fImage = img;
    //if (img->IsAnimated()) {
    //  pImage->timer = new TTimer(this, img->GetAnimDelay());
    //}
      ImageChanged(pImage, img->GetWidth(), img->GetHeight());
      pImage->fPNext = fImageList;
      fImageList = pImage;
   } else {
      pImage = 0;
   }

   delete [] zSrc;

   return pImage;
}

////////////////////////////////////////////////////////////////////////////////
/// Temporary function to read remote pictures

static TImage *ReadRemoteImage(const char *url)
{
   TImage *image = 0;
   FILE *tmp;
   char *buf;
   TUrl fUrl(url);

   TString msg = "GET ";
   msg += fUrl.GetProtocol();
   msg += "://";
   msg += fUrl.GetHost();
   msg += ":";
   msg += fUrl.GetPort();
   msg += "/";
   msg += fUrl.GetFile();
   msg += "\r\n";

   TString uri(url);
   if ((!uri.BeginsWith("http://") && !uri.BeginsWith("https://")) ||
       uri.EndsWith(".html"))
      return 0;
   TSocket *s;
   if (uri.BeginsWith("https://")) {
#ifdef R__SSL
      s = new TSSLSocket(fUrl.GetHost(), fUrl.GetPort());
#else
      ::Error("ReadRemoteImage", "library compiled without SSL, https not supported");
      return 0;
#endif
   }
   else {
      s = new TSocket(fUrl.GetHost(), fUrl.GetPort());
   }
   if (!s->IsValid()) {
      delete s;
      return 0;
   }
   if (s->SendRaw(msg.Data(), msg.Length()) == -1) {
      delete s;
      return 0;
   }
   Int_t size = 1024*1024;
   buf = (char *)calloc(size, sizeof(char));
   if (!buf) {
      delete s;
      return 0;
   }
   if (s->RecvRaw(buf, size) == -1) {
      free(buf);
      delete s;
      return 0;
   }
   TString pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(),
                                     gSystem->BaseName(url));
   tmp = fopen(pathtmp.Data(), "wb");
   if (!tmp) {
      free(buf);
      delete s;
      return 0;
   }
   fwrite(buf, sizeof(char), size, tmp);
   fclose(tmp);
   free(buf);
   image = TImage::Open(pathtmp.Data());
   if (image && !image->IsValid()) {
      delete image;
      image = 0;
   }
   gSystem->Unlink(pathtmp.Data());
   delete s;
   return image;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the default LoadImage() procedure. It just tries to load the
/// image from a file in the local filesystem.

TImage *TGHtml::LoadImage(const char *url, int w, int h)
{
   TImage *image = 0;

   //TGHtmlUri uri(url);

   TString uri(url);
   if ((uri.BeginsWith("http://") || uri.BeginsWith("https://")) &&
       !uri.EndsWith(".html"))
      image = ReadRemoteImage(url);
   else
      image = TImage::Open(url);
   if (image) {
      if (!image->IsValid()) {
         delete image;
         image = 0;
         return 0;
      }
      if ((w > 0 && h > 0) && ((w != (int)image->GetWidth()) ||
          (h != (int)image->GetHeight()))) {
         image->Scale(w, h);
      }
   }
   return image;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the height and width, converting to percent if required
/// ret must be at least 16 characters long

const char *TGHtml::GetPctWidth(TGHtmlElement *p, char *opt, char *ret)
{
   int n, m, val;
   const char *tz, *z;
   TGHtmlElement *pElem = p;

   z = pElem->MarkupArg(opt, "");
   if (!z) return z;
   if (!strchr(z, '%')) return z;
   // coverity[secure_coding]
   if (!sscanf(z, "%d", &n)) return z;
   if (n <= 0 || n > 100) return z;
   if (opt[0] == 'h') {
      val = fCanvas->GetHeight() * 100;
   } else {
      val = fCanvas->GetWidth() * 100;
   }
   if (!fInTd) {
      snprintf(ret, 15, "%d", val / n);
   } else {
      while (pElem && pElem->fType != Html_TD) pElem = pElem->fPPrev;
      if (!pElem) return z;
      tz = pElem->MarkupArg(opt, 0);
      // coverity[secure_coding]
      if (tz && !strchr(tz, '%') && sscanf(tz, "%d", &m)) {
         snprintf(ret, 15, "%d", m * 100 / n);
         return ret;
      }
      pElem = ((TGHtmlCell *)pElem)->fPTable;
      if (!pElem) return z;
      tz = pElem->MarkupArg(opt, 0);
      // coverity[secure_coding]
      if (tz && !strchr(tz, '%') && sscanf(tz, "%d", &m)) {
         snprintf(ret, 15, "%d", m * 100 / n);
         return ret;
      }
      return z;
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine searchs for an image beneath the coordinates x,y
/// and returns the token number of the the image, or -1 if no
/// image found.

int TGHtml::GetImageAt(int x, int y)
{
   TGHtmlBlock *pBlock;
   TGHtmlElement *pElem;
   //int n;

   for (pBlock = fFirstBlock; pBlock; pBlock = pBlock->fBNext) {
      if (pBlock->fTop > y || pBlock->fBottom < y ||
          pBlock->fLeft > x || pBlock->fRight < x) {
         continue;
      }
      for (pElem = pBlock->fPNext; pElem; pElem = pElem->fPNext) {
         if (pBlock->fBNext && pElem == pBlock->fBNext->fPNext) break;
         if (pElem->fType == Html_IMG) {
            return TokenNumber(pElem);
         }
      }
   }

   return -1;
}
