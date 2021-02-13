// $Id: TGHtmlDraw.cxx,v 1.1 2007/05/04 17:07:01 brun Exp $
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

// Routines used to render HTML onto the screen for the TGHtml widget.

#include <cstring>
#include <cstdlib>

#include "TGHtml.h"
#include "TImage.h"
#include "TVirtualX.h"
#include "strlcpy.h"

////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGHtmlBlock::TGHtmlBlock() : TGHtmlElement(Html_Block)
{
   fZ = NULL;
   fTop = fBottom = 0;
   fLeft = fRight = 0;
   fN = 0;
   fPPrev = fPNext = 0;
   fBPrev = fBNext = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// dtor.

TGHtmlBlock::~TGHtmlBlock()
{
   if (fZ) delete[] fZ;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy the given Block after first unlinking it from the element list.
/// Note that this unlinks the block from the element list only -- not from
/// the block list.

void TGHtml::UnlinkAndFreeBlock(TGHtmlBlock *pBlock)
{
   if (pBlock->fPNext) {
      pBlock->fPNext->fPPrev = pBlock->fPPrev;
   } else {
      fPLast = pBlock->fPPrev;
   }
   if (pBlock->fPPrev) {
      pBlock->fPPrev->fPNext = pBlock->fPNext;
   } else {
      fPFirst = pBlock->fPNext;
   }
   pBlock->fPPrev = pBlock->fPNext = 0;
   delete pBlock;
}

////////////////////////////////////////////////////////////////////////////////
/// Append a block to the block list and insert the block into the
/// element list immediately prior to the element given.
///
/// pToken - The token that comes after pBlock
/// pBlock - The block to be appended

void TGHtml::AppendBlock(TGHtmlElement *pToken, TGHtmlBlock *pBlock)
{
   pBlock->fPPrev = pToken->fPPrev;
   pBlock->fPNext = pToken;
   pBlock->fBPrev = fLastBlock;
   pBlock->fBNext = 0;
   if (fLastBlock) {
      fLastBlock->fBNext = pBlock;
   } else {
      fFirstBlock = pBlock;
   }
   fLastBlock = pBlock;
   if (pToken->fPPrev) {
      pToken->fPPrev->fPNext = (TGHtmlElement *) pBlock;
   } else {
      fPFirst = (TGHtmlElement *) pBlock;
   }
   pToken->fPPrev = (TGHtmlElement *) pBlock;
}

////////////////////////////////////////////////////////////////////////////////
/// Print an ordered list index into the given buffer. Use numbering
/// like this:
///
///     A  B  C ... Y Z AA BB CC ... ZZ
///
/// Revert to decimal for indices greater than 52.

static void GetLetterIndex(char *zBuf, int index, int isUpper)
{
   int seed;

   if (index < 1 || index > 52) {
      // coverity[secure_coding]: zBuf is large enough for an integer
      sprintf(zBuf, "%d", index);
      return;
   }

   if (isUpper) {
      seed = 'A';
   } else {
      seed = 'a';
   }

   index--;

   if (index < 26) {
      zBuf[0] = seed + index;
      zBuf[1] = 0;
   } else {
      index -= 26;
      zBuf[0] = seed + index;
      zBuf[1] = seed + index;
      zBuf[2] = 0;
   }

   strcat(zBuf, ".");
}

////////////////////////////////////////////////////////////////////////////////
/// Print an ordered list index into the given buffer.  Use roman
/// numerals.  For indices greater than a few thousand, revert to
/// decimal.

static void GetRomanIndex(char *zBuf, int index, int isUpper)
{
   int i = 0;
   UInt_t j;

   static struct {
      int value;
      const char *name;
   } values[] = {
      { 1000, "m"  },
      {  999, "im" },
      {  990, "xm" },
      {  900, "cm" },
      {  500, "d"  },
      {  499, "id" },
      {  490, "xd" },
      {  400, "cd" },
      {  100, "c"  },
      {   99, "ic" },
      {   90, "xc" },
      {   50, "l"  },
      {   49, "il" },
      {   40, "xl" },
      {   10, "x"  },
      {    9, "ix" },
      {    5, "v"  },
      {    4, "iv" },
      {    1, "i"  },
   };

   if (index < 1 || index >= 5000) {
      // coverity[secure_coding]: zBuf is large enough for an integer
      sprintf(zBuf, "%d", index);
      return;
   }
   for (j = 0; index > 0 && j < sizeof(values)/sizeof(values[0]); j++) {
      int k;
      while (index >= values[j].value) {
         for (k = 0; values[j].name[k]; k++) {
            zBuf[i++] = values[j].name[k];
         }
         index -= values[j].value;
      }
   }
   zBuf[i] = 0;
   if (isUpper) {
      for (i = 0; zBuf[i]; i++) {
         zBuf[i] += 'A' - 'a';
      }
   }

   strcat(zBuf, ".");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the selection background for the given block
///
/// x, y - Virtual coords of top-left of drawable

void TGHtml::DrawSelectionBackground(TGHtmlBlock *pBlock, Drawable_t drawable,
                                     int x, int y)
{
   int xLeft, xRight;        // Left and right bounds of box to draw
   int yTop, yBottom;        // Top and bottom of box
   TGHtmlElement *p = 0;     // First element of the block
   TGFont *font=0;           // Font
   GContext_t gc;            // GC for drawing

   if (pBlock == 0 || (pBlock->fFlags & HTML_Selected) == 0) return;

   xLeft = pBlock->fLeft - x;
   if (pBlock == fPSelStartBlock && fSelStartIndex > 0) {
      if (fSelStartIndex >= pBlock->fN) return;
      p = pBlock->fPNext;
      font = GetFont(p->fStyle.fFont);
      if (font == 0) return;
      if (p->fType == Html_Text) {
         TGHtmlTextElement *tp = (TGHtmlTextElement *) p;
         xLeft = tp->fX - x + font->TextWidth(pBlock->fZ, fSelStartIndex);
      }
   }
   xRight = pBlock->fRight - x;
   if (pBlock == fPSelEndBlock && fSelEndIndex < pBlock->fN) {
      if (p == 0) {
         p = pBlock->fPNext;
         font = GetFont(p->fStyle.fFont);
         if (font == 0) return;
      }
      if (p->fType == Html_Text) {
         TGHtmlTextElement *tp = (TGHtmlTextElement *) p;
         xRight = tp->fX - x + font->TextWidth(pBlock->fZ, fSelEndIndex);
      }
   }
   yTop = pBlock->fTop - y;
   yBottom = pBlock->fBottom - y;
   gc = GetGC(COLOR_Selection, FONT_Any);
   Int_t xx = xLeft;
   Int_t yy = yTop;
   UInt_t width = xRight - xLeft;
   UInt_t height = yBottom - yTop;
   gVirtualX->FillRectangle(drawable, gc, xx, yy, width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a rectangle. The rectangle will have a 3-D appearance if
/// flat is 0 and a flat appearance if flat is 1.
///
/// depth - width of the relief or the flat line

void TGHtml::DrawRect(Drawable_t drawable, TGHtmlElement *src,
                      int x, int y, int w, int h, int depth, int relief)
{
   Int_t xx, yy;
   UInt_t width, height;

   if (depth > 0) {
      int i;
      GContext_t gcLight, gcDark;

      if (relief != HTML_RELIEF_FLAT) {
         int iLight1, iDark1;
         iLight1 = GetLightShadowColor(src->fStyle.fBgcolor);
         gcLight = GetGC(iLight1, FONT_Any);
         iDark1 = GetDarkShadowColor(src->fStyle.fBgcolor);
         gcDark = GetGC(iDark1, FONT_Any);
         if (relief == HTML_RELIEF_SUNKEN) {
            GContext_t gcTemp = gcLight;
            gcLight = gcDark;
            gcDark = gcTemp;
         }
      } else {
         gcLight = GetGC(src->fStyle.fColor, FONT_Any);
         gcDark = gcLight;
      }
      xx = x;
      yy = y;
      width = depth;
      height = h;
      gVirtualX->FillRectangle(drawable, gcLight, xx, yy, width, height);
      xx = x + w - depth;
      gVirtualX->FillRectangle(drawable, gcLight, xx, yy, width, height);
      for (i = 0; i < depth && i < h/2; i++) {
         gVirtualX->DrawLine(drawable, gcLight, x+i, y+i, x+w-i-1, y+i);
         gVirtualX->DrawLine(drawable, gcDark, x+i, y+h-i-1, x+w-i-1, y+h-i-1);
      }
   }
   if (h > depth*2 && w > depth*2) {
      GContext_t gcBg;
      gcBg = GetGC(src->fStyle.fBgcolor, FONT_Any);
      xx = x + depth;
      yy = y + depth;
      width = w - depth*2;
      height = h - depth*2;
      gVirtualX->FillRectangle(drawable, gcBg, xx, yy, width, height);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display a single HtmlBlock. This is where all the drawing happens.

void TGHtml::BlockDraw(TGHtmlBlock *pBlock, Drawable_t drawable,
                       int drawableLeft, int drawableTop,
                       int drawableWidth, int drawableHeight,
                       Pixmap_t pixmap)
{
   TGFont *font;           // Font to use to render text
   GContext_t gc;          // A graphics context
   TGHtmlElement *src;     // TGHtmlElement holding style information
   TGHtmlTable *pTable;    // The table (when drawing part of a table)
   Int_t x, y;             // Where to draw
   UInt_t width, height;

   if (pBlock == 0) return;

   src = pBlock->fPNext;
   while (src && (src->fFlags & HTML_Visible) == 0) src = src->fPNext;

   if (src == 0) return;

   if (pBlock->fN > 0) {
      // We must be dealing with plain old text
      if (src->fType == Html_Text) {
         TGHtmlTextElement *tsrc = (TGHtmlTextElement *) src;
         x = tsrc->fX;
         y = tsrc->fY;
      } else {
         CANT_HAPPEN;
         return;
      }
      if (pBlock->fFlags & HTML_Selected) {
         DrawSelectionBackground(pBlock, drawable, drawableLeft, drawableTop);
      }
      gc = GetGC(src->fStyle.fColor, src->fStyle.fFont);
      font = GetFont(src->fStyle.fFont);
      if (font == 0) return;
      font->DrawChars(drawable, gc, pBlock->fZ, pBlock->fN,
                      x - drawableLeft, y - drawableTop);
      if (src->fStyle.fFlags & STY_Underline) {
         font->UnderlineChars(drawable, gc, pBlock->fZ,
                              x - drawableLeft, y-drawableTop, 0, pBlock->fN);
      }
      if (src->fStyle.fFlags & STY_StrikeThru) {
         x = pBlock->fLeft - drawableLeft;
         y = (pBlock->fTop + pBlock->fBottom) / 2 - drawableTop;
         width = pBlock->fRight - pBlock->fLeft;
         height = 1 + (pBlock->fBottom - pBlock->fTop > 15);
         gVirtualX->FillRectangle(drawable, gc, x, y, width, height);
      }
      if (pBlock == fPInsBlock && fInsStatus > 0) {
         if (fInsIndex < pBlock->fN) {
            TGHtmlTextElement *tsrc = (TGHtmlTextElement *) src;
            x = tsrc->fX - drawableLeft;
            x += font->TextWidth(pBlock->fZ, fInsIndex);
         } else {
            x = pBlock->fRight - drawableLeft;
         }
         if (x > 0) --x;
            gVirtualX->FillRectangle(drawable, gc, x, pBlock->fTop - drawableTop,
                                     2, pBlock->fBottom - pBlock->fTop);
         }
   } else {
      // We are dealing with a single TGHtmlElement which contains something
      // other than plain text.
      int cnt, w;
      char zBuf[30];
      TGHtmlLi *li;
      TGHtmlImageMarkup *image;
      switch (src->fType) {
         case Html_LI:
            li = (TGHtmlLi *) src;
            x = li->fX;
            y = li->fY;
            switch (li->fLtype) {
               case LI_TYPE_Enum_1:
                  // coverity[secure_coding]: zBuf is large enough for an int
                  sprintf(zBuf, "%d.", li->fCnt);
                  break;
               case LI_TYPE_Enum_A:
                  GetLetterIndex(zBuf, li->fCnt, 1);
                  break;
               case LI_TYPE_Enum_a:
                  GetLetterIndex(zBuf, li->fCnt, 0);
                  break;
               case LI_TYPE_Enum_I:
                  GetRomanIndex(zBuf, li->fCnt, 1);
                  break;
               case LI_TYPE_Enum_i:
                  GetRomanIndex(zBuf, li->fCnt, 0);
                  break;
               default:
                  zBuf[0] = 0;
                  break;
            }
            gc = GetGC(src->fStyle.fColor, src->fStyle.fFont);
            switch (li->fLtype) {
               case LI_TYPE_Undefined:
               case LI_TYPE_Bullet1:
                  //gVirtualX->FillArc(drawable, gc,
                  //         x - 7 - drawableLeft, y - 8 - drawableTop, 7, 7,
                  //         0, 360*64);
                  break;

               case LI_TYPE_Bullet2:
                  //gVirtualX->DrawArc(drawable, gc,
                  //         x - 7 - drawableLeft, y - 8 - drawableTop, 7, 7,
                  //         0, 360*64);
                  break;

               case LI_TYPE_Bullet3:
                     gVirtualX->DrawRectangle(drawable, gc, x - 7 - drawableLeft,
                                              y - 8 - drawableTop, 7, 7);
                  break;

               case LI_TYPE_Enum_1:
               case LI_TYPE_Enum_A:
               case LI_TYPE_Enum_a:
               case LI_TYPE_Enum_I:
               case LI_TYPE_Enum_i:
                  cnt = strlen(zBuf);
                  font = GetFont(src->fStyle.fFont);
                  if (font == 0) return;
                  w = font->TextWidth(zBuf, cnt);
                  font->DrawChars(drawable, gc, zBuf, cnt,
                                  x - w - drawableLeft, y - drawableTop);
                  break;
            }
            break;

         case Html_HR: {
            TGHtmlHr *hr = (TGHtmlHr *) src;
            int relief = fRuleRelief;
            switch (relief) {
               case HTML_RELIEF_RAISED:
               case HTML_RELIEF_SUNKEN:
               break;
               default:
                  relief = HTML_RELIEF_FLAT;
                  break;
            }
            DrawRect(drawable, src, hr->fX - drawableLeft, hr->fY - drawableTop,
                     hr->fW, hr->fH, 1, relief);
            break;
         }

         case Html_TABLE: {
            TGHtmlTable *table = (TGHtmlTable *) src;
            int relief = fTableRelief;
            if ((!fBgImage || src->fStyle.fExpbg) && !table->fHasbg) {
               switch (relief) {
                  case HTML_RELIEF_RAISED:
                  case HTML_RELIEF_SUNKEN:
                     break;
                  default:
                     relief = HTML_RELIEF_FLAT;
                     break;
               }

               DrawRect(drawable, src, table->fX - drawableLeft,
                        table->fY - drawableTop, table->fW, table->fH,
                        table->fBorderWidth, relief);
            }

            if (table->fBgImage) {
               DrawTableBgnd(table->fX, table->fY, table->fW, table->fH, pixmap,
                             table->fBgImage);
            }
            break;
         }

         case Html_TH:
         case Html_TD: {
            TGHtmlCell *cell = (TGHtmlCell *) src;
            int depth, relief;
            TImage *bgImg;
            pTable = cell->fPTable;
            if ((!fBgImage || src->fStyle.fExpbg) && !(pTable && pTable->fHasbg)) {
               depth = pTable && (pTable->fBorderWidth > 0);
               switch (fTableRelief) {
                  case HTML_RELIEF_RAISED:  relief = HTML_RELIEF_SUNKEN; break;
                  case HTML_RELIEF_SUNKEN:  relief = HTML_RELIEF_RAISED; break;
                  default:                  relief = HTML_RELIEF_FLAT;   break;
               }
               DrawRect(drawable, src,
                        cell->fX - drawableLeft, cell->fY - drawableTop,
                        cell->fW, cell->fH, depth, relief);
            }
            // See if row has an image
            if (cell->fBgImage) {
               DrawTableBgnd(cell->fX, cell->fY, cell->fW, cell->fH, pixmap,
                             cell->fBgImage);
            } else if (cell->fPRow && (bgImg = ((TGHtmlRef *)cell->fPRow)->fBgImage)) {
               DrawTableBgnd(cell->fX, cell->fY, cell->fW, cell->fH, pixmap, bgImg);
            }
            break;
         }

         case Html_IMG:
            image = (TGHtmlImageMarkup *) src;
            if (image->fPImage) {
               DrawImage(image, drawable, drawableLeft, drawableTop,
                         drawableLeft + drawableWidth,
                         drawableTop + drawableHeight);
            } else if (image->fZAlt) {
               gc = GetGC(src->fStyle.fColor, src->fStyle.fFont);
               font = GetFont(src->fStyle.fFont);
               if (font == 0) return;
               font->DrawChars(drawable, gc,
                               image->fZAlt, strlen(image->fZAlt),
                               image->fX - drawableLeft,
                               image->fY - drawableTop);
            }
            break;

         default:
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw all or part of an image.

void TGHtml::DrawImage(TGHtmlImageMarkup *image, Drawable_t drawable,
                       int drawableLeft, int drawableTop,
                       int drawableRight, int drawableBottom)
{
   int imageTop;          // virtual canvas coordinate for top of image
   int x, y;              // where to place image on the drawable
   int imageX, imageY;    // \__  Subset of image that fits
   int imageW, imageH;    // /    on the drawable

   imageTop = image->fY - image->fAscent;
   y = imageTop - drawableTop;
   if (imageTop + image->fH > drawableBottom) {
      imageH = drawableBottom - imageTop;
   } else {
      imageH = image->fH;
   }
   if (y < 0) {
      imageY = -y;
      imageH += y;
      y = 0;
   } else {
      imageY = 0;
   }
   x = image->fX - drawableLeft;
   if (image->fX + image->fW > drawableRight) {
      imageW = drawableRight - image->fX;
   } else {
      imageW = image->fW;
   }
   if (x < 0) {
      imageX = -x;
      imageW += x;
      x = 0;
   } else {
      imageX = 0;
   }

   TImage *img = image->fPImage->fImage;

   imageH = imageH < 0 ? -imageH : imageH;
   imageW = imageW < 0 ? -imageW : imageW;

   img->PaintImage(drawable, x, y, imageX, imageY, imageW, imageH);
   //gVirtualX->Update(kFALSE);

   image->fRedrawNeeded = 0;
}

////////////////////////////////////////////////////////////////////////////////
///
///TGImage *img = image->image;

void TGHtml::AnimateImage(TGHtmlImage * /*image*/)
{
  //if (!img->IsAnimated()) return;
  //img->NextFrame();
  //delete image->timer;
  //image->timer = new TTimer(this, img->GetAnimDelay());
  //ImageChanged(image, image->fW, image->fH);
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute the following fields of the given block structure:
///
///    base.count         The number of elements described by this
///                       block structure.
///
///    n                  The number of characters of text output
///                       associated with this block.  If the block
///                       renders something other than text (ex: <IMG>)
///                       then set n to 0.
///
///    z                  Pointer to malloced memory containing the
///                       text associated with this block.  NULL if
///                       n is 0.
///
/// Return a pointer to the first TGHtmlElement not covered by the block.

TGHtmlElement *TGHtml::FillOutBlock(TGHtmlBlock *p)
{

   TGHtmlElement *pElem;
   int go, i, n, x, y;
   SHtmlStyle_t style;
   char zBuf[2000];

   // Reset n and z

   if (p->fN) p->fN = 0;

   if (p->fZ) delete[] p->fZ;
   p->fZ = 0;

   // Skip over TGHtmlElements that aren't directly displayed.

   pElem = p->fPNext;
   p->fCount = 0;
   while (pElem && (pElem->fFlags & HTML_Visible) == 0) {
      TGHtmlElement *fPNext = pElem->fPNext;
      if (pElem->fType == Html_Block) {
         UnlinkAndFreeBlock((TGHtmlBlock *) pElem);
      } else {
         p->fCount++;
      }
      pElem = fPNext;
   }
   if (pElem == 0) return 0;

   // Handle "special" elements.

   if (pElem->fType != Html_Text) {
      switch (pElem->fType) {
         case Html_HR: {
            TGHtmlHr *hr = (TGHtmlHr *) pElem;
            p->fTop    = hr->fY - hr->fH;
            p->fBottom = hr->fY;
            p->fLeft   = hr->fX;
            p->fRight  = hr->fX + hr->fW;
            break;
         }

         case Html_LI: {
            TGHtmlLi *li = (TGHtmlLi *) pElem;
            p->fTop    = li->fY - li->fAscent;
            p->fBottom = li->fY + li->fDescent;
            p->fLeft   = li->fX - 10;
            p->fRight  = li->fX + 10;
            break;
         }

         case Html_TD:
         case Html_TH: {
            TGHtmlCell *cell = (TGHtmlCell *) pElem;
            p->fTop    = cell->fY;
            p->fBottom = cell->fY + cell->fH;
            p->fLeft   = cell->fX;
            p->fRight  = cell->fX + cell->fW;
            break;
         }

         case Html_TABLE: {
            TGHtmlTable *table = (TGHtmlTable *) pElem;
            p->fTop    = table->fY;
            p->fBottom = table->fY + table->fH;
            p->fLeft   = table->fX;
            p->fRight  = table->fX + table->fW;
            break;
         }

         case Html_IMG: {
            TGHtmlImageMarkup *image = (TGHtmlImageMarkup *) pElem;
            p->fTop    = image->fY - image->fAscent;
            p->fBottom = image->fY + image->fDescent;
            p->fLeft   = image->fX;
            p->fRight  = image->fX + image->fW;
            break;
         }
      }
      p->fCount++;

      return pElem->fPNext;
   }

   // If we get this far, we must be dealing with text.

   TGHtmlTextElement *text = (TGHtmlTextElement *) pElem;
   n = 0;
   x = text->fX;
   y = text->fY;
   p->fTop = y - text->fAscent;
   p->fBottom = y + text->fDescent;
   p->fLeft = x;
   style = pElem->fStyle;
   go = 1;
   while (pElem) {
      TGHtmlElement *fPNext = pElem->fPNext;
      switch (pElem->fType) {
         case Html_Text: {
            TGHtmlTextElement *txt = (TGHtmlTextElement *) pElem;
            if (pElem->fFlags & STY_Invisible) {
               break;
            }
            if (txt->fSpaceWidth <= 0) {
               //CANT_HAPPEN;
               break;
            }
            if (y != txt->fY
                  ||  style.fFont != pElem->fStyle.fFont
                  ||  style.fColor != pElem->fStyle.fColor
                  ||  (style.fFlags & STY_FontMask)
                  != (pElem->fStyle.fFlags & STY_FontMask)) {
               go = 0;
            } else {
               int sw = txt->fSpaceWidth;
               int nSpace = (txt->fX - x) / sw;
               if (nSpace * sw + x != txt->fX) {
                  go = 0;
               } else if ((n + nSpace + pElem->fCount) >= (int)sizeof(zBuf)) {
                  // go = 0; - this caused a hang, instead lets do what we can
                  for (i = 0; i < nSpace && (n+1) < (int)sizeof(zBuf); ++i) {
                     zBuf[n++] = ' ';
                  }
                  strncpy(&zBuf[n], txt->fZText, sizeof(zBuf) - n - 1);
                  zBuf[sizeof(zBuf)-1] = 0;
                  n += i;
                  x = txt->fX + txt->fW;
               } else {
                  for (i = 0; i < nSpace && (n+1) < (int)sizeof(zBuf); ++i) {
                     zBuf[n++] = ' ';
                  }
                  strncpy(&zBuf[n], txt->fZText, sizeof(zBuf) - n - 1);
                  zBuf[sizeof(zBuf)-1] = 0;
                  n += pElem->fCount;
                  x = txt->fX + txt->fW;
               }
            }
            break;
         }

         case Html_Space:
            if (pElem->fStyle.fFont != style.fFont) {
               pElem = pElem->fPNext;
               go = 0;
            } else if ((style.fFlags & STY_Preformatted) != 0 &&
                       (pElem->fFlags & HTML_NewLine) != 0) {
               pElem = pElem->fPNext;
               go = 0;
            }
            break;

         case Html_Block:
            UnlinkAndFreeBlock((TGHtmlBlock *) pElem);
            break;

         case Html_A:
         case Html_EndA:
            go = 0;
            break;

         default:
            if (pElem->fFlags & HTML_Visible) go = 0;
            break;
      }
      if (go == 0) break;
      p->fCount++;
      pElem = fPNext;
   }
   p->fRight = x;

   while (n > 0 && zBuf[n-1] == ' ') n--;
   p->fZ = new char[n+1];
   strlcpy(p->fZ, zBuf, n+1);
   p->fZ[n] = 0;
   p->fN = n;

   return pElem;
}

////////////////////////////////////////////////////////////////////////////////
/// Scan ahead looking for a place to put a block.  Return a pointer
/// to the element which should come immediately after the block.
///
/// if pCnt != 0, then put the number of elements skipped in *pCnt.
///
/// p    - First candidate for the start of a block
/// pCnt - Write number of elements skipped here

TGHtmlElement *TGHtml::FindStartOfNextBlock(TGHtmlElement *p, int *pCnt)
{
   int cnt = 0;

   while (p && (p->fFlags & HTML_Visible) == 0) {
      TGHtmlElement *fPNext = p->fPNext;
      if (p->fType == Html_Block) {
         UnlinkAndFreeBlock((TGHtmlBlock *) p);
      } else {
         cnt++;
      }
      p = fPNext;
   }
   if (pCnt) *pCnt = cnt;

   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// Add additional blocks to the block list in order to cover
/// all elements on the element list.
///
/// If any old blocks are found on the element list, they must
/// be left over from a prior rendering.  Unlink and delete them.

void TGHtml::FormBlocks()
{
   TGHtmlElement *pElem;

   if (fLastBlock) {
      pElem = FillOutBlock(fLastBlock);
   } else {
      pElem = fPFirst;
   }
   while (pElem) {
      int cnt;
      pElem = FindStartOfNextBlock(pElem, &cnt);
      if (pElem) {
         TGHtmlBlock *pNew = new TGHtmlBlock();
         if (fLastBlock) {
            fLastBlock->fCount += cnt;
         }
         AppendBlock(pElem, pNew);
         pElem = FillOutBlock(pNew);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw table background

void TGHtml::DrawTableBgnd(int l, int t, int w, int h,
                           Drawable_t pixmap, TImage *image)
{
   //int  mx, my, sh, sw, sx, sy, hd;
   int dl, dt, dr, db,  left, top, right, bottom;

   left = l - fVisible.fX;
   top  = t - fVisible.fY;

   dl = fDirtyLeft;
   dt = fDirtyTop;
   dr = fDirtyRight;
   db = fDirtyBottom;

   right = left + w - 1;
   bottom = top + h - 1;
   if (dr == 0 && db == 0) { dr = right; db = bottom; }
   if (left > dr || right < dl || top > db || bottom < dt) return;

#if 0
   int iw = image->GetWidth();
   int ih = image->GetHeight();
   if (iw < 4 && ih < 4) return;  // CPU burners we ignore.
   sx = (left + _visibleStart.x) % iw;   // X offset within image to start from
   sw = iw - sx;                         // Width of section of image to draw.
   for (mx = left - dl; w > 0; mx += sw, sw = iw, sx = 0) {
      if (sw > w) sw = w;
      sy = (top + _visibleStart.y) % ih;  // Y offset within image to start from
      sh = ih - sy;                       // Height of section of image to draw.
      for (my = top - dt, hd = h; hd > 0; my += sh, sh = ih, sy = 0) {
         if (sh > hd) sh = hd;
         // printf("image: %d %d %d %d %d %d\n", sx, sy, sw, sh, mx,my);
         image->Draw(pixmap, GetAnyGC(), sx, sy, sw, sh, mx, my);
         hd -= sh;
      }
      w -= sw;
   }
#else
   if (!image->GetPixmap()) return;
   GContext_t gc = GetAnyGC();
   GCValues_t gcv;
   // unsigned int mask = kGCTile | kGCFillStyle |
   //                     kGCTileStipXOrigin | kGCTileStipYOrigin;
   gcv.fTile      = image->GetPixmap();
   gcv.fFillStyle = kFillTiled;
   gcv.fTsXOrigin = -fVisible.fX - fDirtyLeft;
   gcv.fTsYOrigin = -fVisible.fY - fDirtyTop;
   gVirtualX->ChangeGC(gc, &gcv);

   gVirtualX->FillRectangle(pixmap, gc, left - dl, top - dt, w, h);

   // mask = kGCFillStyle;
   gcv.fFillStyle = kFillSolid;
   gVirtualX->ChangeGC(gc, &gcv);
#endif
}
