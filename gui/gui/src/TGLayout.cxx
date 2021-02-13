// @(#)root/gui:$Id$
// Author: Fons Rademakers   02/01/98

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
// A number of different layout classes (TGLayoutManager,               //
// TGVerticalLayout, TGHorizontalLayout, TGLayoutHints, etc.).          //
//                                                                      //
//                                                                      //
// Concerning the TGMatrixLayout class:                                 //
//                                                                      //
// It arranges frames in a matrix-like way.                             //
// This manager provides :                                              //
// - a column number (0 means unlimited)                                //
// - a row number (0 means unlimited)                                   //
// - horizontal & vertical separators                                   //
//                                                                      //
// Notes : If both column and row are fixed values, any remaining       //
//         frames outside the count won't be managed.                   //
//         Unlimited rows means the frame can expand downward           //
//         (the default behaviour in most UI).                          //
//         Both unlimited rows and columns is undefined (read: will     //
//         crash the algorithm ;-).                                     //
//         With fixed dimensions, frames are always arranged in rows.   //
//         That is: 1st frame is at position (0,0), next one is at      //
//         row(0), column(1) and so on...                               //
//         When specifying one dimension as unlimited (i.e. row=0 or    //
//         column=0) the frames are arranged according to the direction //
//         of the fixed dimension. This layout manager does not make    //
//         use of TGLayoutHints.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGLayout.h"
#include "TGFrame.h"
#include "TList.h"

#include <iostream>


ClassImp(TGLayoutHints);
ClassImp(TGLayoutManager);
ClassImp(TGVerticalLayout);
ClassImp(TGHorizontalLayout);
ClassImp(TGRowLayout);
ClassImp(TGColumnLayout);
ClassImp(TGMatrixLayout);
ClassImp(TGTileLayout);
ClassImp(TGListLayout);
ClassImp(TGListDetailsLayout);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGFrameElement::TGFrameElement(TGFrame *f, TGLayoutHints *l)
{
   fLayout = 0;
   fFrame  = f;
   if (f) f->SetFrameElement(this);

   if (l) {
      l->AddReference();
      fLayout = l;
      l->fPrev = l->fFE;
      l->fFE = this;
   }
   fState = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Decrease ref. count of fLayout.

TGFrameElement::~TGFrameElement()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Print this frame element.

void TGFrameElement::Print(Option_t *option) const
{
   TObject::Print(option);

   std::cout << "\t";
   if (fFrame) {
      std::cout << fFrame->ClassName() << "::" << fFrame->GetName();
   }
   if (fLayout) {
      fLayout->Print(option);
   }
   std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLayoutHints::TGLayoutHints(const TGLayoutHints &lh) : TObject(lh), TRefCnt(lh)
{
   fPadleft = lh.fPadleft; fPadright = lh.fPadright;
   fPadtop  = lh.fPadtop;  fPadbottom = lh.fPadbottom;
   fLayoutHints = lh.fLayoutHints;
   SetRefCount(0);
   fFE = lh.fFE; fPrev = lh.fPrev;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLayoutHints::~TGLayoutHints()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update layout hints of frame elements.

void TGLayoutHints::UpdateFrameElements(TGLayoutHints *l)
{
   if (fFE) fFE->fLayout = l;
   else return;

   TGFrameElement *p = fPrev;

   while (p && p->fLayout) {
      p->fLayout = l;
      p = p->fLayout->fPrev;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Printing.

void TGLayoutHints::Print(Option_t *) const
{
   Bool_t bor = kFALSE;

   if (fLayoutHints & kLHintsLeft) {
      std::cout << "kLHintsLeft";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsCenterX) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsCenterX";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsRight) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsRight";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsTop) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsTop";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsCenterY) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsCenterY";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsBottom) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsBottom";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsExpandX) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsExpandX";
      bor = kTRUE;
   }
   if (fLayoutHints & kLHintsExpandY) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsExpandY";
      bor = kTRUE;
   }
   if (fLayoutHints == kLHintsNoHints) {
      if (bor) std::cout << " | ";
      std::cout << "kLHintsNoHints";
   }
   std::cout << ", fPadtop="    << fPadtop;
   std::cout << ", fPadbottom=" << fPadbottom;
   std::cout << ", fPadleft="   << fPadleft;
   std::cout << ", fPadright="  << fPadright;
   std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create vertical layout manager.

TGVerticalLayout::TGVerticalLayout(TGCompositeFrame *main)
{
   fMain = main;
   fList = fMain->GetList();
}

////////////////////////////////////////////////////////////////////////////////
/// Make a vertical layout of all frames in the list.

void TGVerticalLayout::Layout()
{
   TGFrameElement *ptr;
   TGLayoutHints  *layout;
   Int_t    nb_expand = 0;
   Int_t    top, bottom;
   ULong_t  hints;
   UInt_t   extra_space = 0;
   Int_t    exp = 0;
   Int_t    exp_max = 0;
   Int_t    remain;
   Int_t    x = 0, y = 0;
   Int_t    bw = fMain->GetBorderWidth();
   TGDimension size(0,0), csize(0,0);
   TGDimension msize = fMain->GetSize();
   UInt_t pad_left, pad_top, pad_right, pad_bottom;
   Int_t size_expand=0, esize_expand=0, rem_expand=0, tmp_expand = 0;

   if (!fList) return;

   fModified = kFALSE;

   bottom = msize.fHeight - (top = bw);
   remain = msize.fHeight - (bw << 1);

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         layout = ptr->fLayout;
         size = ptr->fFrame->GetDefaultSize();
         size.fHeight += layout->GetPadTop() + layout->GetPadBottom();
         hints = layout->GetLayoutHints();
         if ((hints & kLHintsExpandY) || (hints & kLHintsCenterY)) {
            nb_expand++;
            exp += size.fHeight;
            if (hints & kLHintsExpandY) exp_max = 0;
            else exp_max = TMath::Max(exp_max, (Int_t)size.fHeight);
         } else {
            remain -= size.fHeight;
            if (remain < 0)
               remain = 0;
         }
      }
   }

   if (nb_expand) {
      size_expand = remain/nb_expand;

      if (size_expand < exp_max)
         esize_expand = (remain - exp)/nb_expand;
      rem_expand = remain % nb_expand;
   }

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         hints = (layout = ptr->fLayout)->GetLayoutHints();
         csize      = ptr->fFrame->GetDefaultSize();
         pad_left   = layout->GetPadLeft();
         pad_top    = layout->GetPadTop();
         pad_right  = layout->GetPadRight();
         pad_bottom = layout->GetPadBottom();

         if (hints & kLHintsRight) {
            x = msize.fWidth - bw - csize.fWidth - pad_right;
         } else if (hints & kLHintsCenterX) {
            x = (msize.fWidth - (bw << 1) - csize.fWidth) >> 1;
         } else { // defaults to kLHintsLeft
            x = pad_left + bw;
         }

         if (hints & kLHintsExpandX) {
            size.fWidth = msize.fWidth - (bw << 1) - pad_left - pad_right;
            x = pad_left + bw;
         } else {
            size.fWidth = csize.fWidth;
         }

         if (hints & kLHintsExpandY) {
            if (size_expand >= exp_max)
               size.fHeight = size_expand - pad_top - pad_bottom;
            else
               size.fHeight = csize.fHeight + esize_expand;

            tmp_expand += rem_expand;
            if (tmp_expand >= nb_expand) {
               size.fHeight++;
               tmp_expand -= nb_expand;
            }
         } else {
            size.fHeight = csize.fHeight;
            if (hints & kLHintsCenterY) {
               if (size_expand >= exp_max) {
                  extra_space = (size_expand - pad_top - pad_bottom - size.fHeight) >> 1;
               } else {
                  extra_space = esize_expand >> 1;
               }
               y += extra_space;
               top += extra_space;
            }
         }

         if (hints & kLHintsBottom) {
            y = bottom - size.fHeight - pad_bottom;
            bottom -= size.fHeight + pad_top + pad_bottom;
         } else { // kLHintsTop by default
            y = top + pad_top;
            top += size.fHeight + pad_top + pad_bottom;
         }

         if (hints & kLHintsCenterY)
            top += extra_space;

         if (x > 32768) x = bw + 1;
         //if (y > 32768) y = bw + 1;
         if (size.fWidth > 32768)
            size.fWidth = 1;
         if (size.fHeight > 32768)
            size.fHeight = 1;
         ptr->fFrame->MoveResize(x, y, size.fWidth, size.fHeight);

         fModified = fModified || (ptr->fFrame->GetX() != x) ||
                    (ptr->fFrame->GetY() != y) ||
                    (ptr->fFrame->GetWidth() != size.fWidth) ||
                    (ptr->fFrame->GetHeight() != size.fHeight);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the vertical layout.

TGDimension TGVerticalLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension     size(0,0), msize = fMain->GetSize(), csize;
   UInt_t options = fMain->GetOptions();

   if ((options & kFixedWidth) && (options & kFixedHeight))
      return msize;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         csize = ptr->fFrame->GetDefaultSize();
         size.fWidth = TMath::Max(size.fWidth, csize.fWidth + ptr->fLayout->GetPadLeft() +
                                  ptr->fLayout->GetPadRight());
         size.fHeight += csize.fHeight + ptr->fLayout->GetPadTop() +
                         ptr->fLayout->GetPadBottom();
      }
   }

   size.fWidth  += fMain->GetBorderWidth() << 1;
   size.fHeight += fMain->GetBorderWidth() << 1;

   if (options & kFixedWidth)  size.fWidth  = msize.fWidth;
   if (options & kFixedHeight) size.fHeight = msize.fHeight;

   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a horizontal layout of all frames in the list.

void TGHorizontalLayout::Layout()
{
   TGFrameElement *ptr;
   TGLayoutHints  *layout;
   Int_t    nb_expand = 0;
   Int_t    left, right;
   ULong_t  hints;
   UInt_t   extra_space = 0;
   Int_t    exp = 0;
   Int_t    exp_max = 0;
   Int_t    remain;
   Int_t    x = 0, y = 0;
   Int_t    bw = fMain->GetBorderWidth();
   TGDimension size, csize;
   TGDimension msize = fMain->GetSize();
   UInt_t pad_left, pad_top, pad_right, pad_bottom;
   Int_t size_expand=0, esize_expand=0, rem_expand=0, tmp_expand = 0;

   if (!fList) return;

   fModified = kFALSE;
   right  = msize.fWidth - (left = bw);
   remain = msize.fWidth - (bw << 1);

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         layout  = ptr->fLayout;
         size    = ptr->fFrame->GetDefaultSize();
         size.fWidth += layout->GetPadLeft() + layout->GetPadRight();
         hints = layout->GetLayoutHints();
         if ((hints & kLHintsExpandX) || (hints & kLHintsCenterX)) {
            nb_expand++;
            exp += size.fWidth;
            if (hints & kLHintsExpandX) exp_max = 0;
            else exp_max = TMath::Max(exp_max, (Int_t)size.fWidth);
         } else {
            remain -= size.fWidth;
            if (remain < 0)
               remain = 0;
         }
      }
   }
   if (nb_expand) {
      size_expand = remain/nb_expand;

      if (size_expand < exp_max) {
         esize_expand = (remain - exp)/nb_expand;
      }
      rem_expand = remain % nb_expand;
   }

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         hints = (layout = ptr->fLayout)->GetLayoutHints();
         csize      = ptr->fFrame->GetDefaultSize();
         pad_left   = layout->GetPadLeft();
         pad_top    = layout->GetPadTop();
         pad_right  = layout->GetPadRight();
         pad_bottom = layout->GetPadBottom();

         if (hints & kLHintsBottom) {
            y = msize.fHeight - bw - csize.fHeight - pad_bottom;
         } else if (hints & kLHintsCenterY) {
            y = (msize.fHeight - (bw << 1) - csize.fHeight) >> 1;
         } else { // kLHintsTop by default
            y = pad_top + bw;
         }

         if (hints & kLHintsExpandY) {
            size.fHeight = msize.fHeight - (bw << 1) - pad_top - pad_bottom;
            y = pad_top + bw;
         } else {
            size.fHeight = csize.fHeight;
         }

         if (hints & kLHintsExpandX) {
            if (size_expand >= exp_max)
               size.fWidth = size_expand - pad_left - pad_right;
            else
               size.fWidth = csize.fWidth +  esize_expand;

            tmp_expand += rem_expand;

            if (tmp_expand >= nb_expand) {
               size.fWidth++;
               tmp_expand -= nb_expand;
            }
         } else {
            size.fWidth = csize.fWidth;
            if (hints & kLHintsCenterX) {
               if (size_expand >= exp_max) {
                  extra_space = (size_expand - pad_left - pad_right - size.fWidth)>> 1;
               } else {
                  extra_space = esize_expand >> 1;
               }
               x += extra_space;
               left += extra_space;
            }
         }

         if (hints & kLHintsRight) {
            x = right - size.fWidth - pad_right;
            right -= size.fWidth + pad_left + pad_right;
         } else { // defaults to kLHintsLeft
            x = left + pad_left;
            left += size.fWidth + pad_left + pad_right;
         }

         if (hints & kLHintsCenterX)
            left += extra_space;

         ptr->fFrame->MoveResize(x, y, size.fWidth, size.fHeight);

         fModified = fModified || (ptr->fFrame->GetX() != x) ||
                    (ptr->fFrame->GetY() != y) ||
                    (ptr->fFrame->GetWidth() != size.fWidth) ||
                    (ptr->fFrame->GetHeight() != size.fHeight);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the horizontal layout.

TGDimension TGHorizontalLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension     size(0,0), msize = fMain->GetSize(), csize;
   UInt_t options = fMain->GetOptions();

   if ((options & kFixedWidth) && (options & kFixedHeight))
      return msize;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         csize = ptr->fFrame->GetDefaultSize();
         size.fWidth += csize.fWidth + ptr->fLayout->GetPadLeft() +
                        ptr->fLayout->GetPadRight();

         size.fHeight = TMath::Max(size.fHeight, csize.fHeight + ptr->fLayout->GetPadTop() +
                                   ptr->fLayout->GetPadBottom());
      }
   }
   size.fWidth  += fMain->GetBorderWidth() << 1;
   size.fHeight += fMain->GetBorderWidth() << 1;

   if (options & kFixedWidth)  size.fWidth = msize.fWidth;
   if (options & kFixedHeight) size.fHeight = msize.fHeight;

   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a row layout of all frames in the list.

void TGRowLayout::Layout()
{
   TGFrameElement *ptr;
   TGDimension     size;
   Int_t  bw = fMain->GetBorderWidth();
   Int_t  x = bw, y = bw;
   fModified = kFALSE;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         size = ptr->fFrame->GetDefaultSize();
         ptr->fFrame->Move(x, y);

         fModified = fModified || (ptr->fFrame->GetX() != x) ||
                    (ptr->fFrame->GetY() != y);

         ptr->fFrame->Layout();
         x += size.fWidth + fSep;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the row layout.

TGDimension TGRowLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension size(0,0), dsize, msize = fMain->GetSize();
   UInt_t options = fMain->GetOptions();

   if ((options & kFixedHeight) && (options & kFixedWidth))
      return msize;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         dsize   = ptr->fFrame->GetDefaultSize();
         size.fHeight  = TMath::Max(size.fHeight, dsize.fHeight);
         size.fWidth  += dsize.fWidth + fSep;
      }
   }

   size.fHeight += fMain->GetBorderWidth() << 1;
   size.fWidth  += fMain->GetBorderWidth() << 1;
   size.fWidth  -= fSep;

   if (options & kFixedHeight) size.fHeight = msize.fHeight;
   if (options & kFixedWidth)  size.fWidth  = msize.fWidth;

   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a column layout of all frames in the list.

void TGColumnLayout::Layout()
{
   TGFrameElement *ptr;
   TGDimension     size;
   Int_t  bw = fMain->GetBorderWidth();
   Int_t  x = bw, y = bw;
   fModified = kFALSE;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         size = ptr->fFrame->GetDefaultSize();
         ptr->fFrame->Move(x, y);
         fModified = fModified || (ptr->fFrame->GetX() != x) ||
                    (ptr->fFrame->GetY() != y);
         ptr->fFrame->Layout();
         y += size.fHeight + fSep;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the column layout.

TGDimension TGColumnLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension     size(0,0), dsize, msize = fMain->GetSize();
   UInt_t options = fMain->GetOptions();

   if (options & kFixedHeight && options & kFixedWidth)
      return msize;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         dsize   = ptr->fFrame->GetDefaultSize();
         size.fHeight += dsize.fHeight + fSep;
         size.fWidth   = TMath::Max(size.fWidth, dsize.fWidth);
      }
   }

   size.fHeight += fMain->GetBorderWidth() << 1;
   size.fHeight -= fSep;
   size.fWidth  += fMain->GetBorderWidth() << 1;

   if (options & kFixedHeight) size.fHeight = msize.fHeight;
   if (options & kFixedWidth)  size.fWidth  = msize.fWidth;

   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// TGMatrixLayout constructor.

TGMatrixLayout::TGMatrixLayout(TGCompositeFrame *main, UInt_t r, UInt_t c,
                               Int_t s, Int_t h)
{
   fMain    = main;
   fList    = fMain->GetList();
   fSep     = s;
   fHints   = h;
   fRows    = r;
   fColumns = c;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a matrix layout of all frames in the list.

void TGMatrixLayout::Layout()
{
   TGFrameElement *ptr;
   TGDimension csize, maxsize(0,0);
   Int_t bw = fMain->GetBorderWidth();
   Int_t x = fSep, y = fSep + bw;
   UInt_t rowcount = fRows, colcount = fColumns;
   fModified = kFALSE;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      csize = ptr->fFrame->GetDefaultSize();
      maxsize.fWidth  = TMath::Max(maxsize.fWidth, csize.fWidth);
      maxsize.fHeight = TMath::Max(maxsize.fHeight, csize.fHeight);
   }

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      ptr->fFrame->Move(x, y);
      fModified = fModified || (ptr->fFrame->GetX() != x) ||
                   (ptr->fFrame->GetY() != y);

      ptr->fFrame->Layout();

      if (fColumns == 0) {
         y += maxsize.fHeight + fSep;
         rowcount--;
         if (rowcount <= 0) {
            rowcount = fRows;
            y = fSep + bw; x += maxsize.fWidth + fSep;
         }
      } else if (fRows == 0) {
         x += maxsize.fWidth + fSep;
         colcount--;
         if (colcount <= 0) {
            colcount = fColumns;
            x = fSep; y += maxsize.fHeight + fSep;
         }
      } else {
         x += maxsize.fWidth + fSep;
         colcount--;
         if (colcount <= 0) {
            rowcount--;
            if (rowcount <= 0) return;
            else {
               colcount = fColumns;
               x = fSep; y += maxsize.fHeight + fSep;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the matrix layout.

TGDimension TGMatrixLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension     size, csize, maxsize(0,0);
   Int_t           count = 0;
   Int_t           bw = fMain->GetBorderWidth();

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      count++;
      csize = ptr->fFrame->GetDefaultSize();
      maxsize.fWidth  = TMath::Max(maxsize.fWidth, csize.fWidth);
      maxsize.fHeight = TMath::Max(maxsize.fHeight, csize.fHeight);
   }

   if (fRows == 0) {
      Int_t rows = (count % fColumns) ? (count / fColumns + 1) : (count / fColumns);
      size.fWidth  = fColumns * (maxsize.fWidth + fSep) + fSep;
      size.fHeight = rows * (maxsize.fHeight + fSep) + fSep + bw;
   } else if (fColumns == 0) {
      Int_t cols = (count % fRows) ? (count / fRows + 1) : (count / fRows);
      size.fWidth  = cols * (maxsize.fWidth + fSep) + fSep;
      size.fHeight = fRows * (maxsize.fHeight + fSep) + fSep + bw;
   } else {
      size.fWidth  = fColumns * (maxsize.fWidth + fSep) + fSep;
      size.fHeight = fRows * (maxsize.fHeight + fSep) + fSep + bw;
   }
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a tile layout manager.

TGTileLayout::TGTileLayout(TGCompositeFrame *main, Int_t sep)
{
   fMain = main;
   fSep  = sep;
   fList = fMain->GetList();
   fModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a tile layout of all frames in the list.

void TGTileLayout::Layout()
{
   TGFrameElement *ptr;
   Int_t   x, y, xw, yw;
   UInt_t  max_width;
   ULong_t hints;
   TGDimension csize, max_osize(0,0), msize = fMain->GetSize();
   fModified = kFALSE;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      csize = ptr->fFrame->GetDefaultSize();
      max_osize.fWidth  = TMath::Max(max_osize.fWidth, csize.fWidth);
      max_osize.fHeight = TMath::Max(max_osize.fHeight, csize.fHeight);
   }

   max_width = TMath::Max(msize.fWidth, max_osize.fWidth + (fSep << 1));
   x = fSep; y = fSep << 1;

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      hints = ptr->fLayout->GetLayoutHints();
      csize = ptr->fFrame->GetDefaultSize();

      if (hints & kLHintsCenterX)
         xw = x + (Int_t)((max_osize.fWidth - csize.fWidth) >> 1);
      else if (hints & kLHintsRight)
         xw = x + (Int_t)max_osize.fWidth - (Int_t)csize.fWidth;
      else // defaults to kLHintsLeft
         xw = x;

      if (hints & kLHintsCenterY)
         yw = y + (Int_t)((max_osize.fHeight - csize.fHeight) >> 1);
      else if (hints & kLHintsBottom)
         yw = y + (Int_t)max_osize.fHeight - (Int_t)csize.fHeight;
      else // defaults to kLHintsTop
         yw = y;

      fModified = fModified || (ptr->fFrame->GetX() != xw) ||
                 (ptr->fFrame->GetY() != yw);
      ptr->fFrame->Move(xw, yw);
      if (hints & kLHintsExpandX)
         ptr->fFrame->Resize(max_osize.fWidth, ptr->fFrame->GetDefaultHeight());
      x += (Int_t)max_osize.fWidth + fSep;

      if (x + max_osize.fWidth > max_width) {
         x = fSep;
         y += (Int_t)max_osize.fHeight + fSep + (fSep >> 1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the tile layout.

TGDimension TGTileLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   Int_t x, y;
   TGDimension max_size, max_osize(0,0), msize = fMain->GetSize();

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      max_size = ptr->fFrame->GetDefaultSize();
      max_osize.fWidth  = TMath::Max(max_osize.fWidth, max_size.fWidth);
      max_osize.fHeight = TMath::Max(max_osize.fHeight, max_size.fHeight);
   }

   max_size.fWidth = TMath::Max(msize.fWidth, max_osize.fWidth + (fSep << 1));

   x = fSep; y = fSep << 1;

   next.Reset();
   // coverity[returned_pointer]
   while ((ptr = (TGFrameElement *) next())) {
      x += max_osize.fWidth + fSep;
      if (x + max_osize.fWidth > max_size.fWidth) {
         x = fSep;
         y += (Int_t)max_osize.fHeight + fSep + (fSep >> 1);  // 3/2
      }
   }
   if (x != fSep) y += max_osize.fHeight + fSep;
   max_size.fHeight = TMath::Max(y, (Int_t)msize.fHeight);

   return max_size;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a tile layout of all frames in the list.

void TGListLayout::Layout()
{
   TGFrameElement *ptr;
   Int_t   x, y, xw, yw;
   UInt_t  max_height;
   ULong_t hints;
   TGDimension csize, max_osize(0,0), msize = fMain->GetSize();
   fModified = kFALSE;

   TIter next(fList);
   // coverity[returned_pointer]
   while ((ptr = (TGFrameElement *) next())) {
      csize = ptr->fFrame->GetDefaultSize();
      max_osize.fWidth  = TMath::Max(max_osize.fWidth, csize.fWidth);
      max_osize.fHeight = TMath::Max(max_osize.fHeight, csize.fHeight);
   }

   max_height = TMath::Max(msize.fHeight, max_osize.fHeight + (fSep << 1));

   x = fSep; y = fSep << 1;

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {

      hints = ptr->fLayout->GetLayoutHints();
      csize = ptr->fFrame->GetDefaultSize();

      if (hints & kLHintsCenterX)
         xw = x + (Int_t)((max_osize.fWidth - csize.fWidth) >> 1);
      else if (hints & kLHintsRight)
         xw = x + (Int_t)max_osize.fWidth - (Int_t)csize.fWidth;
      else // defaults to kLHintsLeft
         xw = x;

      if (hints & kLHintsCenterY)
         yw = y + (Int_t)((max_osize.fHeight - csize.fHeight) >> 1);
      else if (hints & kLHintsBottom)
         yw = y + (Int_t)max_osize.fHeight - (Int_t)csize.fHeight;
      else // defaults to kLHintsTop
         yw = y;

      fModified = fModified || (ptr->fFrame->GetX() != xw) ||
                 (ptr->fFrame->GetY() != yw);
      ptr->fFrame->Move(xw, yw);
      if (hints & kLHintsExpandX)
         ptr->fFrame->Resize(max_osize.fWidth, ptr->fFrame->GetDefaultHeight());
      y += (Int_t)max_osize.fHeight + fSep + (fSep >> 1);

      if (y + max_osize.fHeight > max_height) {
         y = fSep << 1;
         x += (Int_t)max_osize.fWidth + fSep;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the list layout.

TGDimension TGListLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   Int_t x, y;
   TGDimension max_size, max_osize(0,0), msize = fMain->GetSize();

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      max_size = ptr->fFrame->GetDefaultSize();
      max_osize.fWidth  = TMath::Max(max_osize.fWidth, max_size.fWidth);
      max_osize.fHeight = TMath::Max(max_osize.fHeight, max_size.fHeight);
   }

   max_size.fHeight = TMath::Max(msize.fHeight, max_osize.fHeight + (fSep << 1));

   x = fSep; y = fSep << 1;

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      y += (Int_t)max_osize.fHeight + fSep + (fSep >> 1);
      if (y + max_osize.fHeight > max_size.fHeight) {
         y = fSep << 1;
         x += (Int_t)max_osize.fWidth + fSep;
      }
   }
   if (y != (fSep << 1)) x += (Int_t)max_osize.fWidth + fSep;
   max_size.fWidth = TMath::Max(x, (Int_t)msize.fWidth);

   return max_size;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a list details layout of all frames in the list.

void TGListDetailsLayout::Layout()
{
   TGFrameElement *ptr;
   TGDimension     csize, msize = fMain->GetSize();
   Int_t max_oh = 0, x = fSep, y = fSep << 1;
   fModified = kFALSE;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      csize = ptr->fFrame->GetDefaultSize();
      max_oh = TMath::Max(max_oh, (Int_t)csize.fHeight);
   }

   next.Reset();

   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         csize = ptr->fFrame->GetDefaultSize();

         fModified = fModified || (ptr->fFrame->GetX() != x) ||
                     (ptr->fFrame->GetY() != y);

         ptr->fFrame->MoveResize(x, y, msize.fWidth, csize.fHeight);
         ptr->fFrame->Layout();
         y += max_oh + fSep + (fSep >> 1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default dimension of the list details layout.

TGDimension TGListDetailsLayout::GetDefaultSize() const
{
   TGFrameElement *ptr;
   TGDimension csize, max_osize(0,0);
   Int_t y = fSep << 1;

   TIter next(fList);
   while ((ptr = (TGFrameElement *) next())) {
      csize = ptr->fFrame->GetDefaultSize();
      max_osize.fWidth  = TMath::Max(max_osize.fWidth, csize.fWidth);
      max_osize.fHeight = TMath::Max(max_osize.fHeight, csize.fHeight);
   }

   next.Reset();
   while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fState & kIsVisible) {
         y += max_osize.fHeight + fSep + (fSep >> 1);
      }
   }

   return TGDimension( fWidth ? fWidth : max_osize.fWidth, y);
}

// ________________________________________________________________________
void TGLayoutHints::SavePrimitive(std::ostream &out, Option_t * option/*= ""*/)
{

   // Save layout hints as a C++ statement(s) on output stream out

   TString hints;
   UInt_t pad = GetPadLeft()+GetPadRight()+GetPadTop()+GetPadBottom();

   if (!GetLayoutHints()) return;

   if ((option == 0) || strcmp(option, "nocoma"))
      out << ", ";

   if ((fLayoutHints == kLHintsNormal) && (pad == 0)) {
      out << "new TGLayoutHints(kLHintsNormal)";
      return;
   }
   if (fLayoutHints & kLHintsLeft) {
      if (hints.Length() == 0) hints  = "kLHintsLeft";
      else                     hints += " | kLHintsLeft";
   }
   if (fLayoutHints & kLHintsCenterX) {
      if  (hints.Length() == 0) hints  = "kLHintsCenterX";
      else                      hints += " | kLHintsCenterX";
   }
   if (fLayoutHints & kLHintsRight) {
      if (hints.Length() == 0) hints  = "kLHintsRight";
      else                     hints += " | kLHintsRight";
   }
   if (fLayoutHints & kLHintsTop) {
      if (hints.Length() == 0) hints  = "kLHintsTop";
      else                     hints += " | kLHintsTop";
   }
   if (fLayoutHints & kLHintsCenterY) {
      if (hints.Length() == 0) hints  = "kLHintsCenterY";
      else                     hints += " | kLHintsCenterY";
   }
   if (fLayoutHints & kLHintsBottom) {
      if (hints.Length() == 0) hints  = "kLHintsBottom";
      else                     hints += " | kLHintsBottom";
   }
   if (fLayoutHints & kLHintsExpandX) {
      if (hints.Length() == 0) hints  = "kLHintsExpandX";
      else                     hints += " | kLHintsExpandX";
   }
   if (fLayoutHints & kLHintsExpandY) {
      if (hints.Length() == 0) hints  = "kLHintsExpandY";
      else                     hints += " | kLHintsExpandY";
   }

   out << "new TGLayoutHints(" << hints;

   if (pad) {
      out << "," << GetPadLeft() << "," << GetPadRight()
          << "," << GetPadTop()  << "," << GetPadBottom();
   }
   out<< ")";
}

// __________________________________________________________________________
void TGVerticalLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save vertical layout manager as a C++ statement(s) on output stream

   out << "new TGVerticalLayout(" << fMain->GetName() << ")";

}

// __________________________________________________________________________
void TGHorizontalLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save horizontal layout manager as a C++ statement(s) on output stream

   out << "new TGHorizontalLayout(" << fMain->GetName() << ")";
}

// __________________________________________________________________________
void TGRowLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save row layout manager as a C++ statement(s) on output stream

   out << "new TGRowLayout(" << fMain->GetName() << ","
                             << fSep << ")";
}

// __________________________________________________________________________
void TGColumnLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save column layout manager as a C++ statement(s) on output stream

   out << "new TGColumnLayout(" << fMain->GetName() << ","
                                << fSep << ")";

}

// __________________________________________________________________________
void TGMatrixLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save matrix layout manager as a C++ statement(s) on output stream

   out << "new TGMatrixLayout(" << fMain->GetName() << ","
                                << fRows << ","
                                << fColumns << ","
                                << fSep << ","
                                << fHints <<")";

}

// __________________________________________________________________________
void TGTileLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save tile layout manager as a C++ statement(s) on output stream

   out << "new TGTileLayout(" << fMain->GetName() << ","
                              << fSep << ")";

}

// __________________________________________________________________________
void TGListLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save list layout manager as a C++ statement(s) on output stream

   out << "new TGListLayout(" << fMain->GetName() << ","
                              << fSep << ")";

}

// __________________________________________________________________________
void TGListDetailsLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{

   // Save list details layout manager as a C++ statement(s) on out stream

   out << "new TGListDetailsLayout(" << fMain->GetName() << ","
                                     << fSep << "," << fWidth << ")";

}
