// Author: Roel Aaij 21/07/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGFrame.h"
#include "TGWidget.h"
#include "TGWindow.h"
#include "TGResourcePool.h"
#include "TGToolTip.h"
#include "TGPicture.h"
#include "TGTable.h"
#include "TVirtualTableInterface.h"
#include "TVirtualX.h"

ClassImp(TGTableCell);

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGTableCell                                                          //
//                                                                      //
// TGTableCell is the class that represents a single cell in a TGTable. //
//                                                                      //
// This class is for internal use in TGTable only.                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

const TGGC *TGTableCell::fgDefaultGC = 0;
const TGFont *TGTableCell::fgDefaultFont = 0;

////////////////////////////////////////////////////////////////////////////////
/// TGTableCell constructor.

TGTableCell::TGTableCell(const TGWindow *p, TGTable *table, TGString *label,
                         UInt_t row, UInt_t column, UInt_t width, UInt_t height,
                         GContext_t norm, FontStruct_t font, UInt_t option,
                         Bool_t resize)
   : TGFrame(p, width, height, option), fTip(0), fReadOnly(kFALSE),
     fEnabled(kTRUE), fTMode(kTextRight | kTextCenterY), fImage(0),
     fFontStruct(font), fHasOwnFont(kFALSE), fColumn(column), fRow(row),
     fTable(table)
{
   if (label) {
      fLabel = label;
   } else {
      fLabel = new TGString("0");
   }

   fNormGC = norm;
   Init(resize);
}

////////////////////////////////////////////////////////////////////////////////
/// TGTableCell constructor

TGTableCell::TGTableCell(const TGWindow *p, TGTable *table, const char *label,
                         UInt_t row, UInt_t column, UInt_t width, UInt_t height,
                         GContext_t norm, FontStruct_t font, UInt_t option,
                         Bool_t resize)
   : TGFrame(p, width, height, option), fTip(0), fReadOnly(kFALSE),
     fEnabled(kTRUE), fTMode(kTextRight | kTextCenterY), fImage(0),
     fFontStruct(font), fHasOwnFont(kFALSE), fColumn(column), fRow(row),
     fTable(table)
{
   if (label) {
      fLabel = new TGString(label);
   } else {
      fLabel = new TGString("0");
   }

   fNormGC = norm;
   Init(resize);
}

// //______________________________________________________________________________
// TGTableCell::TGTableCell(const TGWindow *p, TGTable *table, TGPicture *image,
//                          TGString *label, UInt_t row, UInt_t column,
//                          GContext_t norm, FontStruct_t font, UInt_t option,
//                          Bool_t resize)
//    : TGFrame(p, 80, 25, option), fTip(0), fReadOnly(kFALSE), fEnabled(kTRUE),
//      fTMode(kTextRight | kTextCenterY), fImage(image), fFontStruct(font),
//      fHasOwnFont(kFALSE), fColumn(column), fRow(row), fTable(table)
// {
//    if (label) {
//       fLabel = label;
//    } else {
//       fLabel = new TGString("0");
//    }
//    fNormGC = norm;

//    Init(resize);
// }

// //______________________________________________________________________________
// TGTableCell::TGTableCell(const TGWindow *p, TGTable *table, TGPicture *image,
//                          const char *label, UInt_t row, UInt_t column,
//                          GContext_t norm, FontStruct_t font, UInt_t option,
//                          Bool_t resize)
//    : TGFrame(p, 80, 25, option), fTip(0), fReadOnly(kFALSE), fEnabled(kTRUE),
//      fTMode(kTextRight | kTextCenterY), fImage(image), fFontStruct(font),
//      fHasOwnFont(kFALSE), fColumn(column), fRow(row), fTable(table)
// {

//    if (label) {
//       fLabel = new TGString(label);
//    } else {
//       fLabel = new TGString("0");
//    }

//    fNormGC = norm;

//    Init(resize);
// }


////////////////////////////////////////////////////////////////////////////////
/// TGTableCell destructor.

TGTableCell::~TGTableCell()
{
   if (fImage) delete fImage;
   if (fLabel) delete fLabel;
   if (fTip) delete fTip;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialise the TGTableCell.

void TGTableCell::Init(Bool_t resize)
{
   Int_t max_ascent = 0, max_descent = 0;

   fTWidth  = gVirtualX->TextWidth(fFontStruct, fLabel->GetString(), fLabel->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;

   // Modifications for individual cell drawing test, original block is marked

   if (fTable) {
      // Original from here
      TGTableHeader *chdr = 0;
      TGTableHeader *rhdr = 0;
      if(resize) {
         chdr = fTable->GetColumnHeader(fColumn);
         rhdr = fTable->GetRowHeader(fRow);
         if (rhdr) {
            SetBackgroundColor(rhdr->GetBackground());
            if (chdr) Resize(chdr->GetWidth(), rhdr->GetHeight());
         }
      }
      SetBackgroundColor(fTable->GetRowBackground(fRow));
      // Upto here
   } else {
      fWidth = 80;
      fHeight = 25;
      Resize(fWidth, fHeight);
      SetBackgroundColor(fgWhitePixel);
   }
   // End of modifications

}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the TGTableCell.

void TGTableCell::DoRedraw()
{
   TGFrame::DoRedraw();

   Int_t x = 0, y = 0;

   // To be done: Add a tooltip with the complete label when it
   // doesn't fit in the cell.
   if (fTWidth > fWidth - 4) fTMode = kTextLeft;

   if (fTMode & kTextLeft) {
      x = 4;
   } else if (fTMode & kTextRight) {
      x = fWidth - fTWidth - 4;
   } else {
      x = (fWidth - fTWidth) / 2;
   }

   if (fTMode & kTextTop) {
      y = 3;
   } else if (fTMode & kTextBottom) {
      y = fHeight - fTHeight - 3;
   } else {
      y = (fHeight - fTHeight - 4) / 2;
   }

   y += fTHeight;

   fLabel->Draw(fId, fNormGC, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Move the TGTableCell and redraw it.

void TGTableCell::MoveDraw(Int_t x, Int_t y)
{
   // Note, this method is unused.

   TGFrame::Move(x, y);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the TGTableCell.

void TGTableCell::Resize(UInt_t width, UInt_t height)
{
   fWidth = width;
   fHeight = height;
   TGFrame::Resize(width, height);
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the TGTableCell.

void TGTableCell::Resize(TGDimension newsize)
{
   Resize(newsize.fWidth, newsize.fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure.

FontStruct_t TGTableCell::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGTableCell::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the cell border.

void TGTableCell::DrawBorder()
{
   // FIXME Borders are drawn very crudely. There is much room for improvement.
   gVirtualX->DrawRectangle(fId, fNormGC, 0, 0, fWidth - 1, fHeight - 1);
}

////////////////////////////////////////////////////////////////////////////////
/// DrawBorder called from DrawCopy.

void TGTableCell::DrawBorder(Handle_t id, Int_t x, Int_t y)
{
   gVirtualX->DrawRectangle(id, fNormGC, x, y, x + fWidth - 1, y +fHeight - 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Highlight the this TGTableCell.

void TGTableCell::Highlight()
{
   // Currently not implemented.
}

// //______________________________________________________________________________
// void TGTableCell::SetRow(UInt_t row)
// {
//    fRow = row;
// }

// //______________________________________________________________________________
// void TGTableCell::SetColumn(UInt_t column)
// {
//    fColumn = column;
// }

////////////////////////////////////////////////////////////////////////////////
/// Draw list view item in other window.
/// List view item is placed and layout in the container frame,
/// but is drawn in viewport.

void TGTableCell::DrawCopy(Handle_t id, Int_t x, Int_t y)
{
   // FIXME this method is only needed if the table container frame is a
   // TGContainer. It is broken and not used in the current implementation.

   Int_t lx = 0, ly = 0;

   if (fTMode & kTextLeft) {
      lx = 4;
   } else if (fTMode & kTextRight) {
      lx = fWidth - fTWidth - 4;
   } else {
      lx = (fWidth - fTWidth) / 2;
   }

   if (fTMode & kTextTop) {
      ly = 3;
   } else if (fTMode & kTextBottom) {
      ly = fHeight - fTHeight - 3;
   } else {
      ly = (fHeight - fTHeight - 4) / 2;
   }

   ly += fTHeight;

   //    if (fActive) {
   //       gVirtualX->SetForeground(fNormGC, fgDefaultSelectedBackground);
   //       gVirtualX->FillRectangle(id, fNormGC, x + lx, y + ly, fTWidth, fTHeight + 1);
   //       gVirtualX->SetForeground(fNormGC, fClient->GetResourcePool()->GetSelectedFgndColor());
   //    } else {

   gVirtualX->SetForeground(fNormGC, fgWhitePixel);
   gVirtualX->FillRectangle(id, fNormGC, x, y, fWidth, fHeight);
   gVirtualX->SetForeground(fNormGC, fgBlackPixel);
   DrawBorder(id, x, y);

   //    }

   fLabel->Draw(id, fNormGC, x + lx, y + ly);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the label of this cell to label.

void TGTableCell::SetLabel(const char *label)
{
   fLabel->SetString(label);

   Int_t max_ascent = 0, max_descent = 0;

   fTWidth  = gVirtualX->TextWidth(fFontStruct, fLabel->GetString(), fLabel->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;

}

////////////////////////////////////////////////////////////////////////////////
/// Set the image that this cell contains to image.

void TGTableCell::SetImage(TGPicture *image)
{
   // Note: currently not used.
   if (fImage) delete fImage;
   fImage = image;
}

// //______________________________________________________________________________
// void TGTableCell::SetBckgndGC(TGGC *gc)
// {
// }

////////////////////////////////////////////////////////////////////////////////
/// Changes text font.
/// If global is kTRUE font is changed globally, otherwise - locally.

void TGTableCell::SetFont(FontStruct_t font)
{
   if (font != fFontStruct) {
      FontH_t v = gVirtualX->GetFontHandle(font);
      if (!v) return;

      fFontStruct = font;
      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);

      gc = pool->GetGC((GCValues_t*)gc->GetAttributes(), kTRUE); // copy
      fHasOwnFont = kTRUE;

      gc->SetFont(v);

      fNormGC = gc->GetGC();
      gClient->NeedRedraw(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font specified by name.
/// If global is true color is changed globally, otherwise - locally.

void TGTableCell::SetFont(const char *fontName)
{
   TGFont *font = fClient->GetFont(fontName);
   if (font) {
      SetFont(font->GetFontStruct());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the text justify mode of the cell to mode.

void TGTableCell::SetTextJustify(Int_t tmode)
{
   fTMode = tmode;
}

////////////////////////////////////////////////////////////////////////////////
/// Select this TGTableCell.

void TGTableCell::Select()
{
   // Note: currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////
/// Select the row that this TGTableCell belongs to.

void TGTableCell::SelectRow()
{
   // Note: currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////
/// Select the column that this TGTableCell belongs to.

void TGTableCell::SelectColumn()
{
   // Note: currently not implemented.
}
