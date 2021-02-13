// Author: Roel Aaij   14/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGTableContainer.h"
#include "TGTableCell.h"
#include "TGWindow.h"
#include "TGTable.h"
#include "TVirtualX.h"

ClassImp(TGTableFrame);
ClassImp(TGTableHeaderFrame);

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGTableFrame and TGTableHeaderFrame                                  //
//                                                                      //
// TGTableFrame contains a composite frame that uses a TGMatrixLayout   //
// to Layout the frames it contains.                                    //
//                                                                      //
// TGTableHeaderFrame implements a frame used to display TGTableHeaders //
// in a TGTable.                                                        //
//                                                                      //
// Both classes are for internal use in TGTable only.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create the container used to view TGTableCells. p.

TGTableFrame::TGTableFrame(const TGWindow *p, UInt_t nrows, UInt_t ncolumns)
   : TQObject(), fFrame(0), fCanvas(0)
{
   fFrame = new TGCompositeFrame(p, 10, 10, kHorizontalFrame,
                                 TGFrame::GetWhitePixel());
   fFrame->Connect("ProcessedEvent(Event_t*)", "TGTableFrame", this,
                   "HandleMouseWheel(Event_t*)");
   fCanvas = 0;
   fFrame->SetLayoutManager(new TGMatrixLayout(fFrame, nrows, ncolumns));

   gVirtualX->GrabButton(fFrame->GetId(), kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse wheel to scroll.

void TGTableFrame::HandleMouseWheel(Event_t *event)
{
   if (event->fType != kButtonPress && event->fType != kButtonRelease)
      return;

   Int_t page = 0;
   if (event->fCode == kButton4 || event->fCode == kButton5) {
      if (!fCanvas) return;
      if (fCanvas->GetContainer()->GetHeight())
         page = Int_t(Float_t(fCanvas->GetViewPort()->GetHeight() *
                              fCanvas->GetViewPort()->GetHeight()) /
                              fCanvas->GetContainer()->GetHeight());
   }

   if (event->fCode == kButton4) {
      //scroll up
      Int_t newpos = fCanvas->GetVsbPosition() - page;
      if (newpos < 0) newpos = 0;
      fCanvas->SetVsbPosition(newpos);
   }
   if (event->fCode == kButton5) {
      // scroll down
      Int_t newpos = fCanvas->GetVsbPosition() + page;
      fCanvas->SetVsbPosition(newpos);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a region of container in viewport.

void TGTableFrame::DrawRegion(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   TGFrameElement *el;
   //   Handle_t id = fId;

   Int_t xx = fCanvas->GetX() + fCanvas->GetHsbPosition() + x; // translate coordinates to current page position
   Int_t yy = fCanvas->GetY() + fCanvas->GetVsbPosition() + y;

   TIter next(fFrame->GetList());

   while ((el = (TGFrameElement *) next())) {
      if ((Int_t(el->fFrame->GetY()) >= yy - (Int_t)el->fFrame->GetHeight()) &&
          (Int_t(el->fFrame->GetX()) >= xx - (Int_t)el->fFrame->GetWidth()) &&
          (Int_t(el->fFrame->GetY()) <= yy + Int_t(h + el->fFrame->GetHeight())) &&
          (Int_t(el->fFrame->GetX()) <= xx + Int_t(w + el->fFrame->GetWidth()))) {

         // draw either in container window or in double-buffer
         //          if (!fMapSubwindows) {
         //             el->fFrame->DrawCopy(id, el->fFrame->GetX() - pos.fX, el->fFrame->GetY() - pos.fY);
         //          } else {
         gClient->NeedRedraw(el->fFrame);
         //          }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TGTableHeaderFrame constuctor.

TGTableHeaderFrame::TGTableHeaderFrame(const TGWindow *p, TGTable *table,
                                       UInt_t w, UInt_t h, EHeaderType type,
                                       UInt_t options) :
   TGCompositeFrame(p, w, h, options), fX0(0), fY0(0), fTable(table)
{
   if (type == kRowHeader) {
      ChangeOptions(GetOptions() | kVerticalFrame);
      fY0 = fTable->GetTableHeader()->GetHeight();
   } else if (type == kColumnHeader) {
      ChangeOptions(GetOptions() | kHorizontalFrame);
      fX0 = fTable->GetTableHeader()->GetWidth();
   } else {
      Error("TGTableHeaderFrame::TGTableHeaderFrame",
            "specify correct header type");
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Draw a region of container in viewport.

void TGTableHeaderFrame::DrawRegion(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   TGFrameElement *el;
   //   Handle_t id = fId;

   Int_t xx = fX0 + x; // translate coordinates to current page position
   Int_t yy = fY0 + y;

   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      if ((Int_t(el->fFrame->GetY()) >= yy - (Int_t)el->fFrame->GetHeight()) &&
          (Int_t(el->fFrame->GetX()) >= xx - (Int_t)el->fFrame->GetWidth()) &&
          (Int_t(el->fFrame->GetY()) <= yy + Int_t(h + el->fFrame->GetHeight())) &&
          (Int_t(el->fFrame->GetX()) <= xx + Int_t(w + el->fFrame->GetWidth()))) {

         fClient->NeedRedraw(el->fFrame);
      }
   }
}
