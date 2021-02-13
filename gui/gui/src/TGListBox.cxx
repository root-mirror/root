// @(#)root/gui:$Id$
// Author: Fons Rademakers   12/01/98

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
// TGListBox, TGLBContainer, TGLBEntry and TGTextLBEntry                //
//                                                                      //
// A listbox is a box, possibly with scrollbar, containing entries.     //
// Currently entries are simple text strings (TGTextLBEntry).           //
// A TGListBox looks a lot like a TGCanvas. It has a TGViewPort         //
// containing a TGLBContainer which contains the entries and it also    //
// has a vertical scrollbar which becomes visible if there are more     //
// items than fit in the visible part of the container.                 //
//                                                                      //
// The TGListBox is user callable. The other classes are service        //
// classes of the listbox.                                              //
//                                                                      //
// Selecting an item in the listbox will generate the event:            //
// kC_COMMAND, kCM_LISTBOX, listbox id, item id.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGPicture.h"
#include "TGListBox.h"
#include "TGResourcePool.h"
#include "TSystem.h"
#include "TMath.h"
#include "TVirtualX.h"

#include <cstdlib>
#include <iostream>


const TGFont *TGTextLBEntry::fgDefaultFont = nullptr;
TGGC         *TGTextLBEntry::fgDefaultGC = nullptr;

ClassImp(TGLBEntry);
ClassImp(TGTextLBEntry);
ClassImp(TGLineLBEntry);
ClassImp(TGLBContainer);
ClassImp(TGListBox);

////////////////////////////////////////////////////////////////////////////////
/// Base class entry constructor.

TGLBEntry::TGLBEntry(const TGWindow *p, Int_t id, UInt_t options, Pixel_t back) :
             TGFrame(p, 10, 10, options | kOwnBackground, back)
{
   fActive = kFALSE;
   fEntryId = id;
   fBkcolor = back;
   fEditDisabled = kEditDisable | kEditDisableGrab;

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle active state of listbox entry.

void TGLBEntry::Activate(Bool_t a)
{
   if (fActive == a) return;
   fActive = a;
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle active state of listbox entry.

void TGLBEntry::Toggle()
{
   fActive = !fActive;
   DoRedraw();
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGTextLBEntry                                                        //
//                                                                      //
// Text string listbox entries.                                         //
// A TGTextLBEntry is for TGListBox internal use.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a text listbox entry. The TGString is adopted.

TGTextLBEntry::TGTextLBEntry(const TGWindow *p, TGString *s, Int_t id,
      GContext_t norm, FontStruct_t font, UInt_t options, ULong_t back) :
   TGLBEntry(p, id, options, back)
{
   fText        = s;
   fTextChanged = kTRUE;
   fFontStruct  = font;
   fNormGC      = norm;
   fTWidth      = 0;

   int max_ascent, max_descent;

   if (fText) fTWidth  = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;
   Resize(fTWidth, fTHeight + 1);
   fEditDisabled = kEditDisable | kEditDisableGrab;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete text listbox entry.

TGTextLBEntry::~TGTextLBEntry()
{
   if (fText) delete fText;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text listbox entry on window/pixmap.

void TGTextLBEntry::DrawCopy(Handle_t id, Int_t x, Int_t y)
{
   int max_ascent, max_descent;

   y += (fHeight - fTHeight) >> 1;

   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   if (fActive) {
      gVirtualX->SetForeground(fNormGC, fgDefaultSelectedBackground );
      gVirtualX->FillRectangle(id,fNormGC, x, y, fWidth, fHeight);
      gVirtualX->SetForeground(fNormGC, fClient->GetResourcePool()->GetSelectedFgndColor());
      fText->Draw(id, fNormGC, x + 3, y + max_ascent);
   } else {
      gVirtualX->SetForeground(fNormGC, fBkcolor);
      gVirtualX->FillRectangle(id, fNormGC, x, y, fWidth, fHeight);
      gVirtualX->SetForeground(fNormGC, GetForeground());
      fText->Draw(id, fNormGC, x + 3, y + max_ascent);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw text listbox entry.

void TGTextLBEntry::DoRedraw()
{
   if (fId) DrawCopy(fId, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change text in text entry.

void TGTextLBEntry::SetText(TGString *new_text)
{
   if (fText) delete fText;
   fText = new_text;
   fTextChanged = kTRUE;

   int max_ascent, max_descent;
   fTWidth = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;

   Resize(fTWidth + 3, fTHeight + 1);

   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure in use for a text listbox entry.

FontStruct_t TGTextLBEntry::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context in use for a text listbox entry.

const TGGC &TGTextLBEntry::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = new TGGC(*gClient->GetResourcePool()->GetFrameGC());
   return *fgDefaultGC;
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLineLBEntry                                                        //
//                                                                      //
// Line style and width listbox entries.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create the line style listbox entry.

TGLineLBEntry::TGLineLBEntry(const TGWindow *p, Int_t id, const char *str,
                             UInt_t w, Style_t style, UInt_t options, ULong_t back) :
   TGTextLBEntry(p, new TGString(str), id, GetDefaultGC()(),
                 GetDefaultFontStruct(), options, back)
{
   GCValues_t gcv;

   gcv.fMask =  kGCLineStyle | kGCLineWidth | kGCFillStyle | kGCDashList;
   fLineWidth = gcv.fLineWidth  = w;
   gcv.fFillStyle  = kFillSolid;
   gcv.fDashLen = 2;
   gcv.fDashOffset = 0;
   memcpy(gcv.fDashes, "\x5\x5", 3);
   gcv.fLineStyle = kLineOnOffDash;
   fLineGC = fClient->GetGC(&gcv, kTRUE);
   SetLineStyle(style);

   int max_ascent, max_descent;

   fTWidth  = gVirtualX->TextWidth(GetDefaultFontStruct(), "8", 1);
   fTWidth += 15;                     // for drawing
   gVirtualX->GetFontProperties(GetDefaultFontStruct(),
                                max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;
   fLineLength = 0;

   Resize(fTWidth, fTHeight + 1);
   fEditDisabled = kEditDisable | kEditDisableGrab;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete line style listbox entry.

TGLineLBEntry::~TGLineLBEntry()
{
   fClient->FreeGC(fLineGC);
}

////////////////////////////////////////////////////////////////////////////////
/// Update line style listbox entry.

void  TGLineLBEntry::Update(TGLBEntry *e)
{
   TGTextLBEntry::Update(e);

   fClient->FreeGC(fLineGC);
   fLineGC = ((TGLineLBEntry *)e)->GetLineGC();
   fLineGC->AddReference();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line style corresponding to the TPad line styles.

void TGLineLBEntry::SetLineStyle(Style_t linestyle)
{
   static const char* dashed = "\x3\x3";
   static const char* dotted= "\x1\x2";
   static const char* dasheddotted = "\x3\x4\x1\x4";
   static const char* ls05 = "\x5\x3\x1\x3";
   static const char* ls06 = "\x5\x3\x1\x3\x1\x3\x1\x3";
   static const char* ls07 = "\x5\x5";
   static const char* ls08 = "\x5\x3\x1\x3\x1\x3";
   static const char* ls09 = "\x20\x5";
   static const char* ls10 = "\x20\x10\x1\x10";


   if (linestyle <= 1)  {
      fLineGC->SetLineStyle(kLineSolid);
   } else {
      switch (linestyle) {
         case 2:
            fLineGC->SetDashList(dashed, 2);
            break;
         case 3:
            fLineGC->SetDashList(dotted, 2);
            break;
         case 4:
            fLineGC->SetDashList(dasheddotted, 4);
            break;
         case 5:
            fLineGC->SetDashList(ls05, 4);
            break;
         case 6:
            fLineGC->SetDashList(ls06, 8);
            break;
         case 7:
            fLineGC->SetDashList(ls07, 2);
            break;
         case 8:
            fLineGC->SetDashList(ls08, 6);
            break;
         case 9:
            fLineGC->SetDashList(ls09, 2);
            break;
         case 10:
            fLineGC->SetDashList(ls10, 4);
            break;
      }
   }
   fLineGC->SetCapStyle(0); // flat cap
   fLineStyle = linestyle;
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change line witdh in an entry.

void TGLineLBEntry::SetLineWidth(Int_t width)
{
   fLineWidth = width;
   fLineGC->SetLineWidth(fLineWidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw copy on window/pixmap.

void TGLineLBEntry::DrawCopy(Handle_t id, Int_t x, Int_t y)
{
   TGTextLBEntry::DrawCopy(id, x, y);
   if (!strcmp(TGTextLBEntry::GetTitle(),"None")) return;
   if (fActive) {
      gVirtualX->SetForeground(fLineGC->GetGC(),
                               fClient->GetResourcePool()->GetSelectedFgndColor());
   } else {
      gVirtualX->SetForeground(fLineGC->GetGC(),
                               fClient->GetResourcePool()->GetBlackColor());
   }
   gVirtualX->DrawLine(id, fLineGC->GetGC(), x + fTWidth + 5, y + fHeight/2,
                       x + fWidth - 5, y + fHeight/2);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw line style listbox entry.

void TGLineLBEntry::DoRedraw()
{
   if (fId) DrawCopy(fId, 0, 0);
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGIconLBEntry                                                        //
//                                                                      //
// Icon + text listbox entry.                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create the icon & text listbox entry.

TGIconLBEntry::TGIconLBEntry(const TGWindow *p, Int_t id, const char *str,
                             const TGPicture *pic,
                             UInt_t /*w*/, Style_t /*style*/, UInt_t options, ULong_t back) :
   TGTextLBEntry(p, new TGString(str), id, GetDefaultGC()(),
                 GetDefaultFontStruct(), options, back)
{
   int max_ascent, max_descent;

   fPicture = pic;
   if (fPicture) {
      fTWidth += fPicture->GetWidth() + 4;
      ((TGPicture *)fPicture)->AddReference();
   }
   else
      fTWidth += 20;
   gVirtualX->GetFontProperties(GetDefaultFontStruct(),
                                max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;
   if (fPicture && fPicture->GetHeight() > fTHeight)
      fTHeight = fPicture->GetHeight();

   Resize(fTWidth, fTHeight + 1);
   fEditDisabled = kEditDisable | kEditDisableGrab;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete icon & text listbox entry.

TGIconLBEntry::~TGIconLBEntry()
{
   fClient->FreePicture(fPicture);
}

////////////////////////////////////////////////////////////////////////////////
/// Update icon & text listbox entry.

void  TGIconLBEntry::Update(TGLBEntry *e)
{
   TGTextLBEntry::Update(e);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw copy on window/pixmap.

void TGIconLBEntry::DrawCopy(Handle_t id, Int_t x, Int_t y)
{
   Int_t off_x = 0;
   if (fPicture) {
      fPicture->Draw(id, fNormGC, x + 2, y);
      off_x = fPicture->GetWidth() + 4;
   }
   TGTextLBEntry::DrawCopy(id, x + off_x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw icon & text listbox entry.

void TGIconLBEntry::DoRedraw()
{
   if (fId) DrawCopy(fId, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the icon of listbox entry containing icon & text.

void TGIconLBEntry::SetPicture(const TGPicture *pic)
{
   fClient->FreePicture(fPicture);

   if (pic) ((TGPicture *)pic)->AddReference();

   fPicture   = pic;
}

/////////////////////////////////////////////////////////////////////////////////
class TGLBFrameElement : public TGFrameElement {
public:
   TGLBFrameElement(TGFrame *f, TGLayoutHints *l) : TGFrameElement(f, l) {}
   virtual ~TGLBFrameElement() {}

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *obj) const {
      if (!fFrame->InheritsFrom(TGTextLBEntry::Class())) {
         return 0;
      }
      TGTextLBEntry *f1 = (TGTextLBEntry*)fFrame;
      TGTextLBEntry *f2 = (TGTextLBEntry *) ((TGFrameElement *) obj)->fFrame;


      double d1, d2;
      const char *t1 = f1->GetText()->Data();
      const char *t2 = f2->GetText()->Data();

      if ((d1 = atof(t1)) && (d2 = atof(t2))) {
         return (d1 > d2);
      } else {
         return strcmp(t1, t2);
      }
      return 0;
   }
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLBContainer                                                        //
//                                                                      //
// A Composite frame that contains a list of TGLBEnties.                //
// A TGLBContainer is for TGListBox internal use.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a listbox container.

TGLBContainer::TGLBContainer(const TGWindow *p, UInt_t w, UInt_t h,
                             UInt_t options, ULong_t back) :
   TGContainer(p, w, h, options, back)
{
   fLastActive = 0;
   fMsgWindow  = p;
   fMultiSelect = kFALSE;
   fChangeStatus = kFALSE;
   fListBox = 0;

   SetWindowName();
   fEditDisabled = kEditDisableGrab | kEditDisableBtnEnable | kEditDisableKeyEnable;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the listbox container.

TGLBContainer::~TGLBContainer()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Layout container

void TGLBContainer::Layout()
{
   TGContainer::Layout();
   TGFrame::Resize(fListBox->GetViewPort()->GetWidth(), fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// redraw

void TGLBContainer::DoRedraw()
{
   return TGContainer::DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Add listbox entry with hints to container. To show entry call
/// MapSubwindows() and Layout().

void TGLBContainer::AddEntry(TGLBEntry *lbe, TGLayoutHints *lhints)
{
   // DEPRECATED: the color should always be set in the TGLBEntry ctor
   //lbe->SetBackgroundColor(fgWhitePixel);

   TGLBFrameElement *nw = new TGLBFrameElement(lbe, lhints ? lhints : fgDefaultHints);
   fList->Add(nw);
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert listbox entry after specified entry with id afterID. If afterID = -1
/// then add entry at head of list. To show entry call MapSubwindows() and
/// Layout().

void TGLBContainer::InsertEntry(TGLBEntry *lbe, TGLayoutHints *lhints, Int_t afterID)
{
   // DEPRECATED: the color should always be set in the TGLBEntry ctor
   //lbe->SetBackgroundColor(fgWhitePixel);

   TGLBEntry      *e;
   TGFrameElement *el, *nw;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      e = (TGLBEntry *) el->fFrame;
      if (e->EntryId() == afterID) break;
   }

   if (!el && afterID != -1) {
      nw = new TGLBFrameElement(lbe, lhints ? lhints : fgDefaultHints);
      fList->Add(nw);
   } else {
      nw = new TGLBFrameElement(lbe, lhints);
      nw->fFrame  = lbe;
      nw->fLayout = lhints;
      nw->fState  = 1;
      //lbe->SetFrameElement(nw);

      if (afterID == -1)
         fList->AddFirst(nw);
      else
         fList->AddAfter(el, nw);
   }
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert listbox entry before the list box entry with a higher id.
/// To show entry call MapSubwindows() and Layout().

void TGLBContainer::AddEntrySort(TGLBEntry *lbe, TGLayoutHints *lhints)
{
   // DEPRECATED: the color should always be set in the TGLBEntry ctor
   //lbe->SetBackgroundColor(fgWhitePixel);

   TGLBEntry      *e;
   TGFrameElement *el, *nw;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      e = (TGLBEntry *) el->fFrame;
      if (e->EntryId() > lbe->EntryId()) break;
   }

   if (!el) {
      nw = new TGLBFrameElement(lbe, lhints ? lhints : fgDefaultHints);
      fList->Add(nw);
   } else {
      nw = new TGLBFrameElement(lbe, lhints);
      nw->fFrame  = lbe;
      nw->fLayout = lhints;
      nw->fState  = 1;
      //lbe->SetFrameElement(nw);

      fList->AddBefore(el, nw);
   }
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the entry with specified id from the listbox container.
/// To update the listbox call Layout().

void TGLBContainer::RemoveEntry(Int_t id)
{
   TGLBEntry      *e;
   TGFrameElement *el;
   TGLayoutHints  *l;

   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      e = (TGLBEntry *) el->fFrame;
      l = el->fLayout;
      if (e->EntryId() == id) {
         if (fLastActive == e) fLastActive = 0;
         e->DestroyWindow();
         fList->Remove(el);  // avoid calling RemoveFrame(e)
         delete el;          // item
         delete e;
         delete l;
         break;
      }
   }
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove entries from from_ID to to_ID (including).
/// To update the listbox call Layout().

void TGLBContainer::RemoveEntries(Int_t from_ID, Int_t to_ID)
{
   TGLBEntry      *e;
   TGFrameElement *el;
   TGLayoutHints  *l;

   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      e = (TGLBEntry *) el->fFrame;
      l = el->fLayout;
      if ((e->EntryId() >= from_ID) && (e->EntryId() <= to_ID)) {
         if (fLastActive == e) fLastActive = 0;
         e->DestroyWindow();
         fList->Remove(el);  // avoid calling RemoveFrame(e)
         delete el;          // idem
         delete e;
         delete l;
      }
   }
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all entries in this container.

void TGLBContainer::RemoveAll()
{
   TGLBEntry      *e;
   TGFrameElement *el;
   TGLayoutHints  *l;

   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      e = (TGLBEntry *) el->fFrame;
      l = el->fLayout;
      if (fLastActive == e) fLastActive = 0;
      e->DestroyWindow();
      fList->Remove(el);  // avoid calling RemoveFrame(e)
      delete el;          // item
      delete e;
      delete l;
   }
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Select the entry with the specified id.
/// Returns the selected TGLBEntry.

TGLBEntry *TGLBContainer::Select(Int_t id)
{
   return Select(id, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Select / deselect the entry with the specified id.
/// Returns the selected TGLBEntry.

TGLBEntry *TGLBContainer::Select(Int_t id, Bool_t sel)
{
   TGLBEntry      *f;
   TGFrameElement *el;

   if (!fMultiSelect && fLastActive) {
      fLastActive->Activate(kFALSE);
      fLastActive = 0;
   }

   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      f = (TGLBEntry *) el->fFrame;
      if (f->EntryId() == id) {
         f->Activate(sel);
         if (fMultiSelect == kFALSE && sel == kTRUE) {
            fLastActive = f;
            fLastActiveEl = el;
         }
         ClearViewPort();
         return f;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns id of selected entry. In case of no selected entry or
/// if multi selection is switched on returns -1.

Int_t TGLBContainer::GetSelected() const
{
   if (fLastActive == 0) return -1;
   return fLastActive->EntryId();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTrue if entry id is selected.

Bool_t TGLBContainer::GetSelection(Int_t id)
{
   TGLBEntry     *f;
   TGFrameElement *el;

   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      f = (TGLBEntry *) el->fFrame;
      if (f->EntryId() == id)
         return f->IsActive();
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds all selected entries (TGLBEntry) of the list box into
/// the list selected.

void TGLBContainer::GetSelectedEntries(TList *selected)
{
   TGLBEntry      *f;
   TGFrameElement *el;

   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      f = (TGLBEntry *) el->fFrame;
      if (f->IsActive()) {
         selected->Add(f);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enables and disables multiple selections of entries.

void TGLBContainer::SetMultipleSelections(Bool_t multi)
{
   TGFrameElement *el;

   fMultiSelect = multi;
   if (!fMultiSelect) {
      // deselect all entries
      TIter next(fList);
      while ((el = (TGFrameElement *) next())) {
         ((TGLBEntry *)(el->fFrame))->Activate(kFALSE);
      }
   }
   fLastActive = 0;
   fLastActiveEl = 0;
   ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to vertical scroll bar.

TGVScrollBar *TGLBContainer::GetVScrollbar() const
{
   return fListBox ? fListBox->GetVScrollbar() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set new vertical scroll bar position.

void TGLBContainer::SetVsbPosition(Int_t newPos)
{
   TGVScrollBar *vb = GetVScrollbar();

   if (vb && vb->IsMapped()) {
      vb->SetPosition(newPos);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in the listbox container.

Bool_t TGLBContainer::HandleButton(Event_t *event)
{
   int xf0, yf0, xff, yff;

   TGLBEntry *f;
   TGFrameElement *el;
   TGLBEntry *last = fLastActive;

   TGPosition pos = GetPagePosition();
   Int_t x = pos.fX + event->fX;
   Int_t y = pos.fY + event->fY;
   Bool_t activate = kFALSE;

   // do not handle "context menu button" during guibuilding
   if (fClient->IsEditable() && (event->fCode == kButton3)) {
      return kTRUE;
   }

   TGVScrollBar *vb = GetVScrollbar();

   if ((event->fCode == kButton4) && vb){
      // scroll 2 lines up (a button down is always followed by a button up)
      Int_t newpos = vb->GetPosition() - 1;
      if (newpos < 0) newpos = 0;
      vb->SetPosition(newpos);
      ClearViewPort();
      return kTRUE;
   }
   if ((event->fCode == kButton5) && vb) {
      // scroll 2 lines down (a button down is always followed by a button up)
      Int_t newpos = vb->GetPosition() + 1;
      vb->SetPosition(newpos);
      ClearViewPort();
      return kTRUE;
   }

   gVirtualX->SetInputFocus(fId);

   if (fMultiSelect) {
      if (event->fType == kButtonPress) {
         TIter next(fList);
         while ((el = (TGFrameElement *) next())) {
            f = (TGLBEntry *) el->fFrame;
            xf0 = f->GetX();
            yf0 = f->GetY();
            xff = xf0 + f->GetWidth();
            yff = yf0 + f->GetHeight();

            activate = fMapSubwindows ? (f->GetId() == (Window_t)event->fUser[0]) :
                        (x > xf0) && (x < xff) && (y > yf0) &&  (y < yff);

            if (activate)  {
               fLastActive = f;
               fLastActiveEl = el;
               f->Toggle();
               fChangeStatus = f->IsActive() ? 1 : 0;
               SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_ITEMCLICK),
                           f->EntryId(), 0);
               break;
            }
         }
      } else {
         fChangeStatus = -1;
      }
   } else {
      if (event->fType == kButtonPress) {
         if (fLastActive) {
            fLastActive->Activate(kFALSE);
            fLastActive = 0;
         }
         TIter next(fList);
         while ((el = (TGFrameElement *) next())) {
            f = (TGLBEntry *) el->fFrame;
            xf0 = f->GetX();
            yf0 = f->GetY();
            xff = xf0 + f->GetWidth();
            yff = yf0 + f->GetHeight();

            activate = fMapSubwindows ? (f->GetId() == (Window_t)event->fUser[0]) :
                        (x > xf0) && (x < xff) && (y > yf0) &&  (y < yff) && !f->IsActive();

            if (activate)  {
               f->Activate(kTRUE);
               fLastActive = f;
               fLastActiveEl = el;
            } else {
               f->Activate(kFALSE);
            }
         }
      } else {
         if (fLastActive) {
            f = fLastActive;
            SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_ITEMCLICK),
                        f->EntryId(), 0);
         }
      }
   }
   if (event->fType == kButtonRelease) {
      fScrolling = kFALSE;
      gSystem->RemoveTimer(fScrollTimer);
   }
   if (fChangeStatus || (last != fLastActive))
      ClearViewPort();
   // trick to avoid mouse move events between the mouse click
   // and the unmapping...
   if (fListBox->GetParent()->InheritsFrom("TGComboBoxPopup"))
      fListBox->GetContainer()->RemoveInput(kPointerMotionMask);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click mouse event in the listbox container.

Bool_t TGLBContainer::HandleDoubleClick(Event_t *ev)
{
   if (!fMultiSelect) {
      if (fLastActive) {
         TGLBEntry *f = fLastActive;
         SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_ITEMDBLCLICK),
                     f->EntryId(), 0);
         DoubleClicked(f, ev->fCode);
         DoubleClicked(f, ev->fCode, ev->fXRoot, ev->fYRoot);
      }
      return kTRUE;
   }
   return TGContainer::HandleDoubleClick(ev);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in listbox container.

Bool_t TGLBContainer::HandleMotion(Event_t *event)
{
   int xf0, yf0, xff, yff;

   static Long64_t was = gSystem->Now();
   Long64_t now = gSystem->Now();

   if ((now-was) < 50) return kFALSE;
   was = now;

   TGLBEntry *f;
   TGFrameElement *el;
   TGPosition pos = GetPagePosition();
   TGDimension dim = GetPageDimension();
   Int_t x = pos.fX + event->fX;
   Int_t y = pos.fY + event->fY;
   Bool_t activate = kFALSE;
   TGLBEntry *last = fLastActive;

   if (fMultiSelect) {

      if ((event->fY < 10) || (event->fY > Int_t(dim.fHeight) - 10)) {
         if (!fScrolling) {
            fScrollTimer->Reset();
            gSystem->AddTimer(fScrollTimer);
         }
         fScrolling = kTRUE;
      }
      else if (fChangeStatus >= 0) {
         TIter next(fList);
         while ((el = (TGFrameElement *) next())) {
            f = (TGLBEntry *) el->fFrame;
            xf0 = f->GetX();
            yf0 = f->GetY();
            xff = xf0 + f->GetWidth();
            yff = yf0 + f->GetHeight();
            activate = fMapSubwindows ? (f->GetId() == (Window_t)event->fUser[0]) :
                        (x > xf0) && (x < xff) && (y > yf0) &&  (y < yff);

            if (activate) {
               if (fChangeStatus != (f->IsActive() ? 1 : 0)) {
                  f->Toggle();
                  ClearViewPort();
                  SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_ITEMCLICK),
                              f->EntryId(), 0);
               }
               break;
            }
         }
      }
   } else if (fListBox->GetParent()->InheritsFrom("TGComboBoxPopup")) {
      TIter next(fList);
      while ((el = (TGFrameElement *) next())) {
         f = (TGLBEntry *) el->fFrame;
         xf0 = f->GetX();
         yf0 = f->GetY();
         xff = xf0 + f->GetWidth();
         yff = yf0 + f->GetHeight();

         activate = fMapSubwindows ? (f->GetId() == (Window_t)event->fUser[0]) :
                        (x > xf0) && (x < xff) && (y > yf0) &&  (y < yff)  && !f->IsActive();

         if (activate)  {
            f->Activate(kTRUE);
            fLastActive = f;
            fLastActiveEl = el;
         } else {
            f->Activate(kFALSE);
         }
         if (last != fLastActive) {
            ClearViewPort();
         }
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Autoscroll while close to & beyond  The Wall

void TGLBContainer::OnAutoScroll()
{
   TGFrameElement* el = 0;
   TGLBEntry *f = 0;
   Int_t yf0, yff;
   Bool_t changed = kFALSE;

   TGDimension dim = GetPageDimension();
   TGPosition pos = GetPagePosition();

   Window_t  dum1, dum2;
   Event_t   ev;
   ev.fType  = kButtonPress;
   Int_t     x, y;

   // Where's the cursor?
   gVirtualX->QueryPointer(fId,dum1,dum2,ev.fXRoot,ev.fYRoot,x,y,ev.fState);
   TGVScrollBar *vb = GetVScrollbar();
   if (vb && y > 0 && y < 10) {
      // scroll 1 line up
      Int_t newpos = vb->GetPosition() - 1;
      if (newpos < 0) newpos = 0;
      vb->SetPosition(newpos);
      changed = kTRUE;
   }
   else if (vb && y > (Int_t)dim.fHeight - 10 && y < (Int_t)dim.fHeight) {
      // scroll 1 line down
      Int_t newpos = vb->GetPosition() + 1;
      vb->SetPosition(newpos);
      changed = kTRUE;
   }
   if (changed && fChangeStatus >= 0) {
      pos = GetPagePosition();
      TIter next(fList);
      while ((el = (TGFrameElement *) next())) {
         f = (TGLBEntry *) el->fFrame;
         yf0 = f->GetY();
         yff = yf0 + f->GetHeight();
         if ((y + pos.fY > yf0) && (y + pos.fY < yff)) {
            if (fChangeStatus != (f->IsActive() ? 1 : 0)) {
               f->Toggle();
               ClearViewPort();
               SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_ITEMCLICK),
                           f->EntryId(), 0);
            }
            break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Activate item.

void TGLBContainer::ActivateItem(TGFrameElement *el)
{
   TGContainer::ActivateItem(el);
   fLastActive = (TGLBEntry *)el->fFrame;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the position in the list box of the entry id.
/// The first position has position no 0. Returns -1 if entry id
/// is not in the list of entries.

Int_t TGLBContainer::GetPos(Int_t id)
{
   Int_t          pos = 0;
   TGLBEntry      *f;
   TGFrameElement *el;

   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      f = (TGLBEntry *) el->fFrame;
      if (f->EntryId() == id)
         return pos;
      pos++;
   }

   return -1;
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGListBox                                                            //
//                                                                      //
// A listbox contains a container frame which is viewed through a       //
// viewport. If the container is larger than the viewport than a        //
// vertical scrollbar is added.                                         //
//                                                                      //
// Selecting an item in the listbox will generate the event:            //
// kC_COMMAND, kCM_LISTBOX, listbox id, item id.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a listbox.

TGListBox::TGListBox(const TGWindow *p, Int_t id,
                     UInt_t options, ULong_t back) :
   TGCompositeFrame(p, 10, 10, options, back)
{
   fMsgWindow = p;
   fWidgetId  = id;

   fItemVsize = 1;
   fIntegralHeight = kTRUE;

   InitListBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a listbox widget.

TGListBox::~TGListBox()
{
   if (!MustCleanup()) {
      delete fVScrollbar;
      delete fVport;
      delete fLbc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initiate the internal classes of a list box.

void TGListBox::InitListBox()
{
   fVport = new TGViewPort(this, 6, 6, kChildFrame | kOwnBackground, fgWhitePixel);
   fVScrollbar = new TGVScrollBar(this, kDefaultScrollBarWidth, 6);
   fLbc = new TGLBContainer(fVport, 10, 10, kVerticalFrame, fgWhitePixel);
   fLbc->fViewPort = fVport;
   fLbc->Associate(this);
   fLbc->SetListBox(this);
   SetContainer(fLbc);

   AddFrame(fVport, 0);
   AddFrame(fVScrollbar, 0);

   fVScrollbar->Associate(this);

   fVScrollbar->AddInput(kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask);
   fLbc->RemoveInput(kPointerMotionMask);
   fLbc->AddInput(kButtonPressMask | kButtonReleaseMask | kButtonMotionMask);

   fVport->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fVScrollbar->SetEditDisabled(kEditDisable | kEditDisableGrab | kEditDisableBtnEnable);
   fLbc->SetEditDisabled(kEditDisableGrab | kEditDisableBtnEnable | kEditDisableKeyEnable);
   fEditDisabled = kEditDisableLayout;

   // layout manager is not used
   delete fLayoutManager;
   fLayoutManager = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw borders of the list box widget.

void TGListBox::DrawBorder()
{
   switch (fOptions & (kSunkenFrame | kRaisedFrame | kDoubleBorder)) {

      case kSunkenFrame | kDoubleBorder:
         gVirtualX->DrawLine(fId, GetShadowGC()(), 0, 0, fWidth-2, 0);
         gVirtualX->DrawLine(fId, GetShadowGC()(), 0, 0, 0, fHeight-2);
         gVirtualX->DrawLine(fId, GetBlackGC()(), 1, 1, fWidth-3, 1);
         gVirtualX->DrawLine(fId, GetBlackGC()(), 1, 1, 1, fHeight-3);
         if (gClient->GetStyle() > 1) break;
         gVirtualX->DrawLine(fId, GetHilightGC()(), 0, fHeight-1, fWidth-1, fHeight-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth-1, fHeight-1, fWidth-1, 0);
         gVirtualX->DrawLine(fId, GetBckgndGC()(),  1, fHeight-2, fWidth-2, fHeight-2);
         gVirtualX->DrawLine(fId, GetBckgndGC()(),  fWidth-2, 1, fWidth-2, fHeight-2);
         break;

      default:
         TGCompositeFrame::DrawBorder();
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add entry with specified string and id to listbox. The id will be
/// used in the event processing routine when the item is selected.
/// The string will be adopted by the listbox.

void TGListBox::AddEntry(TGString *s, Int_t id)
{
   TGTextLBEntry *lbe;
   TGLayoutHints *lhints;

   lbe = new TGTextLBEntry(fLbc, s, id);
   lhints = new TGLayoutHints(kLHintsExpandX | kLHintsTop);
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->AddEntry(lbe, lhints);
}

////////////////////////////////////////////////////////////////////////////////
/// Add entry with specified string and id to listbox. The id will be
/// used in the event processing routine when the item is selected.

void TGListBox::AddEntry(const char *s, Int_t id)
{
   AddEntry(new TGString(s), id);
}

////////////////////////////////////////////////////////////////////////////////
/// Add specified TGLBEntry and TGLayoutHints to listbox. The
/// entry and layout will be adopted and later deleted by the listbox.

void TGListBox::AddEntry(TGLBEntry *lbe, TGLayoutHints *lhints)
{
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->AddEntry(lbe, lhints);
}

////////////////////////////////////////////////////////////////////////////////
/// Add entry with specified string and id to listbox sorted by increasing id.
/// This sorting works proberly only if EntrySort functions are used to add
/// entries without mixing them with other add or insert functions.  The id will be
/// used in the event processing routine when the item is selected.
/// The string will be adopted by the listbox.

void TGListBox::AddEntrySort(TGString *s, Int_t id)
{
   TGTextLBEntry *lbe;
   TGLayoutHints *lhints;

   lbe = new TGTextLBEntry(fLbc, s, id);
   lhints = new TGLayoutHints(kLHintsExpandX | kLHintsTop);
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->AddEntrySort(lbe, lhints);
}

////////////////////////////////////////////////////////////////////////////////
/// Add entry with specified string and id to listbox sorted by increasing id.
/// This sorting works proberly only if EntrySort functions are used to add
/// entries without mixing them with other add or insert functions. The id will be
/// used in the event processing routine when the item is selected.

void TGListBox::AddEntrySort(const char *s, Int_t id)
{
   AddEntrySort(new TGString(s), id);
}

////////////////////////////////////////////////////////////////////////////////
/// Add specified TGLBEntry and TGLayoutHints to listbox sorted by increasing id.
/// This sorting works proberly only if EntrySort functions are used to add
/// entries without mixing them with other add or insert functions. The
/// entry and layout will be adopted and later deleted by the listbox.

void TGListBox::AddEntrySort(TGLBEntry *lbe, TGLayoutHints *lhints)
{
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->AddEntrySort(lbe, lhints);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert entry with specified string and id behind the entry with afterID.
/// The string will be adopted and later deleted by the listbox.

void TGListBox::InsertEntry(TGString *s, Int_t id, Int_t afterID)
{
   TGTextLBEntry *lbe;
   TGLayoutHints *lhints;

   lbe = new TGTextLBEntry(fLbc, s, id);
   lhints = new TGLayoutHints(kLHintsExpandX | kLHintsTop);
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->InsertEntry(lbe, lhints, afterID);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert entry with specified string and id behind the entry with afterID.

void TGListBox::InsertEntry(const char *s, Int_t id, Int_t afterID)
{
   InsertEntry(new TGString(s), id, afterID);
}

////////////////////////////////////////////////////////////////////////////////
/// method used to add entry via context menu

void TGListBox::NewEntry(const char *s)
{
   Int_t selected = fLbc->GetSelected();

   // no selected entry or the last entry
   if ((selected < 0) || (selected == GetNumberOfEntries())) {
      AddEntry(s, GetNumberOfEntries()+1);
   } else {
      InsertEntry(s, GetNumberOfEntries()+1, selected);
   }
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// remove entry with id.
/// If id = -1 - the selected entry/entries is/are removed.
///

void TGListBox:: RemoveEntry(Int_t id)
{
   if (id >= 0) {
      fLbc->RemoveEntry(id);
      Layout();
      return;
   }
   if (!fLbc->GetMultipleSelections()) {
      fLbc->RemoveEntry(fLbc->GetSelected());
      Layout();
      return;
   }
   TList li;
   fLbc->GetSelectedEntries(&li);
   TGLBEntry *e;
   TIter next(&li);

   while ((e = (TGLBEntry*)next())) {
      fLbc->RemoveEntry(e->EntryId());
   }
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all entries.

void TGListBox::RemoveAll()
{
   fLbc->RemoveAll();
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a range of entries defined by from_ID and to_ID.

void TGListBox::RemoveEntries(Int_t from_ID, Int_t to_ID)
{
   fLbc->RemoveEntries(from_ID, to_ID);
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert the specified TGLBEntry and layout hints behind afterID.
/// The entry and layout will be adopted and later deleted by the listbox.

void TGListBox::InsertEntry(TGLBEntry *lbe, TGLayoutHints *lhints, int afterID)
{
   fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
   fLbc->InsertEntry(lbe, lhints, afterID);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns list box entry with specified id.

TGLBEntry *TGListBox::GetEntry(Int_t id) const
{
   TIter next(fLbc->GetList());
   TGFrameElement *el;

   while ((el = (TGFrameElement *)next())) {
      TGLBEntry *lbe = (TGLBEntry *)el->fFrame;
      if (lbe->EntryId() == id) return lbe;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the entry with id to the top of the listbox.

void TGListBox::SetTopEntry(Int_t id)
{
   Int_t idPos;

   idPos = fLbc->GetPos(id);

   // id not found in list of entries
   if (idPos < 0)
      return;

   // call layout to define the range of the scroll bars
   Layout();

   // SetPosition will send a message which will handled by
   // the function TGListBox::ProcessMessage. Now ProcessMessage will
   // set the viewport. SetPosition also will check that the idPos is
   // not out of range.
   fVScrollbar->SetPosition(idPos);
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the listbox widget. If fIntegralHeight is true make the height
/// an integral number of the maximum height of a single entry.

void TGListBox::Resize(UInt_t w, UInt_t h)
{
   if (fIntegralHeight)
      h = TMath::Max(fItemVsize, ((h - (fBorderWidth << 1)) / fItemVsize) * fItemVsize)
                     + (fBorderWidth << 1);

   TGCompositeFrame::Resize(w, h);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Move and resize the listbox widget.

void TGListBox::MoveResize(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   if (fIntegralHeight)
      h = TMath::Max(fItemVsize, ((h - (fBorderWidth << 1)) / fItemVsize) * fItemVsize)
                     + (fBorderWidth << 1);
   TGCompositeFrame::MoveResize(x, y, w, h);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default size of listbox widget.

TGDimension TGListBox::GetDefaultSize() const
{
   UInt_t h;

   if (fIntegralHeight)
      h = TMath::Max(fItemVsize, ((fHeight - (fBorderWidth << 1)) / fItemVsize) * fItemVsize)
                     + (fBorderWidth << 1);
   else
      h = fHeight;

   return TGDimension(fWidth, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Layout the listbox components.

void TGListBox::Layout()
{
   TGFrame *container;
   UInt_t   cw, ch, tch;
   Bool_t   need_vsb;

   need_vsb = kFALSE;

   container = fVport->GetContainer();

   // test whether we need vertical scrollbar or not

   cw = fWidth - (fBorderWidth << 1);
   ch = fHeight - (fBorderWidth << 1);

   container->SetWidth(cw);
   container->SetHeight(ch);

   if (container->GetDefaultHeight() > ch) {
      need_vsb = kTRUE;
      cw -= fVScrollbar->GetDefaultWidth();
      if ((Int_t) cw < 0) {
         Warning("Layout", "width would become too small, setting to 10");
         cw = 10;
      }
      container->SetWidth(cw);
   }

   fVport->MoveResize(fBorderWidth, fBorderWidth, cw, ch);
   container->Layout();

   tch = TMath::Max(container->GetDefaultHeight(), ch);
   container->SetHeight(0); // force a resize in TGFrame::Resize
   container->Resize(cw, tch);
   //fVport->SetPos(0, 0);

   if (need_vsb) {
      fVScrollbar->MoveResize(cw+fBorderWidth, fBorderWidth, fVScrollbar->GetDefaultWidth(), ch);
      fVScrollbar->MapWindow();
   } else {
      fVScrollbar->UnmapWindow();
      fVScrollbar->SetPosition(0);
   }

   fVScrollbar->SetRange((Int_t)TMath::Ceil((Double_t)container->GetHeight()/(Double_t)fItemVsize),
                         fVport->GetHeight()/fItemVsize);
   //fClient->NeedRedraw(container);
   ((TGContainer *)container)->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Sort entries by name

void TGListBox::SortByName(Bool_t ascend)
{
   fLbc->GetList()->Sort(ascend);
   Layout();
   fLbc->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Return id of selected listbox item.

Int_t TGListBox::GetSelected() const
{
   TGLBContainer *ct = (TGLBContainer *) fVport->GetContainer();
   return ct->GetSelected();
}

////////////////////////////////////////////////////////////////////////////////
/// Adds all selected entries (TGLBEntry) of the list box into
/// the list selected.

void TGListBox::GetSelectedEntries(TList *selected)
{
   fLbc->GetSelectedEntries(selected);
}

////////////////////////////////////////////////////////////////////////////////
/// Change background to all entries

void TGListBox::ChangeBackground(Pixel_t back)
{
   fBackground = back;

   TIter next(fLbc->GetList());
   TGFrameElement *el;

   while ((el = (TGFrameElement *)next())) {
      TGLBEntry *lbe = (TGLBEntry *)el->fFrame;
      lbe->SetBackgroundColor(back);
   }
   fLbc->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages generated by the listbox container and forward
/// messages to the listbox message handling window.

Bool_t TGListBox::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch (GET_MSG(msg)) {
      case kC_VSCROLL:
         switch (GET_SUBMSG(msg)) {
            case kSB_SLIDERTRACK:
            case kSB_SLIDERPOS:
               fVport->SetVPos(Int_t(-parm1 * fItemVsize));
               break;
         }
         break;

      case kC_CONTAINER:
         switch (GET_SUBMSG(msg)) {
            case kCT_ITEMCLICK:
               {
                  SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_LISTBOX),
                              fWidgetId, parm1);
                  if (GetMultipleSelections()) SelectionChanged();
                  TGLBEntry *entry = GetSelectedEntry();
                  if (entry) {
                     if (entry->InheritsFrom(TGTextLBEntry::Class())) {
                        const char *text;
                        text = ((TGTextLBEntry*)entry)->GetText()->GetString();
                        Selected(text);
                     }
                     Selected(fWidgetId, (Int_t) parm1);
                     Selected((Int_t) parm1);
                  }
               }
               break;
            case kCT_ITEMDBLCLICK:
               if (!GetMultipleSelections()) {
                  TGLBEntry *entry = GetSelectedEntry();
                  if (entry) {
                     if (entry->InheritsFrom(TGTextLBEntry::Class())) {
                        const char *text;
                        text = ((TGTextLBEntry*)entry)->GetText()->GetString();
                        DoubleClicked(text);
                     }
                     DoubleClicked(fWidgetId, (Int_t) parm1);
                     DoubleClicked((Int_t) parm1);
                  }
               }
               break;
         }
         break;

      default:
         break;

   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Selected signal with list box id and entry id.

void TGListBox::Selected(Int_t widgetId, Int_t id)
{
   Long_t args[2];

   args[0] = widgetId;
   args[1] = id;

   Emit("Selected(Int_t,Int_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit DoubleClicked signal with list box id and entry id.

void TGListBox::DoubleClicked(Int_t widgetId, Int_t id)
{
   Long_t args[2];

   args[0] = widgetId;
   args[1] = id;

   Emit("DoubleClicked(Int_t,Int_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Find entry by name.

TGLBEntry *TGListBox::FindEntry(const char *name) const
{
   TList *list = fLbc->GetList();
   TGFrameElement *el = (TGFrameElement *)list->First();
   while (el) {
      if (el->fFrame->GetTitle() == TString(name))
         return (TGLBEntry *)el->fFrame;
      el = (TGFrameElement *)list->After(el);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a list box widget as a C++ statement(s) on output stream out.

void TGListBox::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetWhitePixel()) SaveUserColor(out, option);

   out << std::endl << "   // list box" << std::endl;

   out<<"   TGListBox *";
   out << GetName() << " = new TGListBox(" << fParent->GetName();

   if (fBackground == GetWhitePixel()) {
      if (GetOptions() == (kSunkenFrame | kDoubleBorder)) {
         if (fWidgetId == -1) {
            out <<");" << std::endl;
         } else {
            out << "," << fWidgetId << ");" << std::endl;
         }
      } else {
         out << "," << fWidgetId << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << fWidgetId << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (!fLbc->GetList()) return;

   TGFrameElement *el;
   TIter next(fLbc->GetList());

   while ((el = (TGFrameElement *) next())) {
      out << "   " << GetName() << "->AddEntry(";
      el->fFrame->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }
   out << "   " << GetName() << "->Resize(" << GetWidth() << "," << GetHeight()
       << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a list box entry widget as a C++ statement(s) on output stream out.

void TGTextLBEntry::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   TString content = GetText()->GetString();
   content.ReplaceAll('\\', "\\\\");
   content.ReplaceAll("\"", "\\\"");
   char quote = '"';
   out << quote << content << quote << "," << EntryId();
}
