// @(#)root/gui:$Id$
// Author: Fons Rademakers   13/01/98

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
// TGTab, TGTabElement, TGTabLayout                                     //
//                                                                      //
// A tab widget contains a set of composite frames each with a little   //
// tab with a name (like a set of folders with tabs).                   //
//                                                                      //
// The TGTab is user callable. The TGTabElement and TGTabLayout are     //
// service classes of the tab widget.                                   //
//                                                                      //
// Clicking on a tab will bring the associated composite frame to the   //
// front and generate the following event:                              //
// kC_COMMAND, kCM_TAB, tab id, 0.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGTab.h"
#include "TGResourcePool.h"
#include "TList.h"
#include "TClass.h"
#include "TGPicture.h"
#include "TVirtualX.h"

#include <iostream>

const TGFont *TGTab::fgDefaultFont = nullptr;
const TGGC   *TGTab::fgDefaultGC = nullptr;

ClassImp(TGTabElement);
ClassImp(TGTabLayout);
ClassImp(TGTab);


////////////////////////////////////////////////////////////////////////////////
/// Create a tab element. Text is adopted by tab element.

TGTabElement::TGTabElement(const TGWindow *p, TGString *text, UInt_t w, UInt_t h,
                           GContext_t norm, FontStruct_t font,
                           UInt_t options, ULong_t back) :
   TGFrame(p, w, h, options, back)
{
   fClosePic     = 0;
   fClosePicD    = 0;
   fShowClose    = kFALSE;
   fActive       = kFALSE;
   fText         = text;
   fBorderWidth  = 0;
   fTWidth       = 0;
   fNormGC       = norm;
   fFontStruct   = font;
   fEditDisabled = kEditDisableGrab | kEditDisableBtnEnable;

   fClosePic      = fClient->GetPicture("closetab.png");
   fClosePicD     = fClient->GetPicture("closetab_d.png");
   int max_ascent, max_descent;
   if (fText)
      fTWidth = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;
   Resize(TMath::Max(fTWidth+12, (UInt_t)45), fTHeight+6);
   fEnabled = kTRUE;
   gVirtualX->GrabButton(fId, kButton1, kAnyModifier, kButtonPressMask, kNone, kNone);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete tab element.

TGTabElement::~TGTabElement()
{
   if (fClosePic) gClient->FreePicture(fClosePic);
   if (fClosePicD) gClient->FreePicture(fClosePicD);
   if (fText) delete fText;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw little tab element.

void TGTabElement::DrawBorder()
{
   gVirtualX->DrawLine(fId, GetHilightGC()(), 0, fHeight-1, 0, 2);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 0, 2, 2, 0);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 2, 0, fWidth-3, 0);
   gVirtualX->DrawLine(fId, GetShadowGC()(),  fWidth-2, 1, fWidth-2, fHeight-1);
   if (gClient->GetStyle() < 2) {
      gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth-2, 1, fWidth-1, 2);
      gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth-1, 2, fWidth-1, fHeight-2);
   }
   gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth-1, fHeight-1, fWidth-1, fHeight-1);

   if (fText) {
      int max_ascent, max_descent;
      gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
      if (fEnabled) {
         fText->Draw(fId, fNormGC, 6, max_ascent+3);
      } else {
         fText->Draw(fId, GetHilightGC()(), 7, max_ascent + 1);
         fText->Draw(fId, GetShadowGC()(), 6, max_ascent);
      }
   }
   if (fShowClose && fClosePic && fClosePicD) {
      if (fEnabled && fActive)
         fClosePic->Draw(fId, fNormGC, fTWidth+12, fHeight/2-7);
      else
         fClosePicD->Draw(fId, fNormGC, fTWidth+12, fHeight/2-7);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button event in the tab widget. Basically we only handle
/// button events in the small tabs.

Bool_t TGTabElement::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress) {
      TGTab* main = (TGTab*)fParent;
      if (main) {
         if (fShowClose && event->fWindow == GetId() &&
             (UInt_t)event->fX > fTWidth+12 && (UInt_t)event->fX < fTWidth+26 &&
             (UInt_t)event->fY > fHeight/2-7 && (UInt_t)event->fY < fHeight/2+7) {
            if (main->GetTabTab(main->GetCurrent()) == this) {
               main->CloseTab(main->GetCurrent()); // emit signal
               //main->RemoveTab(main->GetCurrent());
               return kTRUE;
            }
         }
         TGFrameElement *el;
         TIter next(main->GetList());

         next();   // skip first container

         Int_t i = 0;
         Int_t c = main->GetCurrent();
         while ((el = (TGFrameElement *) next())) {
            if (el->fFrame->GetId() == (Window_t)event->fWindow)
               c = i;
            next(); i++;
         }

         // change tab and generate event
         main->SetTab(c);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default size of tab element.

TGDimension TGTabElement::GetDefaultSize() const
{
   if (fShowClose && fClosePic && fClosePicD)
      return TGDimension(TMath::Max(fTWidth+30, (UInt_t)45), fTHeight+6);
   else
      return TGDimension(TMath::Max(fTWidth+12, (UInt_t)45), fTHeight+6);
}

////////////////////////////////////////////////////////////////////////////////
/// Set new tab text.

void TGTabElement::SetText(TGString *text)
{
   if (fText) delete fText;
   fText = text;

   int max_ascent, max_descent;
   fTWidth = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Show/hide close icon on the tab element, then apply layout
/// to compute correct elements size.

void TGTabElement::ShowClose(Bool_t show)
{
   TGTab* main = (TGTab*)fParent;
   fShowClose = show;
   if (fShowClose && fClosePic && fClosePicD)
      Resize(TMath::Max(fTWidth+30, (UInt_t)45), fTHeight+6);
   else
      Resize(TMath::Max(fTWidth+12, (UInt_t)45), fTHeight+6);
   if (main)
      main->GetLayoutManager()->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a tab layout manager.

TGTabLayout::TGTabLayout(TGTab *main)
{
   fMain = main;
   fList = fMain->GetList();
}

////////////////////////////////////////////////////////////////////////////////
/// Layout the tab widget.

void TGTabLayout::Layout()
{
   Int_t  i, xtab;
   UInt_t tw;
   UInt_t tabh = fMain->GetTabHeight(), bw = fMain->GetBorderWidth();
   UInt_t w = fMain->GetWidth();
   UInt_t h = fMain->GetHeight();

   xtab = 2;

   fMain->GetContainer()->MoveResize(0, tabh, w, h - tabh);

   // first frame is the container, so take next...
   TGFrameElement *el, *elnxt;
   TIter next(fList);
   i = 0;
   next();   // skip first
   while ((el = (TGFrameElement *) next())) {
      elnxt = (TGFrameElement *) next();
      tw = el->fFrame->GetDefaultWidth();
      if (i == fMain->GetCurrent()) {
         el->fFrame->MoveResize(xtab-2, 0, tw+3, tabh+1);
         if (elnxt) elnxt->fFrame->RaiseWindow();
         el->fFrame->RaiseWindow();
      } else {
         el->fFrame->MoveResize(xtab, 2, tw, tabh-1);
         el->fFrame->LowerWindow();
      }
      UInt_t nw = (w - (bw << 1));
      if (nw > 32768) nw = 1;
      UInt_t nh = (h - tabh - (bw << 1));
      if (nh > 32768) nh = 1;
      if (elnxt) {
         elnxt->fFrame->MoveResize(bw, tabh + bw, nw, nh);
         elnxt->fFrame->Layout();
      }
      xtab += (Int_t)tw;
      i++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get default size of tab widget.

TGDimension TGTabLayout::GetDefaultSize() const
{
   TGDimension dsize, dsize_te;
   TGDimension size(0,0), size_te(0,0);

   TGFrameElement *el, *elnxt;
   TIter next(fList);
   next();   // skip first container
   while ((el = (TGFrameElement *)next())) {
      dsize_te = el->fFrame->GetDefaultSize();
      size_te.fWidth += dsize_te.fWidth;
      elnxt = (TGFrameElement *) next();
      if (elnxt) {
         dsize = elnxt->fFrame->GetDefaultSize();
         if (size.fWidth < dsize.fWidth) size.fWidth = dsize.fWidth;
         if (size.fHeight < dsize.fHeight) size.fHeight = dsize.fHeight;
      }
   }

   // check if tab elements make a larger width than the containers
   if (size.fWidth < size_te.fWidth) size.fWidth = size_te.fWidth;

   size.fWidth += fMain->GetBorderWidth() << 1;
   size.fHeight += fMain->GetTabHeight() + (fMain->GetBorderWidth() << 1);

   return size;
}


////////////////////////////////////////////////////////////////////////////////
/// Create tab widget.

TGTab::TGTab(const TGWindow *p, UInt_t w, UInt_t h,
             GContext_t norm, FontStruct_t font,
             UInt_t options, ULong_t back) :
   TGCompositeFrame(p, w, h, options, back)
{
   fMsgWindow  = p;

   fBorderWidth = 2;
   fCurrent     = 0;
   fRemoved     = new TList;

   fNormGC     = norm;
   fFontStruct = font;

   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTabh = max_ascent + max_descent + 6;

   SetLayoutManager(new TGTabLayout(this));

   // we need this in order to avoid border blinking when switching tabs...
   fContainer = new TGCompositeFrame(this, fWidth, fHeight - fTabh,
                       kVerticalFrame | kRaisedFrame | kDoubleBorder);
   AddFrame(fContainer, 0);

   fEditDisabled = kEditDisable | kEditDisableLayout;
   fContainer->SetEditDisabled(kEditDisable | kEditDisableGrab);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete tab widget. This deletes the tab windows and the containers.
/// The tab string is deleted by the TGTabElement dtor.

TGTab::~TGTab()
{
   Cleanup();
   fRemoved->Delete();
   delete fRemoved;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a tab to the tab widget. Returns the new container, which
/// is owned by the tab widget. The text is adopted by the tab widget.

TGCompositeFrame *TGTab::AddTab(TGString *text)
{
   TGTabElement *te = new TGTabElement(this, text, 50, 20, fNormGC, fFontStruct);
   AddFrame(te, 0);

   TGCompositeFrame* cf = new TGCompositeFrame(this, fWidth, fHeight-21);
   AddFrame(cf, 0);
   cf->SetEditDisabled(kEditDisableResize);

   te->MapWindow();
   cf->MapWindow();

   return cf;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a tab to the tab widget. Returns the new container. The container
/// is owned by the tab widget.

TGCompositeFrame *TGTab::AddTab(const char *text)
{
   return AddTab(new TGString(text));
}

////////////////////////////////////////////////////////////////////////////////
/// Add a tab to the tab widget and fill it with given TGCompositeFrame.

void TGTab::AddTab(const char *text, TGCompositeFrame *cf)
{
   AddTab(new TGString(text), cf);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a tab to the tab widget and fill it with given TGCompositeFrame.

void TGTab::AddTab(TGString *text, TGCompositeFrame *cf)
{
   TGTabElement *te = new TGTabElement(this, text, 50, 20, fNormGC, fFontStruct);
   AddFrame(te, 0);

   AddFrame(cf, 0);
   cf->SetEditDisabled(kEditDisableResize);

   te->MapWindow();
   cf->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove container and tab of tab with index tabIndex.
/// Does NOT remove the container contents!

void TGTab::RemoveTab(Int_t tabIndex, Bool_t storeRemoved)
{
   if (tabIndex < 0) {
      tabIndex = fCurrent;
   }

   TGFrameElement *elTab, *elCont;
   Int_t  count = 0;

   // Notify (signal) for removed tab "tabIndex"
   Removed(tabIndex);

   TIter next(fList) ;
   next() ; // skip first container

   while ((elTab = (TGFrameElement *) next())) {
      elCont = (TGFrameElement *) next();

      if (count == tabIndex) {
         elCont->fFrame->UnmapWindow();   // will be destroyed later
         TGFrame *frame = elTab->fFrame;
         RemoveFrame(elTab->fFrame);
         frame->DestroyWindow();
         delete frame;
         if (storeRemoved)
            fRemoved->Add(elCont->fFrame);   // delete only in dtor
         RemoveFrame(elCont->fFrame);
         if (tabIndex == fCurrent) {
            // select another tab only if the current is the one we delete
            SetTab(0);
         } else
            fCurrent--;
         break;
      }
      count++;
   }

   GetLayoutManager()->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Enable or disable tab.

void TGTab::SetEnabled(Int_t tabIndex, Bool_t on)
{
   TGTabElement *te = GetTabTab(tabIndex);
   if (te) {
      te->SetEnabled(on);
      fClient->NeedRedraw(te);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if tab is enabled.

Bool_t TGTab::IsEnabled(Int_t tabIndex) const
{
   TGTabElement *te = GetTabTab(tabIndex);

   return te ? te->IsEnabled() : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make tabIdx the current tab. Utility method called by SetTab and
/// HandleButton().

void TGTab::ChangeTab(Int_t tabIndex, Bool_t emit)
{
   TGTabElement *te = GetTabTab(tabIndex);
   if (!te || !te->IsEnabled()) return;

   if (tabIndex != fCurrent) {
      if (GetTabTab(fCurrent)) {
         GetTabTab(fCurrent)->SetActive(kFALSE);
         fClient->NeedRedraw(GetTabTab(fCurrent));
      }
      TGFrameElement *el, *elnxt;
      UInt_t tw;
      Int_t  xtab  = 2;
      Int_t  count = 0;

      TIter next(fList);
      next();           // skip first container

      fCurrent = tabIndex;
      while ((el = (TGFrameElement *) next())) {
         elnxt = (TGFrameElement *) next();
         tw = el->fFrame->GetDefaultWidth();
         if (count == fCurrent) {
            el->fFrame->MoveResize(xtab-2, 0, tw+3, fTabh+1);
            if (elnxt) elnxt->fFrame->RaiseWindow();
            el->fFrame->RaiseWindow();
         } else {
            el->fFrame->MoveResize(xtab, 2, tw, fTabh-1);
            el->fFrame->LowerWindow();
         }
         xtab += tw;
         count++;
      }
      if (emit) {
         SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_TAB), fCurrent, 0);
         fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_TAB), fCurrent, 0);
         Selected(fCurrent);
      }
      GetTabTab(fCurrent)->SetActive(kTRUE);
      fClient->NeedRedraw(GetTabTab(fCurrent));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Brings the composite frame with the index tabIndex to the
/// front and generate the following event if the front tab has changed:
/// kC_COMMAND, kCM_TAB, tab id, 0.
/// Returns kFALSE if tabIndex is a not valid index

Bool_t TGTab::SetTab(Int_t tabIndex, Bool_t emit)
{
   // check if tabIndex is a valid index
   if (tabIndex < 0)
      return kFALSE;

   // count the tabs
   TIter next(fList);
   Int_t count = 0;
   while (next())
      count++;

   count = count / 2 - 1;
   if (tabIndex > count)
      return kFALSE;

   // change tab and generate event
   ChangeTab(tabIndex, emit);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Brings the composite frame with the name to the
/// front and generate the following event if the front tab has changed:
/// kC_COMMAND, kCM_TAB, tab id, 0.
/// Returns kFALSE if tab with name does not exist.

Bool_t TGTab::SetTab(const char *name, Bool_t emit)
{
   TGFrameElement *el;
   Int_t  count = 0;
   TGTabElement *tab = 0;

   TIter next(fList);
   next();           // skip first container

   while ((el = (TGFrameElement *) next())) {
      next();        // skip tab containter
      tab = (TGTabElement *)el->fFrame;

      if (*(tab->GetText()) == name) {
         // change tab and generate event
         ChangeTab(count, emit);
         return kTRUE;
      }
      count++;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return container of tab with index tabIndex.
/// Return 0 in case tabIndex is out of range.

TGCompositeFrame *TGTab::GetTabContainer(Int_t tabIndex) const
{
   if (tabIndex < 0) return 0;

   TGFrameElement *el;
   Int_t  count = 0;

   TIter next(fList);
   next();           // skip first container

   while (next()) {
      el = (TGFrameElement *) next();
      if (el && count == tabIndex)
         return (TGCompositeFrame *) el->fFrame;
      count++;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the tab container of tab with string name.
/// Returns 0 in case name is not found.

TGCompositeFrame *TGTab::GetTabContainer(const char *name) const
{
   TGFrameElement *el;
   TGTabElement *tab = 0;
   TGCompositeFrame *comp = 0;

   TIter next(fList);
   next();

   while ((el = (TGFrameElement *) next())) {
      tab  = (TGTabElement *) el->fFrame;
      el   = (TGFrameElement *) next();
      comp = (TGCompositeFrame *) el->fFrame;
      if (*tab->GetText() == name){
         return comp;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the tab element of tab with index tabIndex.
/// Returns 0 in case tabIndex is out of range.

TGTabElement *TGTab::GetTabTab(Int_t tabIndex) const
{
   if (tabIndex < 0) return 0;

   TGFrameElement *el;
   Int_t  count = 0;

   TIter next(fList);
   next();           // skip first container

   while ((el = (TGFrameElement *) next())) {
      next();
      if (count == tabIndex)
         return (TGTabElement *) el->fFrame;
      count++;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the tab element of tab with string name.
/// Returns 0 in case name is not found.

TGTabElement *TGTab::GetTabTab(const char *name) const
{
   TGFrameElement *el;
   TGTabElement *tab = 0;

   TIter next(fList);
   next();

   while ((el = (TGFrameElement *) next())) {
      tab = (TGTabElement *)el->fFrame;
      if (name == *(tab->GetText())) {
         return tab;
      }
      next();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of tabs.

Int_t TGTab::GetNumberOfTabs() const
{
   Int_t count = 0;

   TIter next(fList);
   next();           // skip first container

   while (next()) {
      next();
      count++;
   }

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure in use.

FontStruct_t TGTab::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context in use.

const TGGC &TGTab::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Create new tab. Used in context menu.

void TGTab::NewTab(const char *text)
{
   TString name;
   if (text)
      name = text;
   else
      name = TString::Format("tab%d", GetNumberOfTabs()+1);
   AddTab(name.Data());

   GetLayoutManager()->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set text to current tab.

void TGTab::SetText(const char *text)
{
   if (GetCurrentTab()) GetCurrentTab()->SetText(new TGString(text));
   GetLayoutManager()->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Return layout manager.

TGLayoutManager *TGTab::GetLayoutManager() const
{
   TGTab *tab = (TGTab*)this;

   if (tab->fLayoutManager->IsA() != TGTabLayout::Class()) {
      tab->SetLayoutManager(new TGTabLayout(tab));
   }

   return tab->fLayoutManager;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a tab widget as a C++ statement(s) on output stream out.

void TGTab::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   // font + GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC, parFont;
   parFont.Form("%s::GetDefaultFontStruct()",IsA()->GetName());
   parGC.Form("%s::GetDefaultGC()()",IsA()->GetName());

   if ((GetDefaultFontStruct() != fFontStruct) || (GetDefaultGC()() != fNormGC)) {
      TGFont *ufont = gClient->GetResourcePool()->GetFontPool()->FindFont(fFontStruct);
      if (ufont) {
         ufont->SavePrimitive(out, option);
         parFont.Form("ufont->GetFontStruct()");
      }

      TGGC *userGC = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);
      if (userGC) {
         userGC->SavePrimitive(out, option);
         parGC.Form("uGC->GetGC()");
      }
   }

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // tab widget" << std::endl;

   out << "   TGTab *";
   out << GetName() << " = new TGTab(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (GetOptions() == kChildFrame) {
         if (fFontStruct == GetDefaultFontStruct()) {
            if (fNormGC == GetDefaultGC()()) {
               out <<");" << std::endl;
            } else {
               out << "," << parGC.Data() <<");" << std::endl;
            }
         } else {
            out << "," << parGC.Data() << "," << parFont.Data() <<");" << std::endl;
         }
      } else {
         out << "," << parGC.Data() << "," << parFont.Data() << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << parGC.Data() << "," << parFont.Data() << "," << GetOptionString()  << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   TGCompositeFrame *cf;
   TGLayoutManager * lm;
   for (Int_t i=0; i<GetNumberOfTabs(); i++) {
      cf = GetTabContainer(i);
      if (!cf || !GetTabTab(i)) continue;
      out << std::endl << "   // container of " << quote
          << GetTabTab(i)->GetString() << quote << std::endl;
      out << "   TGCompositeFrame *" << cf->GetName() << ";" << std::endl;
      out << "   " << cf->GetName() << " = " << GetName()
                   << "->AddTab(" << quote << GetTabTab(i)->GetString()
                   << quote << ");" << std::endl;
      lm = cf->GetLayoutManager();
      if (lm) {
         if ((cf->GetOptions() & kHorizontalFrame) &&
            (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
            ;
         } else if ((GetOptions() & kVerticalFrame) &&
            (lm->InheritsFrom(TGVerticalLayout::Class()))) {
            ;
         } else {
            out << "   " << cf->GetName() <<"->SetLayoutManager(";
            lm->SavePrimitive(out, option);
            out << ");" << std::endl;
         }
         if (!IsEnabled(i)) {
            out << "   " << GetName() << "->SetEnabled(" << i << ", kFALSE);" << std::endl;
         }
      }
      cf->SavePrimitiveSubframes(out, option);

      if (GetTabTab(i)->IsCloseShown()) {
         out << "   TGTabElement *tab" << i << " = "
             << GetName() << "->GetTabTab(" << i << ");" << std::endl;
         out << "   tab" << i << "->ShowClose(kTRUE);" << std::endl;
      }
      if (GetTabTab(i)->GetBackground() != GetTabTab(i)->GetDefaultFrameBackground()) {
         GetTabTab(i)->SaveUserColor(out, option);
         out << "   TGTabElement *tab" << i << " = "
             << GetName() << "->GetTabTab(" << i << ");" << std::endl;
         out << "   tab" << i << "->ChangeBackground(ucolor);" << std::endl;
      }

   }
   out << std::endl << "   " << GetName() << "->SetTab(" << GetCurrent() << ");" << std::endl;
   out << std::endl << "   " << GetName() << "->Resize(" << GetName()
       << "->GetDefaultSize());" << std::endl;
}

// __________________________________________________________________________
void TGTabLayout::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   // Save tab layout manager as a C++ statement(s) on out stream.

   out << "new TGTabLayout(" << fMain->GetName() << ")";

}
