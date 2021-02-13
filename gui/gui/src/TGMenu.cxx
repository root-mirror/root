// @(#)root/gui:$Id$
// Author: Fons Rademakers   09/01/98

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
// TGMenuBar, TGPopupMenu, TGMenuTitle and TGMenuEntry                  //
//                                                                      //
// The TGMenu.h header contains all different menu classes.             //
//                                                                      //
// The TGMenuBar class implements a menu bar widget. It is used to      //
// specify and provide access to common and frequently used application //
// actions grouped under menu titles (TGMenuTitle class). The menu bar  //
// takes the highest-level of the menu system and it is a starting      //
// point for many interactions. It is always visible and allows using   //
// the keyboard equivalents. The geometry of the menu bar is            //
// automatically set to the parent widget, i.e. the menu bar            //
// automatically resizes itself so that it has the same width as its    //
// parent (typically TGMainFrame). A menu bar contains one or more      //
// popup menus and usually is placed along the top of the application   //
// window. Any popup menu is invisible until the user invokes it by     //
// using the mouse pointer or the keyboard.                             //
//                                                                      //
// Popup menus implemented by TGPopupMenu class are unique in that,     //
// by convention, they are not placed with the other GUI components in  //
// the user interfaces. Instead, a popup menu usually appears either in //
// a menu bar or as a context menu on the TOP of the GUI. For that      //
// reason it needs gClient->GetDefaultRoot() as a parent to get the     //
// pointer to the root (i.e. desktop) window. This way a popup menu     //
// will never be embedded.                                              //
// NOTE: Using gClient->GetRoot() as a parent of TGPopupMenu will not   //
// avoid the possibility of embedding the corresponding popup menu      //
// because the current window hierarchy can be changed by using         //
// gClient->SetRoot() method.                                           //
//                                                                      //
// As a context menus TGPopupMenu shows up after pressing the right     //
// mouse button, over a popup-enabled component. The popup menu then    //
// appears under the mouse pointer.                                     //
//                                                                      //
// Selecting a menu item will generate the event:                       //
// kC_COMMAND, kCM_MENU, menu id, user data.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGMenu.h"
#include "TGResourcePool.h"
#include "TTimer.h"
#include "TMath.h"
#include "TSystem.h"
#include "TList.h"
#include "KeySymbols.h"
#include "TGButton.h"
#include "TQConnection.h"
#include "TParameter.h"
#include "RConfigure.h"
#include "TVirtualX.h"
#include "snprintf.h"

#include <iostream>

const TGGC   *TGPopupMenu::fgDefaultGC = nullptr;
const TGGC   *TGPopupMenu::fgDefaultSelectedGC = nullptr;
const TGGC   *TGPopupMenu::fgDefaultSelectedBackgroundGC = nullptr;
const TGFont *TGPopupMenu::fgDefaultFont = nullptr;
const TGFont *TGPopupMenu::fgHilightFont = nullptr;

const TGGC   *TGMenuTitle::fgDefaultGC = nullptr;
const TGGC   *TGMenuTitle::fgDefaultSelectedGC = nullptr;
const TGFont *TGMenuTitle::fgDefaultFont = nullptr;


ClassImp(TGMenuBar);
ClassImp(TGMenuTitle);
ClassImpQ(TGPopupMenu)


////////////////////////////////////////////////////////////////////////////////

class TPopupDelayTimer : public TTimer {
private:
   TGPopupMenu   *fPopup;   // popup menu
public:
   TPopupDelayTimer(TGPopupMenu *p, Long_t ms) : TTimer(ms, kTRUE) { fPopup = p; }
   Bool_t Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Notify when timer times out and reset the timer.

Bool_t TPopupDelayTimer::Notify()
{
   fPopup->HandleTimer(0);
   Reset();
   return kFALSE;
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGMenuBar member functions.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a menu bar object.

TGMenuBar::TGMenuBar(const TGWindow *p, UInt_t w, UInt_t h, UInt_t options)
   : TGHorizontalFrame(p, w, h, options | kHorizontalFrame)
{
   fCurrent       = 0;
   fTitles        = new TList;
   fStick         = kTRUE;
   fDefaultCursor = fClient->GetResourcePool()->GetGrabCursor();
   fTrash         = new TList();

   gVirtualX->GrabButton(fId, kButton1, kAnyModifier,
                       kButtonPressMask | kButtonReleaseMask | kEnterWindowMask,
                       kNone, kNone);

   fKeyNavigate = kFALSE;

   fMenuMore = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuMore->AddLabel("Hidden Menus");
   fMenuMore->AddSeparator();
   fMenuBarMoreLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   fWithExt = kFALSE;
   fOutLayouts = new TList();
   fNeededSpace = new TList();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete menu bar object. Removes also the hot keys from the main frame,
/// so hitting them will not cause the menus to popup.

TGMenuBar::~TGMenuBar()
{
   TGFrameElement *el;
   TGMenuTitle    *t;
   Int_t           keycode;

   if (!MustCleanup()) {
      fTrash->Delete();
   }
   delete fTrash;

   const TGMainFrame *main = (TGMainFrame *)GetMainFrame();

   if (!MustCleanup()) {
      TIter next(fList);
      while ((el = (TGFrameElement *) next())) {
         t = (TGMenuTitle *) el->fFrame;
         if ((keycode = t->GetHotKeyCode()) != 0 && main) {
            main->RemoveBind(this, keycode, kKeyMod1Mask);
         }
      }
   }

   // delete TGMenuTitles
   if (fTitles && !MustCleanup()) fTitles->Delete();
   delete fTitles;

   delete fOutLayouts;
   fNeededSpace->Delete();
   delete fNeededSpace;
   delete fMenuMore;
   delete fMenuBarMoreLayout;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates whether the >> menu must be shown or not and
/// which menu titles are hidden.

void TGMenuBar::Layout()
{
   if (GetDefaultWidth() > GetWidth()) {
      while (!(GetDefaultWidth() < GetWidth() ||
               GetList()->GetSize() <= 1)) {
         TGFrameElement* entry = GetLastOnLeft();
         if (!entry) break;
         TGMenuTitle* menuTitle = (TGMenuTitle*) entry->fFrame;
         fNeededSpace->AddLast(new TParameter<Int_t>("", menuTitle->GetWidth() +
                                                         entry->fLayout->GetPadLeft() +
                                                         entry->fLayout->GetPadRight() ) );
         fOutLayouts->AddLast( entry->fLayout );
         fMenuMore->AddPopup( menuTitle->GetName(), menuTitle->GetMenu() );
         menuTitle->GetMenu()->Connect("PoppedUp()", "TGMenuBar", this, "PopupConnection()");
         RemovePopup( menuTitle->GetName() );
      }
   }

   if (fNeededSpace->GetSize() > 0) {
      Int_t neededWidth = ((TParameter<Int_t>*) fNeededSpace->Last())->GetVal();
      Bool_t fit = kFALSE;
      if (fNeededSpace->GetSize() > 1)
         fit = GetDefaultWidth() + neededWidth + 5 < GetWidth();
      else
         fit = GetDefaultWidth() + neededWidth - 7 < GetWidth();
      while (fit) {
         TGMenuEntry* menu = (TGMenuEntry*) fMenuMore->GetListOfEntries()->Last();
         TGLayoutHints* layout = (TGLayoutHints*) fOutLayouts->Last();
         ULong_t hints = (layout) ? layout->GetLayoutHints() : 0;
         TGPopupMenu* beforeMenu = 0;
         if (hints & kLHintsRight) {
            TGFrameElement* entry = GetLastOnLeft();
            if (entry) {
               TGMenuTitle* beforeMenuTitle = (TGMenuTitle*) entry->fFrame;
               beforeMenu = beforeMenuTitle->GetMenu();
            }
         }
         if (menu && menu->GetPopup()) {
            menu->GetPopup()->Disconnect("PoppedUp()", this, "PopupConnection()");
            AddPopup( menu->GetName(), menu->GetPopup(), layout, beforeMenu );
         }
         fOutLayouts->Remove( fOutLayouts->Last() );
         fNeededSpace->Remove( fNeededSpace->Last() );
         fMenuMore->DeleteEntry(menu);

         if (fNeededSpace->GetSize() > 0) {
            neededWidth = ((TParameter<Int_t>*)fNeededSpace->Last())->GetVal();
            if (fNeededSpace->GetSize() > 1)
               fit = GetDefaultWidth() + neededWidth + 5 < GetWidth();
            else
               fit = GetDefaultWidth() + neededWidth - 7 < GetWidth();
         } else
            fit = kFALSE;
      }
   }

   if (fNeededSpace->GetSize() > 0) {
      if (!fWithExt) {
         AddPopup(">>", fMenuMore, fMenuBarMoreLayout,
                  ((TGMenuTitle*)((TGFrameElement*)GetList()->First())->fFrame)->GetMenu());
         fWithExt = kTRUE;
      }
   } else {
      RemovePopup(">>");
      fWithExt = kFALSE;
   }

   MapSubwindows();
   TGHorizontalFrame::Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the last visible menu title on the left of the '>>'
/// in the menu bar.

TGFrameElement* TGMenuBar::GetLastOnLeft()
{
   TIter next(GetList());
   while (TGFrameElement *entry = (TGFrameElement*) next()) {

      TGMenuTitle* menuTitle = (TGMenuTitle*) entry->fFrame;
      TGLayoutHints* tmpLayout = (TGLayoutHints*) entry->fLayout;
      ULong_t  hints = tmpLayout->GetLayoutHints();

      if (hints & kLHintsRight && menuTitle->GetMenu() != fMenuMore) {
         return entry;
      }
   }

   return ((TGFrameElement*)GetList()->Last());
}

////////////////////////////////////////////////////////////////////////////////
/// Connects the corresponding cascaded menu to the proper slots,
/// according to the highlighted menu entry in '>>' menu.

void TGMenuBar::PopupConnection()
{
   // Disconnect all previous signals
   TList* slots = fMenuMore->GetListOfSignals();
   TIter next (slots);
   while (TList* connlist = (TList*) next()) {

      const char* signal_name = connlist->GetName();
      TIter next2(connlist);
      while (TQConnection* conn = (TQConnection*) next2()) {
         const char* slot_name = conn->GetName();
         void* receiver = conn->GetReceiver();
         fMenuMore->Disconnect(signal_name, receiver, slot_name);
      }
   }
   fMenuMore->fMsgWindow = 0;

   // Check wheter the current entry is a menu or not (just in case)
   TGMenuEntry* currentEntry = fMenuMore->GetCurrent();
   if (currentEntry->GetType() != kMenuPopup) return;

   // Connect the corresponding active signals to the >> menu
   TGPopupMenu* currentMenu = currentEntry->GetPopup();

   slots = currentMenu->GetListOfSignals();
   TIter next3 (slots);
   while (TList* connlist = (TList*) next3()) {

      const char* signal_name = connlist->GetName();
      if (strcmp(signal_name, "Activated(int)") == 0) {
         TIter next2(connlist);
         while (TQConnection* conn = (TQConnection*) next2()) {

            const char* slot_name = conn->GetName();
            const char* class_name = conn->GetClassName();
            void* receiver = conn->GetReceiver();
            fMenuMore->Connect(signal_name, class_name, receiver, slot_name);
         }
      }
   }

   fMenuMore->fMsgWindow = currentMenu->fMsgWindow;
}

////////////////////////////////////////////////////////////////////////////////
/// If on kTRUE bind arrow, popup menu hot keys, otherwise
/// remove key bindings.

void TGMenuBar::BindKeys(Bool_t on)
{
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Left), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Right), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Up), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Down), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Enter), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Return), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Escape), kAnyModifier, on);

   if (fCurrent && fCurrent->GetMenu()) {
      BindMenu(fCurrent->GetMenu(), on);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If on kTRUE bind subMenu hot keys, otherwise remove key bindings.

void TGMenuBar::BindMenu(TGPopupMenu* subMenu, Bool_t on)
{
   TGMenuEntry *e;
   TIter next(subMenu->GetListOfEntries());

   while ((e = (TGMenuEntry*)next())) {
      Int_t hot = 0;
      if ( e->GetType() == kMenuPopup )
         BindMenu(e->GetPopup(), on);
      if (e->GetLabel()) {
         hot = e->GetLabel()->GetHotChar();
      }
      if (!hot) continue;
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), 0, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyShiftMask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyLockMask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyMod2Mask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyShiftMask | kKeyLockMask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyShiftMask | kKeyMod2Mask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyLockMask  | kKeyMod2Mask, on);
      gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(hot), kKeyShiftMask | kKeyLockMask | kKeyMod2Mask, on);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// If on kTRUE bind hot keys, otherwise remove key binding.

void TGMenuBar::BindHotKey(Int_t keycode, Bool_t on)
{
   const TGMainFrame *main = (TGMainFrame *) GetMainFrame();

   if (!main || !main->InheritsFrom("TGMainFrame")) return;

   if (on) {
      // case unsensitive bindings
      main->BindKey(this, keycode, kKeyMod1Mask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyShiftMask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyLockMask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

      main->BindKey(this, keycode, kKeyMod1Mask | kKeyMod2Mask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
      main->BindKey(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
   } else {
      main->RemoveBind(this, keycode, kKeyMod1Mask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyShiftMask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyLockMask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyMod2Mask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
      main->RemoveBind(this, keycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add popup menu to menu bar. The hot string will be adopted by the
/// menubar (actually the menu title) and deleted when possible.
/// If before is not 0 the menu will be added before it.

void TGMenuBar::AddPopup(TGHotString *s, TGPopupMenu *menu, TGLayoutHints *l,
                         TGPopupMenu *before)
{
   TGMenuTitle *t;
   Int_t keycode;

   AddFrameBefore(t = new TGMenuTitle(this, s, menu), l, before);
   fTitles->Add(t);  // keep track of menu titles for later cleanup in dtor

   if ((keycode = t->GetHotKeyCode()) != 0) {
      BindHotKey(keycode, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add popup via created before menu title.

void TGMenuBar::AddTitle(TGMenuTitle *title, TGLayoutHints *l, TGPopupMenu *before)
{
   Int_t keycode;

   AddFrameBefore(title, l, before);
   fTitles->Add(title);  // keep track of menu titles for later cleanup in dtor

   if ((keycode = title->GetHotKeyCode()) != 0) {
      BindHotKey(keycode, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add popup menu to menu bar. If before is not 0 the menu will be
/// added before it.

void TGMenuBar::AddPopup(const char *s, TGPopupMenu *menu, TGLayoutHints *l,
                         TGPopupMenu *before)
{
   AddPopup(new TGHotString(s), menu, l, before);
}

////////////////////////////////////////////////////////////////////////////////
/// Add popup menu to menu bar.
///
/// Comment:
///    This method is valid  only for horizontal menu bars.
///    The most common case is menu bar containing equidistant titles padding left side.
///       TGMenuBar *bar;
///       bar->AddPopup("title1", 10);
///       bar->AddPopup("title2", 10);
///       ...
///
///    To add equidistant titles  padding right side padleft must be 0.
///       TGMenuBar *bar;
///       bar->AddPopup("title1", 0, 10);
///       bar->AddPopup("title2", 0, 10);
///       ...
///
///    This method guarantee automatic cleanup when menu bar is destroyed.
///    Do not delete returned popup-menu

TGPopupMenu *TGMenuBar::AddPopup(const TString &s, Int_t padleft, Int_t padright,
                                 Int_t padtop, Int_t padbottom)
{
   ULong_t hints = kLHintsTop;

   if (padleft)  {
      hints |= kLHintsLeft;
   } else {
      hints |= kLHintsRight;
   }

   TGLayoutHints *l = new TGLayoutHints(hints, padleft, padright,
                                               padtop, padbottom);
   fTrash->Add(l);

   TGPopupMenu *menu = new TGPopupMenu(fClient->GetDefaultRoot());
   AddPopup(new TGHotString(s), menu, l, 0);
   fTrash->Add(menu);
   return menu;
}

////////////////////////////////////////////////////////////////////////////////
/// Private version of AddFrame for menubar, to make sure that we
/// indeed only add TGMenuTitle objects to it. If before is not 0
/// the menu will be added before it.

void TGMenuBar::AddFrameBefore(TGFrame *f, TGLayoutHints *l,
                               TGPopupMenu *before)
{
   if (!f->InheritsFrom("TGMenuTitle")) {
      Error("AddFrameBefore", "may only add TGMenuTitle objects to a menu bar");
      return;
   }

   if (!before) {
      AddFrame(f, l);
      return;
   }

   TGFrameElement *nw;

   nw = new TGFrameElement;
   nw->fFrame  = f;
   nw->fLayout = l ? l : fgDefaultHints;
   nw->fState  = 1;

   TGFrameElement *el;
   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      TGMenuTitle *t = (TGMenuTitle *) el->fFrame;
      if (t->GetMenu() == before) {
         fList->AddBefore(el, nw);
         return;
      }
   }
   fList->Add(nw);
}

////////////////////////////////////////////////////////////////////////////////
/// Return popup menu with the specified name. Returns 0 if menu is
/// not found. Returnes menu can be used as "before" in AddPopup().
/// Don't use hot key (&) in name.

TGPopupMenu *TGMenuBar::GetPopup(const char *s)
{
   if (!GetList()) return 0;

   TGFrameElement *el;
   TIter next(GetList());
   TString str = s;

   while ((el = (TGFrameElement *) next())) {
      TGMenuTitle *t = (TGMenuTitle *) el->fFrame;
      if (str == t->GetName())
         return t->GetMenu();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove popup menu from menu bar. Returned menu has to be deleted by
/// the user, or can be re-used in another AddPopup(). Returns 0 if
/// menu is not found. Don't use hot key (&) in name.

TGPopupMenu *TGMenuBar::RemovePopup(const char *s)
{
   if (!GetList()) return 0;

   TGFrameElement *el;
   TIter next(GetList());
   TString str = s;

   while ((el = (TGFrameElement *) next())) {
      TGMenuTitle *t = (TGMenuTitle *) el->fFrame;
      if (str == t->GetName()) {
         Int_t keycode;
         if ((keycode = t->GetHotKeyCode())) {
            BindHotKey(keycode, kFALSE);  // remove bind
         }
         TGPopupMenu *m = t->GetMenu();
         fTitles->Remove(t);
         t->DestroyWindow();
         RemoveFrame(t);
         delete t;
         return m;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a mouse motion event in a menu bar.

Bool_t TGMenuBar::HandleMotion(Event_t *event)
{
   if (fKeyNavigate) return kTRUE;

   Int_t        dummy;
   Window_t     wtarget;
   TGMenuTitle *target = 0;

   if (!(event->fState & kButton1Mask))
      fStick = kFALSE; // use some threshold!

   gVirtualX->TranslateCoordinates(fId, fId, event->fX, event->fY,
                                   dummy, dummy, wtarget);
   if (wtarget) target = (TGMenuTitle*) fClient->GetWindowById(wtarget);

   if (fCurrent && target && (target != fCurrent)) {
      // deactivate all others
      TGFrameElement *el;
      TIter next(fList);
      while ((el = (TGFrameElement *) next()))
         ((TGMenuTitle*)el->fFrame)->SetState(kFALSE);

      fStick   = kTRUE;
      fCurrent = target;
      target->SetState(kTRUE);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a mouse button event in a menubar.

Bool_t TGMenuBar::HandleButton(Event_t *event)
{
   Int_t        dummy;
   Window_t     wtarget;
   TGMenuTitle *target;

   // We don't need to check the button number as GrabButton will
   // only allow button1 events

   if (event->fType == kButtonPress) {

      gVirtualX->TranslateCoordinates(fId, fId, event->fX, event->fY,
                                      dummy, dummy, wtarget);
      target = (TGMenuTitle*) fClient->GetWindowById(wtarget);

      if (target != 0) {
         fStick = kTRUE;

         if (target != fCurrent) {
            // deactivate all others
            TGFrameElement *el;
            TIter next(fList);
            while ((el = (TGFrameElement *) next()))
               ((TGMenuTitle*)el->fFrame)->SetState(kFALSE);

            fStick   = kTRUE;
            fCurrent = target;
            target->SetState(kTRUE);

            gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                   kPointerMotionMask, kNone, fDefaultCursor);
         }
      }
   }

   if (event->fType == kButtonRelease) {
      if (fStick) {
         fStick = kFALSE;
         return kTRUE;
      }

      TGFrameElement *el;
      TIter next(fList);
      while ((el = (TGFrameElement *) next()))
         ((TGMenuTitle*)el->fFrame)->SetState(kFALSE);

      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer

      if (fCurrent != 0) {
         target   = fCurrent; // tricky, because WaitFor
         fCurrent = 0;
         if (!fKeyNavigate)
            target->DoSendMessage();
      }
      fKeyNavigate = kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events in a menu bar.

Bool_t TGMenuBar::HandleKey(Event_t *event)
{
   TGMenuTitle *target = 0;
   TGFrameElement *el;
   void *dummy;
   Int_t    ax, ay;
   Window_t wdummy;
   TIter next(fList);

   if (event->fType == kGKeyPress) {
      UInt_t keysym;
      char tmp[2];

      gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);

      if (event->fState & kKeyMod1Mask) {
         while ((el = (TGFrameElement *) next())) {
            target = (TGMenuTitle *) el->fFrame;
            if ((Int_t)event->fCode == target->GetHotKeyCode()) {
               RequestFocus();
               fKeyNavigate = kTRUE;
               break;
            }
         }
         if (el == 0) target = 0;
      } else {
         fKeyNavigate = kTRUE;

         if (fCurrent) {
            TGFrameElement *cur  = 0;
            TGPopupMenu    *menu = 0;
            next.Reset();
            el = 0;
            while ((el = (TGFrameElement *) next())) {
               if (el->fFrame == fCurrent) {
                  cur = el;
                  menu = ((TGMenuTitle*)el->fFrame)->GetMenu();
                  break;
               }
            }

            if (!menu || !menu->fPoppedUp) return kFALSE;

            TGMenuEntry *ce = 0;

            TGPopupMenu* currentMenu = fCurrent->GetMenu();
            TGMenuEntry* currentEntry = currentMenu->GetCurrent();
            while ( currentEntry ) {
               if ( currentEntry->GetType() == kMenuPopup )
                  currentMenu = currentEntry->GetPopup();
               if ( currentEntry != currentMenu->GetCurrent() )
                  currentEntry = currentMenu->GetCurrent();
               else
                  currentEntry = 0;
            }

            TIter next2(currentMenu->GetListOfEntries());

            while ((ce = (TGMenuEntry*)next2())) {
               UInt_t hot = 0;
               if (ce->GetLabel()) hot = ce->GetLabel()->GetHotChar();
               if (!hot || (hot != keysym)) continue;

               currentMenu->Activate(ce);
               if (ce->GetType() != kMenuPopup) {
                  gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
                  fCurrent->SetState(kFALSE);
                  currentMenu->fStick = kFALSE;
                  Event_t ev;
                  ev.fType = kButtonRelease;
                  ev.fWindow = currentMenu->GetId();
                  fCurrent = 0;
                  return currentMenu->HandleButton(&ev);
               }
               else {
                  gVirtualX->TranslateCoordinates(currentMenu->fId,
                                 (ce->fPopup->GetParent())->GetId(),
                                  ce->fEx+currentMenu->fMenuWidth, ce->fEy,
                                  ax, ay, wdummy);
#ifdef R__HAS_COCOA
                  gVirtualX->SetWMTransientHint(ce->fPopup->GetId(), GetId());
#endif
                  ce->fPopup->PlaceMenu(ax-5, ay-1, kFALSE, kFALSE);
               }
            }

            ce = menu->GetCurrent();
            TGPopupMenu *submenu = 0;

            while (ce && (ce->GetType() == kMenuPopup)) {
               submenu = ce->GetPopup();
               if (!submenu->fPoppedUp) break;
               ce =  submenu->GetCurrent();
               menu = submenu;
            }
            switch ((EKeySym)keysym) {
               case kKey_Left:
                  if ((submenu) && (submenu->fPoppedUp)) {
                     submenu->EndMenu(dummy);
                     break;
                  }
                  el = (TGFrameElement*)fList->Before(cur);
                  if (!el) el = (TGFrameElement*)fList->Last();
                  break;
               case kKey_Right:
                  if (submenu) {
                     if (submenu->fPoppedUp) {
                        if (!submenu->GetCurrent()) {
                           ce = (TGMenuEntry*)submenu->GetListOfEntries()->First();
                        } else {
                           submenu->EndMenu(dummy);
                        }
                     }
                     else {
                        gVirtualX->TranslateCoordinates(menu->fId,
                                       (submenu->GetParent())->GetId(),
                                       ce->fEx+menu->fMenuWidth, ce->fEy,
                                       ax, ay, wdummy);
#ifdef R__HAS_COCOA
                        gVirtualX->SetWMTransientHint(submenu->GetId(), GetId());
#endif
                        submenu->PlaceMenu(ax-5, ay-1, kFALSE, kFALSE);
                     }
                     break;
                  }
                  el = (TGFrameElement*)fList->After(cur);
                  if (!el) el = (TGFrameElement*)fList->First();
                  break;
               case kKey_Up:
                  if (ce) ce = (TGMenuEntry*)menu->GetListOfEntries()->Before(ce);
                  while (ce && ((ce->GetType() == kMenuSeparator) ||
                         (ce->GetType() == kMenuLabel) ||
                         !(ce->GetStatus() & kMenuEnableMask))) {
                     ce = (TGMenuEntry*)menu->GetListOfEntries()->Before(ce);
                  }
                  if (!ce) ce = (TGMenuEntry*)menu->GetListOfEntries()->Last();
                  break;
               case kKey_Down:
                  if (ce) ce = (TGMenuEntry*)menu->GetListOfEntries()->After(ce);
                  while (ce && ((ce->GetType() == kMenuSeparator) ||
                         (ce->GetType() == kMenuLabel) ||
                         !(ce->GetStatus() & kMenuEnableMask))) {
                     ce = (TGMenuEntry*)menu->GetListOfEntries()->After(ce);
                  }
                  if (!ce) ce = (TGMenuEntry*)menu->GetListOfEntries()->First();
                  break;
               case kKey_Enter:
               case kKey_Return: {
                  gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
                  fCurrent->SetState(kFALSE);
                  menu->fStick = kFALSE;
                  Event_t ev;
                  ev.fType = kButtonRelease;
                  ev.fWindow = menu->GetId();
                  fCurrent = 0;
                  return menu->HandleButton(&ev);
               }
               case kKey_Escape:
                  gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
                  fCurrent->SetState(kFALSE);
                  fStick = kFALSE;
                  fCurrent = 0;
                  return menu->EndMenu(dummy);
               default:
                  break;
            }
            if (ce) menu->Activate(ce);

            el = el ? el : cur;
            if (el) target = (TGMenuTitle*)el->fFrame;
         } else {
            return kFALSE;
         }
      }

      if (target != 0) {
         fStick = kTRUE;

         if (target != fCurrent) {
            // deactivate all others
            next.Reset();
            while ((el = (TGFrameElement *) next()))
               ((TGMenuTitle*)el->fFrame)->SetState(kFALSE);

            fCurrent = target;
            target->SetState(kTRUE);
            fStick   = kTRUE;

            gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                   kPointerMotionMask, kNone, fDefaultCursor);

            TGMenuEntry *ptr;
            TIter nexte(target->GetMenu()->GetListOfEntries());

            while ((ptr = (TGMenuEntry *) nexte())) {
               if ((ptr->GetStatus() & kMenuEnableMask) &&
                  !(ptr->GetStatus() & kMenuHideMask) &&
                   (ptr->GetType() != kMenuSeparator) &&
                   (ptr->GetType() != kMenuLabel)) break;
            }
            if (ptr)
               target->GetMenu()->Activate(ptr);

            return kTRUE;
         }
      } else {
         return kFALSE;
      }
   }

   if (event->fType == kKeyRelease) {
      if (fStick) {
         fStick = kFALSE;
         return kTRUE;
      }
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer

      next.Reset();
      while ((el = (TGFrameElement *) next()))
         ((TGMenuTitle*)el->fFrame)->SetState(kFALSE);

      if (fCurrent != 0) {
         target   = fCurrent; // tricky, because WaitFor
         fCurrent = 0;
         target->DoSendMessage();
      }
   }

   return kTRUE;
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGPopupMenu member functions.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a popup menu.

TGPopupMenu::TGPopupMenu(const TGWindow *p, UInt_t w, UInt_t h, UInt_t options)
    : TGFrame(p, w, h, options | kOwnBackground)
{
   fNormGC        = GetDefaultGC()();
   fSelGC         = GetDefaultSelectedGC()();
   fSelbackGC     = GetDefaultSelectedBackgroundGC()();
   fFontStruct    = GetDefaultFontStruct();
   fHifontStruct  = GetHilightFontStruct();
   fDefaultCursor = fClient->GetResourcePool()->GetGrabCursor();

   // We need to change the default context to actually use the
   // Menu Fonts.  [Are we actually changing the global settings?]
   GCValues_t    gcval;
   gcval.fMask = kGCFont;
   gcval.fFont = gVirtualX->GetFontHandle(fFontStruct);
   gVirtualX->ChangeGC(fNormGC, &gcval);
   gVirtualX->ChangeGC(fSelGC, &gcval);

   fDelay     = 0;
   fEntryList = new TList;

   // in case any of these magic values is changes, check also Reposition()
   fBorderWidth = 3;
   fMenuHeight  = 6;
   fMenuWidth   = 8;
   fXl          = 16;
   fMsgWindow   = p;
   fStick       = kTRUE;
   fCurrent     = 0;
   fHasGrab     = kFALSE;
   fPoppedUp    = kFALSE;
   fMenuBar     = 0;
   fSplitButton = 0;
   fEntrySep    = 3;

   SetWindowAttributes_t wattr;
   wattr.fMask             = kWAOverrideRedirect | kWASaveUnder;
   wattr.fOverrideRedirect = kTRUE;
   wattr.fSaveUnder        = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &wattr);

   AddInput(kPointerMotionMask | kEnterWindowMask | kLeaveWindowMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a popup menu.

TGPopupMenu::~TGPopupMenu()
{
   gClient->UnregisterPopup(this);

   if (fEntryList) fEntryList->Delete();
   delete fEntryList;
   delete fDelay;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a menu entry. The hotstring is adopted by the menu (actually by
/// the TGMenuEntry) and deleted when possible. A possible picture is
/// borrowed from the picture pool and therefore not adopted.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddEntry(TGHotString *s, Int_t id, void *ud,
                           const TGPicture *p, TGMenuEntry *before)
{
   if (!s) return;
   TGMenuEntry *nw = new TGMenuEntry;
   Ssiz_t tab = s->Index('\t');
   if (tab > 0) {
      TString ts(s->Data());
      TString shortcut = ts(tab+1, s->Length());
      nw->fShortcut = new TGString(shortcut.Data());
      nw->fLabel = new TGHotString(*s);
      nw->fLabel->Remove(tab);
   }
   else {
      nw->fLabel = s;
   }
   nw->fPic      = p;
   nw->fType     = kMenuEntry;
   nw->fEntryId  = id;
   nw->fUserData = ud;
   nw->fPopup    = 0;
   nw->fStatus   = kMenuEnableMask;
   nw->fEx       = 2;
   nw->fEy       = fMenuHeight-2;

   if (before)
      fEntryList->AddBefore(before, nw);
   else
      fEntryList->Add(nw);

   UInt_t tw, ph = 0, pw = 0;
   tw = gVirtualX->TextWidth(fHifontStruct, s->GetString(), s->GetLength());
   if (p) {
      ph = p->GetHeight();
      pw = p->GetWidth();
      if (pw+12 > fXl) { fMenuWidth += pw+12-fXl; fXl = pw+12; }
   }
   if (nw->fShortcut) {
      tw += 10;
      delete s;
   }

   Int_t max_ascent, max_descent;
   nw->fEw = tw + pw /*+8*/+18+12;
   fMenuWidth = TMath::Max(fMenuWidth, nw->fEw);
   gVirtualX->GetFontProperties(fHifontStruct, max_ascent, max_descent);
   nw->fEh = max_ascent + max_descent + fEntrySep;
   if (nw->fEh < ph+fEntrySep) nw->fEh = ph+fEntrySep;
   fMenuHeight += nw->fEh;

   if (before)
      Reposition();
   else
      Resize(fMenuWidth, fMenuHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a menu entry. The string s in not adopted.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddEntry(const char *s, Int_t id, void *ud,
                           const TGPicture *p, TGMenuEntry *before)
{
   AddEntry(new TGHotString(s), id, ud, p, before);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a menu separator to the menu.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddSeparator(TGMenuEntry *before)
{
   TGMenuEntry *nw = new TGMenuEntry;

   nw->fLabel    = 0;
   nw->fPic      = 0;
   nw->fType     = kMenuSeparator;
   nw->fEntryId  = -1;
   nw->fUserData = 0;
   nw->fPopup    = 0;
   nw->fStatus   = kMenuEnableMask;
   nw->fEx       = 2;
   nw->fEy       = fMenuHeight-2;

   if (before)
      fEntryList->AddBefore(before, nw);
   else
      fEntryList->Add(nw);

   nw->fEw = 0;
   nw->fEh = 4;
   fMenuHeight += nw->fEh;

   if (before)
      Reposition();
   else
      Resize(fMenuWidth, fMenuHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a menu label to the menu. The hotstring is adopted by the menu
/// (actually by the TGMenuEntry) and deleted when possible. A possible
/// picture is borrowed from the picture pool and therefore not adopted.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddLabel(TGHotString *s, const TGPicture *p,
                           TGMenuEntry *before)
{
   TGMenuEntry *nw = new TGMenuEntry;

   nw->fLabel    = s;
   nw->fPic      = p;
   nw->fType     = kMenuLabel;
   nw->fEntryId  = -1;
   nw->fUserData = 0;
   nw->fPopup    = 0;
   nw->fStatus   = kMenuEnableMask | kMenuDefaultMask;
   nw->fEx       = 2;
   nw->fEy       = fMenuHeight-2;

   if (before)
      fEntryList->AddBefore(before, nw);
   else
      fEntryList->Add(nw);

   UInt_t tw, ph = 0, pw = 0;
   tw = gVirtualX->TextWidth(fHifontStruct, s->GetString(), s->GetLength());
   if (p) {
      ph = p->GetHeight();
      pw = p->GetWidth();
      if (pw+12 > fXl) { fMenuWidth += pw+12-fXl; fXl = pw+12; }
   }

   Int_t max_ascent, max_descent;
   nw->fEw = tw + pw /*+8*/+18+12;
   fMenuWidth = TMath::Max(fMenuWidth, nw->fEw);
   gVirtualX->GetFontProperties(fHifontStruct, max_ascent, max_descent);
   nw->fEh = max_ascent + max_descent + fEntrySep;
   if (nw->fEh < ph+fEntrySep) nw->fEh = ph+fEntrySep;
   fMenuHeight += nw->fEh;

   if (before)
      Reposition();
   else
      Resize(fMenuWidth, fMenuHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a menu label to the menu. The string s in not adopted.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddLabel(const char *s, const TGPicture *p,
                           TGMenuEntry *before)
{
   AddLabel(new TGHotString(s), p, before);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a (cascading) popup menu to a popup menu. The hotstring is adopted
/// by the menu (actually by the TGMenuEntry) and deleted when possible.
/// If before is not 0, the entry will be added before it.

void TGPopupMenu::AddPopup(TGHotString *s, TGPopupMenu *popup,
                           TGMenuEntry *before, const TGPicture *p)
{
   TGMenuEntry *nw = new TGMenuEntry;

   nw->fLabel    = s;
   nw->fPic      = p;
   nw->fType     = kMenuPopup;
   nw->fEntryId  = -2;
   nw->fUserData = 0;
   nw->fPopup    = popup;
   nw->fStatus   = kMenuEnableMask;
   nw->fEx       = 2;
   nw->fEy       = fMenuHeight-2;

   if (before)
      fEntryList->AddBefore(before, nw);
   else
      fEntryList->Add(nw);

   UInt_t tw = gVirtualX->TextWidth(fHifontStruct, s->GetString(),
                                    s->GetLength());

   UInt_t ph = 0, pw = 8;
   if (p) {
      ph = p->GetHeight();
      pw = p->GetWidth();
      if (pw+12 > fXl) { fMenuWidth += pw+12-fXl; fXl = pw+12; }
   }
   Int_t max_ascent, max_descent;
   nw->fEw = tw + pw+18+12;
   fMenuWidth = TMath::Max(fMenuWidth, nw->fEw);
   gVirtualX->GetFontProperties(fHifontStruct, max_ascent, max_descent);
   nw->fEh = max_ascent + max_descent + fEntrySep;
   if (nw->fEh < ph+fEntrySep) nw->fEh = ph+fEntrySep;
   fMenuHeight += nw->fEh;

   if (before)
      Reposition();
   else
      Resize(fMenuWidth, fMenuHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a (cascading) popup menu to a popup menu. The string s is not
/// adopted. If before is not 0, the entry will be added before it.

void TGPopupMenu::AddPopup(const char *s, TGPopupMenu *popup,
                           TGMenuEntry *before, const TGPicture *p)
{
   AddPopup(new TGHotString(s), popup, before, p);
}

////////////////////////////////////////////////////////////////////////////////
/// Reposition entries in popup menu. Called after menu item has been
/// hidden or removed or inserted at a specified location.

void TGPopupMenu::Reposition()
{
   // in case any of these magic values is changes, check also the ctor.
   fMenuHeight = 6;
   fMenuWidth  = 8;
   fXl         = 16;

   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next())) {

      if (ptr->fStatus & kMenuHideMask) continue;

      if (ptr->fPic) {
         UInt_t pw = ptr->fPic->GetWidth();
         if (pw+12 > fXl) { fMenuWidth += pw+12-fXl; fXl = pw+12; }
      }
      ptr->fEx     = 2;
      ptr->fEy     = fMenuHeight-2;
      fMenuWidth   = TMath::Max(fMenuWidth, ptr->fEw);
      fMenuHeight += ptr->fEh;
   }
   Resize(fMenuWidth, fMenuHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Popup a popup menu. If stick mode is true keep the menu up. If
/// grab_pointer is true the pointer will be grabbed, which means that
/// all pointer events will go to the popup menu, independent of in
/// which window the pointer is.

void TGPopupMenu::PlaceMenu(Int_t x, Int_t y, Bool_t stick_mode, Bool_t grab_pointer)
{
   void *ud;
   EndMenu(ud);

   Int_t  rx, ry;
   UInt_t rw, rh;

   fStick = stick_mode;
   fCurrent = 0;

   // Parent is root window for a popup menu
   gVirtualX->GetWindowSize(fParent->GetId(), rx, ry, rw, rh);

   if (gVirtualX->InheritsFrom("TGWin32")) {
      if ((x > 0) && ((x + abs(rx) + (Int_t)fMenuWidth) > (Int_t)rw))
         x = rw - abs(rx) - fMenuWidth;
      if ((y > 0) && (y + abs(ry) + (Int_t)fMenuHeight > (Int_t)rh))
         y = rh - fMenuHeight;
   } else {
      if (x < 0) x = 0;
      if (x + fMenuWidth > rw) x = rw - fMenuWidth;
      if (y < 0) y = 0;
      if (y + fMenuHeight > rh) y = rh - fMenuHeight;
   }

   Move(x, y);
   MapRaised();

   if (grab_pointer) {
      gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                             kPointerMotionMask, kNone, fDefaultCursor);
      fHasGrab = kTRUE;
   } else {
      fHasGrab = kFALSE;
   }

   fPoppedUp = kTRUE;
   PoppedUp();
   if (fMenuBar) fMenuBar->BindKeys(kTRUE);

   fClient->RegisterPopup(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Close menu and return ID of selected menu item.
/// In case of cascading menus, recursively close all menus.

Int_t TGPopupMenu::EndMenu(void *&userData)
{
   Int_t id;

   if (fDelay) fDelay->Remove();

   // destroy any cascaded children and get any ID

   if (fCurrent != 0) {

      // deactivate the entry
      fCurrent->fStatus &= ~kMenuActiveMask;

      if ((fCurrent->fType == kMenuPopup) && fCurrent->fPopup) {
         id = fCurrent->fPopup->EndMenu(userData);
      } else {
         // return the ID if the entry is enabled, otherwise -1
         if (fCurrent->fStatus & kMenuEnableMask) {
            id       = fCurrent->fEntryId;
            userData = fCurrent->fUserData;
         } else {
            id       = -1;
            userData = 0;
         }
      }

   } else {
      // if no entry selected...
      id       = -1;
      userData = 0;
   }

   // then unmap itself
   UnmapWindow();

   gClient->UnregisterPopup(this);
   if (fMenuBar) fMenuBar->BindKeys(kFALSE);

   if (fPoppedUp) {
      fPoppedUp = kFALSE;
      PoppedDown();
   }

   return id;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button event in the popup menu.

Bool_t TGPopupMenu::HandleButton(Event_t *event)
{
   int   id;
   void *ud = 0;

   if (event->fType == kButtonRelease) {
      if (fStick) {
         fStick = kFALSE;
         return kTRUE;
      }
      //if (fCurrent != 0)
      //   if (fCurrent->fType == kMenuPopup) return kTRUE;
      id = EndMenu(ud);
      if (fHasGrab) gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab
      if (fCurrent != 0) {
         fCurrent->fStatus &= ~kMenuActiveMask;
         if (fCurrent->fStatus & kMenuEnableMask) {
            SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_MENU), id,
                        (Long_t)ud);
            Activated(id);
         }
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle pointer crossing event in popup menu.

Bool_t TGPopupMenu::HandleCrossing(Event_t *event)
{
   if (event->fType == kEnterNotify) {

      TGMenuEntry *ptr;
      TIter next(fEntryList);

      while ((ptr = (TGMenuEntry *) next())) {
         if (ptr->fStatus & kMenuHideMask) continue;

         if ((event->fX >= ptr->fEx) && (event->fX <= ptr->fEx+(Int_t)fMenuWidth-10) &&
             (event->fY >= ptr->fEy) && (event->fY <= ptr->fEy+(Int_t)ptr->fEh))
            break;
      }
      Activate(ptr);
   } else {
      Activate((TGMenuEntry*)0);
   }
   if (fMenuBar) fMenuBar->fKeyNavigate = kFALSE;
   if (fSplitButton) fSplitButton->fKeyNavigate = kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle pointer motion event in popup menu.

Bool_t TGPopupMenu::HandleMotion(Event_t *event)
{
   TGFrame::HandleMotion(event);
   static Int_t twice = 0;
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   if (twice < 2) {
      // hack to eat mouse move events generated by Windows when
      // pressing/releasing a mouse button
      ++twice;
   }
   else {
      twice = 0;
      fStick = kFALSE;   // be careful with this, use some threshold
   }
   while ((ptr = (TGMenuEntry *) next())) {
      if (ptr->fStatus & kMenuHideMask) continue;

      if ((event->fX >= ptr->fEx) && (event->fX <= ptr->fEx+(Int_t)fMenuWidth-4) &&  //fMenuWidth-10??
          (event->fY >= ptr->fEy) && (event->fY <= ptr->fEy+(Int_t)ptr->fEh))
         break;
   }
   Activate(ptr);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate a menu entry in a popup menu.

void TGPopupMenu::Activate(TGMenuEntry *entry)
{
   if (entry == fCurrent) return;

   //-- Deactivate the current entry

   if (fCurrent != 0) {
      void *ud;
      if (entry == 0 && fCurrent->fType == kMenuPopup) return;
      if ((fCurrent->fType == kMenuPopup) && fCurrent->fPopup)
         fCurrent->fPopup->EndMenu(ud);
      fCurrent->fStatus &= ~kMenuActiveMask;
      DrawEntry(fCurrent);
   }

   if (fDelay) fDelay->Remove();

   //-- Activate the new one

   if (entry) {
      entry->fStatus |= kMenuActiveMask;
      DrawEntry(entry);
      if (entry->fType == kMenuPopup) {
         if (!fDelay) fDelay = new TPopupDelayTimer(this, 350);
         fDelay->Reset();
         gSystem->AddTimer(fDelay);
         // after delay expires it will popup the cascading popup menu
         // iff it is still the current entry
      } else if (entry->fType == kMenuEntry) {
         // test...
         SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_MENUSELECT),
                     entry->fEntryId, (Long_t)entry->fUserData);
         Highlighted(entry->fEntryId);
      }
   }
   fCurrent = entry;
}

////////////////////////////////////////////////////////////////////////////////
/// If TPopupDelayTimer times out popup cascading popup menu (if it is
/// still the current entry).

Bool_t TGPopupMenu::HandleTimer(TTimer *)
{
   if (fCurrent != 0) {
      if ((fCurrent->fType == kMenuPopup) && fCurrent->fPopup) {
         Int_t    ax, ay;
         Window_t wdummy;

         gVirtualX->TranslateCoordinates(fId,
                                       (fCurrent->fPopup->GetParent())->GetId(),
                                       fCurrent->fEx+fMenuWidth, fCurrent->fEy,
                                       ax, ay, wdummy);
#ifdef R__HAS_COCOA
         gVirtualX->SetWMTransientHint(fCurrent->fPopup->GetId(), GetId());
#endif
         fCurrent->fPopup->PlaceMenu(ax-5, ay-1, kFALSE, kFALSE);
      }
   }
   fDelay->Remove();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw popup menu.

void TGPopupMenu::DoRedraw()
{
   TGFrame::DoRedraw();

   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      DrawEntry(ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw popup menu entry.

void TGPopupMenu::DrawEntry(TGMenuEntry *entry)
{
   FontStruct_t  font;
   GCValues_t    gcval;

   if (entry->fStatus & kMenuHideMask)
      return;

   if (entry->fStatus & kMenuDefaultMask) {
      font = fHifontStruct;
      gcval.fMask = kGCFont;
      gcval.fFont = gVirtualX->GetFontHandle(font);
      gVirtualX->ChangeGC(fNormGC, &gcval);
      gVirtualX->ChangeGC(fSelGC, &gcval);
   } else {
      font = fFontStruct;
   }

   UInt_t tw = 0;
   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(font, max_ascent, max_descent);
   int tx = entry->fEx + fXl;
   // center text
   int offset = (entry->fEh - (max_ascent + max_descent)) / 2;
   int ty = entry->fEy + max_ascent + offset - 1;
   if (entry->fShortcut)
      tw = 7 + gVirtualX->TextWidth(fFontStruct, entry->fShortcut->Data(), entry->fShortcut->Length());

   switch (entry->fType) {
      case kMenuPopup:
      case kMenuLabel:
      case kMenuEntry:
         if ((entry->fStatus & kMenuActiveMask) && entry->fType != kMenuLabel) {
            gVirtualX->FillRectangle(fId, fSelbackGC, entry->fEx+1, entry->fEy-1,
                                     fMenuWidth-6, entry->fEh);
            if (gClient->GetStyle() > 1)
               gVirtualX->DrawRectangle(fId, GetShadowGC()(), entry->fEx+1, entry->fEy-2,
                                        fMenuWidth-7, entry->fEh);
            if (entry->fType == kMenuPopup)
               DrawTrianglePattern(fSelGC, fMenuWidth-10, entry->fEy+fEntrySep, fMenuWidth-6, entry->fEy+11);
            if (entry->fStatus & kMenuCheckedMask)
               DrawCheckMark(fSelGC, 6, entry->fEy+fEntrySep, 14, entry->fEy+11);
            if (entry->fStatus & kMenuRadioMask)
               DrawRCheckMark(fSelGC, 6, entry->fEy+fEntrySep, 14, entry->fEy+11);
            if (entry->fPic != 0)
               entry->fPic->Draw(fId, fSelGC, 8, entry->fEy+1);
            entry->fLabel->Draw(fId,
                           (entry->fStatus & kMenuEnableMask) ? fSelGC : GetShadowGC()(),
                           tx, ty);
            if (entry->fShortcut)
               entry->fShortcut->Draw(fId, (entry->fStatus & kMenuEnableMask) ? fSelGC : GetShadowGC()(),
                                      fMenuWidth - tw, ty);
         } else {
            if (gClient->GetStyle() > 1)
               gVirtualX->DrawRectangle(fId, GetBckgndGC()(), entry->fEx+1, entry->fEy-2,
                                        fMenuWidth-7, entry->fEh);
            gVirtualX->FillRectangle(fId, GetBckgndGC()(), entry->fEx+1, entry->fEy-1,
                                     fMenuWidth-6, entry->fEh);
            if (entry->fType == kMenuPopup)
               DrawTrianglePattern(fNormGC, fMenuWidth-10, entry->fEy+fEntrySep, fMenuWidth-6, entry->fEy+11);
            if (entry->fStatus & kMenuCheckedMask)
               DrawCheckMark(fNormGC, 6, entry->fEy+fEntrySep, 14, entry->fEy+11);
            if (entry->fStatus & kMenuRadioMask)
               DrawRCheckMark(fNormGC, 6, entry->fEy+fEntrySep, 14, entry->fEy+11);
            if (entry->fPic != 0)
               entry->fPic->Draw(fId, fNormGC, 8, entry->fEy+1);
            if (entry->fStatus & kMenuEnableMask) {
               entry->fLabel->Draw(fId, fNormGC, tx, ty);
               if (entry->fShortcut)
                  entry->fShortcut->Draw(fId, fNormGC, fMenuWidth - tw, ty);
            } else {
               entry->fLabel->Draw(fId, GetHilightGC()(), tx+1, ty+1);
               entry->fLabel->Draw(fId, GetShadowGC()(), tx, ty);
               if (entry->fShortcut) {
                  entry->fShortcut->Draw(fId, GetHilightGC()(), fMenuWidth - tw+1, ty+1);
                  entry->fShortcut->Draw(fId, GetShadowGC()(), fMenuWidth - tw, ty);
               }
            }
         }
         break;

      case kMenuSeparator:
         gVirtualX->DrawLine(fId, GetShadowGC()(),  2, entry->fEy, fMenuWidth-fEntrySep, entry->fEy);
         gVirtualX->DrawLine(fId, GetHilightGC()(), 2, entry->fEy+1, fMenuWidth-fEntrySep, entry->fEy+1);
         break;
   }

   // restore font
   if (entry->fStatus & kMenuDefaultMask) {
      gcval.fFont = gVirtualX->GetFontHandle(fFontStruct);
      gVirtualX->ChangeGC(fNormGC, &gcval);
      gVirtualX->ChangeGC(fSelGC, &gcval);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw border round popup menu.

void TGPopupMenu::DrawBorder()
{
   if (gClient->GetStyle() > 0) {
      // new modern (flat) version
      gVirtualX->DrawLine(fId, GetShadowGC()(), 0, 0, 0, fMenuHeight-1);
      gVirtualX->DrawLine(fId, GetShadowGC()(), 0, fMenuHeight-1, fMenuWidth-1, fMenuHeight-1);
      gVirtualX->DrawLine(fId, GetShadowGC()(), fMenuWidth-1, fMenuHeight-1, fMenuWidth-1, 0);
      gVirtualX->DrawLine(fId, GetShadowGC()(), fMenuWidth-1, 0, 0, 0);
   }
   else {
      // old (raised frame) version
      gVirtualX->DrawLine(fId, GetBckgndGC()(), 0, 0, fMenuWidth-2, 0);
      gVirtualX->DrawLine(fId, GetBckgndGC()(), 0, 0, 0, fMenuHeight-2);
      gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, fMenuWidth-fEntrySep, 1);
      gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, 1, fMenuHeight-fEntrySep);

      gVirtualX->DrawLine(fId, GetShadowGC()(),  1, fMenuHeight-2, fMenuWidth-2, fMenuHeight-2);
      gVirtualX->DrawLine(fId, GetShadowGC()(),  fMenuWidth-2, fMenuHeight-2, fMenuWidth-2, 1);
      gVirtualX->DrawLine(fId, GetBlackGC()(), 0, fMenuHeight-1, fMenuWidth-1, fMenuHeight-1);
      gVirtualX->DrawLine(fId, GetBlackGC()(), fMenuWidth-1, fMenuHeight-1, fMenuWidth-1, 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw triangle pattern. Used for menu entries that are of type
/// kMenuPopup (i.e. cascading menus).

void TGPopupMenu::DrawTrianglePattern(GContext_t gc, Int_t l, Int_t t,
                                      Int_t r, Int_t b)
{
   Point_t  points[3];

   int m = (t + b) >> 1;

   points[0].fX = l;
   points[0].fY = t;
   points[1].fX = l;
   points[1].fY = b;
   points[2].fX = r;
   points[2].fY = m;

   gVirtualX->FillPolygon(fId, gc, points, 3);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw check mark. Used for checked button type menu entries.

void TGPopupMenu::DrawCheckMark(GContext_t gc, Int_t l, Int_t t, Int_t, Int_t b)
{
   Segment_t seg[6];

   t = (t + b - 8) >> 1; ++t;

   seg[0].fX1 = 1+l; seg[0].fY1 = 3+t; seg[0].fX2 = 3+l; seg[0].fY2 = 5+t;
   seg[1].fX1 = 1+l; seg[1].fY1 = 4+t; seg[1].fX2 = 3+l; seg[1].fY2 = 6+t;
   seg[2].fX1 = 1+l; seg[2].fY1 = 5+t; seg[2].fX2 = 3+l; seg[2].fY2 = 7+t;
   seg[3].fX1 = 3+l; seg[3].fY1 = 5+t; seg[3].fX2 = 7+l; seg[3].fY2 = 1+t;
   seg[4].fX1 = 3+l; seg[4].fY1 = 6+t; seg[4].fX2 = 7+l; seg[4].fY2 = 2+t;
   seg[5].fX1 = 3+l; seg[5].fY1 = 7+t; seg[5].fX2 = 7+l; seg[5].fY2 = 3+t;

   gVirtualX->DrawSegments(fId, gc, seg, 6);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw radio check mark. Used for radio button type menu entries.

void TGPopupMenu::DrawRCheckMark(GContext_t gc, Int_t l, Int_t t, Int_t r, Int_t b)
{
   Segment_t seg[5];

   t = (t + b - 5) >> 1; ++t;
   l = (l + r - 5) >> 1; ++l;

   seg[0].fX1 = 1+l; seg[0].fY1 = 0+t; seg[0].fX2 = 3+l; seg[0].fY2 = 0+t;
   seg[1].fX1 = 0+l; seg[1].fY1 = 1+t; seg[1].fX2 = 4+l; seg[1].fY2 = 1+t;
   seg[2].fX1 = 0+l; seg[2].fY1 = 2+t; seg[2].fX2 = 4+l; seg[2].fY2 = 2+t;
   seg[3].fX1 = 0+l; seg[3].fY1 = 3+t; seg[3].fX2 = 4+l; seg[3].fY2 = 3+t;
   seg[4].fX1 = 1+l; seg[4].fY1 = 4+t; seg[4].fX2 = 3+l; seg[4].fY2 = 4+t;

   gVirtualX->DrawSegments(fId, gc, seg, 5);
}

////////////////////////////////////////////////////////////////////////////////
/// Set default entry (default entries are drawn with bold text).

void TGPopupMenu::DefaultEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next())) {
      if (ptr->fEntryId == id)
         ptr->fStatus |= kMenuDefaultMask;
      else
         ptr->fStatus &= ~kMenuDefaultMask;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable entry. By default entries are enabled.

void TGPopupMenu::EnableEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) {
         ptr->fStatus |= kMenuEnableMask;
         if (ptr->fStatus & kMenuHideMask) {
            ptr->fStatus &= ~kMenuHideMask;
            Reposition();
         }
         break;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Disable entry (disabled entries appear in a sunken relieve).

void TGPopupMenu::DisableEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) { ptr->fStatus &= ~kMenuEnableMask; break; }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if menu entry is enabled.

Bool_t TGPopupMenu::IsEntryEnabled(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         return (ptr->fStatus & kMenuEnableMask) ? kTRUE : kFALSE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Hide entry (hidden entries are not shown in the menu).
/// To enable a hidden entry call EnableEntry().

void TGPopupMenu::HideEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) {
         ptr->fStatus |=  kMenuHideMask;
         ptr->fStatus &= ~kMenuEnableMask;
         Reposition();
         break;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if menu entry is hidden.

Bool_t TGPopupMenu::IsEntryHidden(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         return (ptr->fStatus & kMenuHideMask) ? kTRUE : kFALSE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check a menu entry (i.e. add a check mark in front of it).

void TGPopupMenu::CheckEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) { ptr->fStatus |= kMenuCheckedMask; break; }
}

////////////////////////////////////////////////////////////////////////////////
/// Check a menu entry (i.e. add a check mark in front of it).
/// The input argument is user data associated with entry

void TGPopupMenu::CheckEntryByData(void *user_data)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fUserData == user_data) { ptr->fStatus |= kMenuCheckedMask; break; }
}

////////////////////////////////////////////////////////////////////////////////
/// Uncheck menu entry (i.e. remove check mark).

void TGPopupMenu::UnCheckEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) { ptr->fStatus &= ~kMenuCheckedMask; break; }
}

////////////////////////////////////////////////////////////////////////////////
/// Uncheck all entries.

void TGPopupMenu::UnCheckEntries()
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next())) {
      ptr->fStatus &= ~kMenuCheckedMask;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Uncheck a menu entry (i.e. remove check mark in front of it).
/// The input argument is user data associated with entry

void TGPopupMenu::UnCheckEntryByData(void *user_data)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fUserData == user_data) { ptr->fStatus  &= ~kMenuCheckedMask; break; }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if menu item is checked.

Bool_t TGPopupMenu::IsEntryChecked(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         return (ptr->fStatus & kMenuCheckedMask) ? kTRUE : kFALSE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Radio-select entry (note that they cannot be unselected,
/// the selection must be moved to another entry instead).

void TGPopupMenu::RCheckEntry(Int_t id, Int_t IDfirst, Int_t IDlast)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         ptr->fStatus |= kMenuRadioMask | kMenuRadioEntryMask;
      else
         if (ptr->fEntryId >= IDfirst && ptr->fEntryId <= IDlast) {
            ptr->fStatus &= ~kMenuRadioMask;
            ptr->fStatus |=  kMenuRadioEntryMask;
         }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if menu item has radio check mark.

Bool_t TGPopupMenu::IsEntryRChecked(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         return (ptr->fStatus & kMenuRadioMask) ? kTRUE : kFALSE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find entry with specified id. Use the returned entry in DeleteEntry()
/// or as the "before" item in the AddXXXX() methods. Returns 0 if entry
/// is not found. To find entries that don't have an id like the separators,
/// use the GetListOfEntries() method to get the complete entry
/// list and iterate over it and check the type of each entry
/// to find the separators.

TGMenuEntry *TGPopupMenu::GetEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id)
         return ptr;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find entry with specified name. Name must match the original
/// name without hot key symbol, like "Print" and not "&Print".
/// Use the returned entry in DeleteEntry() or as the "before" item
/// in the AddXXXX() methods. Returns 0 if entry is not found.
/// To find entries that don't have a name like the separators,
/// use the GetListOfEntries() method to get the complete entry
/// list and iterate over it and check the type of each entry
/// to find the separators.

TGMenuEntry *TGPopupMenu::GetEntry(const char *s)
{
   return (TGMenuEntry*) fEntryList->FindObject(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete entry with specified id from menu.

void TGPopupMenu::DeleteEntry(Int_t id)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr->fEntryId == id) {
         fEntryList->Remove(ptr);
         delete ptr;
         Reposition();
         if (fCurrent == ptr)
            fCurrent = 0;
         return;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete specified entry from menu.

void TGPopupMenu::DeleteEntry(TGMenuEntry *entry)
{
   TGMenuEntry *ptr;
   TIter next(fEntryList);

   while ((ptr = (TGMenuEntry *) next()))
      if (ptr == entry) {
         fEntryList->Remove(ptr);
         delete ptr;
         Reposition();
         if (fCurrent == ptr)
            fCurrent = 0;
         return;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGPopupMenu::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the selection graphics context in use.

const TGGC &TGPopupMenu::GetDefaultSelectedGC()
{
   if (!fgDefaultSelectedGC)
      fgDefaultSelectedGC = gClient->GetResourcePool()->GetSelectedGC();
   return *fgDefaultSelectedGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the selection background graphics context in use.

const TGGC &TGPopupMenu::GetDefaultSelectedBackgroundGC()
{
   if (!fgDefaultSelectedBackgroundGC)
      fgDefaultSelectedBackgroundGC = gClient->GetResourcePool()->GetSelectedBckgndGC();
   return *fgDefaultSelectedBackgroundGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the default font structure in use.

FontStruct_t TGPopupMenu::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetMenuFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the font structure in use for highlighted menu entries.

FontStruct_t TGPopupMenu::GetHilightFontStruct()
{
   if (!fgHilightFont)
      fgHilightFont = gClient->GetResourcePool()->GetMenuHiliteFont();
   return fgHilightFont->GetFontStruct();
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGMenuTitle member functions.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a menu title. This object is created by a menu bar when adding
/// a popup menu. The menu title adopts the hotstring.

TGMenuTitle::TGMenuTitle(const TGWindow *p, TGHotString *s, TGPopupMenu *menu,
                         GContext_t norm, FontStruct_t font, UInt_t options)
    : TGFrame(p, 1, 1, options)
{
   fLabel      = s;
   fMenu       = menu;
   fFontStruct = font;
   fSelGC      = GetDefaultSelectedGC()();
   fNormGC     = norm;
   fState      = kFALSE;
   fTitleId    = -1;
   fTextColor  = GetForeground();
   fTitleData  = 0;

   Int_t hotchar;
   if (s && (hotchar = s->GetHotChar()) != 0)
      fHkeycode = gVirtualX->KeysymToKeycode(hotchar);
   else
      fHkeycode = 0;

   UInt_t tw = 0;
   Int_t  max_ascent, max_descent;
   if (fLabel)
      tw = gVirtualX->TextWidth(fFontStruct, fLabel->GetString(), fLabel->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   Resize(tw + 8, max_ascent + max_descent + 7);

   if (p && p->InheritsFrom(TGMenuBar::Class())) {
      TGMenuBar *bar = (TGMenuBar*)p;
      fMenu->SetMenuBar(bar);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set state of menu title.

void TGMenuTitle::SetState(Bool_t state)
{
   fState = state;
   if (state) {
      if (fMenu != 0) {
         Int_t    ax, ay;
         Window_t wdummy;

         gVirtualX->TranslateCoordinates(fId, (fMenu->GetParent())->GetId(),
                                         0, 0, ax, ay, wdummy);

         // place the menu just under the window:
#ifdef R__HAS_COCOA
         gVirtualX->SetWMTransientHint(fMenu->GetId(), GetId());
#endif
         fMenu->PlaceMenu(ax-1, ay+fHeight, kTRUE, kFALSE); //kTRUE);
      }
   } else {
      if (fMenu != 0) {
         fTitleId = fMenu->EndMenu(fTitleData);
      }
   }
   fOptions &= ~(kSunkenFrame | kRaisedFrame);
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a menu title.

void TGMenuTitle::DoRedraw()
{
   TGFrame::DoRedraw();

   int x, y, max_ascent, max_descent;
   x = y = 4;

   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   if (fState) {
      gVirtualX->SetForeground(fNormGC, GetDefaultSelectedBackground());
      if (gClient->GetStyle() > 1) {
         gVirtualX->FillRectangle(fId, fNormGC, 1, 2, fWidth-3, fHeight-4);
         gVirtualX->DrawRectangle(fId, GetShadowGC()(), 1, 1, fWidth-3, fHeight-3);
      }
      else {
         gVirtualX->FillRectangle(fId, fNormGC, 0, 0, fWidth, fHeight);
      }
      gVirtualX->SetForeground(fNormGC, GetForeground());
      fLabel->Draw(fId, fSelGC, x, y + max_ascent);
   } else {
      // Use same background color than the menu bar
      Pixel_t back = GetDefaultFrameBackground();
      if (fMenu && fMenu->fMenuBar && fMenu->fMenuBar->GetBackground() != back)
         back = fMenu->fMenuBar->GetBackground();
      gVirtualX->SetForeground(fNormGC, back);
      if (gClient->GetStyle() > 1) {
         gVirtualX->DrawRectangle(fId, fNormGC, 1, 1, fWidth-3, fHeight-3);
         gVirtualX->FillRectangle(fId, fNormGC, 1, 2, fWidth-3, fHeight-4);
      }
      else {
         gVirtualX->FillRectangle(fId, fNormGC, 0, 0, fWidth, fHeight);
      }
      gVirtualX->SetForeground(fNormGC, fTextColor);
      fLabel->Draw(fId, fNormGC, x, y + max_ascent);
      if (fTextColor != GetForeground())
         gVirtualX->SetForeground(fNormGC, GetForeground());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send final selected menu item to be processed.

void TGMenuTitle::DoSendMessage()
{
   if (fMenu)
      if (fTitleId != -1) {
         SendMessage(fMenu->fMsgWindow, MK_MSG(kC_COMMAND, kCM_MENU), fTitleId,
                     (Long_t)fTitleData);
         fMenu->Activated(fTitleId);
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure in use.

FontStruct_t TGMenuTitle::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetMenuFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context in use.

const TGGC &TGMenuTitle::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default selection graphics context in use.

const TGGC &TGMenuTitle::GetDefaultSelectedGC()
{
   if (!fgDefaultSelectedGC)
      fgDefaultSelectedGC = gClient->GetResourcePool()->GetSelectedGC();
   return *fgDefaultSelectedGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a popup menu widget as a C++ statement(s) on output stream out.

void TGPopupMenu::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   out << "   TGPopupMenu *";
   out << GetName() << " = new TGPopupMenu(gClient->GetDefaultRoot()"
       << "," << GetWidth() << "," << GetHeight() << "," << GetOptionString() << ");" << std::endl;

   Bool_t hasradio = kFALSE;
   Int_t r_first, r_last, r_active;
   r_active = r_first = r_last = -1;

   TGMenuEntry *mentry;
   TIter next(GetListOfEntries());

   while ((mentry = (TGMenuEntry *) next())) {
      const char *text;
      Int_t i, lentext, hotpos;
      char shortcut[80];
      char *outext;

      switch (mentry->GetType()) {
         case kMenuEntry:
            text = mentry->GetName();
            lentext = mentry->fLabel->GetLength();
            hotpos = mentry->fLabel->GetHotPos();
            outext = new char[lentext+2];
            i=0;
            while (text && lentext) {
               if (i == hotpos-1) {
                  outext[i] = '&';
                  i++;
               }
               outext[i] = *text;
               i++; text++; lentext--;
            }
            outext[i]=0;
            if (mentry->fShortcut) {
               snprintf(shortcut, 80, "\\t%s", mentry->GetShortcutText());
            } else {
               memset(shortcut, 0, 80);
            }

            {
               TString entrytext = gSystem->UnixPathName(outext);
               gSystem->ExpandPathName(entrytext);
               out << "   " << GetName() << "->AddEntry(" << quote
                   << entrytext // can be a file name
                   << shortcut
                   << quote << "," << mentry->GetEntryId();
            }
            if (mentry->fUserData) {
               out << "," << mentry->fUserData;
            }
            if (mentry->fPic) {
               TString picname = gSystem->UnixPathName(mentry->fPic->GetName());
               gSystem->ExpandPathName(picname);
               out << ",gClient->GetPicture(" << quote << picname << quote << ")";
            }
            out << ");" << std::endl;
            delete [] outext;
            break;
         case kMenuPopup:
            out << std::endl;
            out << "   // cascaded menu " << quote << mentry->GetName() << quote <<std::endl;
            mentry->fPopup->SavePrimitive(out, option);
            text = mentry->GetName();
            lentext = mentry->fLabel->GetLength();
            hotpos = mentry->fLabel->GetHotPos();
            outext = new char[lentext+2];
            i=0;
            while (text && lentext) {
               if (i == hotpos-1) {
                  outext[i] = '&';
                  i++;
               }
               outext[i] = *text;
               i++; text++; lentext--;
            }
            outext[i]=0;

            out << "   " << GetName() << "->AddPopup(" << quote
                << outext << quote << "," << mentry->fPopup->GetName()
                << ");" << std::endl;
            delete [] outext;
            break;
         case kMenuLabel:
            out << "   " << GetName() << "->AddLabel(" << quote
                << mentry->GetName() << quote;
            if (mentry->fPic) {
               out << ",gClient->GetPicture(" << quote
                   << mentry->fPic->GetName()
                   << quote << ")";
            }
            out << ");" << std::endl;
            break;
         case kMenuSeparator:
            out << "   " << GetName() << "->AddSeparator();" << std::endl;
            break;
      }

      if (!(mentry->GetStatus() & kMenuEnableMask)) {
         out<< "   " << GetName() << "->DisableEntry(" << mentry->GetEntryId()
            << ");" << std::endl;
      }
      if (mentry->GetStatus() & kMenuHideMask) {
         out<< "   " << GetName() << "->HideEntry(" << mentry->GetEntryId()
            << ");" << std::endl;
      }
      if (mentry->GetStatus() & kMenuCheckedMask) {
         out<< "   " << GetName() << "->CheckEntry(" << mentry->GetEntryId()
            << ");" << std::endl;
      }
      if (mentry->GetStatus() & kMenuDefaultMask) {
         out<< "   "<< GetName() << "->DefaultEntry(" << mentry->GetEntryId()
            << ");" << std::endl;
      }
      if (mentry->GetStatus() & kMenuRadioEntryMask) {
         if (hasradio) {
            r_last = mentry->GetEntryId();
            if (IsEntryRChecked(mentry->GetEntryId())) r_active = mentry->GetEntryId();
         }
         else {
            r_first = mentry->GetEntryId();
            hasradio = kTRUE;
            if (IsEntryRChecked(mentry->GetEntryId())) r_active = mentry->GetEntryId();
         }
      } else if (hasradio) {
         out << "   " << GetName() << "->RCheckEntry(" << r_active << "," << r_first
             << "," << r_last << ");" << std::endl;
         hasradio = kFALSE;
         r_active = r_first = r_last = -1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a title menu widget as a C++ statement(s) on output stream out.

void TGMenuTitle::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   out << std::endl;
   out << "   // " << quote << fLabel->GetString() << quote <<" menu" << std::endl;

   fMenu->SavePrimitive(out, option);

   const char *text = fLabel->GetString();
   Int_t lentext = fLabel->GetLength();
   Int_t hotpos = fLabel->GetHotPos();
   char *outext = new char[lentext+2];
   Int_t i=0;
   while (lentext) {
      if (i == hotpos-1) {
         outext[i] = '&';
         i++;
      }
      outext[i] = *text;
      i++; text++; lentext--;
   }
   outext[i]=0;
   out << "   " << fParent->GetName() << "->AddPopup(" << quote << outext
       << quote << "," << fMenu->GetName();

   delete [] outext;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a menu bar widget as a C++ statement(s) on output stream out.

void TGMenuBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   out << std::endl;
   out << "   // menu bar" << std::endl;

   out << "   TGMenuBar *";
   out << GetName() << " = new TGMenuBar(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight() << "," << GetOptionString() << ");" << std::endl;
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (!fList) return;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *)next())) {
      el->fFrame->SavePrimitive(out, option);
      el->fLayout->SavePrimitive(out, option);
      out << ");" << std::endl;
   }
}
