// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   20/08/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************

    This file is part of TGMdi an extension to the xclass toolkit.
    Copyright (C) 1998-2002 by Harald Radke, Hector Peraza.

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This application is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

**************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGMdiDecorFrame, TGMdiTitleBar, TGMdiButtons, TGMdiTitleIcon,        //
// TGMdiWinResizer, TGMdiVerticalWinResizer, TGMdiHorizontalWinResizer, //
// and TGMdiCornerWinResizer.                                           //
//                                                                      //
// This file contains all different MDI frame decoration classes.       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include "KeySymbols.h"
#include "TGResourcePool.h"
#include "TGDimension.h"
#include "WidgetMessageTypes.h"
#include "TGMdiFrame.h"
#include "TGMdiDecorFrame.h"
#include "TGMdiMainFrame.h"
#include "TVirtualX.h"


ClassImp(TGMdiDecorFrame);
ClassImp(TGMdiTitleBar);
ClassImp(TGMdiButtons);
ClassImp(TGMdiTitleIcon);
ClassImp(TGMdiWinResizer);
ClassImp(TGMdiVerticalWinResizer);
ClassImp(TGMdiHorizontalWinResizer);
ClassImp(TGMdiCornerWinResizer);

////////////////////////////////////////////////////////////////////////////////
/// TGMdiDecorFrame constructor.
/// The TGMdiDecorFrame is the frame containing MDI decorations like
/// title bar, minimize, maximize, restore and close buttons, and resizers.

TGMdiDecorFrame::TGMdiDecorFrame(TGMdiMainFrame *main, TGMdiFrame *frame,
                                 Int_t w, Int_t h, const TGGC *boxGC,
                                 UInt_t options, Pixel_t back) :
  TGCompositeFrame(main->GetContainer(), w, h,
                   options | kOwnBackground | kVerticalFrame | kFixedSize, back)
{

   fMdiMainFrame = main;
   fEditDisabled = 1;
   fFrame = frame;
   fIsMinimized = fIsMaximized = kFALSE;
   fMinimizedX = fMinimizedY = 0;
   fMinimizedUserPlacement = kFALSE;
   fButtonMask = kMdiDefaultHints;
   SetCleanup(kDeepCleanup);

   SetDecorBorderWidth(kMdiBorderWidth);

   fTitlebar = new TGMdiTitleBar(this, fMdiMainFrame);

   fLHint = new TGLayoutHints(kLHintsExpandX);
   fExpandHint = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);

   AddFrame(fTitlebar, fLHint);
   fTitlebar->LayoutButtons(fButtonMask, kFALSE, kFALSE);

   fUpperHR = new TGMdiVerticalWinResizer(this, main, kMdiResizerTop,
                                         boxGC, kMdiBorderWidth);
   fLowerHR = new TGMdiVerticalWinResizer(this, main, kMdiResizerBottom,
                                         boxGC, kMdiBorderWidth);
   fLeftVR = new TGMdiHorizontalWinResizer(this, main, kMdiResizerLeft,
                                          boxGC, kMdiBorderWidth);
   fRightVR = new TGMdiHorizontalWinResizer(this, main, kMdiResizerRight,
                                           boxGC, kMdiBorderWidth);

   fUpperLeftCR = new TGMdiCornerWinResizer(this, main,
                                kMdiResizerTop | kMdiResizerLeft,
                                boxGC, kMdiBorderWidth);
   fLowerLeftCR = new TGMdiCornerWinResizer(this, main,
                                kMdiResizerBottom | kMdiResizerLeft,
                                boxGC, kMdiBorderWidth);
   fUpperRightCR = new TGMdiCornerWinResizer(this, main,
                                 kMdiResizerTop | kMdiResizerRight,
                                 boxGC, kMdiBorderWidth);
   fLowerRightCR = new TGMdiCornerWinResizer(this, main,
                                 kMdiResizerBottom | kMdiResizerRight,
                                 boxGC, kMdiBorderWidth);

   fUpperHR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fLowerHR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fLeftVR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fRightVR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fUpperLeftCR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fLowerLeftCR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fUpperRightCR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);
   fLowerRightCR->SetMinSize(50, fTitlebar->GetDefaultHeight() + 2 * fBorderWidth);

   AddInput(kStructureNotifyMask | kButtonPressMask);

   fFrame->ReparentWindow(this, fBorderWidth, fTitlebar->GetDefaultHeight());
   fFrame->fParent = this;
   AddFrame(fFrame, fExpandHint);

   MapSubwindows();
   Resize(GetDefaultSize());
   Layout();

   MapWindow();
   TGFrame::SetWindowName();

   fFrame->RaiseWindow();
   fTitlebar->RaiseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiDecorFrame destructor.

TGMdiDecorFrame::~TGMdiDecorFrame()
{
   if (!MustCleanup()) {
      delete fUpperHR;
      delete fLowerHR;
      delete fLeftVR;
      delete fRightVR;
      delete fUpperLeftCR;
      delete fLowerLeftCR;
      delete fUpperRightCR;
      delete fLowerRightCR;
   }
   DestroyWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Set border width of the decor.

void TGMdiDecorFrame::SetDecorBorderWidth(Int_t bw)
{
   fBorderWidth = bw;
}

////////////////////////////////////////////////////////////////////////////////
/// Set-up MDI buttons.

void TGMdiDecorFrame::SetMdiButtons(ULong_t buttons)
{
   fButtonMask = buttons;
   fTitlebar->LayoutButtons(fButtonMask, fIsMinimized, fIsMaximized);
   if (fButtonMask & kMdiSize) {
      fUpperHR->Activate(kTRUE);
      fLowerHR->Activate(kTRUE);
      fLeftVR->Activate(kTRUE);
      fRightVR->Activate(kTRUE);
      fUpperLeftCR->Activate(kTRUE);
      fLowerLeftCR->Activate(kTRUE);
      fUpperRightCR->Activate(kTRUE);
      fLowerRightCR->Activate(kTRUE);
   } else {
      fUpperHR->Activate(kFALSE);
      fLowerHR->Activate(kFALSE);
      fLeftVR->Activate(kFALSE);
      fRightVR->Activate(kFALSE);
      fUpperLeftCR->Activate(kFALSE);
      fLowerLeftCR->Activate(kFALSE);
      fUpperRightCR->Activate(kFALSE);
      fLowerRightCR->Activate(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set resize mode (opaque or transparent)

void TGMdiDecorFrame::SetResizeMode(Int_t mode)
{
   fUpperHR->SetResizeMode(mode);
   fLowerHR->SetResizeMode(mode);
   fLeftVR->SetResizeMode(mode);
   fRightVR->SetResizeMode(mode);
   fUpperLeftCR->SetResizeMode(mode);
   fLowerLeftCR->SetResizeMode(mode);
   fUpperRightCR->SetResizeMode(mode);
   fLowerRightCR->SetResizeMode(mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculates the postion and the size of all decor frame components.

void TGMdiDecorFrame::Layout()
{
   RemoveInput(kStructureNotifyMask);
   TGCompositeFrame::Layout();
   AddInput(kStructureNotifyMask);

   if (fIsMaximized == kFALSE) {
      fUpperLeftCR->Move(0, 0);
      fUpperRightCR->Move(fWidth - fUpperRightCR->GetWidth(), 0);
      fLowerLeftCR->Move(0, fHeight - fLowerLeftCR->GetHeight());
      fLowerRightCR->Move(fWidth - fLowerRightCR->GetWidth(),
                          fHeight - fLowerRightCR->GetHeight());

      fLeftVR->MoveResize(0, fUpperLeftCR->GetHeight(), fLeftVR->GetWidth(),
                          fHeight - fUpperLeftCR->GetHeight() -
                          fLowerLeftCR->GetHeight());
      fUpperHR->MoveResize(fUpperLeftCR->GetWidth(), 0,
                           fWidth - fUpperRightCR->GetWidth() -
                           fUpperLeftCR->GetWidth(), fUpperHR->GetHeight());
      fRightVR->MoveResize(fWidth - fRightVR->GetWidth(),
                           fUpperRightCR->GetHeight(), fRightVR->GetWidth(),
                           fHeight - fUpperLeftCR->GetHeight() -
                           fLowerLeftCR->GetHeight());
      fLowerHR->MoveResize(fLowerLeftCR->GetWidth(), fHeight -
                           fLowerHR->GetHeight(),
                           fWidth - fLowerRightCR->GetWidth() -
                           fLowerLeftCR->GetWidth(), fLowerHR->GetHeight());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set MDI Window name (appearing in the title bar)

void TGMdiDecorFrame::SetWindowName(const char *name)
{
   fTitlebar->GetWinName()->SetText(new TGString(name));
   fTitlebar->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set Window icon (appearing in the title bar)

void TGMdiDecorFrame::SetWindowIcon(const TGPicture *icon)
{
   fTitlebar->GetWinIcon()->SetPicture(icon);
   fClient->NeedRedraw(fTitlebar->GetWinIcon());
}

////////////////////////////////////////////////////////////////////////////////
/// Move the MDI window at position x, y.

void TGMdiDecorFrame::Move(Int_t x, Int_t y)
{
   if (x < 0) {
      fMdiMainFrame->SetHsbPosition(fMdiMainFrame->GetViewPort()->GetWidth());
   }
   if (y < 0) {
      fMdiMainFrame->SetVsbPosition(fMdiMainFrame->GetViewPort()->GetHeight());
   }
   TGCompositeFrame::Move(x, y);
   if (IsMinimized()) fMinimizedUserPlacement = kTRUE;
   if (IsMapped() && !IsMaximized()) fMdiMainFrame->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Move the MDI window at position x, y and set size to w, h.

void TGMdiDecorFrame::MoveResize(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   if (x < 0) {
      fMdiMainFrame->SetHsbPosition(fMdiMainFrame->GetViewPort()->GetWidth());
   }
   if (y < 0) {
      fMdiMainFrame->SetVsbPosition(fMdiMainFrame->GetViewPort()->GetHeight());
   }
   TGCompositeFrame::MoveResize(x, y, w, h);
   if (IsMinimized()) fMinimizedUserPlacement = kTRUE;
   if (IsMapped() && !IsMaximized()) fMdiMainFrame->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle configure notify event.

Bool_t TGMdiDecorFrame::HandleConfigureNotify(Event_t *event)
{
   if ((event->fX < 0) || (event->fY < 0) ||
       (event->fX + event->fWidth > fMdiMainFrame->GetViewPort()->GetWidth()) ||
       (event->fY + event->fHeight > fMdiMainFrame->GetViewPort()->GetHeight())) {
      fMdiMainFrame->Resize();
   }

   if (event->fWindow == fFrame->GetId()) {
      UInt_t newW = event->fWidth + 2 * fBorderWidth;
      UInt_t newH = event->fHeight + 2 * fBorderWidth +
                    fTitlebar->GetDefaultHeight();

      if ((fWidth != newW) || (fHeight != newH)) {
         Resize(newW, newH);
      }
      fMdiMainFrame->Layout();
      return kTRUE;
   }
   fMdiMainFrame->Layout();
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button events.

Bool_t TGMdiDecorFrame::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress) {
      void *ud;
      fTitlebar->GetWinIcon()->GetPopup()->EndMenu(ud);
      SendMessage(fMdiMainFrame, MK_MSG(kC_MDI, kMDI_CURRENT), fId, 0);
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiTitleBar constructor.
/// the TGMdiTitleBar is the frame containing a title (window name)
/// an icon and MDI picture buttons as minimize, maximize, restore,
/// close and help.

TGMdiTitleBar::TGMdiTitleBar(const TGWindow *p, const TGWindow *mdiwin,
                             const char *name) :
   TGCompositeFrame(p, 10, 10, kOwnBackground | kHorizontalFrame)
{
   fMdiWin = mdiwin;
   fEditDisabled = kTRUE;
   fWinName = 0;
   fMidButPressed = fLeftButPressed = fRightButPressed = kFALSE;

   AddInput(kButtonPressMask | kButtonReleaseMask | kButtonMotionMask);

   fLHint = new TGLayoutHints(kLHintsNormal);
   fLeftHint = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 1, 1, 1);
   //fMiddleHint = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 1, 1, 1);
   fMiddleHint = new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 1, 1, 1, 1);
   fRightHint = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 1, 2, 1, 1);

   fLFrame = new TGCompositeFrame(this, 10, 10, kHorizontalFrame);
   fMFrame = new TGCompositeFrame(this, 10, 10, kHorizontalFrame);
   fRFrame = new TGCompositeFrame(this, 10, 10, kHorizontalFrame);
   AddFrame(fLFrame, fLeftHint);
   AddFrame(fMFrame, fMiddleHint);
   AddFrame(fRFrame, fRightHint);

   fWinIcon = new TGMdiTitleIcon(fLFrame, this,
                 fClient->GetPicture("mdi_default.xpm"),
                 16, 16);
   fLFrame->AddFrame(fWinIcon, fLHint);

   fWinName = new TGLabel(fMFrame, new TGString(name));
   fWinName->SetTextJustify(kTextLeft);
   fMFrame->AddFrame(fWinName, fLHint);

   fButtons = new TGMdiButtons(fRFrame, this);
   fRFrame->AddFrame(fButtons, fLHint);

   MapWindow();
   MapSubwindows();
   Layout();
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiTitleBar destructor.

TGMdiTitleBar::~TGMdiTitleBar()
{
   if (!MustCleanup()) {
      delete fLHint;
      delete fLeftHint;
      delete fMiddleHint;
      delete fRightHint;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculates the position of every enabled (displayed) buttons.

void TGMdiTitleBar::LayoutButtons(UInt_t buttonmask,Bool_t isMinimized,
                                  Bool_t isMaximized)
{
   fWinIcon->GetPopup()->EnableEntry(kMdiMove);

   if (buttonmask & kMdiSize) {
      fWinIcon->GetPopup()->EnableEntry(kMdiSize);
   } else {
      fWinIcon->GetPopup()->DisableEntry(kMdiSize);
   }

   if (buttonmask & kMdiMenu) {
      fLFrame->ShowFrame(fWinIcon);
   } else {
      fLFrame->HideFrame(fWinIcon);
   }

   if (buttonmask & kMdiClose) {
      fButtons->ShowFrame(fButtons->GetButton(4));
      fWinIcon->GetPopup()->EnableEntry(kMdiClose);
   } else {
      fButtons->HideFrame(fButtons->GetButton(4));
      fWinIcon->GetPopup()->DisableEntry(kMdiClose);
   }

   if (buttonmask & kMdiHelp) {
      fButtons->ShowFrame(fButtons->GetButton(3));
   } else {
      fButtons->HideFrame(fButtons->GetButton(3));
   }

   if ((buttonmask & kMdiMaximize) && (!isMaximized)) {
      fButtons->ShowFrame(fButtons->GetButton(2));
      fWinIcon->GetPopup()->EnableEntry(kMdiMaximize);
   } else {
      fButtons->HideFrame(fButtons->GetButton(2));
      fWinIcon->GetPopup()->DisableEntry(kMdiMaximize);
   }

   if (isMinimized | isMaximized) {
      fButtons->ShowFrame(fButtons->GetButton(1));
      fWinIcon->GetPopup()->EnableEntry(kMdiRestore);
      fWinIcon->GetPopup()->DisableEntry(kMdiSize);
      if (isMaximized) fWinIcon->GetPopup()->DisableEntry(kMdiMove);
   } else {
      fButtons->HideFrame(fButtons->GetButton(1));
      fWinIcon->GetPopup()->DisableEntry(kMdiRestore);
   }

   if ((buttonmask & kMdiMinimize) && (!isMinimized)) {
      fButtons->ShowFrame(fButtons->GetButton(0));
      fWinIcon->GetPopup()->EnableEntry(kMdiMinimize);
   } else {
      fButtons->HideFrame(fButtons->GetButton(0));
      fWinIcon->GetPopup()->DisableEntry(kMdiMinimize);
   }

   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set title bar color (blue or grey, depends on active state).

void TGMdiTitleBar::SetTitleBarColors(UInt_t fore, UInt_t back, TGFont *font)
{
   SetBackgroundColor(back);

   fClient->GetFont(font->GetName());
   fWinName->SetTextFont(font);
   fWinName->SetTextColor(fore, kFALSE);
   fMFrame->SetBackgroundColor(back);
   fWinName->SetBackgroundColor(back);
   fWinIcon->SetBackgroundColor(back);
   fButtons->SetBackgroundColor(back);
   Layout();
   fClient->NeedRedraw(this);
   fClient->NeedRedraw(fWinName);
   fClient->NeedRedraw(fMFrame);
   fClient->NeedRedraw(fButtons);
   fClient->NeedRedraw(fWinIcon);
   fWinIcon->DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click in title bar (maximize window)

Bool_t TGMdiTitleBar::HandleDoubleClick(Event_t *event)
{
   if ((event->fType == kButtonPress) && (event->fCode == kButton1)) {
      SendMessage(fMdiWin, MK_MSG(kC_MDI, kMDI_MAXIMIZE), fParent->GetId(), 0);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse click on title bar.

Bool_t TGMdiTitleBar::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress) {
      void *ud;
      GetWinIcon()->GetPopup()->EndMenu(ud);
      gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kMove));
      switch (event->fCode) {
         case kButton1:
            fX0 = event->fX;
            fY0 = event->fY;
            fLeftButPressed = kTRUE;
            SendMessage(fMdiWin, MK_MSG(kC_MDI, kMDI_CURRENT), fParent->GetId(), 0);
            break;

         case kButton2:
            fMidButPressed = kTRUE;
            break;

         case kButton3:
            gVirtualX->LowerWindow(fParent->GetId());
            fRightButPressed = kTRUE;
            break;
      }
   } else {
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
      gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kPointer));
      switch (event->fCode) {

         case kButton1:
            //if (!fClient->IsEditable()) ((TGMainFrame *)((TGMdiMainFrame *)fMdiWin)->GetMainFrame())->Layout();
            fLeftButPressed = kFALSE;
            break;

         case kButton2:
            fMidButPressed = kFALSE;
            break;

         case kButton3:
            fRightButPressed = kFALSE;
            break;
      }

      TGFrame *f = GetFrameFromPoint(event->fX, event->fY);
      if (f && (f != this)) {
         TranslateCoordinates(f, event->fX, event->fY, event->fX, event->fY);
         f->HandleButton(event);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages for title bar.

Bool_t TGMdiTitleBar::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
            case kCM_MENU:
               SendMessage(fMdiWin, MK_MSG(kC_MDI, (EWidgetMessageTypes)parm1),
                           fParent->GetId(), parm2);
               break;
         }
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events in title bar (used to move MDI window).

Bool_t TGMdiTitleBar::HandleMotion(Event_t *event)
{
   if (event->fWindow != fId) return kTRUE;
   if (!fLeftButPressed) return kTRUE;

   Int_t x = ((TGFrame *)fParent)->GetX();
   Int_t y = ((TGFrame *)fParent)->GetY();
   ((TGFrame *)fParent)->Move(x + event->fX - fX0, y + event->fY - fY0);

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// This is called from TGMdiMainFrame on Restore().

void TGMdiTitleBar::AddFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   icon->ReparentWindow(fLFrame);
   buttons->ReparentWindow(fRFrame);
   fLFrame->AddFrame(icon, fLHint);
   fLFrame->ShowFrame(icon);
   fRFrame->AddFrame(buttons, fLHint);
   fRFrame->ShowFrame(buttons);
}


////////////////////////////////////////////////////////////////////////////////
/// This is called from TGMdiMainFrame on Maximize().

void TGMdiTitleBar::RemoveFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   fLFrame->RemoveFrame(icon);
   fRFrame->RemoveFrame(buttons);
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiButtons constructor.
/// the TGMdiButtons is the frame containing MDI picture buttons like
/// minimize, maximize, restore, close and help.

TGMdiButtons::TGMdiButtons(const TGWindow *p, const TGWindow *titlebar) :
   TGCompositeFrame(p, 10, 10, kHorizontalFrame)
{
   fDefaultHint = new TGLayoutHints(kLHintsNormal, 0, 0, 1, 0);
   fCloseHint = new TGLayoutHints(kLHintsNormal, 2, 0, 1, 0);
   fEditDisabled = kTRUE;
   fMsgWindow = 0;

   //--- Minimize button

   fButton[0] = new TGPictureButton(this,
                     fClient->GetPicture("mdi_minimize.xpm"),
                     kMdiMinimize);
   fButton[0]->SetToolTipText("Minimize");
   AddFrame(fButton[0], fDefaultHint);
   fButton[0]->SetBackgroundColor(GetDefaultFrameBackground());
   fButton[0]->Associate(titlebar);

   //--- Restore button

   fButton[1] = new TGPictureButton(this,
                     fClient->GetPicture("mdi_restore.xpm"),
                     kMdiRestore);
   fButton[1]->SetToolTipText("Restore");
   AddFrame(fButton[1], fDefaultHint);
   fButton[1]->SetBackgroundColor(GetDefaultFrameBackground());
   fButton[1]->Associate(titlebar);

   //--- Maximize button

   fButton[2] = new TGPictureButton(this,
                     fClient->GetPicture("mdi_maximize.xpm"),
                     kMdiMaximize);
   fButton[2]->SetToolTipText("Maximize");
   AddFrame(fButton[2], fDefaultHint);
   fButton[2]->SetBackgroundColor(GetDefaultFrameBackground());
   fButton[2]->Associate(titlebar);

   //--- Help button

   fButton[3] = new TGPictureButton(this,
                     fClient->GetPicture("mdi_help.xpm"),
                     kMdiHelp);
   fButton[3]->SetToolTipText("Help");
   AddFrame(fButton[3], fDefaultHint);
   fButton[3]->SetBackgroundColor(GetDefaultFrameBackground());
   fButton[3]->Associate(titlebar);

   //--- Close button

   fButton[4] = new TGPictureButton(this,
                     fClient->GetPicture("mdi_close.xpm"),
                     kMdiClose);
   fButton[4]->SetToolTipText("Close");
   AddFrame(fButton[4], fCloseHint);
   fButton[4]->SetBackgroundColor(GetDefaultFrameBackground());
   fButton[4]->Associate(titlebar);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiButtons destructor.

TGMdiButtons::~TGMdiButtons()
{
   if (!MustCleanup()) {
      delete fDefaultHint;
      delete fCloseHint;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiTitleIcon constructor.
/// the TGMdiTitleIcon is the left icon used also for the MDI
/// popup menu allowing access to MDI commands as : restore,
/// move, size, minimize and close.

TGMdiTitleIcon::TGMdiTitleIcon(const TGWindow *p, const TGWindow *titlebar,
                               const TGPicture *pic, Int_t w, Int_t h) :
   TGIcon(p, pic, w, h)
{
   fMsgWindow = titlebar;
   fEditDisabled = kTRUE;

   //--- MDI system menu

   fPopup = new TGPopupMenu(fClient->GetDefaultRoot());// fClient->GetRoot());
   fPopup->AddEntry(new TGHotString("&Restore"), kMdiRestore);
   fPopup->AddEntry(new TGHotString("&Move"), kMdiMove);
   fPopup->AddEntry(new TGHotString("&Size"), kMdiSize);
   fPopup->AddEntry(new TGHotString("Mi&nimize"), kMdiMinimize);
   fPopup->AddEntry(new TGHotString("Ma&ximize"), kMdiMaximize);
   fPopup->AddSeparator();
   fPopup->AddEntry(new TGHotString("&Close  Ctrl+F4"), kMdiClose);
   fPopup->DisableEntry(kMdiRestore);
   fPopup->Associate(titlebar);

   AddInput(kButtonPressMask | kButtonReleaseMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiTitleIcon destructor.

TGMdiTitleIcon::~TGMdiTitleIcon()
{
   delete fPopup;
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw icon.

void TGMdiTitleIcon::DoRedraw()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);
   TGIcon::DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click event on MDI icon (close the window)

Bool_t TGMdiTitleIcon::HandleDoubleClick(Event_t *event)
{
   if (event->fCode == kButton1) {
      void *ud;
      fPopup->EndMenu(ud);
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_MENU), kMDI_CLOSE, 0);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button event on MDI icon (popup menu)

Bool_t TGMdiTitleIcon::HandleButton(Event_t *event)
{
   SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_MENU), kMDI_CURRENT, 0);

   if (event->fType == kButtonPress) {
      fPopup->PlaceMenu(event->fXRoot - event->fX,
                        event->fYRoot + (fHeight - event->fY),
                        kTRUE, kTRUE);
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiWinResizer constructor.
/// The TGMdiWinResizer is a frame allowing to resize MDI window.
/// Could be horizontal, vertical or corner resizer (see derived classes
/// TGMdiVerticalWinResizer, TGMdiHorizontalWinResizer, and
/// TGMdiCornerWinResizer).

TGMdiWinResizer::TGMdiWinResizer(const TGWindow *p, const TGWindow *mdiwin,
                   Int_t pos, const TGGC *boxgc, Int_t linew,
                   Int_t mdioptions, Int_t w, Int_t h, UInt_t options) :
   TGFrame(p, w, h, options)
{
   fWinX = fWinY = fWinW = fWinH = fOldX = fOldY = fOldW = fOldH = 0;
   fNewX = fNewY = fNewW = fNewH = fX0 = fY0 = 0;

   fWidgetFlags = kWidgetIsEnabled;

   fMdiWin = mdiwin;
   fMdiOptions = mdioptions;
   fPos = pos;

   fBoxGC = boxgc;
   fLineW = linew;

   fMinW = 50;
   fMinH = 20;

   fMidButPressed = fLeftButPressed = fRightButPressed = kFALSE;

   gVirtualX->GrabButton(fId, kButton1, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask | kButtonMotionMask,
                         kNone, kNone);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button events in resizer (grab button and resize).

Bool_t TGMdiWinResizer::HandleButton(Event_t *event)
{
   if (!IsEnabled()) return kTRUE;

   if (event->fType == kButtonPress) {
      void *ud;
      ((TGMdiDecorFrame *)fParent)->GetTitleBar()->GetWinIcon()->GetPopup()->EndMenu(ud);
      switch (event->fCode) {
         case kButton1:
            SendMessage(fMdiWin, MK_MSG(kC_MDI, kMDI_CURRENT), fParent->GetId(), 0);
            fNewX = fOldX = fWinX = ((TGFrame *)fParent)->GetX();
            fNewY = fOldY = fWinY = ((TGFrame *)fParent)->GetY();
            fWinW = ((TGFrame *)fParent)->GetWidth();
            fWinH = ((TGFrame *)fParent)->GetHeight();
            fX0 = event->fXRoot;
            fY0 = event->fYRoot;
            fNewW = fWinW;
            fNewH = fWinH;
            if (fMdiOptions != kMdiOpaque) {
               DrawBox(fNewX, fNewY, fNewW, fNewH);
            }
            fLeftButPressed = kTRUE;
            gVirtualX->GrabPointer(fId, kButtonReleaseMask | kPointerMotionMask,
                                   kNone, kNone, kTRUE, kFALSE);
            break;

         case kButton2:
            fMidButPressed = kTRUE;
            break;

         case kButton3:
            fRightButPressed = kTRUE;
            break;
      }
   } else {
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
      switch (event->fCode) {
         case kButton1:
            if (fMdiOptions != kMdiOpaque) {
               DrawBox(fNewX, fNewY, fNewW, fNewH);
               ((TGFrame *)fParent)->MoveResize(fNewX, fNewY, fNewW, fNewH);
            }
            //if (!fClient->IsEditable()) ((TGMainFrame *)((TGMdiMainFrame *)fMdiWin)->GetMainFrame())->Layout();
            fLeftButPressed = kFALSE;
            break;

         case kButton2:
            fMidButPressed = kFALSE;
            break;

         case kButton3:
            fRightButPressed = kFALSE;
            break;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw resize box (rectangle).

void TGMdiWinResizer::DrawBox(Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   TGMdiMainFrame *m = (TGMdiMainFrame *) fMdiWin;

   gVirtualX->DrawRectangle(m->GetContainer()->GetId(), fBoxGC->GetGC(),
       x + fLineW / 2, y + fLineW / 2, width - fLineW, height - fLineW);
}

////////////////////////////////////////////////////////////////////////////////
/// Move (resize) parent MDI window.

void TGMdiWinResizer::MoveResizeIt()
{
   if (fMdiOptions == kMdiOpaque) {
      ((TGFrame *)fParent)->MoveResize(fNewX, fNewY, fNewW, fNewH);
   } else {
      DrawBox(fOldX, fOldY, fOldW, fOldH);
      DrawBox(fNewX, fNewY, fNewW, fNewH);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiVerticalWinResizer constructor.

TGMdiVerticalWinResizer::TGMdiVerticalWinResizer(const TGWindow *p,
             const TGWindow *mdiwin, Int_t pos, const TGGC *boxGC, Int_t linew,
             Int_t mdioptions, Int_t w, Int_t h) :
   TGMdiWinResizer(p, mdiwin, pos, boxGC, linew, mdioptions,
                   w, h, kFixedHeight | kOwnBackground)
{
   gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kArrowVer));
}

////////////////////////////////////////////////////////////////////////////////
/// Handle motion events in resizer (resize associated MDI window).

Bool_t TGMdiVerticalWinResizer::HandleMotion(Event_t *event)
{
   if (((TGMdiDecorFrame *)fParent)->IsMinimized()) return kTRUE;

   fOldX = fNewX;
   fOldY = fNewY;
   fOldW = fNewW;
   fOldH = fNewH;

   Int_t dy = event->fYRoot - fY0;

   if (!fLeftButPressed) return kTRUE;

   switch (fPos) {
      case kMdiResizerTop:
         if (fWinH - dy < fMinH) dy = fWinH - fMinH;
         fNewY = fWinY + dy;
         fNewH = fWinH - dy;
         break;

      case kMdiResizerBottom:
         if (fWinH + dy < fMinH) dy = fMinH - fWinH;
         fNewY = fWinY;
         fNewH = fWinH + dy;
         break;
   }

   MoveResizeIt();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw vertical resizer frame border.

void TGMdiVerticalWinResizer::DrawBorder()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);
   if (fPos == kMdiResizerTop) {
      gVirtualX->DrawLine(fId, GetHilightGC()(), 0, 1, fWidth - 1, 1);
   } else {
      gVirtualX->DrawLine(fId, GetShadowGC()(), 0, fHeight - 2, fWidth - 1,
          fHeight - 2);
      gVirtualX->DrawLine(fId, GetBlackGC()(), 0, fHeight - 1, fWidth - 1,
          fHeight - 1);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiCornerWinResizer constructor.

TGMdiCornerWinResizer::TGMdiCornerWinResizer(const TGWindow *p,
           const TGWindow *mdiwin, Int_t pos, const TGGC *boxGC, Int_t linew,
           Int_t mdioptions, Int_t w, Int_t h) :
   TGMdiWinResizer(p, mdiwin, pos, boxGC, linew, mdioptions,
                   w, h, kFixedSize | kOwnBackground)
{
   Cursor_t defaultCursor = kNone;
   fEditDisabled = kTRUE;

   switch (fPos) {
      case (kMdiResizerTop | kMdiResizerLeft):
         defaultCursor = gVirtualX->CreateCursor(kTopLeft);
         break;

      case (kMdiResizerBottom | kMdiResizerLeft):
         defaultCursor = gVirtualX->CreateCursor(kBottomLeft);
         break;

      case (kMdiResizerTop | kMdiResizerRight):
         defaultCursor = gVirtualX->CreateCursor(kTopRight);
         break;

      case (kMdiResizerBottom | kMdiResizerRight):
         defaultCursor = gVirtualX->CreateCursor(kBottomRight);
         break;
   }
   gVirtualX->SetCursor(fId, defaultCursor);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle motion events in resizer (resize associated MDI window).

Bool_t TGMdiCornerWinResizer::HandleMotion(Event_t *event)
{
   if (((TGMdiDecorFrame *)fParent)->IsMinimized()) return kTRUE;

   fOldX = fNewX;
   fOldY = fNewY;
   fOldW = fNewW;
   fOldH = fNewH;

   Int_t dx = event->fXRoot - fX0;
   Int_t dy = event->fYRoot - fY0;

   if (!fLeftButPressed) return kTRUE;

   switch (fPos) {
      case (kMdiResizerTop | kMdiResizerLeft):
         if (fWinW - dx < fMinW) dx = fWinW - fMinW;
         if (fWinH - dy < fMinH) dy = fWinH - fMinH;
         fNewX = fWinX + dx;
         fNewW = fWinW - dx;
         fNewY = fWinY + dy;
         fNewH = fWinH - dy;
         break;

      case (kMdiResizerBottom | kMdiResizerLeft):
         if (fWinW - dx < fMinW) dx = fWinW - fMinW;
         if (fWinH + dy < fMinH) dy = fMinH - fWinH;
         fNewX = fWinX + dx;
         fNewW = fWinW - dx;
         fNewY = fWinY;
         fNewH = fWinH + dy;
         break;

      case (kMdiResizerTop | kMdiResizerRight):
         if (fWinW + dx < fMinW) dx = fMinW - fWinW;
         if (fWinH - dy < fMinH) dy = fWinH - fMinH;
         fNewX = fWinX;
         fNewW = fWinW + dx;
         fNewY = fWinY + dy;
         fNewH = fWinH - dy;
         break;

      case (kMdiResizerBottom | kMdiResizerRight):
         if (fWinW + dx < fMinW) dx = fMinW - fWinW;
         if (fWinH + dy < fMinH) dy = fMinH - fWinH;
         fNewX = fWinX;
         fNewW = fWinW + dx;
         fNewY = fWinY;
         fNewH = fWinH + dy;
         break;
   }

   MoveResizeIt();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw corner resizer frame border.

void TGMdiCornerWinResizer::DrawBorder()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);

   switch (fPos) {
      case (kMdiResizerTop | kMdiResizerLeft):
         gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, fWidth - 1, 1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, 1, fHeight - 1);
         break;

      case (kMdiResizerBottom | kMdiResizerLeft):
         gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 0, 1, fHeight - 1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), 1, fHeight - 2,
             fWidth - 1, fHeight - 2);
         gVirtualX->DrawLine(fId, GetBlackGC()(), 0, fHeight - 1,
             fWidth - 1, fHeight - 1);
         break;

      case (kMdiResizerTop | kMdiResizerRight):
         gVirtualX->DrawLine(fId, GetHilightGC()(), 0, 1, fWidth - 1, 1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), fWidth - 2, 1,
             fWidth - 2, fHeight - 1);
         gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth - 1, 0,
             fWidth - 1, fHeight - 1);
         break;

      case (kMdiResizerBottom | kMdiResizerRight):
         gVirtualX->DrawLine(fId, GetShadowGC()(), fWidth - 2, 0,
             fWidth - 2, fHeight - 2);
         gVirtualX->DrawLine(fId, GetShadowGC()(), 0, fHeight - 2,
             fWidth - 1, fHeight - 2);
         gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth - 1, 0,
             fWidth - 1, fHeight - 1);
         gVirtualX->DrawLine(fId, GetBlackGC()(), 0, fHeight - 1,
             fWidth - 1, fHeight - 1);
         break;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGMdiHorizontalWinResizer constructor.

TGMdiHorizontalWinResizer::TGMdiHorizontalWinResizer(const TGWindow *p,
               const TGWindow *mdiwin, Int_t pos, const TGGC *boxGC, Int_t linew,
               Int_t mdioptions, Int_t w, Int_t h) :
   TGMdiWinResizer(p, mdiwin, pos, boxGC, linew, mdioptions,
                   w, h, kFixedWidth | kOwnBackground)
{
   fEditDisabled = kTRUE;
   gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kArrowHor));
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle motion events in resizer (resize associated MDI window).

Bool_t TGMdiHorizontalWinResizer::HandleMotion(Event_t *event)
{
   if (((TGMdiDecorFrame *)fParent)->IsMinimized()) return kTRUE;

   fOldX = fNewX;
   fOldY = fNewY;
   fOldW = fNewW;
   fOldH = fNewH;

   Int_t dx = event->fXRoot - fX0;

   if (!fLeftButPressed) return kTRUE;

   switch (fPos) {
      case (kMdiResizerLeft):
         if (fWinW - dx < fMinW) dx = fWinW - fMinW;
         fNewX = fWinX + dx;
         fNewW = fWinW - dx;
         break;

      case (kMdiResizerRight):
         if (fWinW + dx < fMinW) dx = fMinW - fWinW;
         fNewX = fWinX;
         fNewW = fWinW + dx;
         break;
   }

   MoveResizeIt();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw horizontal resizer frame border.

void TGMdiHorizontalWinResizer::DrawBorder()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);
   if (fPos == kMdiResizerLeft) {
      gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 0, 1, fHeight - 1);
   } else {
      gVirtualX->DrawLine(fId, GetShadowGC()(), fWidth - 2, 0, fWidth - 2, fHeight - 1);
      gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth - 1, 0, fWidth - 1, fHeight - 1);
   }
}
