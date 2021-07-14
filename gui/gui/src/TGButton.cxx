// @(#)root/gui:$Id: ee86415852b0e43b57190b9645717cf508b7920e $
// Author: Fons Rademakers   06/01/98

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


/** \class TGButton
    \ingroup guiwidgets

A button abstract base class. It defines general button behaviour.

*/


/** \class TGTextButton
    \ingroup guiwidgets

Yield an action as soon as it is clicked. This buttons usually provides fast access to
frequently used or critical commands. It may appear alone or placed in a group.

The action it performs can be inscribed with a meaningful tooltip
set by `SetToolTipText(const char* text, Long_t delayms=400).

The text button has a label indicating the action to be taken when
the button is pressed. The text can be a hot string ("&Exit") that
defines the label "Exit" and keyboard mnemonics Alt+E for button
selection. A button label can be changed by SetText(new_label).

Selecting a text or picture button will generate the event:
  - kC_COMMAND, kCM_BUTTON, button id, user data.
*/


/** \class TGPictureButton
    \ingroup guiwidgets

Yield an action as soon as it is clicked. This buttons usually provides fast access to
frequently used or critical commands. It may appear alone or placed in a group.

The action it performs can be inscribed with a meaningful tooltip
set by `SetToolTipText(const char* text, Long_t delayms=400).

The text button has a label indicating the action to be taken when
the button is pressed. The text can be a hot string ("&Exit") that
defines the label "Exit" and keyboard mnemonics Alt+E for button
selection. A button label can be changed by SetText(new_label).

Selecting a text or picture button will generate the event:
  - kC_COMMAND, kCM_BUTTON, button id, user data.
*/


/** \class TGCheckButton
    \ingroup guiwidgets

Selects different options. Like text buttons, they have text or hot string as a label.

Selecting a check button will generate the event:
  - kC_COMMAND, kCM_CHECKBUTTON, button id, user data.

If a command string has been specified (via SetCommand()) then this
command string will be executed via the interpreter whenever a
button is selected. A command string can contain the macros:

  - $MSG   -- kC_COMMAND, kCMCHECKBUTTON packed message
              (use GET_MSG() and GET_SUBMSG() to unpack)
  - $PARM1 -- button id
  - $PARM2 -- user data pointer

Before executing these macros are expanded into the respective Long_t's
*/


/** \class TGRadioButton
    \ingroup guiwidgets

Selects different options. Like text buttons, they have text or hot string as a label.

Radio buttons are grouped usually in logical sets of two or more
buttons to present mutually exclusive choices.

Selecting a radio button will generate the event:
  - kC_COMMAND, kCM_RADIOBUTTON, button id, user data.

If a command string has been specified (via SetCommand()) then this
command string will be executed via the interpreter whenever a
button is selected. A command string can contain the macros:

  - $MSG   -- kC_COMMAND, kCMRADIOBUTTON packed message
              (use GET_MSG() and GET_SUBMSG() to unpack)
  - $PARM1 -- button id
  - $PARM2 -- user data pointer

Before executing these macros are expanded into the respective Long_t's
*/


/** \class TGSplitButton
    \ingroup guiwidgets

Implements a button with added menu functionality.
There are 2 modes of operation available.

If the button is split, a menu will popup when the menu area of the
button is clicked. Activating a menu item changes the functionality
of the button by having it emit a additional signal when it is
clicked. The signal emitted when the button is clicked, is the
ItemClicked(Int_t) signal with a different fixed value for the
Int_t that corresponds to the id of the activated menu entry.

If the button is not split, clicking it will popup the menu and the
ItemClicked(Int_t) signal will be emitted when a menu entry is
activated. The value of the Int_t is again equal to the value of
the id of the activated menu entry.

The mode of operation of a SplitButton can be changed on the fly
by calling the SetSplit(Bool_t) method.
*/


#include "TGButton.h"
#include "TGWidget.h"
#include "TGPicture.h"
#include "TGToolTip.h"
#include "TGButtonGroup.h"
#include "TGResourcePool.h"
#include "TSystem.h"
#include "TImage.h"
#include "TEnv.h"
#include "TClass.h"
#include "TGMenu.h"
#include "KeySymbols.h"
#include "TVirtualX.h"

#include <iostream>


const TGGC *TGButton::fgHibckgndGC = nullptr;
const TGGC *TGButton::fgDefaultGC = nullptr;

const TGFont *TGTextButton::fgDefaultFont = nullptr;

const TGFont *TGCheckButton::fgDefaultFont = nullptr;
const TGGC   *TGCheckButton::fgDefaultGC = nullptr;

const TGFont *TGRadioButton::fgDefaultFont = nullptr;
const TGGC   *TGRadioButton::fgDefaultGC = nullptr;

Window_t TGButton::fgReleaseBtn = 0;

ClassImp(TGButton);
ClassImp(TGTextButton);
ClassImp(TGPictureButton);
ClassImp(TGCheckButton);
ClassImp(TGRadioButton);
ClassImp(TGSplitButton);

////////////////////////////////////////////////////////////////////////////////
/// Create button base class part.

TGButton::TGButton(const TGWindow *p, Int_t id, GContext_t norm, UInt_t options)
    : TGFrame(p, 1, 1, options)
{
   fWidgetId    = id;
   fWidgetFlags = kWidgetWantFocus;
   fMsgWindow   = p;
   fUserData    = 0;
   fTip         = 0;
   fGroup       = 0;
   fStyle       = 0;
   fTWidth = fTHeight = 0;

   fNormGC   = norm;
   fState    = kButtonUp;
   fStayDown = kFALSE;
   fWidgetFlags = kWidgetIsEnabled;

//   fStyle = gClient->GetStyle();
//   if (fStyle > 0) {
//      fOptions &= ~(kRaisedFrame | kDoubleBorder);
//   }

   // coverity[returned_null]
   // coverity[dereference]
   if (p && p->IsA()->InheritsFrom(TGButtonGroup::Class())) {
      TGButtonGroup *bg = (TGButtonGroup*) p;
      bg->Insert(this, id);
   }

   fBgndColor = fBackground;
   fHighColor = gClient->GetResourcePool()->GetHighLightColor();

   gVirtualX->GrabButton(fId, kButton1, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask,
                         kNone, kNone);

   AddInput(kEnterWindowMask | kLeaveWindowMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete button.

TGButton::~TGButton()
{
   // remove from button group
   if (fGroup) {
      fGroup->Remove(this);
      fGroup = 0;
   }

   if (fTip) delete fTip;
}

////////////////////////////////////////////////////////////////////////////////
/// Set button state.

void TGButton::SetState(EButtonState state, Bool_t emit)
{
   Bool_t was = !IsDown();   // kTRUE if button was off

   if (state == kButtonDisabled)
      fWidgetFlags &= ~kWidgetIsEnabled;
   else
      fWidgetFlags |= kWidgetIsEnabled;
   if (state != fState) {
      switch (state) {
         case kButtonEngaged:
         case kButtonDown:
            fOptions &= ~kRaisedFrame;
            fOptions |= kSunkenFrame;
            break;
         case kButtonDisabled:
         case kButtonUp:
            if (fStyle > 0) {
               fOptions &= ~kRaisedFrame;
               fOptions &= ~kSunkenFrame;
            }
            else {
               fOptions &= ~kSunkenFrame;
               fOptions |= kRaisedFrame;
            }
            break;
      }
      fState = state;
      DoRedraw();
      if (emit || fGroup) EmitSignals(was);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the button style (modern or classic).

void TGButton::SetStyle(UInt_t newstyle)
{
   fStyle = newstyle;
   fBgndColor = fBackground;
   if (fStyle > 0) {
      ChangeOptions(GetOptions() & ~kRaisedFrame);
   }
   else {
      ChangeOptions(GetOptions() | kRaisedFrame);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the button style (modern or classic).

void TGButton::SetStyle(const char *style)
{
   fBgndColor = fBackground;
   if (style && strstr(style, "modern")) {
      fStyle = 1;
      ChangeOptions(GetOptions() & ~kRaisedFrame);
   }
   else {
      fStyle = 0;
      ChangeOptions(GetOptions() | kRaisedFrame);
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGButton::IsDown() const
{
   if (fStyle > 0)
      return (fOptions & kSunkenFrame);
   return !(fOptions & kRaisedFrame);
}

////////////////////////////////////////////////////////////////////////////////

void TGButton::SetDown(Bool_t on, Bool_t emit)
{
   // Set button state down according to the parameter 'on'.

   if (GetState() == kButtonDisabled) return;

   SetState(on ? kButtonDown : kButtonUp, emit);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets new button-group for this button.

void TGButton::SetGroup(TGButtonGroup *group)
{
   fGroup = group;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event.

Bool_t TGButton::HandleButton(Event_t *event)
{
   Bool_t click = kFALSE;

   if (fTip) fTip->Hide();

   if (fState == kButtonDisabled) return kTRUE;

   Bool_t in = (event->fX >= 0) && (event->fY >= 0) &&
               (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   // We don't need to check the button number as GrabButton will
   // only allow button1 events
   if (event->fType == kButtonPress) {
      fgReleaseBtn = 0;

      if (fState == kButtonEngaged) {
         return kTRUE;
      }
      if (in) SetState(kButtonDown, kTRUE);
   } else { // ButtonRelease
      if (fState == kButtonEngaged) {
         if (in) SetState(kButtonUp, kTRUE);
         click = kTRUE;
      } else {
         click = (fState == kButtonDown) && in;
         if (click && fStayDown) {
            if (in) {
               SetState(kButtonEngaged, kTRUE);
               fgReleaseBtn = 0;
            }
         } else {
            if (in) {
               SetState(kButtonUp, kTRUE);
               fgReleaseBtn = fId;
            }
         }
      }
   }
   if (click) {
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }
   if ((fStyle > 0) && (event->fType == kButtonRelease)) {
      fBgndColor = fBackground;
   }
   DoRedraw();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit button signals.

void TGButton::EmitSignals(Bool_t was)
{
   Bool_t now = !IsDown();       // kTRUE if button now is off

   // emit signals
   if (was && !now) {
      Pressed();                 // emit Pressed  = was off , now on
      if (fStayDown) Clicked();  // emit Clicked
   }
   if (!was && now) {
      Released();                // emit Released = was on , now off
      Clicked();                 // emit Clicked
   }
   if ((was != now) && IsToggleButton()) Toggled(!now); // emit Toggled  = was != now
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGButton::HandleCrossing(Event_t *event)
{
   if (fTip) {
      if (event->fType == kEnterNotify)
         fTip->Reset();
      else
         fTip->Hide();
   }

   if (fStyle > 0) {
      if ((event->fType == kEnterNotify) && (fState != kButtonDisabled)) {
         fBgndColor = fHighColor;
      } else {
         fBgndColor = fBackground;
      }
      if (event->fType == kLeaveNotify) {
         fBgndColor = fBackground;
      }
      DoRedraw();
   }

   if ((fgDbw != event->fWindow) || (fgReleaseBtn == event->fWindow)) return kTRUE;

   if (!(event->fState & (kButton1Mask | kButton2Mask | kButton3Mask)))
      return kTRUE;

   if (fState == kButtonEngaged || fState == kButtonDisabled) return kTRUE;

   if (event->fType == kLeaveNotify) {
      fgReleaseBtn = fId;
      SetState(kButtonUp, kFALSE);
   }
   DoRedraw();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set tool tip text associated with this button. The delay is in
/// milliseconds (minimum 250). To remove tool tip call method with
/// text = 0.

void TGButton::SetToolTipText(const char *text, Long_t delayms)
{
   if (fTip) {
      delete fTip;
      fTip = 0;
   }

   if (text && strlen(text))
      fTip = new TGToolTip(fClient->GetDefaultRoot(), this, text, delayms);
}

////////////////////////////////////////////////////////////////////////////////
/// Set enabled or disabled state of button

void TGButton::SetEnabled(Bool_t e)
{
   SetState(e ? kButtonUp : kButtonDisabled);

   if (e) fWidgetFlags |= kWidgetIsEnabled;
   else   fWidgetFlags &= ~kWidgetIsEnabled;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGButton::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return graphics context for highlighted frame background.

const TGGC &TGButton::GetHibckgndGC()
{
   if (!fgHibckgndGC) {
      GCValues_t gval;
      gval.fMask = kGCForeground | kGCBackground | kGCTile |
                   kGCFillStyle  | kGCGraphicsExposures;
      gval.fForeground = gClient->GetResourcePool()->GetFrameHiliteColor();
      gval.fBackground = gClient->GetResourcePool()->GetFrameBgndColor();
      gval.fFillStyle  = kFillTiled;
      gval.fTile       = gClient->GetResourcePool()->GetCheckeredPixmap();
      gval.fGraphicsExposures = kFALSE;
      fgHibckgndGC = gClient->GetGC(&gval, kTRUE);
   }
   return *fgHibckgndGC;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a text button widget. The hotstring will be adopted and deleted
/// by the text button.

TGTextButton::TGTextButton(const TGWindow *p, TGHotString *s, Int_t id,
                           GContext_t norm, FontStruct_t font,
                           UInt_t options) : TGButton(p, id, norm, options)
{
   fLabel = s;
   fFontStruct = font;

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a text button widget.

TGTextButton::TGTextButton(const TGWindow *p, const char *s, Int_t id,
                           GContext_t norm, FontStruct_t font,
                           UInt_t options) : TGButton(p, id, norm, options)
{
   fLabel = new TGHotString(!p && !s ? GetName() : s);
   fFontStruct = font;

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a text button widget and set cmd string at same time.

TGTextButton::TGTextButton(const TGWindow *p, const char *s, const char *cmd,
                           Int_t id, GContext_t norm, FontStruct_t font,
                           UInt_t options) : TGButton(p, id, norm, options)
{
   fLabel = new TGHotString(s);
   fFontStruct = font;
   fCommand = cmd;

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialization used by the different ctors.

void TGTextButton::Init()
{
   int hotchar;

   fTMode       = kTextCenterX | kTextCenterY;
   fHKeycode    = 0;
   fHasOwnFont  = kFALSE;
   fPrevStateOn =
   fStateOn     = kFALSE;
   fWrapLength  = -1;
   fMLeft = fMRight = fMTop = fMBottom = 0;

   TGFont *font = fClient->GetFontPool()->FindFont(fFontStruct);
   if (!font) {
      font = fClient->GetFontPool()->GetFont(fgDefaultFont);
      if (font) fFontStruct = font->GetFontStruct();
   }
   if (font) {
      fTLayout = font->ComputeTextLayout(fLabel->GetString(), fLabel->GetLength(),
                                         fWrapLength, kTextLeft, 0,
                                         &fTWidth, &fTHeight);
   }
   Resize();
   fWidth = fTWidth;
   fHeight = fTHeight;

   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0) {
         const TGMainFrame *main = (TGMainFrame *) GetMainFrame();
         if (main) {
            main->BindKey(this, fHKeycode, kKeyMod1Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
         }
      }
   }
   SetWindowAttributes_t wattr;
   wattr.fMask = kWAWinGravity | kWABitGravity;
   wattr.fBitGravity = 5; // center
   wattr.fWinGravity = 1;
   gVirtualX->ChangeWindowAttributes(fId, &wattr);

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a text button widget.

TGTextButton::~TGTextButton()
{
   if (fHKeycode && (fParent->MustCleanup() != kDeepCleanup)) {
      const TGMainFrame *main = (TGMainFrame *) GetMainFrame();
      if (main) {
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
      }
   }
   if (fLabel) delete fLabel;
   if (fHasOwnFont) {
      TGGCPool *pool = fClient->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      pool->FreeGC(gc);
   }

   delete fTLayout;
}

////////////////////////////////////////////////////////////////////////////////
/// layout text button

void TGTextButton::Layout()
{
   SafeDelete(fTLayout);

   TGFont *font = fClient->GetFontPool()->FindFont(fFontStruct);
   if (!font) {
      font = fClient->GetFontPool()->GetFont(fgDefaultFont);
      if (font) fFontStruct = font->GetFontStruct();
   }
   if (font) {
      fTLayout = font->ComputeTextLayout(fLabel->GetString(), fLabel->GetLength(),
                                         fWrapLength, kTextLeft, 0,
                                         &fTWidth, &fTHeight);
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set new button text.

void TGTextButton::SetText(TGHotString *new_label)
{
   int hotchar;
   const TGMainFrame *main = (TGMainFrame *) GetMainFrame();

   if (fLabel) {
      if (main && fHKeycode) {
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
      }
      delete fLabel;
   }

   fLabel = new_label;
   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if (main && ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0)) {
         main->BindKey(this, fHKeycode, kKeyMod1Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
      }
   }

   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set new button text.

void TGTextButton::SetText(const TString &new_label)
{
   SetText(new TGHotString(new_label));
}

////////////////////////////////////////////////////////////////////////////////
/// Set text justification. Mode is an OR of the bits:
/// kTextTop, kTextBottom, kTextLeft, kTextRight, kTextCenterX and
/// kTextCenterY.

void TGTextButton::SetTextJustify(Int_t mode)
{
   fTMode = mode;

   SetWindowAttributes_t wattr;
   wattr.fMask = kWAWinGravity | kWABitGravity;
   wattr.fWinGravity = 1;

   switch (mode) {
      case kTextTop | kTextLeft:
         wattr.fBitGravity = 1; //NorthWestGravity
         break;
      case kTextTop | kTextCenterX:
      case kTextTop:
         wattr.fBitGravity = 2; //NorthGravity
         break;
      case kTextTop | kTextRight:
         wattr.fBitGravity = 3; //NorthEastGravity
         break;
      case kTextLeft | kTextCenterY:
      case kTextLeft:
         wattr.fBitGravity = 4; //WestGravity
         break;
      case kTextCenterY | kTextCenterX:
         wattr.fBitGravity = 5; //CenterGravity
         break;
      case kTextRight | kTextCenterY:
      case kTextRight:
         wattr.fBitGravity = 6; //EastGravity
         break;
      case kTextBottom | kTextLeft:
         wattr.fBitGravity = 7; //SouthWestGravity
         break;
      case kTextBottom | kTextCenterX:
      case kTextBottom:
         wattr.fBitGravity = 8; //SouthGravity
         break;
      case kTextBottom | kTextRight:
         wattr.fBitGravity = 9; //SouthEastGravity
         break;
      default:
         wattr.fBitGravity = 5; //CenterGravity
         break;
   }

   gVirtualX->ChangeWindowAttributes(fId, &wattr);
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the text button.

void TGTextButton::DoRedraw()
{
   int x, y;
   UInt_t w = GetWidth() - 1;
   UInt_t h = GetHeight()- 1;

   if ((fStyle > 0) && !(fOptions & kOwnBackground))
      gVirtualX->SetWindowBackground(fId, fBgndColor);
   TGFrame::DoRedraw();

   if (fTMode & kTextLeft) {
      x = fMLeft + 4;
   } else if (fTMode & kTextRight) {
      x = fWidth - fTWidth - fMRight - 4;
   } else {
      x = (fWidth - fTWidth + fMLeft - fMRight) >> 1;
   }

   if (fTMode & kTextTop) {
      y = fMTop + 3;
   } else if (fTMode & kTextBottom) {
      y = fHeight - fTHeight - fMBottom - 3;
   } else {
      y = (fHeight - fTHeight + fMTop - fMBottom) >> 1;
   }

   if (fState == kButtonDown || fState == kButtonEngaged) {
      ++x; ++y;
      w--; h--;
   }
   if (fStyle == 0) {
      if (fState == kButtonEngaged) {
         gVirtualX->FillRectangle(fId, GetHibckgndGC()(), 2, 2, fWidth-4, fHeight-4);
         gVirtualX->DrawLine(fId, GetHilightGC()(), 2, 2, fWidth-3, 2);
      }
   }

   Int_t hotpos = fLabel->GetHotPos();

   if (fStyle > 0) {
      gVirtualX->DrawRectangle(fId, TGFrame::GetShadowGC()(), 0, 0, w, h);
   }
   if (fState == kButtonDisabled) {
      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      if (gc) {
         Pixel_t fore = gc->GetForeground();
         Pixel_t hi = GetHilightGC().GetForeground();
         Pixel_t sh = GetShadowGC().GetForeground();

         gc->SetForeground(hi);
         fTLayout->DrawText(fId, gc->GetGC(), x + 1, y + 1, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x + 1, y + 1, hotpos - 1);

         gc->SetForeground(sh);
         fTLayout->DrawText(fId, gc->GetGC(), x, y, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x, y, hotpos - 1);
         gc->SetForeground(fore);
      }
   } else {
      fTLayout->DrawText(fId, fNormGC, x, y, 0, -1);
      if (hotpos) fTLayout->UnderlineChar(fId, fNormGC, x, y, hotpos - 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event. This function will be called when the hotkey is hit.

Bool_t TGTextButton::HandleKey(Event_t *event)
{
   if (fState == kButtonDisabled || !(event->fState & kKeyMod1Mask)) return kFALSE;
   
   Bool_t click = kFALSE;
   Bool_t was = !IsDown();   // kTRUE if button was off

   if (event->fType == kGKeyPress) {
      gVirtualX->SetKeyAutoRepeat(kFALSE);
   } else {
      gVirtualX->SetKeyAutoRepeat(kTRUE);
   }

   if (fTip && event->fType == kGKeyPress) fTip->Hide();

   // We don't need to check the key number as GrabKey will only
   // allow fHotchar events if Alt button is pressed (kKeyMod1Mask)

   if ((event->fType == kGKeyPress) && (event->fState & kKeyMod1Mask)) {
      if (fState == kButtonEngaged) return kTRUE;
      SetState(kButtonDown);
   } else if ((event->fType == kKeyRelease) && (event->fState & kKeyMod1Mask)) {
      if (fState == kButtonEngaged /*&& !allowRelease*/) return kTRUE;
      click = (fState == kButtonDown);
      if (click && fStayDown) {
         SetState(kButtonEngaged);
      } else {
         SetState(kButtonUp);
      }
   }
   if (click) {
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }
   EmitSignals(was);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// returns default size

TGDimension TGTextButton::GetDefaultSize() const
{
   UInt_t w = GetOptions() & kFixedWidth ? fWidth : fTWidth + fMLeft + fMRight + 8;
   UInt_t h = GetOptions() & kFixedHeight ? fHeight : fTHeight + fMTop + fMBottom + 7;
   return TGDimension(w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure.

FontStruct_t TGTextButton::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font.
/// If global is kTRUE font is changed globally, otherwise - locally.

void TGTextButton::SetFont(FontStruct_t font, Bool_t global)
{
   if (font != fFontStruct) {
      FontH_t v = gVirtualX->GetFontHandle(font);
      if (!v) return;

      fFontStruct = font;
      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);

      if (gc && !global) {
         gc = pool->GetGC((GCValues_t*)gc->GetAttributes(), kTRUE); // copy
         fHasOwnFont = kTRUE;
      }
      if (gc) {
         gc->SetFont(v);
         fNormGC = gc->GetGC();
      }
      Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font specified by name.
/// If global is true color is changed globally, otherwise - locally.

void TGTextButton::SetFont(const char *fontName, Bool_t global)
{
   TGFont *font = fClient->GetFont(fontName);
   if (font) {
      SetFont(font->GetFontStruct(), global);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text color.
/// If global is true color is changed globally, otherwise - locally.

void TGTextButton::SetTextColor(Pixel_t color, Bool_t global)
{
   TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
   TGGC *gc = pool->FindGC(fNormGC);

   if (gc && !global) {
      gc = pool->GetGC((GCValues_t*)gc->GetAttributes(), kTRUE); // copy
      fHasOwnFont = kTRUE;
   }
   if (gc) {
      gc->SetForeground(color);
      fNormGC = gc->GetGC();
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if text attributes are unique,
/// returns kFALSE if text attributes are shared (global).

Bool_t TGTextButton::HasOwnFont() const
{
   return fHasOwnFont;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a picture button widget. The picture is not adopted and must
/// later be freed by the user once the picture button is deleted (a single
/// picture reference might be used by other buttons).

TGPictureButton::TGPictureButton(const TGWindow *p, const TGPicture *pic,
      Int_t id, GContext_t norm, UInt_t option) : TGButton(p, id, norm, option)
{
   if (!pic) {
      Error("TGPictureButton", "pixmap not found or the file format is not supported for button %d", id);
      fPic = fClient->GetPicture("mb_question_s.xpm");
   } else {
      fPic = pic;
   }

   if (fPic) {
      fTWidth  = fPic->GetWidth();
      fTHeight = fPic->GetHeight();

      Resize(fTWidth  + (fBorderWidth << 1) + fBorderWidth + 1,
             fTHeight + (fBorderWidth << 1) + fBorderWidth); // *3
   }
   fPicD = 0;
   fOwnDisabledPic = kFALSE;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a picture button widget and set action command. The picture is
/// not adopted and must later be freed by the user once the picture button
/// is deleted (a single picture reference might be used by other
/// buttons).

TGPictureButton::TGPictureButton(const TGWindow *p, const TGPicture *pic,
      const char *cmd, Int_t id, GContext_t norm, UInt_t option)
   : TGButton(p, id, norm, option)
{
   if (!pic) {
      Error("TGPictureButton", "pixmap not found or the file format is not supported for button\n%s",
            cmd ? cmd : "");
      fPic = fClient->GetPicture("mb_question_s.xpm");
   } else {
      fPic = pic;
   }

   fCommand = cmd;

   if (fPic) {
      fTWidth  = fPic->GetWidth();
      fTHeight = fPic->GetHeight();

      Resize(fTWidth  + (fBorderWidth << 1) + fBorderWidth + 1,
             fTHeight + (fBorderWidth << 1) + fBorderWidth); // *3
   }
   fPicD = 0;
   fOwnDisabledPic = kFALSE;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a picture button. Where pic is the file name of the picture.

TGPictureButton::TGPictureButton(const TGWindow *p, const char *pic,
   Int_t id, GContext_t norm, UInt_t option) : TGButton(p, id, norm, option)
{
   if (!pic || !pic[0]) {
      if (p) Error("TGPictureButton", "pixmap not found or the file format is not supported for button");
      fPic = fClient->GetPicture("mb_question_s.xpm");
   } else {
      fPic = fClient->GetPicture(pic);
   }

   if (fPic) {
      fTWidth  = fPic->GetWidth();
      fTHeight = fPic->GetHeight();

      Resize(fTWidth  + (fBorderWidth << 1) + fBorderWidth + 1,
             fTHeight + (fBorderWidth << 1) + fBorderWidth); // *3
   }
   fPicD = 0;
   fOwnDisabledPic = kFALSE;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGPictureButton::~TGPictureButton()
{
   if (fOwnDisabledPic) fClient->FreePicture(fPicD);
}

////////////////////////////////////////////////////////////////////////////////
/// Change a picture in a picture button. The picture is not adopted and
/// must later be freed by the user once the picture button is deleted
/// (a single picture reference might be used by other buttons).

void TGPictureButton::SetPicture(const TGPicture *new_pic)
{
   if (!new_pic) {
      Error("SetPicture", "pixmap not found or the file format is not supported for button %d\n%s",
            fWidgetId, fCommand.Data());
      return;
   }

   fPic = new_pic;

   if (fState == kButtonDisabled) {
      fClient->FreePicture(fPicD);
      fPicD = 0;
   }

   fTWidth  = fPic->GetWidth();
   fTHeight = fPic->GetHeight();

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw picture button.

void TGPictureButton::DoRedraw()
{
   if (!fPic) {
      TGFrame::DoRedraw();
      return;
   }

   int x = (fWidth - fTWidth) >> 1;
   int y = (fHeight - fTHeight) >> 1;
   UInt_t w = GetWidth() - 1;
   UInt_t h = GetHeight()- 1;

   if ((fStyle > 0) && !(fOptions & kOwnBackground))
      gVirtualX->SetWindowBackground(fId, fBgndColor);
   TGFrame::DoRedraw();

   if (fState == kButtonDown || fState == kButtonEngaged) {
      ++x; ++y;
      w--; h--;
   }
   if (fStyle == 0) {
      if (fState == kButtonEngaged) {
         gVirtualX->FillRectangle(fId, GetHibckgndGC()(), 2, 2, fWidth-4, fHeight-4);
         gVirtualX->DrawLine(fId, GetHilightGC()(), 2, 2, fWidth-3, 2);
      }
   }

   const TGPicture *pic = fPic;
   if (fState == kButtonDisabled) {
      if (!fPicD) CreateDisabledPicture();
      pic = fPicD ? fPicD : fPic;
   }
   if (fStyle > 0) {
      if (fBgndColor == fHighColor) {
         gVirtualX->DrawRectangle(fId, TGFrame::GetShadowGC()(), 0, 0, w, h);
      }
   }

   pic->Draw(fId, fNormGC, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates disabled picture.

void TGPictureButton::CreateDisabledPicture()
{
   TImage *img = TImage::Create();
   if (!img) return;
   TImage *img2 = TImage::Create();
   if (!img2) {
      if (img) delete img;
      return;
   }
   TString back = gEnv->GetValue("Gui.BackgroundColor", "#c0c0c0");
   img2->FillRectangle(back.Data(), 0, 0, fPic->GetWidth(), fPic->GetHeight());
   img->SetImage(fPic->GetPicture(), fPic->GetMask());
   Pixmap_t mask = img->GetMask();
   img2->Merge(img, "overlay");

   TString name = "disbl_";
   name += fPic->GetName();
   fPicD = fClient->GetPicturePool()->GetPicture(name.Data(), img2->GetPixmap(),
                                                 mask);
   fOwnDisabledPic = kTRUE;
   delete img;
   delete img2;
}

////////////////////////////////////////////////////////////////////////////////
/// Changes disabled picture.

void TGPictureButton::SetDisabledPicture(const TGPicture *pic)
{
   if (!pic) return;

   if (fOwnDisabledPic && fPicD) fClient->FreePicture(fPicD);
   fPicD = pic;
   ((TGPicture*)pic)->AddReference();
   fOwnDisabledPic = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a check button widget. The hotstring will be adopted and deleted
/// by the check button.

TGCheckButton::TGCheckButton(const TGWindow *p, TGHotString *s, Int_t id,
                             GContext_t norm, FontStruct_t font, UInt_t option)
   : TGTextButton(p, s, id, norm, font, option)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a check button widget.

TGCheckButton::TGCheckButton(const TGWindow *p, const char *s, Int_t id,
                             GContext_t norm, FontStruct_t font, UInt_t option)
   : TGTextButton(p, s, id, norm, font, option)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a check button widget.

TGCheckButton::TGCheckButton(const TGWindow *p, const char *s, const char *cmd,
                             Int_t id, GContext_t norm, FontStruct_t font,
                             UInt_t option) : TGTextButton(p, s, cmd, id, norm, font, option)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Common check button initialization.

void TGCheckButton::Init()
{
   fPrevState =
   fState     = kButtonUp;
   fHKeycode = 0;

   fOn  = fClient->GetPicture("checked_t.xpm");
   fOff = fClient->GetPicture("unchecked_t.xpm");
   fDisOn  = fClient->GetPicture("checked_dis_t.xpm");
   fDisOff = fClient->GetPicture("unchecked_dis_t.xpm");

   Resize();

   if (!fOn) {
      Error("TGCheckButton", "checked_t.xpm not found or the file format is not supported.");
   } else if (!fOff) {
      Error("TGCheckButton", "unchecked_t.xpm not found or the file format is not supported.");
   } else if (!fDisOn) {
      Error("TGCheckButton", "checked_dis_t.xpm not found or the file format is not supported.");
   } else if (!fDisOff) {
      Error("TGCheckButton", "unchecked_dis_t.xpm not found or the file format is not supported.");
   }
   int hotchar;

   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0) {
         const TGMainFrame *main = (TGMainFrame *) GetMainFrame();
         if (main) {
            main->BindKey(this, fHKeycode, kKeyMod1Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
         }
      }
   }
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a check button.

TGCheckButton::~TGCheckButton()
{
   if (fOn)  fClient->FreePicture(fOn);
   if (fOff) fClient->FreePicture(fOff);
   if (fDisOn)  fClient->FreePicture(fDisOn);
   if (fDisOff) fClient->FreePicture(fDisOff);
}

////////////////////////////////////////////////////////////////////////////////
/// default size

TGDimension TGCheckButton::GetDefaultSize() const
{
   UInt_t w = !fTWidth ? fOff->GetWidth() : fTWidth + fOff->GetWidth() + 9;
   UInt_t h = !fTHeight ? fOff->GetHeight() : fTHeight + 2;

   w = GetOptions() & kFixedWidth ? fWidth : w;
   h = GetOptions() & kFixedHeight ? fHeight : h;

   return TGDimension(w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Set check button state.

void TGCheckButton::SetState(EButtonState state, Bool_t emit)
{
   if (state == kButtonDisabled)
      fWidgetFlags &= ~kWidgetIsEnabled;
   else
      fWidgetFlags |= kWidgetIsEnabled;
   PSetState(state, emit);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signals.

void TGCheckButton::EmitSignals(Bool_t /*wasUp*/)
{
   if (fState == kButtonUp)   Released();            // emit Released
   if (fState == kButtonDown) Pressed();             // emit Pressed
   Clicked();                                        // emit Clicked
   Toggled(fStateOn);                                // emit Toggled
}

////////////////////////////////////////////////////////////////////////////////
/// Set check button state.

void TGCheckButton::PSetState(EButtonState state, Bool_t emit)
{
   if (state != fState) {
      if (state == kButtonUp) {
         if (fPrevState == kButtonDisabled) {
            if (fStateOn) {
               fState = kButtonDown;
               fPrevState = kButtonDown;
            } else {
               fState = state;
               fPrevState = state;
            }
         } else if (fPrevState == kButtonDown) {
            fStateOn = kFALSE;
            fState = state;
            fPrevState = state;
         }
      } else if (state == kButtonDown) {
         fStateOn = kTRUE;
         fState = state;
         fPrevState = state;
      } else {
         fState = state;
         fPrevState = state;
      }
      if (emit) {
         // button signals
         EmitSignals();
      }
      DoRedraw();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the state of a check button to disabled and either on or
/// off.

void TGCheckButton::SetDisabledAndSelected(Bool_t enable)
{
   if (!enable) {
      if (fState == kButtonDisabled && fStateOn) {
         PSetState(kButtonUp, kFALSE);         // enable button
         PSetState(kButtonUp, kFALSE);         // set button up
         PSetState(kButtonDisabled, kFALSE);   // disable button
      } else {
         PSetState(kButtonUp, kFALSE);
         PSetState(kButtonDisabled, kFALSE);
      }
   } else {
      PSetState(kButtonDown, kFALSE);          // set button down
      PSetState(kButtonDisabled, kFALSE);      // disable button
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event.

Bool_t TGCheckButton::HandleButton(Event_t *event)
{
   Bool_t click = kFALSE;

   if (fTip) fTip->Hide();

   if (fState == kButtonDisabled) return kTRUE;

   Bool_t in = (event->fX >= 0) && (event->fY >= 0) &&
               (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   // We don't need to check the button number as GrabButton will
   // only allow button1 events
   if (event->fType == kButtonPress) {
      fgReleaseBtn = 0;
      if (in) {
         fOptions |= kSunkenFrame;
         Pressed();
      }
   } else { // ButtonRelease
      if (in) {
         PSetState((fPrevState == kButtonUp) ? kButtonDown : kButtonUp, kFALSE);
         click = kTRUE;
         fPrevStateOn = fStateOn;
         Released();
      }
      fgReleaseBtn = fId;
      fOptions &= ~kSunkenFrame;
   }
   if (click) {
      Clicked();
      Toggled(fStateOn);
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_CHECKBUTTON),
                  fWidgetId, (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_CHECKBUTTON),
                           fWidgetId, (Long_t) fUserData);
   }
   DoRedraw();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGCheckButton::HandleCrossing(Event_t *event)
{
   if (fTip) {
      if (event->fType == kEnterNotify)
         fTip->Reset();
      else
         fTip->Hide();
   }

   if ((fgDbw != event->fWindow) || (fgReleaseBtn == event->fWindow)) return kTRUE;

   if (!(event->fState & (kButton1Mask | kButton2Mask | kButton3Mask)))
      return kTRUE;

   if (fState == kButtonDisabled) return kTRUE;

   if (event->fType == kEnterNotify) {
      fOptions |= kSunkenFrame;
   } else {
      fOptions &= ~kSunkenFrame;
   }
   DoRedraw();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event. This function will be called when the hotkey is hit.

Bool_t TGCheckButton::HandleKey(Event_t *event)
{
   Bool_t click = kFALSE;

   if (event->fType == kGKeyPress)
      gVirtualX->SetKeyAutoRepeat(kFALSE);
   else
      gVirtualX->SetKeyAutoRepeat(kTRUE);

   if (fTip && event->fType == kGKeyPress) fTip->Hide();

   if (fState == kButtonDisabled) return kTRUE;

   // We don't need to check the key number as GrabKey will only
   // allow fHotchar events if Alt button is pressed (kKeyMod1Mask)

   if ((event->fType == kGKeyPress) && (event->fState & kKeyMod1Mask)) {
      PSetState((fPrevState == kButtonUp) ? kButtonDown : kButtonUp, kTRUE);
   } else if ((event->fType == kKeyRelease) && (event->fState & kKeyMod1Mask)) {
      click = (fState != fPrevState);
      fPrevState = fState;
   }
   if (click) {
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_CHECKBUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_CHECKBUTTON), fWidgetId,
                           (Long_t) fUserData);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the check button widget.

void TGCheckButton::DoRedraw()
{
   int x, y, y0;

   TGFrame::DoRedraw();

   x = 20;
   y = (fHeight - fTHeight) >> 1;

   y0 = !fTHeight ? 0 : y + 1;
   if (fOn && fOff) {
      Int_t smax = TMath::Max(fOn->GetHeight(), fOff->GetHeight());
      y0 = ((Int_t)fHeight <= smax) ? 0 : 1 + (((Int_t)fHeight - smax) >> 1);
   }

   if (fStateOn) {
      if (fOn) fOn->Draw(fId, fNormGC, 0, y0);
   } else {
      if (fOff) fOff->Draw(fId, fNormGC, 0, y0);
   }

   Int_t hotpos = fLabel->GetHotPos();

   if (fState == kButtonDisabled) {
      if (fStateOn == kTRUE) {
         if (fDisOn) fDisOn->Draw(fId, fNormGC, 0, y0);
      } else {
         if (fDisOff) fDisOff->Draw(fId, fNormGC, 0, y0);
      }

      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      if (gc) {
         Pixel_t fore = gc->GetForeground();
         Pixel_t hi = GetHilightGC().GetForeground();
         Pixel_t sh = GetShadowGC().GetForeground();

         gc->SetForeground(hi);
         fTLayout->DrawText(fId, gc->GetGC(), x + 1, y + 1, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x, y, hotpos - 1);

         gc->SetForeground(sh);
         fTLayout->DrawText(fId, gc->GetGC(), x, y, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x, y, hotpos - 1);

         gc->SetForeground(fore);
      }
   } else {
      fTLayout->DrawText(fId, fNormGC, x, y, 0, -1);
      if (hotpos) fTLayout->UnderlineChar(fId, fNormGC, x, y, hotpos - 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure.

FontStruct_t TGCheckButton::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGCheckButton::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a radio button widget. The hotstring will be adopted and deleted
/// by the radio button.

TGRadioButton::TGRadioButton(const TGWindow *p, TGHotString *s, Int_t id,
                             GContext_t norm, FontStruct_t font, UInt_t option)
   : TGTextButton(p, s, id, norm, font, option)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a radio button widget.

TGRadioButton::TGRadioButton(const TGWindow *p, const char *s, Int_t id,
                             GContext_t norm, FontStruct_t font, UInt_t option)
   : TGTextButton(p, s, id, norm, font, option)
{
   Init();
}
////////////////////////////////////////////////////////////////////////////////
/// Create a radio button widget.

TGRadioButton::TGRadioButton(const TGWindow *p, const char *s, const char *cmd,
                             Int_t id, GContext_t norm,
                             FontStruct_t font, UInt_t option)
    : TGTextButton(p, s, cmd, id, norm, font, option)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Common radio button initialization.

void TGRadioButton::Init()
{
   fPrevState =
   fState     = kButtonUp;
   fHKeycode  = 0;

   fOn  = fClient->GetPicture("rbutton_on.xpm");
   fOff = fClient->GetPicture("rbutton_off.xpm");
   fDisOn  = fClient->GetPicture("rbutton_dis_on.xpm");
   fDisOff = fClient->GetPicture("rbutton_dis_off.xpm");

   if (!fOn || !fOff || !fDisOn || !fDisOff)
      Error("TGRadioButton", "rbutton_*.xpm not found or the file format is not supported.");

   Resize();

   int hotchar;

   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0) {
         const TGMainFrame *main = (TGMainFrame *) GetMainFrame();
         if (main) {
            main->BindKey(this, fHKeycode, kKeyMod1Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
         }
      }
   }

   if (fParent->IsA()->InheritsFrom(TGButtonGroup::Class())) {
      ((TGButtonGroup*)fParent)->SetRadioButtonExclusive(kTRUE);
   }
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a radio button.

TGRadioButton::~TGRadioButton()
{
   if (fOn)  fClient->FreePicture(fOn);
   if (fOff) fClient->FreePicture(fOff);
   if (fDisOn)  fClient->FreePicture(fDisOn);
   if (fDisOff) fClient->FreePicture(fDisOff);
}

////////////////////////////////////////////////////////////////////////////////
/// default size

TGDimension TGRadioButton::GetDefaultSize() const
{
   UInt_t w = !fTWidth ? ( fOff?fOff->GetWidth():10) : fTWidth + fOff->GetWidth() + 10;
   UInt_t h = !fTHeight ? ( fOff?fOff->GetHeight():2) : fTHeight + 2;

   w = GetOptions() & kFixedWidth ? fWidth : w;
   h = GetOptions() & kFixedHeight ? fHeight : h;

   return TGDimension(w, h);
}
////////////////////////////////////////////////////////////////////////////////
/// Set radio button state.

void TGRadioButton::SetState(EButtonState state, Bool_t emit)
{
   if (state == kButtonDisabled)
      fWidgetFlags &= ~kWidgetIsEnabled;
   else
      fWidgetFlags |= kWidgetIsEnabled;
   PSetState(state, emit);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the state of a radio button to disabled and either on or
/// off.

void TGRadioButton::SetDisabledAndSelected(Bool_t enable)
{
   if (!enable) {
      if (fState == kButtonDisabled && fStateOn) {
         PSetState(kButtonUp, kFALSE);         // enable button
         PSetState(kButtonUp, kFALSE);         // set button up
         PSetState(kButtonDisabled, kFALSE);   // disable button
      } else {
         PSetState(kButtonUp, kFALSE);
         PSetState(kButtonDisabled, kFALSE);
      }
   } else {
      PSetState(kButtonDown, kFALSE);          // set button down
      PSetState(kButtonDisabled, kFALSE);      // disable button
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signals.

void TGRadioButton::EmitSignals(Bool_t /*wasUp*/)
{
   if (fState == kButtonUp) Released();              // emit Released
   if (fState == kButtonDown) Pressed();             // emit Pressed
   Clicked();                                        // emit Clicked
   Toggled(fStateOn);                                // emit Toggled
}

////////////////////////////////////////////////////////////////////////////////
/// Set radio button state.

void TGRadioButton::PSetState(EButtonState state, Bool_t emit)
{
   if (state != fState) {
      //      fPrevState = fState = state;
      if (state == kButtonUp) {
         if (fPrevState == kButtonDisabled) {
            if (fStateOn) {
               fState = kButtonDown;
               fPrevState = kButtonDown;
            } else {
               fState = state;
               fPrevState = state;
            }
         } else if (fPrevState == kButtonDown) {
            fStateOn = kFALSE;
            fState = state;
            fPrevState = state;
         }
      } else if (state == kButtonDown) {
         fStateOn = kTRUE;
         fState = state;
         fPrevState = state;
      } else {
         fState = state;
         fPrevState = state;
      }
      if (emit) {
         // button signals
         EmitSignals();
      }
      DoRedraw();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event.

Bool_t TGRadioButton::HandleButton(Event_t *event)
{
   Bool_t click = kFALSE;
   Bool_t toggled = kFALSE;

   if (fTip) fTip->Hide();

   if (fState == kButtonDisabled) return kFALSE;


   Bool_t in = (event->fX >= 0) && (event->fY >= 0) &&
               (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   if (event->fType == kButtonPress) { // button pressed
      fgReleaseBtn = 0;
      if (in) {
         fOptions |= kSunkenFrame;
         Pressed();
      }
   } else { // ButtonRelease
      if (in) {
         if (!fStateOn) {
            PSetState(kButtonDown, kFALSE);
            toggled = kTRUE;
         }
         fPrevStateOn = fStateOn;
         Released();
         click = kTRUE;
      }
      fOptions &= ~kSunkenFrame;
      fgReleaseBtn = fId;
   }
   if (click) {
      Clicked();
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_RADIOBUTTON),
                  fWidgetId, (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_RADIOBUTTON),
                           fWidgetId, (Long_t) fUserData);
   }
   if (toggled) {
      Toggled(fStateOn);
   }
   DoRedraw();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGRadioButton::HandleCrossing(Event_t *event)
{
   if (fTip) {
      if (event->fType == kEnterNotify)
         fTip->Reset();
      else
         fTip->Hide();
   }

   if ((fgDbw != event->fWindow) || (fgReleaseBtn == event->fWindow)) return kTRUE;

   if (!(event->fState & (kButton1Mask | kButton2Mask | kButton3Mask)))
      return kTRUE;

   if (fState == kButtonDisabled) return kTRUE;

   if (event->fType == kEnterNotify) {
      fOptions |= kSunkenFrame;
   } else {
      fOptions &= ~kSunkenFrame;
   }
   DoRedraw();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event. This function will be called when the hotkey is hit.

Bool_t TGRadioButton::HandleKey(Event_t *event)
{
   if (event->fType == kGKeyPress)
      gVirtualX->SetKeyAutoRepeat(kFALSE);
   else
      gVirtualX->SetKeyAutoRepeat(kTRUE);

   if (fTip && event->fType == kGKeyPress)
      fTip->Hide();

   if (fState == kButtonDisabled) return kTRUE;

   // We don't need to check the key number as GrabKey will only
   // allow fHotchar events if Alt button is pressed (kKeyMod1Mask)

   if ((event->fType == kGKeyPress) && (event->fState & kKeyMod1Mask)) {
      PSetState(kButtonDown, kTRUE);
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_RADIOBUTTON),
                  fWidgetId, (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_RADIOBUTTON),
                           fWidgetId, (Long_t) fUserData);
   } else if ((event->fType == kKeyRelease) && (event->fState & kKeyMod1Mask)) {
      fPrevState = fState;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a radio button.

void TGRadioButton::DoRedraw()
{
   Int_t tx, ty, y0;

   TGFrame::DoRedraw();

   tx = 20;
   ty = (fHeight - fTHeight) >> 1;

//   pw = 12;
   y0 = !fTHeight ? 0 : ty + 1;
   if (fOn && fOff) {
      Int_t smax = TMath::Max(fOn->GetHeight(), fOff->GetHeight());
      y0 = ((Int_t)fHeight <= smax) ? 0 : 1 + (((Int_t)fHeight - smax) >> 1);
   }

   if (fStateOn) {
      if (fOn) fOn->Draw(fId, fNormGC, 0, y0);
   } else {
      if (fOff) fOff->Draw(fId, fNormGC, 0, y0);
   }

   Int_t hotpos = fLabel->GetHotPos();

   if (fState == kButtonDisabled) {
      if (fStateOn == kTRUE) {
         if (fDisOn) fDisOn->Draw(fId, fNormGC, 0, y0);
      } else {
         if (fDisOff) fDisOff->Draw(fId, fNormGC, 0, y0);
      }

      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      if (gc) {
         Pixel_t fore = gc->GetForeground();
         Pixel_t hi = GetHilightGC().GetForeground();
         Pixel_t sh = GetShadowGC().GetForeground();

         gc->SetForeground(hi);
         fTLayout->DrawText(fId, gc->GetGC(), tx + 1, ty + 1, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), tx, ty, hotpos - 1);

         gc->SetForeground(sh);
         fTLayout->DrawText(fId, gc->GetGC(), tx, ty, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), tx, ty, hotpos - 1);

         gc->SetForeground(fore);
      }
   } else {
      fTLayout->DrawText(fId, fNormGC, tx, ty, 0, -1);
      if (hotpos) fTLayout->UnderlineChar(fId, fNormGC, tx, ty, hotpos-1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure.

FontStruct_t TGRadioButton::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGRadioButton::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a button widget as a C++ statement(s) on output stream out.

void TGButton::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (fState == kButtonDown) {
      out << "   " << GetName() << "->SetState(kButtonDown);"  << std::endl;
   }
   if (fState == kButtonDisabled) {
      out << "   " << GetName() << "->SetState(kButtonDisabled);"  << std::endl;
   }
   if (fState == kButtonEngaged) {
      out << "   " << GetName() << "->SetState(kButtonEngaged);"  << std::endl;
   }
   if (fBackground != fgDefaultFrameBackground) {
      SaveUserColor(out, option);
      out << "   " << GetName() << "->ChangeBackground(ucolor);" << std::endl;
   }

   if (fTip) {
      TString tiptext = fTip->GetText()->GetString();
      tiptext.ReplaceAll("\n", "\\n");
      out << "   ";
      out << GetName() << "->SetToolTipText(" << quote
          << tiptext << quote << ");"  << std::endl;
   }
   if (strlen(fCommand)) {
      out << "   " << GetName() << "->SetCommand(" << quote << fCommand
          << quote << ");" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a text button widget as a C++ statement(s) on output stream out.

void TGTextButton::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   TString outext(fLabel->GetString());
   if (fLabel->GetHotPos() > 0)
      outext.Insert(fLabel->GetHotPos()-1, "&");
   if (outext.First('\n') >= 0)
      outext.ReplaceAll("\n", "\\n");

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

   out << "   TGTextButton *";
   out << GetName() << " = new TGTextButton(" << fParent->GetName()
       << "," << quote << outext.Data() << quote;

   if (GetOptions() == (kRaisedFrame | kDoubleBorder)) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            if (fWidgetId == -1) {
               out << ");" << std::endl;
            } else {
               out << "," << fWidgetId <<");" << std::endl;
            }
         } else {
            out << "," << fWidgetId << "," << parGC << ");" << std::endl;
         }
      } else {
         out << "," << fWidgetId << "," << parGC << "," << parFont << ");" << std::endl;
      }
   } else {
      out << "," << fWidgetId << "," << parGC << "," << parFont << "," << GetOptionString() << ");" << std::endl;
   }

   out << "   " << GetName() << "->SetTextJustify(" << fTMode << ");" << std::endl;
   out << "   " << GetName() << "->SetMargins(" << fMLeft << "," << fMRight << ",";
   out << fMTop << "," << fMBottom << ");" << std::endl;
   out << "   " << GetName() << "->SetWrapLength(" << fWrapLength << ");" << std::endl;

   out << "   " << GetName() << "->Resize(" << GetWidth() << "," << GetHeight()
       << ");" << std::endl;

   TGButton::SavePrimitive(out,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a picture button widget as a C++ statement(s) on output stream out.

void TGPictureButton::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (!fPic) {
      Error("SavePrimitive()", "pixmap not found or the file format is not supported for picture button %d ", fWidgetId);
      return;
   }

   // GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC;
   parGC.Form("%s::GetDefaultGC()()",IsA()->GetName());

   if (GetDefaultGC()() != fNormGC) {
      TGGC *userGC = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);
      if (userGC) {
         userGC->SavePrimitive(out, option);
         parGC.Form("uGC->GetGC()");
      }
   }

   char quote = '"';
   TString picname = gSystem->UnixPathName(fPic->GetName());
   gSystem->ExpandPathName(picname);

   out <<"   TGPictureButton *";

   out << GetName() << " = new TGPictureButton(" << fParent->GetName()
       << ",gClient->GetPicture(" << quote
       << picname << quote << ")";

   if (GetOptions() == (kRaisedFrame | kDoubleBorder)) {
      if (fNormGC == GetDefaultGC()()) {
         if (fWidgetId == -1) {
            out << ");" << std::endl;
         } else {
            out << "," << fWidgetId << ");" << std::endl;
         }
      } else {
         out << "," << fWidgetId << "," << parGC.Data() << ");" << std::endl;
      }
   } else {
      out << "," << fWidgetId << "," << parGC.Data() << "," << GetOptionString()
          << ");" << std::endl;
   }

   TGButton::SavePrimitive(out,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a check button widget as a C++ statement(s) on output stream out.

void TGCheckButton::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   TString outext(fLabel->GetString());
   if (fLabel->GetHotPos() > 0)
      outext.Insert(fLabel->GetHotPos()-1, "&");
   if (outext.First('\n') >= 0)
      outext.ReplaceAll("\n", "\\n");

   out <<"   TGCheckButton *";
   out << GetName() << " = new TGCheckButton(" << fParent->GetName()
       << "," << quote << outext.Data() << quote;

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

   if (GetOptions() == kChildFrame) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            if (fWidgetId == -1) {
               out << ");" << std::endl;
            } else {
               out << "," << fWidgetId << ");" << std::endl;
            }
         } else {
            out << "," << fWidgetId << "," << parGC << ");" << std::endl;
         }
      } else {
         out << "," << fWidgetId << "," << parGC << "," << parFont << ");" << std::endl;
      }
   } else {
      out << "," << fWidgetId << "," << parGC << "," << parFont << "," << GetOptionString() << ");" << std::endl;
   }

   TGButton::SavePrimitive(out,option);
   if (fState == kButtonDisabled) {
      if (IsDisabledAndSelected())
         out << "   " << GetName() << "->SetDisabledAndSelected(kTRUE);" << std::endl;
      else
         out << "   " << GetName() << "->SetDisabledAndSelected(kFALSE);" << std::endl;
   }
   out << "   " << GetName() << "->SetTextJustify(" << fTMode << ");" << std::endl;
   out << "   " << GetName() << "->SetMargins(" << fMLeft << "," << fMRight << ",";
   out << fMTop << "," << fMBottom << ");" << std::endl;
   out << "   " << GetName() << "->SetWrapLength(" << fWrapLength << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a radio button widget as a C++ statement(s) on output stream out.

void TGRadioButton::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   TString outext(fLabel->GetString());
   if (fLabel->GetHotPos() > 0)
      outext.Insert(fLabel->GetHotPos()-1, "&");
   if (outext.First('\n') >= 0)
      outext.ReplaceAll("\n", "\\n");

   out << "   TGRadioButton *";
   out << GetName() << " = new TGRadioButton(" << fParent->GetName()
       << "," << quote << outext.Data() << quote;

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

   if (GetOptions() == (kChildFrame)) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            if (fWidgetId == -1) {
               out <<");" << std::endl;
            } else {
               out << "," << fWidgetId << ");" << std::endl;
            }
         } else {
            out << "," << fWidgetId << "," << parGC << ");" << std::endl;
         }
      } else {
         out << "," << fWidgetId << "," << parGC << "," << parFont << ");" << std::endl;
      }
   } else {
      out << "," << fWidgetId << "," << parGC << "," << parFont << "," << GetOptionString() << ");" << std::endl;
   }

   TGButton::SavePrimitive(out,option);
   if (fState == kButtonDisabled) {
      if (IsDisabledAndSelected())
         out << "   " << GetName() << "->SetDisabledAndSelected(kTRUE);" << std::endl;
      else
         out << "   " << GetName() << "->SetDisabledAndSelected(kFALSE);" << std::endl;
   }
   out << "   " << GetName() << "->SetTextJustify(" << fTMode << ");" << std::endl;
   out << "   " << GetName() << "->SetMargins(" << fMLeft << "," << fMRight << ",";
   out << fMTop << "," << fMBottom << ");" << std::endl;
   out << "   " << GetName() << "->SetWrapLength(" << fWrapLength << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a menu button widget. The hotstring will be adopted and
/// deleted by the menu button. This constructor creates a
/// menubutton with a popup menu attached that appears when the
/// button for it is clicked. The popup menu is adopted.

TGSplitButton::TGSplitButton(const TGWindow *p, TGHotString* menulabel,
                           TGPopupMenu *popmenu, Bool_t split, Int_t id,
                           GContext_t norm, FontStruct_t fontstruct, UInt_t options)
                           : TGTextButton(p, menulabel, id, norm, fontstruct, options)
{
   fFontStruct = fontstruct;
   fMBWidth = 16;
   fMenuLabel = new TGHotString(*menulabel);
   fPopMenu = popmenu;
   fPopMenu->fSplitButton = this;
   fSplit = split;
   fTMode = 0;
   fHKeycode = 0;
   fMBState = kButtonUp; fDefaultCursor = fClient->GetResourcePool()->GetGrabCursor();
   fKeyNavigate = kFALSE;
   fWidestLabel = "";
   fHeighestLabel = "";

   // Find and set the correct size for the menu and the button.
   TGMenuEntry *entry = 0;
   TGHotString lstring(*fMenuLabel);
   TGHotString hstring(*fMenuLabel);
   const TList *list = fPopMenu->GetListOfEntries();
   UInt_t lwidth = 0, lheight = 0;
   UInt_t twidth = 0, theight = 0;

   TGFont *font = fClient->GetFontPool()->FindFont(fFontStruct);
   if (!font) {
      font = fClient->GetFontPool()->GetFont(fgDefaultFont);
      if (font) fFontStruct = font->GetFontStruct();
   }

   if (font) font->ComputeTextLayout(lstring, lstring.GetLength(),
                                     fWrapLength, kTextLeft, 0,
                                     &lwidth, &lheight);

   TIter iter(list);
   entry = (TGMenuEntry *)iter.Next();
   while (entry != 0) {
      if (entry->GetType() == kMenuEntry) {
         const TGHotString string(*(entry->GetLabel()));
         if (font) font->ComputeTextLayout(string, string.GetLength(),
                                           fWrapLength, kTextLeft, 0,
                                           &twidth, &theight);
         if(twidth > lwidth) {
            lstring = string;
         }
         if(theight > lheight) {
            hstring = string;
         }
      }
      entry = (TGMenuEntry *)iter.Next();
   }
   fWidestLabel = lstring;
   fHeighestLabel =  hstring;

   if (font) {
      UInt_t dummy = 0;
      font->ComputeTextLayout(lstring, lstring.GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &fTWidth, &dummy);
      font->ComputeTextLayout(hstring, hstring.GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &dummy, &fTHeight);
   }
   fTBWidth = fTWidth + 8;
   fHeight = fTHeight + 7;
   Resize(fTBWidth + fMBWidth, fHeight);

   ChangeOptions(GetOptions() | kFixedSize);

   // Save the id of the 1st item on the menu.
   TIter iter1(list);
   do {
      entry = (TGMenuEntry *)iter1.Next();
      if ((entry) && (entry->GetStatus() & kMenuEnableMask) &&
          !(entry->GetStatus() & kMenuHideMask) &&
          (entry->GetType() != kMenuSeparator) &&
          (entry->GetType() != kMenuLabel)) break;
      entry = (TGMenuEntry *)iter1.Next();
   } while (entry);
   if (entry) fEntryId = entry->GetEntryId();

   // An additional connection that is needed.
   fPopMenu->Connect("Activated(Int_t)", "TGSplitButton", this, "HandleMenu(Int_t)");
   SetSplit(fSplit);

   Init();
}


////////////////////////////////////////////////////////////////////////////////
/// Common initialization used by the different ctors.

void TGSplitButton::Init()
{
   Int_t hotchar;

   fTMode       = kTextCenterX | kTextCenterY;
   fHKeycode    = 0;
   fHasOwnFont  = kFALSE;
   fPrevStateOn =
   fStateOn     = kFALSE;
   fMBState     = kButtonUp;

   SetSize(TGDimension(fWidth, fHeight));

   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0) {
         const TGMainFrame *main = (TGMainFrame *) GetMainFrame();
         if (main) {
            main->BindKey(this, fHKeycode, kKeyMod1Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
            main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
         }
      }
   }
   SetWindowAttributes_t wattr;
   wattr.fMask = kWAWinGravity | kWABitGravity;
   wattr.fBitGravity = 5; // center
   wattr.fWinGravity = 1;
   gVirtualX->ChangeWindowAttributes(fId, &wattr);

   // Make sure motion is detected too.
   AddInput(kPointerMotionMask | kEnterWindowMask | kLeaveWindowMask);

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a split button widget.

TGSplitButton::~TGSplitButton()
{
   if (fPopMenu) delete fPopMenu;
   if (fMenuLabel) delete fMenuLabel;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw triangle (arrow) on which user can click to open Popup.

void TGSplitButton::DrawTriangle(const GContext_t gc, Int_t x, Int_t y)
{
   Point_t points[3];

   points[0].fX = x;
   points[0].fY = y;
   points[1].fX = x + 5;
   points[1].fY = y;
   points[2].fX = x + 2;
   points[2].fY = y + 3;

   gVirtualX->FillPolygon(fId, gc, points, 3);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the size of the button.

void TGSplitButton::CalcSize()
{
   Int_t max_ascent, max_descent;
   fTWidth = gVirtualX->TextWidth(fFontStruct, fLabel->GetString(), fLabel->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;

   fTBWidth = fTWidth + 8;
   fHeight = fTHeight + 7;
   fWidth = fTBWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in case the button is split.

Bool_t TGSplitButton::HandleSButton(Event_t *event)
{
   if (fState == kButtonDisabled) return kFALSE;

   Bool_t activate = kFALSE;
   Bool_t bclick = kFALSE;
   static Bool_t mbpress = kFALSE;
   static Bool_t tbpress = kFALSE;
   static Bool_t outpress = kFALSE;

   Bool_t inTB = (event->fX >= 0) && (event->fY >= 0) &&
                 (event->fX <= (Int_t)fTBWidth) && (event->fY <= (Int_t)fHeight);

   Bool_t inMB = (event->fX >= (Int_t)(fWidth -fMBWidth)) && (event->fY >= 0) &&
      (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   // We don't need to check the button number as GrabButton will
   // only allow button1 events

   if (inTB) {
      if (event->fType == kButtonPress) {
         mbpress = kFALSE;
         tbpress = kTRUE;
         fgReleaseBtn = 0;
         if (fState == kButtonEngaged) {
            return kTRUE;
         }
         SetState(kButtonDown);
         Pressed();
      } else { // ButtonRelease
         if (fMBState == kButtonDown) {
            SetMBState(kButtonUp);
         }
         if (fState == kButtonEngaged && tbpress) {
            SetState(kButtonUp);
            Released();
            bclick = kTRUE;
         } else if (fState == kButtonDown && tbpress) {
            tbpress = kFALSE;
            if (fStayDown) {
               bclick = kTRUE;
               SetState(kButtonEngaged);
               fgReleaseBtn = 0;
            } else {
               bclick = kTRUE;
               SetState(kButtonUp);
               Released();
               fgReleaseBtn = fId;
            }
         }else {
            SetState(kButtonUp);
         }
      }
   } else if (inMB) {
      if (event->fType == kButtonPress) {
         fgReleaseBtn = 0;
         mbpress = kTRUE;
         tbpress = kFALSE;
         if (fMBState == kButtonEngaged) {
            return kTRUE;
         }
         SetMBState(kButtonDown);
         MBPressed();
         gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                kPointerMotionMask, kNone, fDefaultCursor);
      } else { // ButtonRelease
         if (fState == kButtonDown) {
            SetState(kButtonUp);
         }
         if (fMBState == kButtonEngaged && mbpress) {
            mbpress = kFALSE;
            SetMBState(kButtonUp);
            SetMenuState(kFALSE);
            MBReleased();
            MBClicked();
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
         } else if (fMBState == kButtonDown && mbpress) {
            MBClicked();
            SetMBState(kButtonEngaged);
            SetMenuState(kTRUE);
            fgReleaseBtn = 0;
         } else {
            SetMBState(kButtonUp);
         }
      }
   } else {
      if (event->fType == kButtonPress) {
         fgReleaseBtn = 0;
         outpress = kTRUE;
      } else { // ButtonRelease
         if(outpress) {
            outpress = kFALSE;
            SetMBState(kButtonUp);
            SetMenuState(kFALSE);
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
            activate = kTRUE;
         }
      }
   }
   if (bclick) {
      Clicked();
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }
   if (activate) {
      TGMenuEntry *entry =  fPopMenu->GetCurrent();
      if (entry) {
         if ((entry->GetStatus() & kMenuEnableMask) &&
             !(entry->GetStatus() & kMenuHideMask) &&
             (entry->GetType() != kMenuSeparator) &&
             (entry->GetType() != kMenuLabel)) {
            Int_t id = entry->GetEntryId();
            fPopMenu->Activated(id);
         }
      }
   }
   //   if (mbclick) {
   //      MBClicked();
   //      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
   //             (Long_t) fUserData);
   //      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
   //                      (Long_t) fUserData);
   // }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event in case of split menu.

Bool_t TGSplitButton::HandleSCrossing(Event_t *event)
{
   if (fTip) {
      if (event->fType == kEnterNotify)
         fTip->Reset();
      else
         fTip->Hide();
   }

   if ((fgDbw != event->fWindow) || (fgReleaseBtn == event->fWindow)) return kTRUE;

   if (!(event->fState & (kButton1Mask | kButton2Mask | kButton3Mask)))
      return kTRUE;

   if (fState == kButtonEngaged || fState == kButtonDisabled) return kTRUE;

   Bool_t inTB = (event->fX <= (Int_t)fTBWidth);

   //   Bool_t inMB = (event->fX >= (Int_t)(fWidth -fMBWidth)) && (event->fY >= 0) &&
   //      (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   if (event->fType == kEnterNotify) {
      if (inTB) {
         SetState(kButtonDown, kFALSE);
      } else {
         if(fMBState == kButtonEngaged)  return kTRUE;
         SetMBState(kButtonDown);
      }
   } else {
      // kLeaveNotify
      if(fState == kButtonDown) {
         SetState(kButtonUp, kFALSE);
      }
      if (fMBState == kButtonEngaged) return kTRUE;
      SetMBState(kButtonUp);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event. This function will be called when the hotkey is hit.

Bool_t TGSplitButton::HandleSKey(Event_t *event)
{
   if (fState == kButtonDisabled) return kFALSE;

   Bool_t click = kFALSE;

   if (event->fType == kGKeyPress) {
      gVirtualX->SetKeyAutoRepeat(kFALSE);
   } else {
      gVirtualX->SetKeyAutoRepeat(kTRUE);
   }

   if (fTip && event->fType == kGKeyPress) fTip->Hide();

   // We don't need to check the key number as GrabKey will only
   // allow fHotchar events if Alt button is pressed (kKeyMod1Mask)

   if ((event->fType == kGKeyPress) && (event->fState & kKeyMod1Mask)) {
      if (fState == kButtonEngaged) return kTRUE;
      SetState(kButtonDown);
      Pressed();
   } else if ((event->fType == kKeyRelease) && (event->fState & kKeyMod1Mask)) {
      if (fState == kButtonEngaged) {
         SetState(kButtonUp);
         Released();
      }
      if (fStayDown) {
         SetState(kButtonEngaged);
      } else {
         SetState(kButtonUp);
         Released();
      }
      click = kTRUE;
   }
   if (click) {
      Clicked();
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Popup the attached menu.

void TGSplitButton::SetMenuState(Bool_t state)
{
   if (state) {
      Int_t    ax, ay;
      Window_t wdummy;

      if (fSplit) {
         Int_t n_entries = 0;
         TGMenuEntry *entry = 0;
         TIter next(fPopMenu->GetListOfEntries());

         while ((entry = (TGMenuEntry *) next())) {
            if ((entry->GetType() != kMenuSeparator) &&
                (entry->GetType() != kMenuLabel)) {
               n_entries++;
            }
         }
         if (n_entries <= 1) {
            Info("TGSplitButton", "Only one entry in the menu.");
            return;
         }
      }

      gVirtualX->TranslateCoordinates(fId, fPopMenu->GetParent()->GetId(),
                                      0, 0, ax, ay, wdummy);

      // place the menu just under the window:
      fPopMenu->PlaceMenu(ax-1, ay+fHeight, kTRUE, kFALSE); //kTRUE);
      BindKeys(kTRUE);
      BindMenuKeys(kTRUE);
   } else {
      fPopMenu->EndMenu(fUserData);
      BindKeys(kFALSE);
      BindMenuKeys(kFALSE);
      fPopMenu->EndMenu(fUserData);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the text button.

void TGSplitButton::DoRedraw()
{
   int x, y;
   TGFrame::DoRedraw();

   if (fState == kButtonDisabled) fMBState = kButtonDisabled;
   else if (fMBState == kButtonDisabled) fMBState = kButtonUp;

   if (fTMode & kTextLeft) {
      x = fMLeft + 4;
   } else if (fTMode & kTextRight) {
      x = fWidth - fTWidth -fMBWidth - fMRight - 4;
   } else {
      x = (fWidth - fTWidth -fMBWidth + fMLeft - fMRight) >> 1;
   }

   if (fTMode & kTextTop) {
      y = fMTop + 3;
   } else if (fTMode & kTextBottom) {
      y = fHeight - fTHeight - fMBottom - 3;
   } else {
      y = (fHeight - fTHeight + fMTop - fMBottom) >> 1;
   }

   if (fState == kButtonDown || fState == kButtonEngaged) { ++x; ++y; }
   if (fState == kButtonEngaged) {
      gVirtualX->FillRectangle(fId, GetHibckgndGC()(), 2, 2, fWidth-4, fHeight-4);
      gVirtualX->DrawLine(fId, GetHilightGC()(), 2, 2, fWidth-3, 2);
   }

   Int_t hotpos = fLabel->GetHotPos();

   if (fState == kButtonDisabled) {
      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      if (gc) {
         Pixel_t fore = gc->GetForeground();
         Pixel_t hi = GetHilightGC().GetForeground();
         Pixel_t sh = GetShadowGC().GetForeground();

         gc->SetForeground(hi);
         fTLayout->DrawText(fId, gc->GetGC(), x + 1, y + 1, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x + 1, y + 1, hotpos - 1);

         gc->SetForeground(sh);
         fTLayout->DrawText(fId, gc->GetGC(), x, y, 0, -1);
         if (hotpos) fTLayout->UnderlineChar(fId, gc->GetGC(), x, y, hotpos - 1);
         gc->SetForeground(fore);
      }
   } else {
      fTLayout->DrawText(fId, fNormGC, x, y, 0, -1);
      if (hotpos) fTLayout->UnderlineChar(fId, fNormGC, x, y, hotpos - 1);
   }

   // Draw the parts of the button needed when a menu is attached.

   // triangle position
   x = fWidth - 11;
   y = fHeight - 10;

   if (fSplit) {
      // separator position
      Int_t lx = fWidth - fMBWidth;
      Int_t ly = 2;
      Int_t lh = fHeight - 2;

      if(fMBState == kButtonDown || fMBState == kButtonEngaged) {
         x++;
         y++;
      }

      gVirtualX->DrawLine(fId, GetShadowGC()(),  lx, ly + 2, lx, lh - 4);
      gVirtualX->DrawLine(fId, GetHilightGC()(), lx + 1, ly + 2, lx + 1, lh - 3);
      gVirtualX->DrawLine(fId, GetHilightGC()(), lx, lh - 3, lx + 1, lh - 3);

      if (fMBState == kButtonEngaged) {
         gVirtualX->FillRectangle(fId, GetHibckgndGC()(), fTBWidth + 1, 1, fMBWidth - 3, fHeight - 3);
      }

      if (fMBState == kButtonDisabled) {
         DrawTriangle(GetHilightGC()(), x + 1, y + 1);
         DrawTriangle(GetShadowGC()(), x, y);
      } else {
         DrawTriangle(fNormGC, x, y);
      }

   } else {
      x -= 2;
      if(fState == kButtonDown || fState == kButtonEngaged) {
         x++;
         y++;
      }
      if (fState == kButtonDisabled) {
         DrawTriangle(GetHilightGC()(), x + 1, y + 1);
         DrawTriangle(GetShadowGC()(), x, y);
      } else {
         DrawTriangle(fNormGC, x, y);
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// If on kTRUE bind arrow, popup menu hot keys, otherwise
/// remove key bindings.

void TGSplitButton::BindKeys(Bool_t on)
{
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Up), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Down), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Enter), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Return), kAnyModifier, on);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Escape), kAnyModifier, on);
}

////////////////////////////////////////////////////////////////////////////////
/// If on kTRUE bind Menu hot keys, otherwise remove key bindings.

void TGSplitButton::BindMenuKeys(Bool_t on)
{
   TGMenuEntry *e = 0;
   TIter next(fPopMenu->GetListOfEntries());

   while ((e = (TGMenuEntry*)next())) {
      Int_t hot = 0;
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
/// returns default size

TGDimension TGSplitButton::GetDefaultSize() const
{
   UInt_t w = GetOptions() & kFixedWidth ? fWidth + fMBWidth : fTWidth + fMLeft + fMRight + fMBWidth + 8;
   UInt_t h = GetOptions() & kFixedHeight ? fHeight : fTHeight + fMTop + fMBottom + 7;
   return TGDimension(w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Set new button text.

void TGSplitButton::SetText(TGHotString *new_label)
{
   Int_t hotchar;
   static Bool_t longlabeltip = kFALSE;
   const TGMainFrame *main = (TGMainFrame *) GetMainFrame();

   TGFont *font = fClient->GetFontPool()->FindFont(fFontStruct);
   if (!font) {
      font = fClient->GetFontPool()->GetFont(fgDefaultFont);
      if (font) fFontStruct = font->GetFontStruct();
   }

   UInt_t width = 0, bwidth = 0, dummy;
   if (font) {
      font->ComputeTextLayout(new_label->GetString(), new_label->GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &width, &dummy);
      font->ComputeTextLayout(fWidestLabel.GetString(), fWidestLabel.GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &bwidth, &dummy);
   }
   if (width > bwidth) {
      if (!fTip) {
         SetToolTipText(new_label->GetString());
         longlabeltip = kTRUE;
      }
      Info("TGSplitbutton", "Length of new label to long, label truncated.");
      new_label->Resize(fWidestLabel.GetLength());
   } else if (new_label->GetLength() <= fWidestLabel.GetLength() && longlabeltip) {
      if (fTip) delete fTip;
      fTip = 0;
      longlabeltip = kFALSE;
   }

   if (fLabel) {
      if (main && fHKeycode) {
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
         main->RemoveBind(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
      }
      delete fLabel;
   }

   fLabel = new_label;
   if ((hotchar = fLabel->GetHotChar()) != 0) {
      if (main && ((fHKeycode = gVirtualX->KeysymToKeycode(hotchar)) != 0)) {
         main->BindKey(this, fHKeycode, kKeyMod1Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyLockMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);

         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
         main->BindKey(this, fHKeycode, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
      }
   }

   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set new button text.

void TGSplitButton::SetText(const TString &new_label)
{
   SetText(new TGHotString(new_label));
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font.
/// If global is kTRUE font is changed globally, otherwise - locally.

void TGSplitButton::SetFont(FontStruct_t font, Bool_t global)
{
   if (font != fFontStruct) {
      FontH_t v = gVirtualX->GetFontHandle(font);
      if (!v) return;

      fFontStruct = font;
      TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);

      if ((gc) && !global) {
         gc = pool->GetGC((GCValues_t*)gc->GetAttributes(), kTRUE); // copy
         fHasOwnFont = kTRUE;
      }
      if (gc) {
         gc->SetFont(v);
         fNormGC = gc->GetGC();
      }
      fClient->NeedRedraw(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font specified by name.
/// If global is true color is changed globally, otherwise - locally.

void TGSplitButton::SetFont(const char *fontName, Bool_t global)
{
   TGFont *font = fClient->GetFont(fontName);
   if (font) {
      SetFont(font->GetFontStruct(), global);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the state of the Menu Button part

void TGSplitButton::SetMBState(EButtonState state)
{
   if (state != fMBState) {
      fMBState = state;
      DoRedraw();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the split status of a button.

void TGSplitButton::SetSplit(Bool_t split)
{
   if(split) {
      fStayDown = kFALSE;
      Disconnect(fPopMenu, "PoppedDown()");
      fPopMenu->Connect("PoppedDown()", "TGSplitButton", this, "SetMBState(=kButtonUp)");
      fPopMenu->Connect("PoppedDown()", "TGSplitButton", this, "MBReleased()");

      TGMenuEntry *entry = fPopMenu->GetEntry(fEntryId);
      if (entry) {
         TGHotString *tmp = new TGHotString(*(entry->GetLabel()));
         SetText(tmp);

         TString str("ItemClicked(=");
         str += entry->GetEntryId();
         str += ")";
         Connect("Clicked()", "TGSplitButton", this, str);
         fEntryId = entry->GetEntryId();
         fPopMenu->HideEntry(fEntryId);
      }
   } else {
      fStayDown = kTRUE;
      Disconnect(fPopMenu, "PoppedDown()");
      Disconnect(this, "Clicked()", this);
      fPopMenu->Connect("PoppedDown()", "TGSplitButton", this, "SetState(=kButtonUp)");
      fPopMenu->Connect("PoppedDown()", "TGSplitButton", this, "Released()");
      fPopMenu->EnableEntry(fEntryId);
      TGHotString *tmp = new TGHotString(*fMenuLabel);
      SetText(tmp);
   }

   fSplit = split;
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button events.

Bool_t TGSplitButton::HandleButton(Event_t *event)
{
   if (fState == kButtonDisabled) return kFALSE;

   if (fSplit) return HandleSButton(event);

   Bool_t in = (event->fX >= 0) && (event->fY >= 0) &&
               (event->fX <= (Int_t)fWidth) && (event->fY <= (Int_t)fHeight);

   Bool_t activate = kFALSE;
   Bool_t click = kFALSE;

   if (in) {
      if (event->fType == kButtonPress) {
         fgReleaseBtn = 0;
         if (fState == kButtonEngaged) {
            return kTRUE;
         }
         SetState(kButtonDown);
         Pressed();
         gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                kPointerMotionMask, kNone, fDefaultCursor);
      } else { // ButtonRelease
         if (fState == kButtonEngaged) {
            SetState(kButtonUp);
            SetMenuState(kFALSE);
            Released();
            click = kTRUE;
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
         } else {
            click = (fState == kButtonDown);
            if (click && fStayDown) {
               SetState(kButtonEngaged);
               SetMenuState(kTRUE);
               fgReleaseBtn = 0;
            } else {
               SetState(kButtonUp);
               Released();
               fgReleaseBtn = fId;
            }
         }
         fKeyNavigate = kFALSE;
      }
   } else {
      if (event->fType == kButtonPress) {
         fgReleaseBtn = 0;
      } else { // ButtonRelease
         SetState(kButtonUp);
         SetMenuState(kFALSE);
         gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
         activate = kTRUE;
      }
   }
   if (click) {
      Clicked();
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }
   if (activate && !fKeyNavigate) {
      TGMenuEntry *entry =  fPopMenu->GetCurrent();
      if (entry) {
         if ((entry->GetStatus() & kMenuEnableMask) &&
             !(entry->GetStatus() & kMenuHideMask) &&
             (entry->GetType() != kMenuSeparator) &&
             (entry->GetType() != kMenuLabel)) {
            Int_t id = entry->GetEntryId();
            fPopMenu->Activated(id);
         }
      }
   }

   return kTRUE;

}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGSplitButton::HandleCrossing(Event_t *event)
{
   if (fSplit) {
      return HandleSCrossing(event);
   } else {
      return TGButton::HandleCrossing(event);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event. This function will be called when the hotkey is hit.

Bool_t TGSplitButton::HandleKey(Event_t *event)
{
   Bool_t click = kFALSE;

   if (fState == kButtonDisabled) return kTRUE;

   if(fSplit) return HandleSKey(event);

   if (event->fType == kGKeyPress) {
      gVirtualX->SetKeyAutoRepeat(kFALSE);
   } else {
      gVirtualX->SetKeyAutoRepeat(kTRUE);
   }

   if (fTip && event->fType == kGKeyPress) fTip->Hide();

   // We don't need to check the key number as GrabKey will only
   // allow fHotchar events if Alt button is pressed (kKeyMod1Mask)
   if (event->fState & kKeyMod1Mask) {
      RequestFocus();
      fKeyNavigate = kTRUE;
      if (event->fType == kGKeyPress) {
         if (fState == kButtonEngaged) return kTRUE;
         SetState(kButtonDown);
         Pressed();
      } else if (event->fType == kKeyRelease) {
         click = kTRUE;
         if (fState == kButtonEngaged) {
            SetState(kButtonUp);
            SetMenuState(kFALSE);
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
         } else if (fState == kButtonDown && fStayDown) {
            SetState(kButtonEngaged);
            SetMenuState(kTRUE);
            gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                   kPointerMotionMask, kNone, fDefaultCursor);
            TGMenuEntry *entry = 0;
            TIter next(fPopMenu->GetListOfEntries());

            while ((entry = (TGMenuEntry *) next())) {
               if ((entry->GetStatus() & kMenuEnableMask) &&
                   !(entry->GetStatus() & kMenuHideMask) &&
                   (entry->GetType() != kMenuSeparator) &&
                   (entry->GetType() != kMenuLabel)) break;
            }
            if (entry) {
               fPopMenu->Activate(entry);
            }
         } else {
            Released();
            SetState(kButtonUp);
         }
      }
   } else {
      fKeyNavigate = kTRUE;
      if (event->fType == kGKeyPress) {
         Event_t ev;
         ev.fX = ev.fY = 1;
         UInt_t keysym;
         char tmp[2];

         gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);

         TGMenuEntry *ce = 0;
         TIter next(fPopMenu->GetListOfEntries());

         while ((ce = (TGMenuEntry*)next())) {
            UInt_t hot = 0;
            if (ce->GetLabel()) hot = ce->GetLabel()->GetHotChar();
            if (!hot || (hot != keysym)) continue;

            fPopMenu->Activate(ce);
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
            SetMenuState(kFALSE);
            ev.fType = kButtonRelease;
            ev.fWindow = fPopMenu->GetId();
            fKeyNavigate = kFALSE;
            return HandleButton(&ev);
         }

         ce = fPopMenu->GetCurrent();

         switch ((EKeySym)keysym) {
         case kKey_Up:
            if (ce) ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->Before(ce);
            while (ce && ((ce->GetType() == kMenuSeparator) ||
                          (ce->GetType() == kMenuLabel) ||
                          !(ce->GetStatus() & kMenuEnableMask))) {
               ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->Before(ce);
            }
            if (!ce) ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->Last();
            break;
         case kKey_Down:
            if (ce) ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->After(ce);
            while (ce && ((ce->GetType() == kMenuSeparator) ||
                          (ce->GetType() == kMenuLabel) ||
                          !(ce->GetStatus() & kMenuEnableMask))) {
               ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->After(ce);
            }
            if (!ce) ce = (TGMenuEntry*)fPopMenu->GetListOfEntries()->First();
            break;
         case kKey_Enter:
         case kKey_Return:
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
            SetMenuState(kFALSE);
            ev.fType = kButtonRelease;
            ev.fWindow = fPopMenu->GetId();
            fKeyNavigate = kFALSE;
            HandleButton(&ev);
            break;
         case kKey_Escape:
            gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
            SetMenuState(kFALSE);
            break;
         default:
            break;
         }
         if (ce) fPopMenu->Activate(ce);
      }
   }
   if (click) {
      Clicked();
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                  (Long_t) fUserData);
      fClient->ProcessLine(fCommand, MK_MSG(kC_COMMAND, kCM_BUTTON), fWidgetId,
                           (Long_t) fUserData);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a motion event in a TGSplitButton.

Bool_t TGSplitButton::HandleMotion(Event_t *event)
{
   if (fKeyNavigate) return kTRUE;

   if (fSplit) {
      if (fMBState == kButtonDown) {
         if (event->fX < (Int_t)fTBWidth) {
            SetMBState(kButtonUp);
            SetState(kButtonDown);
         }
      } else if (fState == kButtonDown) {
         if (event->fX > (Int_t)fTBWidth) {
            SetState(kButtonUp);
            SetMBState(kButtonDown);
         }

      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// layout text button

void TGSplitButton::Layout()
{
   UInt_t dummya = 0, dummyb = 0;
   SafeDelete(fTLayout);

   TGFont *font = fClient->GetFontPool()->FindFont(fFontStruct);
   if (!font) {
      font = fClient->GetFontPool()->GetFont(fgDefaultFont);
      if (font) fFontStruct = font->GetFontStruct();
   }
   if (font) {
      fTLayout = font->ComputeTextLayout(fLabel->GetString(),
                                         fLabel->GetLength(),
                                         fWrapLength, kTextLeft, 0,
                                         &dummya, &dummyb);
      UInt_t dummy = 0;
      font->ComputeTextLayout(fWidestLabel.GetString(), fWidestLabel.GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &fTWidth, &dummy);
      font->ComputeTextLayout(fHeighestLabel.GetString(), fHeighestLabel.GetLength(),
                              fWrapLength, kTextLeft, 0,
                              &dummy, &fTHeight);
   }
   fTBWidth = fTWidth + 8;
   fWidth = fTBWidth + fMBWidth;
   fHeight = fTHeight + 7;
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a menu item activation.

void TGSplitButton::HandleMenu(Int_t id)
{
   SetMenuState(kFALSE);

   if (fSplit) {
      SetMBState(kButtonUp);
      Disconnect(this, "Clicked()", this);
      // connect clicked to the ItemClicked signal with the correct id
      Connect("Clicked()", "TGSplitButton", this,
              TString::Format("ItemClicked(=%d)", id));

      // reenable hidden entries
      const TList *list = fPopMenu->GetListOfEntries();
      TIter iter(list);
      fPopMenu->EnableEntry(fEntryId);
      TGMenuEntry *entry = fPopMenu->GetEntry(id);
      if (entry) {
         TGHotString *label = entry->GetLabel();
         TGHotString *tmp = new TGHotString(*label);
         SetText(tmp);
      }
      fPopMenu->HideEntry(id);
      if (entry) fEntryId = entry->GetEntryId();
   } else {
      SetState(kButtonUp);
      ItemClicked(id);
   }
   DoRedraw();
}
