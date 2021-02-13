// @(#)root/gui:$Id$
// Author: Fons Rademakers   03/01/98

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
// TGFrame, TGCompositeFrame, TGVerticalFrame, TGHorizontalFrame,       //
// TGMainFrame, TGTransientFrame and TGGroupFrame                       //
//                                                                      //
// The frame classes describe the different "dressed" GUI windows.      //
//                                                                      //
// The TGFrame class is a subclasses of TGWindow, and is used as base   //
// class for some simple widgets (buttons, labels, etc.).               //
// It provides:                                                         //
//  - position & dimension fields                                       //
//  - an 'options' attribute (see constant above)                       //
//  - a generic event handler                                           //
//  - a generic layout mechanism                                        //
//  - a generic border                                                  //
//                                                                      //
// The TGCompositeFrame class is the base class for composite widgets   //
// (menu bars, list boxes, etc.).                                       //
// It provides:                                                         //
//  - a layout manager                                                  //
//  - a frame container (TList *)                                       //
//                                                                      //
// The TGVerticalFrame and TGHorizontalFrame are composite frame that   //
// layout their children in vertical or horizontal way.                 //
//                                                                      //
// The TGMainFrame class defines top level windows that interact with   //
// the system Window Manager.                                           //
//                                                                      //
// The TGTransientFrame class defines transient windows that typically  //
// are used for dialogs windows.                                        //
//                                                                      //
// The TGGroupFrame is a composite frame with a border and a title.     //
// It is typically used to group a number of logically related widgets  //
// visually together.                                                   //
//                                                                      //
//Begin_Html
/*
<img src="gif/tgcompositeframe_classtree.gif">
*/
//End_Html
//////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "TGFrame.h"
#include "TGResourcePool.h"
#include "TGPicture.h"
#include "TList.h"
#include "TApplication.h"
#include "TTimer.h"
#include "TClass.h"

#include "TObjString.h"
#include "TBits.h"
#include "TColor.h"
#include "TROOT.h"
#include "TDatime.h"
#include "KeySymbols.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TSystem.h"
#include "TVirtualDragManager.h"
#include "TGuiBuilder.h"
#include "TQConnection.h"
#include "TGButton.h"
#include "TGSplitter.h"
#include "TGDNDManager.h"
#include "TImage.h"
#include "TObjectSpy.h"
#include "TVirtualX.h"

#include <iostream>
#include <fstream>


Bool_t      TGFrame::fgInit = kFALSE;
Pixel_t     TGFrame::fgDefaultFrameBackground = 0;
Pixel_t     TGFrame::fgDefaultSelectedBackground = 0;
Pixel_t     TGFrame::fgWhitePixel = 0;
Pixel_t     TGFrame::fgBlackPixel = 0;
const TGGC *TGFrame::fgBlackGC = nullptr;
const TGGC *TGFrame::fgWhiteGC = nullptr;
const TGGC *TGFrame::fgHilightGC = nullptr;
const TGGC *TGFrame::fgShadowGC = nullptr;
const TGGC *TGFrame::fgBckgndGC = nullptr;
Time_t      TGFrame::fgLastClick = 0;
UInt_t      TGFrame::fgLastButton = 0;
Int_t       TGFrame::fgDbx = 0;
Int_t       TGFrame::fgDby = 0;
Window_t    TGFrame::fgDbw = 0;
UInt_t      TGFrame::fgUserColor = 0;

const TGFont *TGGroupFrame::fgDefaultFont = nullptr;
const TGGC   *TGGroupFrame::fgDefaultGC = nullptr;

TGLayoutHints *TGCompositeFrame::fgDefaultHints = nullptr;

static const char *gSaveMacroTypes[] = {
   "ROOT macros", "*.C",
   "GIF",         "*.gif",
   "PNG",         "*.png",
   "JPEG",        "*.jpg",
   "TIFF",        "*.tiff",
   "XPM",         "*.xpm",
   "All files",   "*",
   0,             0
};

TList *gListOfHiddenFrames = new TList();

ClassImp(TGFrame);
ClassImp(TGCompositeFrame);
ClassImp(TGVerticalFrame);
ClassImp(TGHorizontalFrame);
ClassImp(TGMainFrame);
ClassImp(TGTransientFrame);
ClassImp(TGGroupFrame);
ClassImp(TGHeaderFrame);


////////////////////////////////////////////////////////////////////////////////
/// Create a TGFrame object. Options is an OR of the EFrameTypes.

TGFrame::TGFrame(const TGWindow *p, UInt_t w, UInt_t h,
                 UInt_t options, Pixel_t back)
   : TGWindow(p, 0, 0, w, h, 0, 0, 0, 0, 0, options)
{
   if (!fgInit && gClient) {
      TGFrame::GetDefaultFrameBackground();
      TGFrame::GetDefaultSelectedBackground();
      TGFrame::GetWhitePixel();
      TGFrame::GetBlackPixel();
      TGFrame::GetBlackGC();
      TGFrame::GetWhiteGC();
      TGFrame::GetHilightGC();
      TGFrame::GetShadowGC();
      TGFrame::GetBckgndGC();
      fgInit = kTRUE;
   }

   SetWindowAttributes_t wattr;

   fDNDState   = 0;
   fBackground = back;
   fOptions    = options;
   fWidth = w; fHeight = h; fX = fY = fBorderWidth = 0;
   fMinWidth    = 0;
   fMinHeight   = 0;
   fMaxWidth    = kMaxUInt;
   fMaxHeight   = kMaxUInt;
   fFE          = 0;

   if (fOptions & (kSunkenFrame | kRaisedFrame))
      fBorderWidth = (gClient->GetStyle() > 1) ? 1 : (fOptions & kDoubleBorder) ? 2 : 1;

   wattr.fMask = kWABackPixel | kWAEventMask;
   wattr.fBackgroundPixel = back;
   wattr.fEventMask = kExposureMask;
   if (fOptions & kMainFrame) {
      wattr.fEventMask |= kStructureNotifyMask;
      gVirtualX->ChangeWindowAttributes(fId, &wattr);
      //if (fgDefaultBackgroundPicture)
      //   SetBackgroundPixmap(fgDefaultBackgroundPicture->GetPicture());
   } else {
      gVirtualX->ChangeWindowAttributes(fId, &wattr);
      //if (!(fOptions & kOwnBackground))
      //   SetBackgroundPixmap(kParentRelative);
   }
   fEventMask = (UInt_t) wattr.fEventMask;

   if ((fOptions & kDoubleBorder) && (gClient->GetStyle() > 1))
      ChangeOptions(fOptions ^ kDoubleBorder);

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a frame using an externally created window. For example
/// to register the root window (called by TGClient), or a window
/// created via TVirtualX::InitWindow() (id is obtained with
/// TVirtualX::GetWindowID()).

TGFrame::TGFrame(TGClient *c, Window_t id, const TGWindow *parent)
   : TGWindow(c, id, parent)
{
   if (!fgInit && gClient) {
      TGFrame::GetDefaultFrameBackground();
      TGFrame::GetDefaultSelectedBackground();
      TGFrame::GetWhitePixel();
      TGFrame::GetBlackPixel();
      TGFrame::GetBlackGC();
      TGFrame::GetWhiteGC();
      TGFrame::GetHilightGC();
      TGFrame::GetShadowGC();
      TGFrame::GetBckgndGC();
      fgInit = kTRUE;
   }

   WindowAttributes_t attributes;
   // Initialize some values - needed for batch mode!
   attributes.fX      = 0;
   attributes.fY      = 0;
   attributes.fWidth  = 100;
   attributes.fHeight = 100;
   attributes.fBorderWidth = 4;
   attributes.fYourEventMask = 0;
   gVirtualX->GetWindowAttributes(id, attributes);

   fDNDState    = 0;
   fX           = attributes.fX;
   fY           = attributes.fY;
   fWidth       = attributes.fWidth;
   fHeight      = attributes.fHeight;
   fBorderWidth = attributes.fBorderWidth;
   fEventMask   = (UInt_t) attributes.fYourEventMask;
   fBackground  = 0;
   fOptions     = 0;
   fMinWidth    = 0;
   fMinHeight   = 0;
   fMaxWidth    = kMaxUInt;
   fMaxHeight   = kMaxUInt;
   fFE          = 0;

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGFrame::~TGFrame()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Delete window. Use single shot timer to call final delete method.
/// We use this indirect way since deleting the window in its own
/// execution "thread" can cause side effects because frame methods
/// can still be called while the window object has already been deleted.

void TGFrame::DeleteWindow()
{
   if (gDNDManager) {
      if (gDNDManager->GetMainFrame() == this)
         gDNDManager->SetMainFrame(0);
   }
   if (!TestBit(kDeleteWindowCalled)) {
      // coverity[returned_null]
      // coverity[dereference]
      TTimer::SingleShot(150, IsA()->GetName(), this, "ReallyDelete()");
   }
   SetBit(kDeleteWindowCalled);
}

////////////////////////////////////////////////////////////////////////////////
/// Change frame background color.

void TGFrame::ChangeBackground(Pixel_t back)
{
   fBackground = back;
   gVirtualX->SetWindowBackground(fId, back);
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return frame foreground color.

Pixel_t TGFrame::GetForeground() const
{
   return fgBlackPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Set background color (override from TGWindow base class).
/// Same effect as ChangeBackground().

void TGFrame::SetBackgroundColor(Pixel_t back)
{
   fBackground = back;
   TGWindow::SetBackgroundColor(back);
}

////////////////////////////////////////////////////////////////////////////////
/// Change frame options. Options is an OR of the EFrameTypes.

void TGFrame::ChangeOptions(UInt_t options)
{
   if ((options & (kDoubleBorder | kSunkenFrame | kRaisedFrame)) !=
      (fOptions & (kDoubleBorder | kSunkenFrame | kRaisedFrame))) {
      if (!InheritsFrom(TGGroupFrame::Class())) {
         if (options & (kSunkenFrame | kRaisedFrame))
            fBorderWidth = (gClient->GetStyle() > 1) ? 1 : (fOptions & kDoubleBorder) ? 2 : 1;
         else
            fBorderWidth = 0;
      }
   }

   fOptions = options;
}

////////////////////////////////////////////////////////////////////////////////
/// Add events specified in the emask to the events the frame should handle.

void TGFrame::AddInput(UInt_t emask)
{
   fEventMask |= emask;
   gVirtualX->SelectInput(fId, fEventMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove events specified in emask from the events the frame should handle.

void TGFrame::RemoveInput(UInt_t emask)
{
   fEventMask &= ~emask;
   gVirtualX->SelectInput(fId, fEventMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw 3D rectangle on the frame border.

void TGFrame::Draw3dRectangle(UInt_t type, Int_t x, Int_t y,
                              UInt_t w, UInt_t h)
{
   switch (type) {
      case kSunkenFrame:
         gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y,     x+w-2, y);
         gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y,     x,     y+h-2);
         gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y+h-1, x+w-1, y+h-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), x+w-1, y+h-1, x+w-1, y);
         break;

      case kSunkenFrame | kDoubleBorder:
         if (gClient->GetStyle() < 2) {
            gVirtualX->DrawLine(fId, GetShadowGC()(), x,     y,     x+w-2, y);
            gVirtualX->DrawLine(fId, GetShadowGC()(), x,     y,     x,     y+h-2);
            gVirtualX->DrawLine(fId, GetBlackGC()(),  x+1,   y+1,   x+w-3, y+1);
            gVirtualX->DrawLine(fId, GetBlackGC()(),  x+1,   y+1,   x+1,   y+h-3);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y+h-1, x+w-1, y+h-1);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x+w-1, y+h-1, x+w-1, y);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  x+1,   y+h-2, x+w-2, y+h-2);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  x+w-2, y+1,   x+w-2, y+h-2);
         }
         else {
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y,     x+w-2, y);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y,     x,     y+h-2);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y+h-1, x+w-1, y+h-1);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x+w-1, y+h-1, x+w-1, y);
         }
         break;

      case kRaisedFrame:
         gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x+w-2, y);
         gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x,     y+h-2);
         gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y+h-1, x+w-1, y+h-1);
         gVirtualX->DrawLine(fId, GetShadowGC()(),  x+w-1, y+h-1, x+w-1, y);
         break;

      case kRaisedFrame | kDoubleBorder:
         if (gClient->GetStyle() < 2) {
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x+w-2, y);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x,     y+h-2);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  x+1,   y+1,   x+w-3, y+1);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  x+1,   y+1,   x+1,   y+h-3);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x+1,   y+h-2, x+w-2, y+h-2);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x+w-2, y+h-2, x+w-2, y+1);
            gVirtualX->DrawLine(fId, GetBlackGC()(),   x,     y+h-1, x+w-1, y+h-1);
            gVirtualX->DrawLine(fId, GetBlackGC()(),   x+w-1, y+h-1, x+w-1, y);
         }
         else {
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x+w-2, y);
            gVirtualX->DrawLine(fId, GetHilightGC()(), x,     y,     x,     y+h-2);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x,     y+h-1, x+w-1, y+h-1);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  x+w-1, y+h-1, x+w-1, y);
         }
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw frame border.

void TGFrame::DrawBorder()
{
   Draw3dRectangle(fOptions & (kSunkenFrame | kRaisedFrame | kDoubleBorder),
                   0, 0, fWidth, fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the frame.

void TGFrame::DoRedraw()
{
   gVirtualX->ClearArea(fId, fBorderWidth, fBorderWidth,
                   fWidth - (fBorderWidth << 1), fHeight - (fBorderWidth << 1));

   // border will only be drawn if we have a 3D option hint
   // (kRaisedFrame or kSunkenFrame)
   DrawBorder();
}

////////////////////////////////////////////////////////////////////////////////
/// This event is generated when the frame is resized.

Bool_t TGFrame::HandleConfigureNotify(Event_t *event)
{
   if ((event->fWidth != fWidth) || (event->fHeight != fHeight)) {
      fWidth  = event->fWidth;
      fHeight = event->fHeight;
      Layout();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle all frame events. Events are dispatched to the specific
/// event handlers.

Bool_t TGFrame::HandleEvent(Event_t *event)
{
   if (gDragManager && !fClient->IsEditDisabled() &&
       gDragManager->HandleEvent(event)) return kTRUE;

   TObjectSpy deleteCheck(this);

   switch (event->fType) {

      case kExpose:
         HandleExpose(event);
         break;

      case kConfigureNotify:
         while (gVirtualX->CheckEvent(fId, kConfigureNotify, *event))
            ;
         // protection
         if ((event->fWidth < 32768) && (event->fHeight  < 32768)){
            ProcessedConfigure(event);  // emit signal
            HandleConfigureNotify(event);
         }
         break;

      case kGKeyPress:
      case kKeyRelease:
         HandleKey(event);
         break;

      case kFocusIn:
      case kFocusOut:
         HandleFocusChange(event);
         break;

      case kButtonPress:
         {
            Int_t dbl_clk = kFALSE;

            if ((event->fTime - fgLastClick < 350) &&
                (event->fCode == fgLastButton) &&
                (TMath::Abs(event->fXRoot - fgDbx) < 6) &&
                (TMath::Abs(event->fYRoot - fgDby) < 6) &&
                (event->fWindow == fgDbw))
               dbl_clk = kTRUE;

            fgLastClick  = event->fTime;
            fgLastButton = event->fCode;
            fgDbx = event->fXRoot;
            fgDby = event->fYRoot;
            fgDbw = event->fWindow;

            if (dbl_clk) {
               if ((event->fState & kKeyControlMask) &&
                    !GetEditDisabled() && gGuiBuilder) {
                  StartGuiBuilding(!IsEditable());
                  return kTRUE;
               }

               if (!HandleDoubleClick(event)) {
                  HandleButton(event);
               }
            } else {
               HandleButton(event);
            }
         }
         break;

      case kButtonDoubleClick:
         {
            fgLastClick  = event->fTime;
            fgLastButton = event->fCode;
            fgDbx = event->fXRoot;
            fgDby = event->fYRoot;
            fgDbw = event->fWindow;

            HandleDoubleClick(event);
         }
         break;

      case kButtonRelease:
         HandleButton(event);
         break;

      case kEnterNotify:
      case kLeaveNotify:
         HandleCrossing(event);
         break;

      case kMotionNotify:
         while (gVirtualX->CheckEvent(fId, kMotionNotify, *event))
            ;
         HandleMotion(event);
         break;

      case kClientMessage:
         HandleClientMessage(event);
         break;

      case kSelectionNotify:
         HandleSelection(event);
         break;

      case kSelectionRequest:
         HandleSelectionRequest(event);
         break;

      case kSelectionClear:
         HandleSelectionClear(event);
         break;

      case kColormapNotify:
         HandleColormapChange(event);
         break;

      default:
         //Warning("HandleEvent", "unknown event (%#x) for (%#x)", event->fType, fId);
         break;
   }

   if (deleteCheck.GetObject())
      ProcessedEvent(event);  // emit signal

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///   std::cout << fWidth << "x" << fHeight << std::endl;

TGDimension TGFrame::GetDefaultSize() const
{
   return TGDimension(fWidth, fHeight);
}


////////////////////////////////////////////////////////////////////////////////
/// Move frame.

void TGFrame::Move(Int_t x, Int_t y)
{
   if (x != fX || y != fY) {
      TGWindow::Move(x, y);
      fX = x; fY = y;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the frame.
/// If w=0 && h=0 - Resize to default size

void TGFrame::Resize(UInt_t w, UInt_t h)
{
   if (w != fWidth || h != fHeight) {
      TGDimension siz(0,0);
      siz = GetDefaultSize();
      fWidth = w ? w : siz.fWidth;
      fHeight = h ? h : siz.fHeight;
      TGWindow::Resize(fWidth, fHeight);
      Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the frame.

void TGFrame::Resize(TGDimension size)
{
   Resize(size.fWidth, size.fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Move and/or resize the frame.
/// If w=0 && h=0 - Resize to default size

void TGFrame::MoveResize(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   // we do it anyway as we don't know if it's only a move or only a resize
   TGDimension siz(0,0);
   siz = GetDefaultSize();
   fWidth = w ? w : siz.fWidth;
   fHeight = h ? h : siz.fHeight;
   fX = x; fY = y;
   TGWindow::MoveResize(x, y, fWidth, fHeight);
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Send message (i.e. event) to window w. Message is encoded in one long
/// as message type and up to two long parameters.

void TGFrame::SendMessage(const TGWindow *w, Long_t msg, Long_t parm1, Long_t parm2)
{
   Event_t event;

   if (w) {
      event.fType   = kClientMessage;
      event.fFormat = 32;
      event.fHandle = gROOT_MESSAGE;

      event.fWindow  = w->GetId();
      event.fUser[0] = msg;
      event.fUser[1] = parm1;
      event.fUser[2] = parm2;
      event.fUser[3] = 0;
      event.fUser[4] = 0;

      gVirtualX->SendEvent(w->GetId(), &event);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a client message. Client messages are the ones sent via
/// TGFrame::SendMessage (typically by widgets).

Bool_t TGFrame::HandleClientMessage(Event_t *event)
{
   if (gDNDManager) {
      gDNDManager->HandleClientMessage(event);
   }
   if (event->fHandle == gROOT_MESSAGE) {
      ProcessMessage(event->fUser[0], event->fUser[1], event->fUser[2]);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get default frame background.

ULong_t TGFrame::GetDefaultFrameBackground()
{
   static Bool_t init = kFALSE;
   if (!init && gClient) {
      fgDefaultFrameBackground = gClient->GetResourcePool()->GetFrameBgndColor();
      init = kTRUE;
   }
   return fgDefaultFrameBackground;
}

////////////////////////////////////////////////////////////////////////////////
/// Get default selected frame background.

ULong_t TGFrame::GetDefaultSelectedBackground()
{
   static Bool_t init = kFALSE;
   if (!init && gClient) {
      fgDefaultSelectedBackground = gClient->GetResourcePool()->GetSelectedBgndColor();
      init = kTRUE;
   }
   return fgDefaultSelectedBackground;
}

////////////////////////////////////////////////////////////////////////////////
/// Get white pixel value.

ULong_t TGFrame::GetWhitePixel()
{
   static Bool_t init = kFALSE;
   if (!init && gClient) {
      fgWhitePixel = gClient->GetResourcePool()->GetWhiteColor();
      init  = kTRUE;
   }
   return fgWhitePixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Get black pixel value.

ULong_t TGFrame::GetBlackPixel()
{
   static Bool_t init = kFALSE;
   if (!init && gClient) {
      fgBlackPixel = gClient->GetResourcePool()->GetBlackColor();
      init = kTRUE;
   }
   return fgBlackPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Get black graphics context.

const TGGC &TGFrame::GetBlackGC()
{
   if (!fgBlackGC && gClient)
      fgBlackGC = gClient->GetResourcePool()->GetBlackGC();
   return *fgBlackGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Get white graphics context.

const TGGC &TGFrame::GetWhiteGC()
{
   if (!fgWhiteGC && gClient)
      fgWhiteGC = gClient->GetResourcePool()->GetWhiteGC();
   return *fgWhiteGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Get highlight color graphics context.

const TGGC &TGFrame::GetHilightGC()
{
   if (!fgHilightGC && gClient)
      fgHilightGC = gClient->GetResourcePool()->GetFrameHiliteGC();
   return *fgHilightGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Get shadow color graphics context.

const TGGC &TGFrame::GetShadowGC()
{
   if (!fgShadowGC && gClient)
      fgShadowGC = gClient->GetResourcePool()->GetFrameShadowGC();
   return *fgShadowGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Get background color graphics context.

const TGGC &TGFrame::GetBckgndGC()
{
   if (!fgBckgndGC && gClient)
      fgBckgndGC = gClient->GetResourcePool()->GetFrameBckgndGC();
   return *fgBckgndGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Get time of last mouse click.

Time_t TGFrame::GetLastClick()
{
   return fgLastClick;
}

////////////////////////////////////////////////////////////////////////////////
/// Print window id.

void TGFrame::Print(Option_t *option) const
{
   TString opt = option;
   if (opt.Contains("tree")) {
      TGWindow::Print(option);
      return;
   }

   std::cout <<  option << ClassName() << ":\tid=" << fId << " parent=" << fParent->GetId();
   std::cout << " x=" << fX << " y=" << fY;
   std::cout << " w=" << fWidth << " h=" << fHeight << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// SetDragType

void TGFrame::SetDragType(Int_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// SetDropType

void TGFrame::SetDropType(Int_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns drag source type.
/// If frame is not "draggable" - return zero

Int_t TGFrame::GetDragType() const
{
   return fClient->IsEditable();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns drop target type.
/// If frame cannot accept drop - return zero

Int_t TGFrame::GetDropType() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Go into GUI building mode.

void TGFrame::StartGuiBuilding(Bool_t on)
{
   if (GetEditDisabled()) return;
   if (!gDragManager) gDragManager = TVirtualDragManager::Instance();
   if (!gDragManager) return;

   TGCompositeFrame *comp = 0;

   if (InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame *)this;
   } else if (fParent->InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame*)fParent;
   }
   if (comp) comp->SetEditable(on);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a composite frame. A composite frame has in addition to a TGFrame
/// also a layout manager and a list of child frames.

TGCompositeFrame::TGCompositeFrame(const TGWindow *p, UInt_t w, UInt_t h,
         UInt_t options, Pixel_t back) : TGFrame(p, w, h, options, back)
{
   fLayoutManager = 0;
   fList          = new TList;
   fLayoutBroken  = kFALSE;
   fMustCleanup   = kNoCleanup;
   fMapSubwindows = fParent->IsMapSubwindows();
   if (!fgDefaultHints)
      fgDefaultHints = new TGLayoutHints;

   if (fOptions & kHorizontalFrame)
      SetLayoutManager(new TGHorizontalLayout(this));
   else
      SetLayoutManager(new TGVerticalLayout(this));

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a frame using an externally created window. For example
/// to register the root window (called by TGClient), or a window
/// created via TVirtualX::InitWindow() (id is obtained with TVirtualX::GetWindowID()).

TGCompositeFrame::TGCompositeFrame(TGClient *c, Window_t id, const TGWindow *parent)
   : TGFrame(c, id, parent)
{
   fLayoutManager = 0;
   fList          = new TList;
   fLayoutBroken  = kFALSE;
   fMustCleanup   = kNoCleanup;
   fMapSubwindows = fParent->IsMapSubwindows();
   if (!fgDefaultHints)
      fgDefaultHints = new TGLayoutHints;

   SetLayoutManager(new TGVerticalLayout(this));

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a composite frame.

TGCompositeFrame::~TGCompositeFrame()
{
   if (fMustCleanup != kNoCleanup) {
      Cleanup();
   } else {
      TGFrameElement *el = 0;
      TIter next(fList);

      while ((el = (TGFrameElement *) next())) {
         fList->Remove(el);
         delete el;
      }
   }

   delete fList;
   delete fLayoutManager;
   fList = 0;
   fLayoutManager = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if frame is being edited.

Bool_t TGCompositeFrame::IsEditable() const
{
   return (fClient->GetRoot() == (TGWindow*)this);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch ON/OFF edit mode.
/// If edit mode is ON it is possible:
///
///  1. embed other ROOT GUI application (a la ActiveX)
///
///  For example:
///    TGMainFrame *m = new TGMainFrame(gClient->GetRoot(), 500, 500);
///    m->SetEditable();
///    gSystem->Load("$ROOTSYS/test/Aclock"); // load Aclock demo
///    Aclock a;
///    gROOT->Macro("$ROOTSYS/tutorials/gui/guitest.C");
///    m->SetEditable(0);
///    m->MapWindow();
///

void TGCompositeFrame::SetEditable(Bool_t on)
{
   if (on && ((fEditDisabled & kEditDisable) ||
              (fEditDisabled & kEditDisableLayout))) return;

   if (on) {
      fClient->SetRoot(this);
   } else {
      fClient->SetRoot(0);
   }
   if (gDragManager) gDragManager->SetEditable(on);
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup and delete all objects contained in this composite frame.
/// This will delete all objects added via AddFrame().
/// CAUTION: all objects (frames and layout hints) must be unique, i.e.
/// cannot be shared.

void TGCompositeFrame::Cleanup()
{
   if (!fList) return;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      if (el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (!gVirtualX->InheritsFrom("TGX11") && !gVirtualX->InheritsFrom("TGCocoa"))
            el->fFrame->DestroyWindow();
         delete el->fFrame;
      }

      if (el->fLayout && (el->fLayout != fgDefaultHints) &&
          (el->fLayout->References() > 0)) {
         el->fLayout->RemoveReference();
         if (!el->fLayout->References()) {
            el->fLayout->fFE = 0;
            delete el->fLayout;
         }
      }
      fList->Remove(el);
      delete el;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the layout manager for the composite frame.
/// The layout manager is adopted by the frame and will be deleted
/// by the frame.

void TGCompositeFrame::SetLayoutManager(TGLayoutManager *l)
{
   if (l) {
      delete fLayoutManager;
      fLayoutManager = l;
   } else
      Error("SetLayoutManager", "no layout manager specified");
}

////////////////////////////////////////////////////////////////////////////////
/// Set broken layout. No Layout method is called.

void TGCompositeFrame::SetLayoutBroken(Bool_t on)
{
   fLayoutBroken = on;
}

////////////////////////////////////////////////////////////////////////////////
/// Set edit disable flag for this frame and subframes
///
///  - if (on & kEditDisable) - disable edit for this frame and all subframes.

void TGCompositeFrame::SetEditDisabled(UInt_t on)
{
   fEditDisabled = on;
   UInt_t set = on & kEditDisable;

   // propagate only kEditDisable
   if (set == kEditDisable) {

      TGFrameElement *el;
      TIter next(fList);
      while ((el = (TGFrameElement *) next())) {
         if (el->fFrame) {
            el->fFrame->SetEditDisabled(set);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change composite frame options. Options is an OR of the EFrameTypes.

void TGCompositeFrame::ChangeOptions(UInt_t options)
{
   TGFrame::ChangeOptions(options);

   if (options & kHorizontalFrame)
      SetLayoutManager(new TGHorizontalLayout(this));
   else if (options & kVerticalFrame)
      SetLayoutManager(new TGVerticalLayout(this));
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on automatic cleanup of child frames in dtor.
///
/// if mode = kNoCleanup    - no automatic cleanup
/// if mode = kLocalCleanup - automatic cleanup in this composite frame only
/// if mode = kDeepCleanup  - automatic deep cleanup in this composite frame
///                           and all child composite frames (hierarchical)
///
/// Attention!
///    Hierarchical cleaning is dangerous and must be used with caution.
///    There are many GUI components (in ROOT and in user code) which do not
///    use Clean method in destructor ("custom deallocation").
///    Adding such component to GUI container which is using hierarchical
///    cleaning will produce seg. violation when container is deleted.
///    The reason is double deletion: first whem Clean method is invoked,
///    then at "custom deallocation".
///    We are going to correct all ROOT code to make it to be
///    consitent with hierarchical cleaning scheeme.

void TGCompositeFrame::SetCleanup(Int_t mode)
{
   if (mode == fMustCleanup)
      return;

   fMustCleanup = mode;

   if (fMustCleanup == kDeepCleanup) {
      TGFrameElement *el;
      TIter next(fList);

      while ((el = (TGFrameElement *) next())) {
         if (el->fFrame->InheritsFrom(TGCompositeFrame::Class())) {
            el->fFrame->SetCleanup(kDeepCleanup);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find frame-element holding frame f.

TGFrameElement* TGCompositeFrame::FindFrameElement(TGFrame *f) const
{
   if (!fList) return 0;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next()))
      if (el->fFrame == f)
         return el;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add frame to the composite frame using the specified layout hints.
/// If no hints are specified default hints TGLayoutHints(kLHintsNormal,0,0,0,0)
/// will be used. Most of the time, however, you will want to provide
/// specific hints. User specified hints can be reused many times
/// and need to be destroyed by the user. The added frames cannot not be
/// added to different composite frames but still need to be deleted by
/// the user.

void TGCompositeFrame::AddFrame(TGFrame *f, TGLayoutHints *l)
{
   TGFrameElement *nw = new TGFrameElement(f, l ? l : fgDefaultHints);
   fList->Add(nw);

   // in case of recusive cleanup, propagate cleanup setting to all
   // child composite frames
   if (fMustCleanup == kDeepCleanup)
      f->SetCleanup(kDeepCleanup);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all frames from composite frame.

void TGCompositeFrame::RemoveAll()
{
   if (!fList) return;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      fList->Remove(el);
      if (el->fLayout) el->fLayout->RemoveReference();
//       el->fFrame->SetFrameElement(0);
      delete el;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove frame from composite frame.

void TGCompositeFrame::RemoveFrame(TGFrame *f)
{
   TGFrameElement *el = FindFrameElement(f);

   if (el) {
      fList->Remove(el);
      if (el->fLayout) el->fLayout->RemoveReference();
      f->SetFrameElement(0);
      delete el;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map all sub windows that are part of the composite frame.

void TGCompositeFrame::MapSubwindows()
{
   if (!fMapSubwindows) {
      //MapWindow();
      return;
   }

   TGWindow::MapSubwindows();

   if (!fList) return;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      if (el->fFrame) {
         el->fFrame->MapSubwindows();
         TGFrameElement *fe = el->fFrame->GetFrameElement();
         if (fe) fe->fState |= kIsVisible;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Hide sub frame.

void TGCompositeFrame::HideFrame(TGFrame *f)
{
   TGFrameElement *el = FindFrameElement(f);

   if (el) {
      el->fState = 0;
      el->fFrame->UnmapWindow();
      Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show sub frame.

void TGCompositeFrame::ShowFrame(TGFrame *f)
{
   TGFrameElement *el = FindFrameElement(f);

   if (el) {
      el->fState = 1;
      el->fFrame->MapWindow();
      Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get state of sub frame.

Int_t TGCompositeFrame::GetState(TGFrame *f) const
{
   TGFrameElement *el = FindFrameElement(f);

   if (el)
      return el->fState;
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get state of sub frame.

Bool_t TGCompositeFrame::IsVisible(TGFrame *f) const
{
   TGFrameElement *el = FindFrameElement(f);

   if (el)
      return (el->fState & kIsVisible);
   else
      return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get state of sub frame.

Bool_t TGCompositeFrame::IsArranged(TGFrame *f) const
{
   TGFrameElement *el = FindFrameElement(f);

   if (el)
      return (el->fState & kIsArranged);
   else
      return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Layout the elements of the composite frame.

void TGCompositeFrame::Layout()
{
   if (IsLayoutBroken()) return;
   fLayoutManager->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Print all frames in this composite frame.

void TGCompositeFrame::Print(Option_t *option) const
{
   TString opt = option;
   if (opt.Contains("tree")) {
      TGWindow::Print(option);
      return;
   }

   TGFrameElement *el;
   TIter next(fList);
   TString tab = option;

   TGFrame::Print(tab.Data());
   tab += "   ";
   while ((el = (TGFrameElement*)next())) {
      el->fFrame->Print(tab.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color for this frame and all subframes.

void TGCompositeFrame::ChangeSubframesBackground(Pixel_t back)
{
   TGFrame::ChangeBackground(back);
   TGFrameElement *el;

   TIter next(fList);

   while ((el = (TGFrameElement*)next())) {
      el->fFrame->SetBackgroundColor(back);
      if (el->fFrame->InheritsFrom(TGCompositeFrame::Class())) {
         ((TGCompositeFrame*)el->fFrame)->ChangeSubframesBackground(back);
      }
      fClient->NeedRedraw(el->fFrame);
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Get frame located at specified point.

TGFrame *TGCompositeFrame::GetFrameFromPoint(Int_t x, Int_t y)
{
   if (!Contains(x, y)) return 0;

   if (!fList) return this;

   TGFrame *f;
   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      //if (el->fFrame->IsVisible()) { //for this need to move IsVisible to TGFrame
      if (el->fState & kIsVisible) {
         f = el->fFrame->GetFrameFromPoint(x - el->fFrame->GetX(),
                                           y - el->fFrame->GetY());
         if (f) return f;
      }
   }
   return this;
}

////////////////////////////////////////////////////////////////////////////////
/// Translate coordinates to child frame.

Bool_t TGCompositeFrame::TranslateCoordinates(TGFrame *child, Int_t x, Int_t y,
                                              Int_t &fx, Int_t &fy)
{
   if (child == this) {
      fx = x;
      fy = y;
      return kTRUE;
   }

   if (!Contains(x, y)) return kFALSE;

   if (!fList) return kFALSE;

   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {
      if (el->fFrame == child) {
         fx = x - el->fFrame->GetX();
         fy = y - el->fFrame->GetY();
         return kTRUE;
      } else if (el->fFrame->IsComposite()) {
         if (((TGCompositeFrame *)el->fFrame)->TranslateCoordinates(child,
              x - el->fFrame->GetX(), y - el->fFrame->GetY(), fx, fy))
            return kTRUE;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag enter event.

Bool_t TGCompositeFrame::HandleDragEnter(TGFrame *)
{
   if (fClient && fClient->IsEditable() &&
       (fId != fClient->GetRoot()->GetId())) {

      // the dragged frame cannot be droppped
      if (fEditDisabled & (kEditDisable | kEditDisableLayout)) return kFALSE;

      //
      if (IsEditable()) {
         return kTRUE;
      }

      Float_t r, g, b;
      TColor::Pixel2RGB(fBackground, r, g, b);
      r *= 1.12;
      g *= 1.13;
      b *= 1.12;
      Pixel_t back = TColor::RGB2Pixel(r, g, b);
      gVirtualX->SetWindowBackground(fId, back);
      DoRedraw();
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag leave event.

Bool_t TGCompositeFrame::HandleDragLeave(TGFrame *)
{
   if (fClient && fClient->IsEditable() &&
       (fId != fClient->GetRoot()->GetId())) {

      if (fEditDisabled & (kEditDisable | kEditDisableLayout)) return kFALSE;

      gVirtualX->SetWindowBackground(fId, fBackground);
      DoRedraw();
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag motion event.

Bool_t TGCompositeFrame::HandleDragMotion(TGFrame *)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drop event.

Bool_t TGCompositeFrame::HandleDragDrop(TGFrame *frame, Int_t x, Int_t y,
                                        TGLayoutHints *lo)
{
   if (fClient && fClient->IsEditable() && frame && (x >= 0) && (y >= 0) &&
       (x + frame->GetWidth() <= fWidth) && (y + frame->GetHeight() <= fHeight)) {

      if (fEditDisabled & (kEditDisable | kEditDisableLayout)) return kFALSE;

      frame->ReparentWindow(this, x, y);
      AddFrame(frame, lo);
      frame->MapWindow();
      SetEditable(kTRUE);
      return kTRUE;
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a top level main frame. A main frame interacts
/// with the window manager.

TGMainFrame::TGMainFrame(const TGWindow *p, UInt_t w, UInt_t h,
        UInt_t options) : TGCompositeFrame(p, w, h, options | kMainFrame)
{
   // WMDeleteNotify causes the system to send a kClientMessage to the
   // window with fFormat=32 and fUser[0]=gWM_DELETE_WINDOW when window
   // closed via WM

   gVirtualX->WMDeleteNotify(fId);

   fBindList = new TList;

   fMWMValue    = 0;
   fMWMFuncs    = 0;
   fMWMInput    = 0;
   fWMX         = -1;
   fWMY         = -1;
   fWMWidth     = (UInt_t) -1;
   fWMHeight    = (UInt_t) -1;
   fWMMinWidth  = (UInt_t) -1;
   fWMMinHeight = (UInt_t) -1;
   fWMMaxWidth  = (UInt_t) -1;
   fWMMaxHeight = (UInt_t) -1;
   fWMWidthInc  = (UInt_t) -1;
   fWMHeightInc = (UInt_t) -1;
   fWMInitState = (EInitialState) 0;

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_s),
                      kKeyControlMask, kTRUE);
   if (p == fClient->GetDefaultRoot()) {
      fMWMValue    = kMWMDecorAll;
      fMWMFuncs    = kMWMFuncAll;
      fMWMInput    = kMWMInputModeless;
      gVirtualX->SetMWMHints(fId, fMWMValue, fMWMFuncs, fMWMInput);
   }
   // if parent is editing/embedable add this frame to the parent
   if (fClient->IsEditable() && (p == fClient->GetRoot())) {
      TGCompositeFrame *frame;
      if (p && p->InheritsFrom(TGCompositeFrame::Class())) {
         frame = (TGCompositeFrame*)p;
         frame->AddFrame(this, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

         // used during paste operation
         if (gDragManager && gDragManager->IsPasting()) {
            gDragManager->SetPasteFrame(this);
         }
      }
   }
   //AddInput(kButtonPressMask); // to allow Drag and Drop
   // Create Drag&Drop Manager and define a few DND types
   fDNDTypeList = new Atom_t[3];
   fDNDTypeList[0] = gVirtualX->InternAtom("application/root", kFALSE);
   fDNDTypeList[1] = gVirtualX->InternAtom("text/uri-list", kFALSE);
   fDNDTypeList[2] = 0;
   if (!gDNDManager)
      new TGDNDManager(this, fDNDTypeList);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// TGMainFrame destructor.

TGMainFrame::~TGMainFrame()
{
   delete [] fDNDTypeList;
   if (fBindList) {
      fBindList->Delete();
      delete fBindList;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Opens dialog window allowing user to save the frame contents
/// as a ROOT macro or as an image.
/// Returns kTRUE if something was saved.
///
/// This is bound to Ctrl-S by default.

Bool_t TGMainFrame::SaveFrameAsCodeOrImage()
{
   static TString dir(".");
   static Bool_t overwr = kFALSE;

   Bool_t repeat_save;
   do {
      repeat_save = kFALSE;

      TGFileInfo fi;
      TGMainFrame *main = (TGMainFrame*)GetMainFrame();
      fi.fFileTypes = gSaveMacroTypes;
      fi.SetIniDir(dir);
      fi.fOverwrite = overwr;
      new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);
      if (!fi.fFilename) return kFALSE;
      dir = fi.fIniDir;
      overwr = fi.fOverwrite;
      TString fname = gSystem->UnixPathName(fi.fFilename);
      if (fname.EndsWith(".C"))
         main->SaveSource(fname.Data(), "");
      else {
         TImage::EImageFileTypes gtype = TImage::kUnknown;
         if (fname.EndsWith("gif")) {
            gtype = TImage::kGif;
         } else if (fname.EndsWith(".png")) {
            gtype = TImage::kPng;
         } else if (fname.EndsWith(".jpg")) {
            gtype = TImage::kJpeg;
         } else if (fname.EndsWith(".tiff")) {
            gtype = TImage::kTiff;
         } else if (fname.EndsWith(".xpm")) {
            gtype = TImage::kXpm;
         }
         if (gtype != TImage::kUnknown) {
            Int_t saver = gErrorIgnoreLevel;
            gErrorIgnoreLevel = kFatal;
            TImage *img = TImage::Create();
            RaiseWindow();
            img->FromWindow(GetId());
            img->WriteImage(fname, gtype);
            gErrorIgnoreLevel = saver;
            delete img;
         }
         else {
            Int_t retval;
            new TGMsgBox(fClient->GetDefaultRoot(), this, "Error...",
                         TString::Format("file (%s) cannot be saved with this extension",
                                         fname.Data()), kMBIconExclamation,
                         kMBRetry | kMBCancel, &retval);
            repeat_save = (retval == kMBRetry);
         }
      }
   } while (repeat_save);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events.

Bool_t TGMainFrame::HandleKey(Event_t *event)
{
   if (fBindList) {

      TIter next(fBindList);
      TGMapKey *m;
      TGFrame  *w = 0;

      while ((m = (TGMapKey *) next())) {
         if (m->fKeyCode == event->fCode) {
            w = (TGFrame *) m->fWindow;
            if (w->HandleKey(event)) return kTRUE;
         }
      }
   }

   if ((event->fType == kGKeyPress) && (event->fState & kKeyControlMask)) {
      UInt_t keysym;
      char str[2];
      gVirtualX->LookupString(event, str, sizeof(str), keysym);

      if ((keysym & ~0x20) == kKey_S) { // case insensitive ctrl-s
         return SaveFrameAsCodeOrImage();
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Bind key to a window.

Bool_t TGMainFrame::BindKey(const TGWindow *w, Int_t keycode, Int_t modifier) const
{
   TList *list = fBindList;
   Handle_t id = fId;

   if (fClient->IsEditable()) {
      TGMainFrame *main = (TGMainFrame*)GetMainFrame();
      list = main->GetBindList();
      id = main->GetId();
   }

   if (list) {
      TGMapKey *m = new TGMapKey(keycode, (TGWindow *)w);
      list->Add(m);
      gVirtualX->GrabKey(id, keycode, modifier, kTRUE);
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove key binding.

void TGMainFrame::RemoveBind(const TGWindow *, Int_t keycode, Int_t modifier) const
{
   if (fBindList) {
      TIter next(fBindList);
      TGMapKey *m;
      while ((m = (TGMapKey *) next())) {
         if (m->fKeyCode == (UInt_t) keycode) {
            fBindList->Remove(m);
            delete m;
            gVirtualX->GrabKey(fId, keycode, modifier, kFALSE);
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button events.

Bool_t TGMainFrame::HandleButton(Event_t *event)
{
   if (event->fType == kButtonRelease) {
      if (gDNDManager->IsDragging()) gDNDManager->Drop();
   }
   return TGCompositeFrame::HandleButton(event);
}


////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events.

Bool_t TGMainFrame::HandleMotion(Event_t *event)
{
   if (gDNDManager && gDNDManager->IsDragging()) {
      gDNDManager->Drag(event->fXRoot, event->fYRoot,
                        TGDNDManager::GetDNDActionCopy(), event->fTime);
   }
   return TGCompositeFrame::HandleMotion(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle primary selection event.

Bool_t TGMainFrame::HandleSelection(Event_t *event)
{
   if ((Atom_t)event->fUser[1] == TGDNDManager::GetDNDSelection()) {
      if (gDNDManager)
         return gDNDManager->HandleSelection(event);
   }
   return TGCompositeFrame::HandleSelection(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle selection request event.

Bool_t TGMainFrame::HandleSelectionRequest(Event_t *event)
{
   if ((Atom_t)event->fUser[1] == TGDNDManager::GetDNDSelection()) {
      if (gDNDManager)
         return gDNDManager->HandleSelectionRequest(event);
   }
   return TGCompositeFrame::HandleSelectionRequest(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle client messages sent to this frame.

Bool_t TGMainFrame::HandleClientMessage(Event_t *event)
{
   TGCompositeFrame::HandleClientMessage(event);

   if ((event->fFormat == 32) && ((Atom_t)event->fUser[0] == gWM_DELETE_WINDOW) &&
       (event->fHandle != gROOT_MESSAGE)) {
      Emit("CloseWindow()");
      if (TestBit(kNotDeleted) && !TestBit(kDontCallClose))
         CloseWindow();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send close message to self. This method should be called from
/// a button to close this window.

void TGMainFrame::SendCloseMessage()
{
   Event_t event;

   event.fType   = kClientMessage;
   event.fFormat = 32;
   event.fHandle = gWM_DELETE_WINDOW;

   event.fWindow  = GetId();
   event.fUser[0] = (Long_t) gWM_DELETE_WINDOW;
   event.fUser[1] = 0;
   event.fUser[2] = 0;
   event.fUser[3] = 0;
   event.fUser[4] = 0;

   gVirtualX->SendEvent(GetId(), &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Close and delete main frame. We get here in response to ALT+F4 or
/// a window manager close command. To terminate the application when this
/// happens override this method and call gApplication->Terminate(0) or
/// make a connection to this signal (if after the slot this method
/// should not be called call DontCallClose() in the slot).
/// By default the window will be deleted.

void TGMainFrame::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Typically call this method in the slot connected to the CloseWindow()
/// signal to prevent the calling of the default or any derived CloseWindow()
/// methods to prevent premature or double deletion of this window.

void TGMainFrame::DontCallClose()
{
   SetBit(kDontCallClose);
}

////////////////////////////////////////////////////////////////////////////////
/// Set window name. This is typically done via the window manager.

void TGMainFrame::SetWindowName(const char *name)
{
   if (!name) {
      TGWindow::SetWindowName();
   } else {
      fWindowName = name;
      gVirtualX->SetWindowName(fId, (char *)name);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set window icon name. This is typically done via the window manager.

void TGMainFrame::SetIconName(const char *name)
{
   fIconName = name;
   gVirtualX->SetIconName(fId, (char *)name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set window icon pixmap by name. This is typically done via the window
/// manager. Icon can be in any image format supported by TImage, e.g.
/// GIF, XPM, PNG, JPG .. or even PS, PDF (see EImageFileTypes in TImage.h
/// for the full list of supported formats).
///
/// For example,
///    main_frame->SetIconPixmap("/home/root/icons/bld_rgb.png");

const TGPicture *TGMainFrame::SetIconPixmap(const char *iconName)
{
   fIconPixmap = iconName;
   const TGPicture *iconPic = fClient->GetPicture(iconName);
   if (iconPic) {
      Pixmap_t pic = iconPic->GetPicture();
      gVirtualX->SetIconPixmap(fId, pic);
      return iconPic;
   } else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set window icon by xpm array. That allows to have icons
/// builtin to the source code.
///
/// For example,
///    #include "/home/root/icons/bld_rgb.xpm"
///    //bld_rgb.xpm contains char *bld_rgb[] array
///    main_frame->SetIconPixmap(bld_rgb);

void TGMainFrame::SetIconPixmap(char **xpm_array)
{
   TImage *img = TImage::Create();
   if (!img) return;
   img->SetImageBuffer(xpm_array, TImage::kXpm);
   Pixmap_t pic = img->GetPixmap();
   if (pic) {
      gVirtualX->SetIconPixmap(fId, pic);
   } else {
      Warning("SetIconPixmap", "Failed to set window icon from xpm array.");
   }
   delete img;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the windows class and resource name. Used to get the right
/// resources from the resource database. However, ROOT applications
/// will typically use the .rootrc file for this.

void TGMainFrame::SetClassHints(const char *className, const char *resourceName)
{
   fClassName    = className;
   fResourceName = resourceName;
   gVirtualX->SetClassHints(fId, (char *)className, (char *)resourceName);
}

////////////////////////////////////////////////////////////////////////////////
/// Set decoration style for MWM-compatible wm (mwm, ncdwm, fvwm?).

void TGMainFrame::SetMWMHints(UInt_t value, UInt_t funcs, UInt_t input)
{
   if (fClient->IsEditable() && (fParent == fClient->GetRoot())) return;

   fMWMValue = value;
   fMWMFuncs = funcs;
   fMWMInput = input;
   gVirtualX->SetMWMHints(fId, value, funcs, input);
}

////////////////////////////////////////////////////////////////////////////////
/// Give the window manager a window position hint.

void TGMainFrame::SetWMPosition(Int_t x, Int_t y)
{
   if (fClient->IsEditable() && (fParent == fClient->GetRoot())) return;

   fWMX = x;
   fWMY = y;
   gVirtualX->SetWMPosition(fId, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Give the window manager a window size hint.

void TGMainFrame::SetWMSize(UInt_t w, UInt_t h)
{
   if (fClient->IsEditable() && (fParent == fClient->GetRoot())) return;

   fWMWidth  = w;
   fWMHeight = h;
   gVirtualX->SetWMSize(fId, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Give the window manager minimum and maximum size hints. Also
/// specify via winc and hinc the resize increments.

void TGMainFrame::SetWMSizeHints(UInt_t wmin, UInt_t hmin,
                                 UInt_t wmax, UInt_t hmax,
                                 UInt_t winc, UInt_t hinc)
{
   if (fClient->IsEditable() && (fParent == fClient->GetRoot())) return;

   fMinWidth    = fWMMinWidth  = wmin;
   fMinHeight   = fWMMinHeight = hmin;
   fMaxWidth    = fWMMaxWidth  = wmax;
   fMaxHeight   = fWMMaxHeight = hmax;
   fWMWidthInc  = winc;
   fWMHeightInc = hinc;
   gVirtualX->SetWMSizeHints(fId, wmin, hmin, wmax, hmax, winc, hinc);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the initial state of the window. Either kNormalState or kIconicState.

void TGMainFrame::SetWMState(EInitialState state)
{
   if (fClient->IsEditable() && (fParent == fClient->GetRoot())) return;

   fWMInitState = state;
   gVirtualX->SetWMState(fId, state);
}


////////////////////////////////////////////////////////////////////////////////
/// Create a transient window. A transient window is typically used for
/// dialog boxes.

TGTransientFrame::TGTransientFrame(const TGWindow *p, const TGWindow *main,
                                   UInt_t w, UInt_t h, UInt_t options)
   : TGMainFrame(p, w, h, options | kTransientFrame)
{
   fMain = main;
   if (!fMain && gClient)
      fMain = gClient->GetRoot();

   if (fMain) {
      gVirtualX->SetWMTransientHint(fId, fMain->GetId());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Position transient frame centered relative to the parent frame.
/// If fMain is 0 (i.e. TGTransientFrame is acting just like a
/// TGMainFrame) and croot is true, the window will be centered on
/// the root window, otherwise no action is taken and the default
/// wm placement will be used.

void TGTransientFrame::CenterOnParent(Bool_t croot, EPlacement pos)
{
   Int_t x=0, y=0, ax, ay;
   Window_t wdummy;

   UInt_t dw = fClient->GetDisplayWidth();
   UInt_t dh = fClient->GetDisplayHeight();

   if (fMain) {

      switch (pos) {
         case kCenter:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - fWidth) >> 1;
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - fHeight) >> 1;
            break;
         case kRight:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - (fWidth >> 1));
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - fHeight) >> 1;
            break;
         case kLeft:
            x = (Int_t)(-1 * (Int_t)(fWidth >> 1));
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - fHeight) >> 1;
            break;
         case kTop:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - fWidth) >> 1;
            y = (Int_t)(-1 * (Int_t)(fHeight >> 1));
            break;
         case kBottom:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - fWidth) >> 1;
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - (fHeight >> 1));
            break;
         case kTopLeft:
            x = (Int_t)(-1 * (Int_t)(fWidth >> 1));
            y = (Int_t)(-1 * (Int_t)(fHeight >> 1));
            break;
         case kTopRight:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - (fWidth >> 1));
            y = (Int_t)(-1 * (Int_t)(fHeight >> 1));
            break;
         case kBottomLeft:
            x = (Int_t)(-1 * (Int_t)(fWidth >> 1));
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - (fHeight >> 1));
            break;
         case kBottomRight:
            x = (Int_t)(((TGFrame *) fMain)->GetWidth() - (fWidth >> 1));
            y = (Int_t)(((TGFrame *) fMain)->GetHeight() - (fHeight >> 1));
            break;
      }

      gVirtualX->TranslateCoordinates(fMain->GetId(), GetParent()->GetId(),
                                      x, y, ax, ay, wdummy);
      if (!gVirtualX->InheritsFrom("TGWin32")) {
         if (ax < 10)
            ax = 10;
         else if (ax + fWidth + 10 > dw)
            ax = dw - fWidth - 10;

         if (ay < 20)
            ay = 20;
         else if (ay + fHeight + 50 > dh)
            ay = dh - fHeight - 50;
      }

   } else if (croot) {

      switch (pos) {
         case kCenter:
            x = (dw - fWidth) >> 1;
            y = (dh - fHeight) >> 1;
            break;
         case kRight:
            x = dw - (fWidth >> 1);
            y = (dh - fHeight) >> 1;
            break;
         case kLeft:
            x = -1 * (Int_t)(fWidth >> 1);
            y = (dh - fHeight) >> 1;
            break;
         case kTop:
            x = (dw - fWidth) >> 1;
            y = -1 * (Int_t)(fHeight >> 1);
            break;
         case kBottom:
            x = (dw - fWidth) >> 1;
            y = dh - (fHeight >> 1);
            break;
         case kTopLeft:
            x = -1 * (Int_t)(fWidth >> 1);
            y = -1 * (Int_t)(fHeight >> 1);
            break;
         case kTopRight:
            x = dw - (fWidth >> 1);
            y = -1 * (Int_t)(fHeight >> 1);
            break;
         case kBottomLeft:
            x = -1 * (Int_t)(fWidth >> 1);
            y = dh - (fHeight >> 1);
            break;
         case kBottomRight:
            x = dw - (fWidth >> 1);
            y = dh - (fHeight >> 1);
            break;
      }

      ax = x;
      ay = y;

   } else {

      return;

   }

   Move(ax, ay);
   SetWMPosition(ax, ay);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a group frame. The title will be adopted and deleted by the
/// group frame.

TGGroupFrame::TGGroupFrame(const TGWindow *p, TGString *title,
                           UInt_t options, GContext_t norm,
                           FontStruct_t font, Pixel_t back) :
   TGCompositeFrame(p, 1, 1, options, back)
{
   fText       = title;
   fFontStruct = font;
   fNormGC     = norm;
   fTitlePos   = kLeft;
   fHasOwnFont = kFALSE;

   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fBorderWidth = max_ascent + max_descent + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a group frame.

TGGroupFrame::TGGroupFrame(const TGWindow *p, const char *title,
                           UInt_t options, GContext_t norm,
                           FontStruct_t font, Pixel_t back) :
   TGCompositeFrame(p, 1, 1, options, back)
{
   fText       = new TGString(!p && !title ? GetName() : title);
   fFontStruct = font;
   fNormGC     = norm;
   fTitlePos   = kLeft;
   fHasOwnFont = kFALSE;

   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fBorderWidth = max_ascent + max_descent + 1;

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a group frame.

TGGroupFrame::~TGGroupFrame()
{
   if (fHasOwnFont) {
      TGGCPool *pool = fClient->GetGCPool();
      TGGC *gc = pool->FindGC(fNormGC);
      pool->FreeGC(gc);
   }
   delete fText;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns default size.

TGDimension TGGroupFrame::GetDefaultSize() const
{
   UInt_t tw = gVirtualX->TextWidth(fFontStruct, fText->GetString(),
                                    fText->GetLength()) + 24;

   TGDimension dim = TGCompositeFrame::GetDefaultSize();

   return  tw>dim.fWidth ? TGDimension(tw, dim.fHeight) : dim;
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the group frame. Need special DoRedraw() since we need to
/// redraw with fBorderWidth=0.

void TGGroupFrame::DoRedraw()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);

   DrawBorder();
}


////////////////////////////////////////////////////////////////////////////////
/// Changes text color.
/// If local is true color is changed locally, otherwise - globally.

void TGGroupFrame::SetTextColor(Pixel_t color, Bool_t local)
{
   TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
   TGGC *gc = pool->FindGC(fNormGC);

   if (gc && local) {
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
/// Changes text font.
/// If local is true font is changed locally - otherwise globally.

void TGGroupFrame::SetTextFont(FontStruct_t font, Bool_t local)
{
   FontH_t v = gVirtualX->GetFontHandle(font);
   if (!v) return;

   fFontStruct = font;

   TGGCPool *pool =  fClient->GetResourcePool()->GetGCPool();
   TGGC *gc = pool->FindGC(fNormGC);

   if (gc && local) {
      gc = pool->GetGC((GCValues_t*)gc->GetAttributes(), kTRUE); // copy
      fHasOwnFont = kTRUE;
   }
   if (gc) {
      gc->SetFont(v);
      fNormGC = gc->GetGC();
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Changes text font specified by name.
/// If local is true font is changed locally - otherwise globally.

void TGGroupFrame::SetTextFont(const char *fontName, Bool_t local)
{
   TGFont *font = fClient->GetFont(fontName);

   if (font) {
      SetTextFont(font->GetFontStruct(), local);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if text attributes are unique,
/// returns kFALSE if text attributes are shared (global).

Bool_t TGGroupFrame::HasOwnFont() const
{
   return fHasOwnFont;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw border of around the group frame.
///
/// if frame is kRaisedFrame  - a frame border is of "wall style",
/// otherwise of "groove style".

void TGGroupFrame::DrawBorder()
{
   Int_t x, y, l, t, r, b, gl, gr, sep, max_ascent, max_descent;

   UInt_t tw = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   l = 0;
   t = (max_ascent + max_descent + 2) >> 1;
   r = fWidth - 1;
   // next three lines are for backward compatibility in case of horizontal layout
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if ((lm->InheritsFrom(TGHorizontalLayout::Class())) ||
       (lm->InheritsFrom(TGMatrixLayout::Class())))
      b = fHeight - 1;
   else
      b = fHeight - t;

   sep = 3;
   UInt_t rr = 5 + (sep << 1) + tw;

   switch (fTitlePos) {
      case kRight:
         gl = fWidth>rr ? Int_t(fWidth - rr) : 5 + sep;
         break;
      case kCenter:
         gl = fWidth>tw ? Int_t((fWidth - tw)>>1) - sep : 5 + sep;
         break;
      case kLeft:
      default:
         gl = 5 + sep;
   }
   gr = gl + tw + (sep << 1);

   switch (fOptions & (kSunkenFrame | kRaisedFrame)) {
      case kRaisedFrame:
         gVirtualX->DrawLine(fId, GetHilightGC()(),  l,   t,   gl,  t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), l+1, t+1, gl,  t+1);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  gr,  t,   r-1, t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), gr,  t+1, r-2, t+1);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  r-1, t,   r-1, b-1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), r,   t,   r,   b);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  r-1, b-1, l,   b-1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), r,   b,   l,   b);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  l,   b-1, l,   t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), l+1, b-2, l+1, t+1);
         break;
      case kSunkenFrame:
      default:
         gVirtualX->DrawLine(fId, GetShadowGC()(),  l,   t,   gl,  t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), l+1, t+1, gl,  t+1);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  gr,  t,   r-1, t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), gr,  t+1, r-2, t+1);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  r-1, t,   r-1, b-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), r,   t,   r,   b);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  r-1, b-1, l,   b-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), r,   b,   l,   b);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  l,   b-1, l,   t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), l+1, b-2, l+1, t+1);
         break;
   }

   x = gl + sep;
   y = 1;

   fText->Draw(fId, fNormGC, x, y + max_ascent);
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change title of the group frame. Titlte TGString is adopted
/// by the TGGroupFrame.

void TGGroupFrame::SetTitle(TGString *title)
{
   if (!title) {
      Warning("SetTitle", "title cannot be 0, try \"\"");
      title = new TGString("");
   }

   delete fText;

   fText = title;
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change title of the group frame.

void TGGroupFrame::SetTitle(const char *title)
{
   if (!title) {
      Error("SetTitle", "title cannot be 0, try \"\"");
      return;
   }

   SetTitle(new TGString(title));
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure in use.

FontStruct_t TGGroupFrame::GetDefaultFontStruct()
{
   if (!fgDefaultFont && gClient)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context in use.

const TGGC &TGGroupFrame::GetDefaultGC()
{
   if (!fgDefaultGC && gClient)
      fgDefaultGC = gClient->GetResourcePool()->GetFrameGC();
   return *fgDefaultGC;
}

////////////////////////////////////////////////////////////////////////////////
/// Header Frame constructor.

TGHeaderFrame::TGHeaderFrame(const TGWindow *p, UInt_t w, UInt_t h,
                 UInt_t options, Pixel_t back) :
  TGHorizontalFrame(p, w, h, options | kVerticalFrame, back)
{
   fSplitCursor = kNone;
   fSplitCursor = gVirtualX->CreateCursor(kArrowHor);
   fOverSplitter = false;
   fOverButton = -1;
   fLastButton = -1;
   fNColumns   = 1;
   fColHeader  = 0;
   fSplitHeader = 0;

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask,
                         kNone, kNone);
   AddInput(kPointerMotionMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Set columns information in the header frame.

void TGHeaderFrame::SetColumnsInfo(Int_t nColumns, TGTextButton  **colHeader,
               TGVFileSplitter  **splitHeader)
{
   fNColumns = nColumns;
   fColHeader = colHeader;
   fSplitHeader = splitHeader;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in header frame.

Bool_t TGHeaderFrame::HandleButton(Event_t* event)
{
   if ( event->fY > 0 &&
        event->fY <= (Int_t) this->GetHeight() ) {
      for (Int_t i = 1; i < fNColumns; ++i ) {
         if ( event->fX < fColHeader[i]->GetX() &&
            event->fX >= fColHeader[i-1]->GetX() ) {
            if ( fOverSplitter ) {
               if ( event->fX <= fColHeader[i-1]->GetX() + 5 )
                  fSplitHeader[i-2]->HandleButton(event);
               else
                  fSplitHeader[i-1]->HandleButton(event);
            } else {
               if ( event->fType == kButtonPress ) {
                  fLastButton = i - 1;
               } else {
                  fLastButton = -1;
               }
               event->fX -= fColHeader[i-1]->GetX();
               fColHeader[i-1]->HandleButton(event);
            }
            break;
         }
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click mouse event in header frame.

Bool_t TGHeaderFrame::HandleDoubleClick(Event_t *event)
{
   if ( event->fY > 0 &&
        event->fY <= (Int_t) this->GetHeight() ) {
      for (Int_t i = 1; i < fNColumns; ++i ) {
         if ( event->fX < fColHeader[i]->GetX() &&
            event->fX >= fColHeader[i-1]->GetX() ) {
            if ( fOverSplitter ) {
               if ( event->fX <= fColHeader[i-1]->GetX() + 5 )
                  fSplitHeader[i-2]->HandleDoubleClick(event);
               else
                  fSplitHeader[i-1]->HandleDoubleClick(event);
            } else {
               event->fX -= fColHeader[i-1]->GetX();
               fColHeader[i-1]->HandleDoubleClick(event);
            }
            break;
         }
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events in header frame.

Bool_t TGHeaderFrame::HandleMotion(Event_t* event)
{
   if ( event->fY > 0 &&
        event->fY <= (Int_t) this->GetHeight() ) {
      Bool_t inMiddle = false;

      for (Int_t i = 1; i < fNColumns; ++i ) {
         if ( event->fX > fColHeader[i]->GetX() - 5 &&
            event->fX < fColHeader[i]->GetX() + 5 ) {
            inMiddle = true;
         }
         if ( event->fX < fColHeader[i]->GetX() &&
            event->fX >= fColHeader[i-1]->GetX() ) {
            fOverButton = i - 1;
         }
      }
      fOverSplitter = inMiddle;
      if ( fOverSplitter ) {
         gVirtualX->SetCursor(fId, fSplitCursor);
      }
      else {
         gVirtualX->SetCursor(fId, kNone);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a user color in a C++ macro file - used in SavePrimitive().

void TGFrame::SaveUserColor(std::ostream &out, Option_t *option)
{
   char quote = '"';

   if (gROOT->ClassSaved(TGFrame::Class())) {
      out << std::endl;
   } else {
      //  declare a color variable to reflect required user changes
      out << std::endl;
      out << "   ULong_t ucolor;        // will reflect user color changes" << std::endl;
   }
   ULong_t ucolor;
   if (option && !strcmp(option, "slider"))
      ucolor = GetDefaultFrameBackground();
   else
      ucolor = GetBackground();
   if ((ucolor != fgUserColor) || (ucolor == GetWhitePixel())) {
      const char *ucolorname = TColor::PixelAsHexString(ucolor);
      out << "   gClient->GetColorByName(" << quote << ucolorname << quote
          << ",ucolor);" << std::endl;
      fgUserColor = ucolor;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a frame option string - used in SavePrimitive().

TString TGFrame::GetOptionString() const
{
   TString options;

   if (!GetOptions()) {
      options = "kChildFrame";
   } else {
      if (fOptions & kMainFrame) {
         if (options.Length() == 0) options  = "kMainFrame";
         else                       options += " | kMainFrame";
      }
      if (fOptions & kVerticalFrame) {
         if (options.Length() == 0) options  = "kVerticalFrame";
         else                       options += " | kVerticalFrame";
      }
      if (fOptions & kHorizontalFrame) {
         if (options.Length() == 0) options  = "kHorizontalFrame";
         else                       options += " | kHorizontalFrame";
      }
      if (fOptions & kSunkenFrame) {
         if (options.Length() == 0) options  = "kSunkenFrame";
         else                       options += " | kSunkenFrame";
      }
      if (fOptions & kRaisedFrame) {
         if (options.Length() == 0) options  = "kRaisedFrame";
         else                       options += " | kRaisedFrame";
      }
      if (fOptions & kDoubleBorder) {
         if (options.Length() == 0) options  = "kDoubleBorder";
         else                       options += " | kDoubleBorder";
      }
      if (fOptions & kFitWidth) {
         if (options.Length() == 0) options  = "kFitWidth";
         else                       options += " | kFitWidth";
      }
      if (fOptions & kFixedWidth) {
         if (options.Length() == 0) options  = "kFixedWidth";
         else                       options += " | kFixedWidth";
      }
      if (fOptions & kFitHeight) {
         if (options.Length() == 0) options  = "kFitHeight";
         else                       options += " | kFitHeight";
      }
      if (fOptions & kFixedHeight) {
         if (options.Length() == 0) options  = "kFixedHeight";
         else                       options += " | kFixedHeight";
      }
      if (fOptions & kOwnBackground) {
         if (options.Length() == 0) options  = "kOwnBackground";
         else                       options += " | kOwnBackground";
      }
      if (fOptions & kTransientFrame) {
         if (options.Length() == 0) options  = "kTransientFrame";
         else                       options += " | kTransientFrame";
      }
      if (fOptions & kTempFrame) {
         if (options.Length() == 0) options  = "kTempFrame";
         else                       options += " | kTempFrame";
      }
   }
   return options;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns MWM decoration hints as a string - used in SavePrimitive().

TString TGMainFrame::GetMWMvalueString() const
{
   TString hints;

   if (fMWMValue) {
      if (fMWMValue & kMWMDecorAll) {
         if (hints.Length() == 0) hints  = "kMWMDecorAll";
         else                     hints += " | kMWMDecorAll";
      }
      if (fMWMValue & kMWMDecorBorder) {
         if (hints.Length() == 0) hints  = "kMWMDecorBorder";
         else                     hints += " | kMWMDecorBorder";
      }
      if (fMWMValue & kMWMDecorResizeH) {
         if (hints.Length() == 0) hints  = "kMWMDecorResizeH";
         else                     hints += " | kMWMDecorResizeH";
      }
      if (fMWMValue & kMWMDecorTitle) {
         if (hints.Length() == 0) hints  = "kMWMDecorTitle";
         else                     hints += " | kMWMDecorTitle";
      }
      if (fMWMValue & kMWMDecorMenu) {
         if (hints.Length() == 0) hints  = "kMWMDecorMenu";
         else                     hints += " | kMWMDecorMenu";
      }
      if (fMWMValue & kMWMDecorMinimize) {
         if (hints.Length() == 0) hints  = "kMWMDecorMinimize";
         else                     hints += " | kMWMDecorMinimize";
      }
      if (fMWMValue & kMWMDecorMaximize) {
         if (hints.Length() == 0) hints  = "kMWMDecorMaximize";
         else                     hints += " | kMWMDecorMaximize";
      }
   }
   return hints;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns MWM function hints as a string - used in SavePrimitive().

TString TGMainFrame::GetMWMfuncString() const
{
   TString hints;

   if (fMWMFuncs) {

      if (fMWMFuncs & kMWMFuncAll) {
         if (hints.Length() == 0) hints  = "kMWMFuncAll";
         else                     hints += " | kMWMFuncAll";
      }
      if (fMWMFuncs & kMWMFuncResize) {
         if (hints.Length() == 0) hints  = "kMWMFuncResize";
         else                     hints += " | kMWMFuncResize";
      }
      if (fMWMFuncs & kMWMFuncMove) {
         if (hints.Length() == 0) hints  = "kMWMFuncMove";
         else                     hints += " | kMWMFuncMove";
      }
      if (fMWMFuncs & kMWMFuncMinimize) {
         if (hints.Length() == 0) hints  = "kMWMFuncMinimize";
         else                     hints += " | kMWMFuncMinimize";
      }
      if (fMWMFuncs & kMWMFuncMaximize) {
         if (hints.Length() == 0) hints  = "kMWMFuncMaximize";
         else                     hints += " | kMWMFuncMaximize";
      }
      if (fMWMFuncs & kMWMFuncClose) {
         if (hints.Length() == 0) hints  = "kMWMFuncClose";
         else                     hints += " | kMWMFuncClose";
      }
   }
   return hints;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns MWM input mode hints as a string - used in SavePrimitive().

TString TGMainFrame::GetMWMinpString() const
{
   TString hints;

   if (fMWMInput == 0) hints = "kMWMInputModeless";

   if (fMWMInput == 1) hints = "kMWMInputPrimaryApplicationModal";

   if (fMWMInput == 2) hints = "kMWMInputSystemModal";

   if (fMWMInput == 3) hints = "kMWMInputFullApplicationModal";

   return hints;
}

////////////////////////////////////////////////////////////////////////////////
/// Auxilary protected method  used to save subframes.

void TGCompositeFrame::SavePrimitiveSubframes(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fLayoutBroken)
      out << "   " << GetName() << "->SetLayoutBroken(kTRUE);" << std::endl;

   if (!fList) return;

   char quote = '"';

   TGFrameElement *el;
   static TGHSplitter *hsplit = 0;
   static TGVSplitter *vsplit = 0;
   TList *signalslist;
   TList *connlist;
   TQConnection *conn;
   TString signal_name, slot_name;

   TIter next(fList);

   while ((el = (TGFrameElement *) next())) {

      // Don't save hidden (unmapped) frames having a parent different
      // than this frame. Solves a problem with shared frames
      // (e.g. shared menus in the new Browser)
      if ((!(el->fState & kIsVisible)) && (el->fFrame->GetParent() != this))
         continue;

      // Remember if the frame to be saved is a TG(H,V)Splitter
      // See comments below and in TG[H/V]Splitter::SavePrimitive()
      if (el->fFrame->InheritsFrom("TGVSplitter")) {
         vsplit = (TGVSplitter *)el->fFrame;
         if (vsplit->GetLeft())
            vsplit = 0;
      }
      else if (el->fFrame->InheritsFrom("TGHSplitter")) {
         hsplit = (TGHSplitter *)el->fFrame;
         if (hsplit->GetAbove())
            hsplit = 0;
      }
      el->fFrame->SavePrimitive(out, option);
      out << "   " << GetName() << "->AddFrame(" << el->fFrame->GetName();
      el->fLayout->SavePrimitive(out, option);
      out << ");"<< std::endl;
      if (IsLayoutBroken()) {
         out << "   " << el->fFrame->GetName() << "->MoveResize(";
         out << el->fFrame->GetX() << "," << el->fFrame->GetY() << ",";
         out << el->fFrame->GetWidth() << ","  << el->fFrame->GetHeight();
         out << ");" << std::endl;
      }
      // TG(H,V)Splitter->SetFrame(theframe) can only be saved _AFTER_
      // having saved "theframe", when "theframe" is either at right
      // or below the splitter (that means after the splitter in the
      // list of frames), otherwise "theframe" would be undefined
      // (aka used before to be created)...
      if (vsplit && el->fFrame == vsplit->GetFrame()) {
         out << "   " << vsplit->GetName() << "->SetFrame(" << vsplit->GetFrame()->GetName();
         if (vsplit->GetLeft()) out << ",kTRUE);" << std::endl;
         else                 out << ",kFALSE);"<< std::endl;
         vsplit = 0;
      }
      if (hsplit && el->fFrame == hsplit->GetFrame()) {
         out << "   " << hsplit->GetName() << "->SetFrame(" << hsplit->GetFrame()->GetName();
         if (hsplit->GetAbove()) out << ",kTRUE);" << std::endl;
         else                  out << ",kFALSE);"<< std::endl;
         hsplit = 0;
      }

      if (!(el->fState & kIsVisible)) {
         gListOfHiddenFrames->Add(el->fFrame);
      }

      // saving signals/slots
      signalslist = (TList*)el->fFrame->GetListOfSignals();
      if (!signalslist)  continue;
      connlist = (TList*)signalslist->Last();
      if (connlist) {
         conn = (TQConnection*)connlist->Last();
         if (conn) {
            signal_name = connlist->GetName();
            slot_name = conn->GetName();
            Int_t eq = slot_name.First('=');
            Int_t rb = slot_name.First(')');
            if (eq != -1)
               slot_name.Remove(eq, rb-eq);
            out << "   " << el->fFrame->GetName() << "->Connect(" << quote << signal_name
                << quote << ", 0, 0, " << quote << slot_name << quote << ");" << std::endl;

            TList *lsl = (TList *)gROOT->GetListOfSpecials()->FindObject("ListOfSlots");
            if (lsl) {
               TObjString *slotel = (TObjString *)lsl->FindObject(slot_name);
               if (!slotel)
                  lsl->Add(new TObjString(slot_name));
            }
         }
      }
   }
   out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a composite frame widget as a C++ statement(s) on output stream out.

void TGCompositeFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   if (!strcmp(GetName(),"")) {
      SetName(Form("fCompositeframe%d",fgCounter));
      fgCounter++;
   }

   out << std::endl << "   // composite frame" << std::endl;
   out << "   TGCompositeFrame *";
   out << GetName() << " = new TGCompositeFrame(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out << ");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager if it differs from the composite frame type
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager *lm = GetLayoutManager();
   if ((GetOptions() & kHorizontalFrame) &&
       (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
      ;
   } else if ((GetOptions() & kVerticalFrame) &&
              (lm->InheritsFrom(TGVerticalLayout::Class()))) {
      ;
   } else {
      out << "   " << GetName() <<"->SetLayoutManager(";
      lm->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }

   SavePrimitiveSubframes(out, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the GUI main frame widget in a C++ macro file.

void TGMainFrame::SaveSource(const char *filename, Option_t *option)
{
   // iteration over all active classes to exclude the base ones
   TString opt = option;
   TBits *bc = new TBits();
   TClass *c1, *c2, *c3;
   UInt_t k = 0;      // will mark k-bit of TBits if the class is a base class

   TIter nextc1(gROOT->GetListOfClasses());
   //gROOT->GetListOfClasses()->ls();    // valid. test
   while((c1 = (TClass *)nextc1())) {

      //   resets bit TClass::kClassSaved for all classes
      c1->ResetBit(TClass::kClassSaved);

      TIter nextc2(gROOT->GetListOfClasses());
      while ((c2 = (TClass *)nextc2())) {
         if (c1==c2) continue;
         else {
            c3 = c2->GetBaseClass(c1);
            if (c3 != 0) {
               bc->SetBitNumber(k, kTRUE);
               break;
            }
         }
      }
      k++;
   }

   TList *ilist = new TList();   // will contain include file names without '.h'
   ilist->SetName("ListOfIncludes");
   gROOT->GetListOfSpecials()->Add(ilist);
   k=0;

   //   completes list of include file names
   TIter nextdo(gROOT->GetListOfClasses());
   while ((c2 = (TClass *)nextdo())) {
      // for used GUI header files
      if (bc->TestBitNumber(k) == 0 && c2->InheritsFrom(TGObject::Class()) == 1) {
         // for any used ROOT header files activate the line below, comment the line above
         //if (bc->TestBitNumber(k) == 0) {
         const char *iname;
         iname = c2->GetDeclFileName();
         if (iname[0] && strstr(iname,".h")) {
            const char *lastsl = strrchr(iname,'/');
            if (lastsl) iname = lastsl + 1;
            char *tname = new char[strlen(iname)+1];
            Int_t i=0;
            while (*iname != '.') {
               tname[i] = *iname;
               i++; iname++;
            }
            tname[i] = 0;    //tname = include file name without '.h'

            TObjString *iel = (TObjString *)ilist->FindObject(tname);
            if (!iel) {
               ilist->Add(new TObjString(tname));
            }
            // Weird, but when saving a canvas, the following two classes
            // may be missing if the toolbar has not been displayed...
            if (strstr(tname, "TRootCanvas")) {
               if (!ilist->FindObject("TGDockableFrame"))
                  ilist->Add(new TObjString("TGDockableFrame"));
               if (!ilist->FindObject("TG3DLine"))
                  ilist->Add(new TObjString("TG3DLine"));
            }
            delete [] tname;
         }
         k++;  continue;
      }
      k++;
   }

   char quote = '"';
   std::ofstream out;

   TString ff = filename && strlen(filename) ? filename : "Rootappl.C";

   // Computes the main method name.
   const char *fname = gSystem->BaseName(ff.Data());
   Int_t lenfname = strlen(fname);
   char *sname = new char[lenfname+1];

   Int_t i = 0;
   while ((*fname != '.') && (i < lenfname)) {
      sname[i] = *fname;
      i++; fname++;
   }
   if (i == lenfname)
      ff += ".C";
   sname[i] = 0;

   out.open(ff.Data(), std::ios::out);
   if (!out.good()) {
      Error("SaveSource", "cannot open file: %s", ff.Data());
      delete [] sname;
      return;
   }

   // writes include files in C++ macro
   TObjString *inc;
   ilist = (TList *)gROOT->GetListOfSpecials()->FindObject("ListOfIncludes");

   if (!ilist) {
      delete [] sname;
      return;
   }

   // write macro header, date/time stamp as string, and the used Root version
   TDatime t;
   out <<"// Mainframe macro generated from application: "<< gApplication->Argv(0) << std::endl;
   out <<"// By ROOT version "<< gROOT->GetVersion() <<" on "<<t.AsSQLString()<< std::endl;
   out << std::endl;

   TIter nexti(ilist);
   while((inc = (TObjString *)nexti())) {
         out << "#ifndef ROOT_" << inc->GetString() << std::endl;
         out << "#include " << quote << inc->GetString() << ".h" << quote << std::endl;
         out << "#endif" << std::endl;
         if (strstr(inc->GetString(),"TRootEmbeddedCanvas")) {
            out << "#ifndef ROOT_TCanvas" << std::endl;
            out << "#include " << quote << "TCanvas.h" << quote << std::endl;
            out << "#endif" << std::endl;
         }
   }
   out << std::endl << "#include " << quote << "Riostream.h" << quote << std::endl;
   // deletes created ListOfIncludes
   gROOT->GetListOfSpecials()->Remove(ilist);
   ilist->Delete();
   delete ilist;
   delete bc;

   // writes the macro entry point equal to the fname
   out << std::endl;
   out << "void " << sname << "()" << std::endl;
   out <<"{"<< std::endl;
   delete [] sname;

   gListOfHiddenFrames->Clear();

   // saivng slots
   TList *lSlots = new TList;
   lSlots->SetName("ListOfSlots");
   gROOT->GetListOfSpecials()->Add(lSlots);

   TGMainFrame::SavePrimitive(out, option);

   if (strlen(fClassName) || strlen(fResourceName)) {
      out << "   " << GetName() << "->SetClassHints(" << quote << fClassName
          << quote << "," << quote << fResourceName << quote << ");" << std::endl;
   }

   GetMWMHints(fMWMValue, fMWMFuncs, fMWMInput);
   if (fMWMValue || fMWMFuncs || fMWMInput) {
      out << "   " << GetName() << "->SetMWMHints(";
      out << GetMWMvalueString() << "," << std::endl;
      out << "                        ";
      out << GetMWMfuncString() << "," << std::endl;
      out << "                        ";
      out << GetMWMinpString() << ");"<< std::endl;
   }

///   GetWMPosition(fWMX, fWMY);
///   if ((fWMX != -1) || (fWMY != -1)) {
///      out <<"   "<<GetName()<<"->SetWMPosition("<<fWMX<<","<<fWMY<<");"<<std::endl;
///   }   // does not work - fixed via Move() below...

   GetWMSize(fWMWidth, fWMHeight);
   if (fWMWidth != UInt_t(-1) || fWMHeight != UInt_t(-1)) {
      out <<"   "<<GetName()<<"->SetWMSize("<<fWMWidth<<","<<fWMHeight<<");"<<std::endl;
   }

   GetWMSizeHints(fWMMinWidth, fWMMinHeight, fWMMaxWidth, fWMMaxHeight, fWMWidthInc, fWMHeightInc);
   if (fWMMinWidth != UInt_t(-1) || fWMMinHeight != UInt_t(-1) ||
      fWMMaxWidth != UInt_t(-1) || fWMMaxHeight != UInt_t(-1) ||
      fWMWidthInc != UInt_t(-1) || fWMHeightInc != UInt_t(-1)) {
      out <<"   "<<GetName()<<"->SetWMSizeHints("<<fWMMinWidth<<","<<fWMMinHeight
          <<","<<fWMMaxWidth<<","<<fWMMaxHeight
          <<","<<fWMWidthInc<<","<<fWMHeightInc <<");"<<std::endl;
   }

   out << "   " <<GetName()<< "->MapSubwindows();" << std::endl;

   TIter nexth(gListOfHiddenFrames);
   TGFrame *fhidden;
   while ((fhidden = (TGFrame*)nexth())) {
      out << "   " <<fhidden->GetName()<< "->UnmapWindow();" << std::endl;
   }

   out << std::endl;
   gListOfHiddenFrames->Clear();

   Bool_t usexy = kFALSE;
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if (lm->InheritsFrom("TGXYLayout"))
      usexy = kTRUE;

   if (!usexy)
      out << "   " <<GetName()<< "->Resize("<< GetName()<< "->GetDefaultSize());" << std::endl;
   else
      out << "   " <<GetName()<< "->Resize("<< GetWidth()<<","<<GetHeight()<<");"<<std::endl;

   out << "   " <<GetName()<< "->MapWindow();" <<std::endl;

   GetWMPosition(fWMX, fWMY);
   if ((fWMX != -1) || (fWMY != -1)) {
      out <<"   "<<GetName()<<"->Move("<<fWMX<<","<<fWMY<<");"<<std::endl;
   }

   // needed in case the frame was resized
   // otherwhice the frame became bigger showing all hidden widgets (layout algorithm)
   if (!usexy) out << "   " <<GetName()<< "->Resize("<< GetWidth()<<","<<GetHeight()<<");"<<std::endl;
   out << "}  " << std::endl;

   // writing slots
   TList *sl = (TList *)gROOT->GetListOfSpecials()->FindObject("ListOfSlots");
   if (sl) {
      TIter nextsl(sl);
      TObjString *slobj;
      Int_t pnumber = 1;

      while ((slobj = (TObjString*) nextsl())) {
         TString s = slobj->GetString();
         TString p = "";
         Int_t lb, rb, eq;
         lb = s.First('(');
         rb = s.First(')');
         eq = s.First('=');
         out << std::endl;

         if (rb - lb > 1 && eq == -1) {
            p = TString::Format(" par%d", pnumber);
            s.Insert(rb, p);
            pnumber++;
            out << "void " << s << std::endl;
            out << "{" << std::endl;
            s = slobj->GetString();
            s[rb] = ' ';
            out << "   std::cout << " << quote << "Slot " << s  << quote
                << " <<" << p << " << " << quote << ")" << quote
                << " << std::endl; " << std::endl;
            } else {
               if (eq != -1) {
                  s.Remove(eq, rb-eq);
                  out << "void " << s << std::endl;
                  out << "{" << std::endl;
                  out << "   std::cout << " << quote << "Slot " << s
                      << quote << " << std::endl; " << std::endl;
               } else {
                  out << "void " << slobj->GetString() << std::endl;
                  out << "{" << std::endl;
                  out << "   std::cout << " << quote << "Slot " << slobj->GetString()
                      << quote << " << std::endl; " << std::endl;
               }
            }
         out << "}" << std::endl;
      }
      gROOT->GetListOfSpecials()->Remove(sl);
      sl->Delete();
      delete sl;
   }
   out.close();

   if (!opt.Contains("quiet"))
      printf(" C++ macro file %s has been generated\n", gSystem->BaseName(ff.Data()));

   // reset bit TClass::kClassSaved for all classes
   nextc1.Reset();
   while((c1=(TClass*)nextc1())) {
      c1->ResetBit(TClass::kClassSaved);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a main frame widget as a C++ statement(s) on output stream out.

void TGMainFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fParent != gClient->GetDefaultRoot()) { // frame is embedded
      fOptions &= ~kMainFrame;
      TGCompositeFrame::SavePrimitive(out, option);
      fOptions |= kMainFrame;
      return;
   }

   char quote = '"';

   out << std::endl << "   // main frame" << std::endl;
   out << "   TGMainFrame *";
   out << GetName() << " = new TGMainFrame(gClient->GetRoot(),10,10,"   // layout alg.
       << GetOptionString() << ");" <<std::endl;
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager if it differs from the main frame type
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if ((GetOptions() & kHorizontalFrame) &&
       (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
      ;
   } else if ((GetOptions() & kVerticalFrame) &&
              (lm->InheritsFrom(TGVerticalLayout::Class()))) {
      ;
   } else {
      out << "   " << GetName() <<"->SetLayoutManager(";
      lm->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }

   SavePrimitiveSubframes(out, option);

   if (strlen(fWindowName)) {
      out << "   " << GetName() << "->SetWindowName(" << quote << GetWindowName()
          << quote << ");" << std::endl;
   }
   if (strlen(fIconName)) {
      out <<"   "<<GetName()<< "->SetIconName("<<quote<<GetIconName()<<quote<<");"<<std::endl;
   }
   if (strlen(fIconPixmap)) {
      out << "   " << GetName() << "->SetIconPixmap(" << quote << GetIconPixmap()
          << quote << ");" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a horizontal frame widget as a C++ statement(s) on output stream out.

void TGHorizontalFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // horizontal frame" << std::endl;
   out << "   TGHorizontalFrame *";
   out << GetName() << " = new TGHorizontalFrame(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out << ");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager if it differs from the main frame type
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if ((GetOptions() & kHorizontalFrame) &&
       (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
      ;
   } else if ((GetOptions() & kVerticalFrame) &&
              (lm->InheritsFrom(TGVerticalLayout::Class()))) {
      ;
   } else {
      out << "   " << GetName() <<"->SetLayoutManager(";
      lm->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }

   SavePrimitiveSubframes(out, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a vertical frame widget as a C++ statement(s) on output stream out.

void TGVerticalFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // vertical frame" << std::endl;
   out << "   TGVerticalFrame *";
   out << GetName() << " = new TGVerticalFrame(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager if it differs from the main frame type
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if ((GetOptions() & kHorizontalFrame) &&
       (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
      ;
   } else if ((GetOptions() & kVerticalFrame) &&
              (lm->InheritsFrom(TGVerticalLayout::Class()))) {
      ;
   } else {
      out << "   " << GetName() <<"->SetLayoutManager(";
      lm->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }

   SavePrimitiveSubframes(out, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a frame widget as a C++ statement(s) on output stream out.

void TGFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << "   TGFrame *";
   out << GetName() << " = new TGFrame("<< fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a group frame widget as a C++ statement(s) on output stream out.

void TGGroupFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   // font + GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC, parFont;
   // coverity[returned_null]
   // coverity[dereference]
   parFont.Form("%s::GetDefaultFontStruct()",IsA()->GetName());
   // coverity[returned_null]
   // coverity[dereference]
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

   out << std::endl << "   // " << quote << GetTitle() << quote << " group frame" << std::endl;
   out << "   TGGroupFrame *";
   out << GetName() <<" = new TGGroupFrame("<<fParent->GetName()
       << "," << quote << GetTitle() << quote;

   if (fBackground == GetDefaultFrameBackground()) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            if (GetOptions() & kVerticalFrame) {
               out <<");" << std::endl;
            } else {
               out << "," << GetOptionString() <<");" << std::endl;
            }
         } else {
            out << "," << GetOptionString() << "," << parGC.Data() <<");" << std::endl;
         }
      } else {
         out << "," << GetOptionString() << "," << parGC.Data() << "," << parFont.Data() << ");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << "," << parGC.Data() << "," << parFont.Data() << ",ucolor);"  << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (GetTitlePos() != -1)
      out << "   " << GetName() <<"->SetTitlePos(";
   if (GetTitlePos() == 0)
      out << "TGGroupFrame::kCenter);" << std::endl;
   if (GetTitlePos() == 1)
      out << "TGGroupFrame::kRight);" << std::endl;

   SavePrimitiveSubframes(out, option);

   // setting layout manager
   out << "   " << GetName() <<"->SetLayoutManager(";
   // coverity[returned_null]
   // coverity[dereference]
   GetLayoutManager()->SavePrimitive(out, option);
   out << ");"<< std::endl;

   out << "   " << GetName() <<"->Resize(" << GetWidth() << ","
       << GetHeight() << ");" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Save the GUI tranzient frame widget in a C++ macro file.

void TGTransientFrame::SaveSource(const char *filename, Option_t *option)
{
   // iterate over all active classes to exclude the base ones

   TString opt = option;
   TBits *bc = new TBits();
   TClass *c1, *c2, *c3;
   UInt_t k = 0;      // will mark k-bit of TBits if the class is a base class

   TIter nextc1(gROOT->GetListOfClasses());
   while((c1 = (TClass *)nextc1())) {

      //   resets bit TClass::kClassSaved for all classes
      c1->ResetBit(TClass::kClassSaved);

      TIter nextc2(gROOT->GetListOfClasses());
      while ((c2 = (TClass *)nextc2())) {
         if (c1==c2) continue;
         else {
            c3 = c2->GetBaseClass(c1);
            if (c3 != 0) {
               bc->SetBitNumber(k, kTRUE);
               break;
            }
         }
      }
      k++;
   }

   TList *ilist = new TList();   // will contain include file names without '.h'
   ilist->SetName("ListOfIncludes");
   gROOT->GetListOfSpecials()->Add(ilist);
   k=0;

   // completes list of include file names
   TIter nextdo(gROOT->GetListOfClasses());
   while ((c2 = (TClass *)nextdo())) {
      // to have only used GUI header files
      if (bc->TestBitNumber(k) == 0 && c2->InheritsFrom(TGObject::Class()) == 1) {
         // for any used ROOT header files activate the line below, comment the line above
         //if (bc->TestBitNumber(k) == 0) {
         const char *iname;
         iname = c2->GetDeclFileName();
         if (iname[0] && strstr(iname,".h")) {
            const char *lastsl = strrchr(iname,'/');
            if (lastsl) iname = lastsl + 1;
            char *tname = new char[strlen(iname)+1];
            Int_t i=0;
            while (*iname != '.') {
               tname[i] = *iname;
               i++; iname++;
            }
            tname[i] = 0;    //tname = include file name without '.h'

            TObjString *iel = (TObjString *)ilist->FindObject(tname);
            if (!iel) {
               ilist->Add(new TObjString(tname));
            }
            delete [] tname;
         }
         k++;  continue;
      }
      k++;
   }

   char quote = '"';
   std::ofstream out;

   TString ff = filename && strlen(filename) ? filename : "Rootdlog.C";

   // Computes the main method name.
   const char *fname = gSystem->BaseName(ff.Data());
   Int_t lenfname = strlen(fname);
   char *sname = new char[lenfname+1];

   Int_t i = 0;
   while ((*fname != '.') && (i < lenfname)) {
      sname[i] = *fname;
      i++; fname++;
   }
   if (i == lenfname)
      ff += ".C";
   sname[i] = 0;

   out.open(ff.Data(), std::ios::out);
   if (!out.good()) {
      Error("SaveSource", "cannot open file: %s", ff.Data());
      delete [] sname;
      return;
   }

   // writes include files in C++ macro
   TObjString *inc;
   ilist = (TList *)gROOT->GetListOfSpecials()->FindObject("ListOfIncludes");

   if (!ilist) {
      delete [] sname;
      return;
   }

   // write macro header, date/time stamp as string, and the used Root version
   TDatime t;
   out <<"// Dialog macro generated from application: "<< gApplication->Argv(0) << std::endl;
   out <<"// By ROOT version "<< gROOT->GetVersion() <<" on "<<t.AsSQLString()<< std::endl;
   out << std::endl;

   out << "#if !defined( __CINT__) || defined (__MAKECINT__)" << std::endl << std::endl;

   TIter nexti(ilist);
   while((inc = (TObjString *)nexti())) {
      out <<"#ifndef ROOT_"<< inc->GetString() << std::endl;
      out <<"#include "<< quote << inc->GetString() <<".h"<< quote << std::endl;
      out <<"#endif" << std::endl;
      if (strstr(inc->GetString(),"TRootEmbeddedCanvas")) {
         out <<"#ifndef ROOT_TCanvas"<< std::endl;
         out <<"#include "<< quote <<"TCanvas.h"<< quote << std::endl;
         out <<"#endif" << std::endl;
      }
   }
   out << std::endl << "#include " << quote << "Riostream.h" << quote << std::endl;
   out << std::endl << "#endif" << std::endl;
   // deletes created ListOfIncludes
   gROOT->GetListOfSpecials()->Remove(ilist);
   ilist->Delete();
   delete ilist;
   delete bc;

   // writes the macro entry point equal to the fname
   out << std::endl;
   out << "void " << sname << "()" << std::endl;
   delete [] sname;

   //  Save GUI widgets as a C++ macro in a file
   out <<"{"<< std::endl;

   gListOfHiddenFrames->Clear();

   // saivng slots
   TList *lSlots = new TList;
   lSlots->SetName("ListOfSlots");
   gROOT->GetListOfSpecials()->Add(lSlots);

   TGTransientFrame::SavePrimitive(out, option);

   if (strlen(fClassName) || strlen(fResourceName)) {
      out<<"   "<<GetName()<< "->SetClassHints("<<quote<<fClassName<<quote
                                            <<"," <<quote<<fResourceName<<quote
                                            <<");"<<std::endl;
   }

   GetMWMHints(fMWMValue, fMWMFuncs, fMWMInput);
   if (fMWMValue || fMWMFuncs || fMWMInput) {
      out << "   " << GetName() << "->SetMWMHints(";
      out << GetMWMvalueString() << "," << std::endl;
      out << "                        ";
      out << GetMWMfuncString() << "," << std::endl;
      out << "                        ";
      out << GetMWMinpString() << ");"<< std::endl;
   }

   GetWMPosition(fWMX, fWMY);
   if ((fWMX != -1) || (fWMY != -1)) {
      out <<"   "<<GetName()<<"->SetWMPosition("<<fWMX<<","<<fWMY<<");"<<std::endl;
   }

   GetWMSize(fWMWidth, fWMHeight);
   if (fWMWidth != UInt_t(-1) || fWMHeight != UInt_t(-1)) {
      out <<"   "<<GetName()<<"->SetWMSize("<<fWMWidth<<","<<fWMHeight<<");"<<std::endl;
   }

   GetWMSizeHints(fWMMinWidth,fWMMinHeight,fWMMaxWidth,fWMMaxHeight,fWMWidthInc,fWMHeightInc);
   if (fWMMinWidth != UInt_t(-1) || fWMMinHeight != UInt_t(-1) ||
       fWMMaxWidth != UInt_t(-1) || fWMMaxHeight != UInt_t(-1) ||
       fWMWidthInc != UInt_t(-1) || fWMHeightInc != UInt_t(-1)) {

      out <<"   "<<GetName()<<"->SetWMSizeHints("<<fWMMinWidth<<","<<fWMMinHeight
          <<","<<fWMMaxWidth<<","<<fWMMaxHeight <<","<<fWMWidthInc<<","<<fWMHeightInc
          <<");"<<std::endl;
   }

   GetWMPosition(fWMX, fWMY);
   if ((fWMX != -1) || (fWMY != -1)) {
      out <<"   "<<GetName()<<"->Move("<<fWMX<<","<<fWMY<<");"<<std::endl;
   }

   out << "   " <<GetName()<< "->MapSubwindows();" << std::endl;

   TIter nexth(gListOfHiddenFrames);
   TGFrame *fhidden;
   while ((fhidden = (TGFrame*)nexth())) {
      out << "   " <<fhidden->GetName()<< "->UnmapWindow();" << std::endl;
   }
   out << std::endl;
   gListOfHiddenFrames->Clear();

   Bool_t usexy = kFALSE;
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if (lm->InheritsFrom("TGXYLayout"))
      usexy = kTRUE;

   if (!usexy)
      out << "   " <<GetName()<< "->Resize("<< GetName()<< "->GetDefaultSize());" << std::endl;
   else
      out << "   " <<GetName()<< "->Resize("<< GetWidth()<<","<<GetHeight()<<");"<<std::endl;

   out << "   " <<GetName()<< "->MapWindow();" <<std::endl;
   if (!usexy) out << "   " <<GetName()<< "->Resize();" << std::endl;
   out << "}  " << std::endl;

   // writing slots
   TList *sl = (TList *)gROOT->GetListOfSpecials()->FindObject("ListOfSlots");
   if (sl) {
      TIter nextsl(sl);
      TObjString *slobj;
      Int_t pnumber = 1;

      while ((slobj = (TObjString*) nextsl())) {
         TString s = slobj->GetString();
         TString p = "";
         Int_t lb, rb, eq;
         lb = s.First('(');
         rb = s.First(')');
         eq = s.First('=');
         out << std::endl;

         if (rb - lb > 1 && eq == -1) {
            p = TString::Format(" par%d", pnumber);
            s.Insert(rb, p);
            pnumber++;
            out << "void " << s << std::endl;
            out << "{" << std::endl;
            s = slobj->GetString();
            s[rb] = ' ';
            out << "   std::cout << " << quote << "Slot " << s  << quote
                << " <<" << p << " << " << quote << ")" << quote
                << " << std::endl; " << std::endl;
            } else {
               if (eq != -1) {
                  s.Remove(eq, rb-eq);
                  out << "void " << s << std::endl;
                  out << "{" << std::endl;
                  out << "   std::cout << " << quote << "Slot " << s
                      << quote << " << std::endl; " << std::endl;
               } else {
                  out << "void " << slobj->GetString() << std::endl;
                  out << "{" << std::endl;
                  out << "   std::cout << " << quote << "Slot " << slobj->GetString()
                      << quote << " << std::endl; " << std::endl;
               }
            }
         out << "}" << std::endl;
      }
      gROOT->GetListOfSpecials()->Remove(sl);
      sl->Delete();
      delete sl;
   }

   out.close();

   if (!opt.Contains("quiet"))
      printf(" C++ macro file %s has been generated\n", gSystem->BaseName(ff.Data()));

   // reset bit TClass::kClassSaved for all classes
   nextc1.Reset();
   while((c1=(TClass*)nextc1())) {
      c1->ResetBit(TClass::kClassSaved);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a transient frame widget as a C++ statement(s) on output stream out.

void TGTransientFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   out << std::endl << "   // transient frame" << std::endl;
   out << "   TGTransientFrame *";
   out << GetName()<<" = new TGTransientFrame(gClient->GetRoot(),0"
       << "," << GetWidth() << "," << GetHeight() << "," << GetOptionString() <<");" << std::endl;

   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager if it differs from transient frame type
   // coverity[returned_null]
   // coverity[dereference]
   TGLayoutManager * lm = GetLayoutManager();
   if ((GetOptions() & kHorizontalFrame) &&
       (lm->InheritsFrom(TGHorizontalLayout::Class()))) {
      ;
   } else if ((GetOptions() & kVerticalFrame) &&
              (lm->InheritsFrom(TGVerticalLayout::Class()))) {
      ;
   } else {
      out << "   " << GetName() <<"->SetLayoutManager(";
      lm->SavePrimitive(out, option);
      out << ");"<< std::endl;
   }

   SavePrimitiveSubframes(out, option);

   if (strlen(fWindowName)) {
      out << "   " << GetName() << "->SetWindowName(" << quote << GetWindowName()
          << quote << ");" << std::endl;
   }
   if (strlen(fIconName)) {
      out <<"   "<<GetName()<< "->SetIconName("<<quote<<GetIconName()<<quote<<");"<<std::endl;
   }
   if (strlen(fIconPixmap)) {
      out << "   " << GetName() << "->SetIconPixmap(" << quote << GetIconPixmap()
          << quote << ");" << std::endl;
   }
}
