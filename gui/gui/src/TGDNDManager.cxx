// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   19/04/07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGFrame.h"
#include "TTimer.h"
#include "TGDNDManager.h"
#include "TRootCanvas.h"
#include "TVirtualX.h"


#define ROOTDND_PROTOCOL_VERSION      4
#define XA_ATOM ((Atom_t) 4)
#define XA_WINDOW ((Atom_t) 33)

Atom_t TGDNDManager::fgDNDAware         = kNone;
Atom_t TGDNDManager::fgDNDSelection     = kNone;
Atom_t TGDNDManager::fgDNDProxy         = kNone;

Atom_t TGDNDManager::fgDNDEnter         = kNone;
Atom_t TGDNDManager::fgDNDLeave         = kNone;
Atom_t TGDNDManager::fgDNDPosition      = kNone;
Atom_t TGDNDManager::fgDNDStatus        = kNone;
Atom_t TGDNDManager::fgDNDDrop          = kNone;
Atom_t TGDNDManager::fgDNDFinished      = kNone;
Atom_t TGDNDManager::fgDNDVersion       = kNone;

Atom_t TGDNDManager::fgDNDActionCopy    = kNone;
Atom_t TGDNDManager::fgDNDActionMove    = kNone;
Atom_t TGDNDManager::fgDNDActionLink    = kNone;
Atom_t TGDNDManager::fgDNDActionAsk     = kNone;
Atom_t TGDNDManager::fgDNDActionPrivate = kNone;

Atom_t TGDNDManager::fgDNDTypeList      = kNone;
Atom_t TGDNDManager::fgDNDActionList    = kNone;
Atom_t TGDNDManager::fgDNDActionDescrip = kNone;

Atom_t TGDNDManager::fgXAWMState     = kNone;
Atom_t TGDNDManager::fgXCDNDData     = kNone;

Bool_t TGDNDManager::fgInit = kFALSE;

// TODO:
// - add an TGFrame::HandleDNDStatus event handler?
// - implement INCR protocol
// - cache several requests?

TGDNDManager *gDNDManager = 0;

Cursor_t TGDragWindow::fgDefaultCursor = kNone;

//_____________________________________________________________________________
//
// TGDragWindow
//
// Window used as drag icon during drag and drop operations.
//_____________________________________________________________________________

ClassImp(TGDragWindow);

////////////////////////////////////////////////////////////////////////////////
/// TGDragWindow constructor.

TGDragWindow::TGDragWindow(const TGWindow *p, Pixmap_t pic, Pixmap_t mask,
                           UInt_t options, Pixel_t back) :
   TGFrame(p, 32, 32, options, back)
{
   if (fgDefaultCursor == kNone) {
      fgDefaultCursor = gVirtualX->CreateCursor(kTopLeft);
   }

   fPic = pic;
   fMask = mask;

   SetWindowAttributes_t wattr;

   wattr.fMask = kWAOverrideRedirect | kWASaveUnder;
   wattr.fSaveUnder = kTRUE;
   wattr.fOverrideRedirect = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &wattr);

   int x, y;

   gVirtualX->GetWindowSize(fPic, x, y, fPw, fPh);

   wattr.fMask = kWAOverrideRedirect;
   wattr.fOverrideRedirect = kTRUE;

   // This input window is used to make the dragging smoother when using
   // highly complicated shapped windows (like labels and semitransparent
   // icons), for some obscure reason most of the motion events get lost
   // while the pointer is over the shaped window.

   //fInput = gVirtualX->CreateWindow(fParent->GetId(), 0, 0, fWidth,
   //                                 fHeight, 0, 0, 0, 0, &wattr, 0);
   fInput = fId;

   Resize(GetDefaultSize());

   gVirtualX->ShapeCombineMask(fId, 0, 0, fMask);

   gVirtualX->SetCursor(fId, fgDefaultCursor);
}

////////////////////////////////////////////////////////////////////////////////
/// TGDragWindow destructor.

TGDragWindow::~TGDragWindow()
{
   //gVirtualX->DestroyWindow(fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Map TGDragWindow.

void TGDragWindow::MapWindow()
{
   TGFrame::MapWindow();
   //gVirtualX->MapWindow(fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Unmap TGDragWindow.

void TGDragWindow::UnmapWindow()
{
   TGFrame::UnmapWindow();
   //gVirtualX->UnmapWindow(fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Raise TGDragWindow.

void TGDragWindow::RaiseWindow()
{
   TGFrame::RaiseWindow();
   //gVirtualX->RaiseWindow(fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Lower TGDragWindow.

void TGDragWindow::LowerWindow()
{
   //gVirtualX->LowerWindow(fInput);
   TGFrame::LowerWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Map and Raise TGDragWindow.

void TGDragWindow::MapRaised()
{
   TGFrame::MapRaised();
   //gVirtualX->MapRaised(fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Layout TGDragWindow.

void TGDragWindow::Layout()
{
   gVirtualX->ShapeCombineMask(fId, 0, 0, fMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw TGDragWindow.

void TGDragWindow::DoRedraw()
{
   gVirtualX->CopyArea(fPic, fId, GetBckgndGC()(), 0, 0, fWidth, fHeight, 0, 0);
}

//_____________________________________________________________________________
//
// TGDNDManager
//
// Central Drag and Drop manager for ROOT.
//_____________________________________________________________________________

ClassImp(TGDNDManager);

////////////////////////////////////////////////////////////////////////////////
/// TGDNDManager constructor.

TGDNDManager::TGDNDManager(TGFrame *toplevel, Atom_t * /*typelist*/)
{
   if (gDNDManager)
      // coverity[uninit_member]: already done
      return;

   fMain = toplevel;
   fVersion = ROOTDND_PROTOCOL_VERSION;
   fUseVersion = kTRUE;
   //fTypelist = typelist;
   fTypelist = new Atom_t[3];
   fTypelist[0] = gVirtualX->InternAtom("application/root", kFALSE);
   fTypelist[1] = gVirtualX->InternAtom("text/uri-list", kFALSE);
   fTypelist[2] = 0;

   if (!fgInit) {
      InitAtoms();
      fgInit = kTRUE;
   }

   //Reset();
   fDropTimeout = 0;

   fSource = kNone;
   fTarget = kNone;
   fTargetIsDNDAware = kFALSE;
   fStatusPending = kFALSE;
   fDropAccepted = kFALSE;  // this would become obsoleted by _acceptedAction
   fAcceptedAction = kNone; // target's accepted action
   fLocalAction = kNone;    // our last specified action when we act as source
   fDragging = kFALSE;
   fDragWin = 0;
   fLocalSource = 0;
   fLocalTarget = 0;
   fPic = fMask = kNone;
   fDraggerTypes = 0;
   fDropType = kNone;
   fHotx = fHoty = 0;

   fGrabEventMask = kButtonPressMask | kButtonReleaseMask | kButtonMotionMask;

   fDNDNoDropCursor = gVirtualX->CreateCursor(kNoDrop); // kNoDrop

   // set the aware prop

   fProxyOurs = kFALSE;
   gDNDManager = this;
}

////////////////////////////////////////////////////////////////////////////////
/// TGDNDManager destructor.

TGDNDManager::~TGDNDManager()
{
   // remove the proxy prop if we own it
   if (fProxyOurs)
      RemoveRootProxy();

   // remove the aware prop ant the types list, if any
   if (fMain) {
      gVirtualX->DeleteProperty(fMain->GetId(), fgDNDAware);
      gVirtualX->DeleteProperty(fMain->GetId(), fgDNDTypeList);
   }
   if (fDropTimeout) delete fDropTimeout;

   // delete the drag pixmap, if any
   if (fDragWin) {
      fDragWin->DeleteWindow();
      fDragWin = 0;
   }
   if (fPic != kNone) gVirtualX->DeletePixmap(fPic);
   if (fMask != kNone) gVirtualX->DeletePixmap(fMask);

   if (fDraggerTypes) delete[] fDraggerTypes;
   if (fTypelist) delete[] fTypelist;
}

Atom_t TGDNDManager::GetDNDAware() { return fgDNDAware; }
Atom_t TGDNDManager::GetDNDSelection() { return fgDNDSelection; }
Atom_t TGDNDManager::GetDNDProxy() { return fgDNDProxy; }
Atom_t TGDNDManager::GetDNDEnter() { return fgDNDEnter; }
Atom_t TGDNDManager::GetDNDLeave() { return fgDNDLeave; }
Atom_t TGDNDManager::GetDNDPosition() { return fgDNDPosition; }
Atom_t TGDNDManager::GetDNDStatus() { return fgDNDStatus; }
Atom_t TGDNDManager::GetDNDDrop() { return fgDNDDrop; }
Atom_t TGDNDManager::GetDNDFinished() { return fgDNDFinished; }
Atom_t TGDNDManager::GetDNDVersion() { return fgDNDVersion; }
Atom_t TGDNDManager::GetDNDActionCopy() { return fgDNDActionCopy; }
Atom_t TGDNDManager::GetDNDActionMove() { return fgDNDActionMove; }
Atom_t TGDNDManager::GetDNDActionLink() { return fgDNDActionLink; }
Atom_t TGDNDManager::GetDNDActionAsk() { return fgDNDActionAsk; }
Atom_t TGDNDManager::GetDNDActionPrivate() { return fgDNDActionPrivate; }
Atom_t TGDNDManager::GetDNDTypeList() { return fgDNDTypeList; }
Atom_t TGDNDManager::GetDNDActionList() { return fgDNDActionList; }
Atom_t TGDNDManager::GetDNDActionDescrip() { return fgDNDActionDescrip; }
Atom_t TGDNDManager::GetXCDNDData() { return fgXCDNDData; }

////////////////////////////////////////////////////////////////////////////////
/// Initialize drag and drop atoms.

void TGDNDManager::InitAtoms()
{
   // awareness
   fgDNDAware = gVirtualX->InternAtom("XdndAware", kFALSE);

   // selection
   fgDNDSelection = gVirtualX->InternAtom("XdndSelection", kFALSE);

   // proxy window
   fgDNDProxy = gVirtualX->InternAtom("XdndProxy", kFALSE);

   // messages
   fgDNDEnter    = gVirtualX->InternAtom("XdndEnter", kFALSE);
   fgDNDLeave    = gVirtualX->InternAtom("XdndLeave", kFALSE);
   fgDNDPosition = gVirtualX->InternAtom("XdndPosition", kFALSE);
   fgDNDStatus   = gVirtualX->InternAtom("XdndStatus", kFALSE);
   fgDNDDrop     = gVirtualX->InternAtom("XdndDrop", kFALSE);
   fgDNDFinished = gVirtualX->InternAtom("XdndFinished", kFALSE);

   // actions
   fgDNDActionCopy    = gVirtualX->InternAtom("XdndActionCopy", kFALSE);
   fgDNDActionMove    = gVirtualX->InternAtom("XdndActionMove", kFALSE);
   fgDNDActionLink    = gVirtualX->InternAtom("XdndActionLink", kFALSE);
   fgDNDActionAsk     = gVirtualX->InternAtom("XdndActionAsk", kFALSE);
   fgDNDActionPrivate = gVirtualX->InternAtom("XdndActionPrivate", kFALSE);

   // types list
   fgDNDTypeList      = gVirtualX->InternAtom("XdndTypeList", kFALSE);
   fgDNDActionList    = gVirtualX->InternAtom("XdndActionList", kFALSE);
   fgDNDActionDescrip = gVirtualX->InternAtom("XdndActionDescription", kFALSE);

   // misc
   fgXAWMState = gVirtualX->InternAtom("WM_STATE", kFALSE);
   fgXCDNDData = gVirtualX->InternAtom("_XC_DND_DATA", kFALSE);
}

static int ArrayLength(Atom_t *a)
{
   // Returns length of array a.

   int n;

   for (n = 0; a[n]; n++) { }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if window win is DND aware.

Bool_t TGDNDManager::IsDNDAware(Window_t win, Atom_t *typelist)
{
   return gVirtualX->IsDNDAware(win, typelist);
}

////////////////////////////////////////////////////////////////////////////////
/// Search for DND aware window at position x,y.

Window_t TGDNDManager::FindWindow(Window_t root, int x, int y, int maxd)
{
   if (maxd <= 0) return kNone;

   if (fDragWin && fDragWin->HasWindow(root)) return kNone;

   return gVirtualX->FindRWindow(root, fDragWin ? fDragWin->GetId() : 0,
                                 fDragWin ? fDragWin->GetInputId() : 0,
                                 x, y, maxd);
}


////////////////////////////////////////////////////////////////////////////////
/// Get root window proxy.

Window_t TGDNDManager::GetRootProxy()
{
   Atom_t actual;
   Int_t format = 32;
   ULong_t count, remaining;
   unsigned char *data = 0;
   Window_t win, proxy = kNone;

   // search for XdndProxy property on the root window...

   // XSync(_dpy, kFALSE);      // get to known state...
   gVirtualX->UpdateWindow(0);

   //oldhandler = XSetErrorHandler(TGDNDManager::CatchXError);
   //target_error = kFALSE;

   gVirtualX->GetProperty(gVirtualX->GetDefaultRootWindow(),
                          fgDNDProxy, 0, 1, kFALSE, XA_WINDOW,
                          &actual, &format, &count, &remaining, &data);

   if ((actual == XA_WINDOW) && (format == 32) && (count > 0) && data) {

      // found the XdndProxy property, now check for the proxy window...
      win = *((Window_t *) data);
      delete[] data;
      data = 0;

      gVirtualX->GetProperty(win, fgDNDProxy, 0, 1, kFALSE, XA_WINDOW,
                             &actual, &format, &count, &remaining, &data);

      // XSync(_dpy, kFALSE);      // force the error...
      gVirtualX->UpdateWindow(0);

      if ((actual == XA_WINDOW) && (format == 32) && (count > 0) && data) {
         if (*((Window_t *) data) == win) {

            // proxy window exists and is correct
            proxy = win;
         }
      }
   }
   if (data) delete[] data;
   //oldhandler = XSetErrorHandler(oldhandler);
   return proxy;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND related client messages.

Bool_t TGDNDManager::HandleClientMessage(Event_t *event)
{
   if (event->fHandle == fgDNDEnter) {
      HandleDNDEnter((Window_t) event->fUser[0], event->fUser[1],
                     (Atom_t *) &event->fUser[2]);

   } else if (event->fHandle == fgDNDLeave) {
      HandleDNDLeave((Window_t) event->fUser[0]);

   } else if (event->fHandle == fgDNDPosition) {
      HandleDNDPosition((Window_t) event->fUser[0],
                       (Int_t) (event->fUser[2] >> 16) & 0xFFFF,  // x_root
                       (Int_t) (event->fUser[2] & 0xFFFF),        // y_root
                       (Atom_t) event->fUser[4],                  // action
                       (Time_t) event->fUser[3]);                 // timestamp

   } else if (event->fHandle == fgDNDStatus) {
      Rectangle_t skip;
      skip.fX      = (event->fUser[2] >> 16) & 0xFFFF;
      skip.fY      = (event->fUser[2] & 0xFFFF);
      skip.fWidth  = (event->fUser[3] >> 16) & 0xFFFF;
      skip.fHeight = (event->fUser[3] & 0xFFFF);

      HandleDNDStatus((Window_t) event->fUser[0],
                      (int) (event->fUser[1] & 0x1),
                       skip, (Atom_t) event->fUser[4]);

   } else if (event->fHandle == fgDNDDrop) {
      HandleDNDDrop((Window_t) event->fUser[0], (Time_t) event->fUser[2]);

   } else if (event->fHandle == fgDNDFinished) {
      HandleDNDFinished((Window_t) event->fUser[0]);

   } else {
      return kFALSE;  // not for us...
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Drop timeout.

Bool_t TGDNDManager::HandleTimer(TTimer *t)
{
   if (t == fDropTimeout) {
      // The drop operation timed out without receiving
      // status confirmation from the target. Send a
      // leave message instead (and notify the user or widget).
      delete fDropTimeout;
      fDropTimeout = 0;

      SendDNDLeave(fTarget);
      fStatusPending = kFALSE;

      if (fLocalSource) fLocalSource->HandleDNDFinished();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND enter message to target window.

void TGDNDManager::SendDNDEnter(Window_t target)
{
   Int_t i, n;
   Event_t event;

   event.fType   = kClientMessage;
   event.fWindow = target;
   event.fHandle = fgDNDEnter;
   event.fFormat = 32;

   event.fUser[0] = fMain->GetId();  // from;

   n = ArrayLength(fTypelist);

   event.fUser[1] = ((n > 3) ? 1L : 0L) | (fUseVersion << 24);

   // set the first 1-3 data types

   for (i = 0; i < 3; ++i)
      event.fUser[2+i] = (i < n) ? fTypelist[i] : kNone;

   if (fLocalSource) {
      TDNDData *dnddata = 0;
      Atom_t dataType;

      // get the data type from the drag source widget
      if (fLocalSource)
         dnddata = fLocalSource->GetDNDData(0);
      dataType = dnddata ? (Atom_t) dnddata->fDataType : (Atom_t) kNone;
      event.fUser[2] = dataType;
      event.fUser[3] = kNone;
      event.fUser[4] = kNone;
   }

   gVirtualX->SendEvent(target, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND leave message to target window.

void TGDNDManager::SendDNDLeave(Window_t target)
{
   Event_t event;

   event.fType    = kClientMessage;
   event.fWindow  = target;
   event.fHandle  = fgDNDLeave;
   event.fFormat  = 32;

   event.fUser[0] = fMain->GetId();  // from;
   event.fUser[1] = 0L;

   event.fUser[2] = 0L;
   event.fUser[3] = 0L;
   event.fUser[4] = 0L;

   gVirtualX->SendEvent(target, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND position message to target window.

void TGDNDManager::SendDNDPosition(Window_t target, int x, int y,
                                  Atom_t action, Time_t timestamp)
{
   Event_t event;

   event.fType    = kClientMessage;
   event.fWindow  = target;
   event.fHandle  = fgDNDPosition;
   event.fFormat  = 32;

   event.fUser[0] = fMain->GetId();  // from;
   event.fUser[1] = 0L;

   event.fUser[2] = (x << 16) | y;   // root coodinates
   event.fUser[3] = timestamp;       // timestamp for retrieving data
   event.fUser[4] = action;          // requested action

   gVirtualX->SendEvent(target, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND status message to source window.

void TGDNDManager::SendDNDStatus(Window_t source, Atom_t action)
{
   Event_t event;

   event.fType    = kClientMessage;
   event.fWindow  = source;
   event.fHandle  = fgDNDStatus;
   event.fFormat  = 32;

   event.fUser[0] = fMain->GetId();    // from;
   event.fUser[1] = (action == kNone) ? 0L : 1L;

   event.fUser[2] = 0L;                // empty rectangle
   event.fUser[3] = 0L;
   event.fUser[4] = action;            // accepted action

   gVirtualX->SendEvent(source, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND drop message to target window.

void TGDNDManager::SendDNDDrop(Window_t target)
{
   Event_t event;

   event.fType    = kClientMessage;
   event.fWindow  = target;
   event.fHandle  = fgDNDDrop;
   event.fFormat  = 32;

   event.fUser[0] = fMain->GetId();    // from;
   event.fUser[1] = 0L;                // reserved
   event.fUser[2] = 0L; //CurrentTime;       // timestamp
   event.fUser[3] = 0L;
   event.fUser[4] = 0L;

   gVirtualX->SendEvent(target, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Send DND finished message to source window.

void TGDNDManager::SendDNDFinished(Window_t source)
{
   Event_t event;

   event.fType    = kClientMessage;
   event.fWindow  = source;
   event.fHandle  = fgDNDFinished;
   event.fFormat  = 32;

   event.fUser[0] = fMain->GetId();    // from;
   event.fUser[1] = 0L;                // reserved
   event.fUser[2] = 0L;
   event.fUser[3] = 0L;
   event.fUser[4] = 0L;

   gVirtualX->SendEvent(source, &event);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND enter event.

Bool_t TGDNDManager::HandleDNDEnter(Window_t src, Long_t vers, Atom_t dataTypes[3])
{
   fSource = src;

   if (fDraggerTypes) delete[] fDraggerTypes;

   if (vers & 1) {  // more than 3 data types?
      Atom_t type, *a;
      Int_t format = 32;
      ULong_t i, count, remaining;
      unsigned char *data = 0;

      gVirtualX->GetProperty(src, fgDNDTypeList,
                             0, 0x8000000L, kFALSE, XA_ATOM,
                             &type, &format, &count, &remaining, &data);

      if (type != XA_ATOM || format != 32 || !data) {
         count = 0;
      }

      fDraggerTypes = new Atom_t[count+4];

      a = (Atom_t *) data;
      for (i = 0; i < count; i++)
         fDraggerTypes[i] = a[i];

      fDraggerTypes[i] = kNone;

      if (data) delete[] data;

   } else {
      fDraggerTypes = new Atom_t[4];

      fDraggerTypes[0] = dataTypes[0];
      fDraggerTypes[1] = dataTypes[1];
      fDraggerTypes[2] = dataTypes[2];

      fDraggerTypes[3] = kNone;
   }

   // the following is not strictly neccessary, unless the previous
   // dragging application crashed without sending XdndLeave
   if (fLocalTarget) fLocalTarget->HandleDNDLeave();
   fLocalTarget = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND leave event.

Bool_t TGDNDManager::HandleDNDLeave(Window_t /*src*/)
{
   fSource = kNone;
   if (fLocalTarget) fLocalTarget->HandleDNDLeave();
   fLocalTarget = 0;

   if (fDraggerTypes) delete[] fDraggerTypes;
   fDraggerTypes = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND position event.

Bool_t TGDNDManager::HandleDNDPosition(Window_t source, Int_t x_root, Int_t y_root,
                                      Atom_t action, Time_t /*timestamp*/)
{
   Int_t x = 0, y = 0;
   Window_t child;
   TGFrame *f = 0, *main = 0;
   TGWindow *w = 0;
   Window_t wtarget = 0;

   wtarget = FindWindow(gVirtualX->GetDefaultRootWindow(), x_root, y_root, 15);

   if (wtarget) {
      gVirtualX->TranslateCoordinates(gVirtualX->GetDefaultRootWindow(),
                                      wtarget, x_root, y_root, x, y, child);
      w = gClient->GetWindowById(wtarget);
      if (w)
         f = dynamic_cast<TGFrame *>(w);
   }

   if (f != fLocalTarget) {
      if (fLocalTarget) fLocalTarget->HandleDNDLeave();
      fLocalTarget = f;
      if (fLocalTarget) {
         main = (TGFrame *)fLocalTarget->GetMainFrame();
         main->RaiseWindow();
         if (fMain == 0)
            fMain = main;
         fDropType = fLocalTarget->HandleDNDEnter(fDraggerTypes);
      }
   }
   // query the target widget to determine whether it accepts the
   // required action
   if (fLocalTarget) {
      action = (fDropType == kNone) ? kNone :
              fLocalTarget->HandleDNDPosition(x, y, action, x_root, y_root);
   } else if (fProxyOurs) {
      action = fMain->HandleDNDPosition(x, y, action, x_root, y_root);
   } else {
      action = kNone;
   }
   SendDNDStatus(source, fLocalAction = action);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND status event.

Bool_t TGDNDManager::HandleDNDStatus(Window_t target, Int_t accepted,
                                    Rectangle_t /*area*/, Atom_t action)
{
   if (target) {
      fStatusPending = kFALSE;
      if (accepted) {
         fDropAccepted = kTRUE;
         fAcceptedAction = action;
         if (fDragWin)
            gVirtualX->ChangeActivePointerGrab(fDragWin->GetId(),
                                               fGrabEventMask, kNone);
      } else {
         fDropAccepted = kFALSE;
         fAcceptedAction = kNone;
         if (fDragWin)
            gVirtualX->ChangeActivePointerGrab(fDragWin->GetId(),
                                               fGrabEventMask,
                                               fDNDNoDropCursor);
      }
      if (fDropTimeout) {   // were we waiting for this to do the drop?
         delete fDropTimeout;
         fDropTimeout = 0;
         SendDNDDrop(fTarget);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND drop event.

Bool_t TGDNDManager::HandleDNDDrop(Window_t source, Time_t timestamp)
{
   // to get the data, we must call XConvertSelection with
   // the timestamp in XdndDrop, wait for SelectionNotify
   // to arrive to retrieve the data, and when we are finished,
   // send a XdndFinished message to the source.

   if (fMain && fDropType != kNone) {
      gVirtualX->ChangeProperties(fMain->GetId(), fgXCDNDData, fDropType,
                                  8, (unsigned char *) 0, 0);

      gVirtualX->ConvertSelection(fMain->GetId(), fgDNDSelection, fDropType,
                                  fgXCDNDData, timestamp);
   }

   fSource = source;
   if (fMain) SendDNDFinished(source);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle DND finished event.

Bool_t TGDNDManager::HandleDNDFinished(Window_t /*target*/)
{
   if (fLocalSource) fLocalSource->HandleDNDFinished();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle selection request event.

Bool_t TGDNDManager::HandleSelectionRequest(Event_t *event)
{
   if ((Atom_t)event->fUser[1] == fgDNDSelection) {
      Event_t xevent;
      TDNDData *dnddata = 0;
      char *data;
      int len;

      // get the data from the drag source widget
      if (fLocalSource)
         dnddata = fLocalSource->GetDNDData(event->fUser[2]);

      data = dnddata ? (char *) dnddata->fData : (char *) "";
      len  = dnddata ? dnddata->fDataLength : 0;

      if ((Atom_t)event->fUser[3] == kNone) {
         //printf("warning: kNone property specified in SelectionRequest\n");
         event->fUser[3] = fgXCDNDData;
      }

      gVirtualX->ChangeProperties(event->fUser[0], event->fUser[3],
                                  event->fUser[2], 8,
                                  (unsigned char *) data, len);

      xevent.fType    = kSelectionNotify;
      xevent.fTime    = event->fTime;
      xevent.fUser[0] = event->fUser[0]; // requestor
      xevent.fUser[1] = event->fUser[1]; // selection
      xevent.fUser[2] = event->fUser[2]; // target;
      xevent.fUser[3] = event->fUser[3]; // property;
      gVirtualX->SendEvent(event->fUser[0], &xevent);

      return kTRUE;
   } else {
      return kFALSE;  // not for us...
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle selection event.

Bool_t TGDNDManager::HandleSelection(Event_t *event)
{
   if ((Atom_t)event->fUser[1] == fgDNDSelection) {
      Atom_t actual = fDropType;
      Int_t format = 8;
      ULong_t count, remaining;
      unsigned char *data = 0;


      gVirtualX->GetProperty(event->fUser[0], event->fUser[3],
                             0, 0x8000000L, kTRUE, event->fUser[2],
                             &actual, &format, &count, &remaining, &data);

      if ((actual != fDropType) || (format != 8) || (count == 0) || !data) {
         if (data) delete[] data;

         return kFALSE;
      }

      if (fSource != kNone) SendDNDFinished(fSource);

      // send the data to the target widget

      if (fLocalTarget) {
         TDNDData dndData(actual, data, count, fLocalAction);
         fLocalTarget->HandleDNDDrop(&dndData);
         if (fDraggerTypes) delete[] fDraggerTypes;
         fDraggerTypes = 0;
      }

      fSource = kNone;
      fLocalAction = kNone;

//      delete[] data;

      return kTRUE;

   } else {
      return kFALSE;  // not for us...
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set drag window pixmaps and hotpoint.

void TGDNDManager::SetDragPixmap(Pixmap_t pic, Pixmap_t mask,
                                int hot_x, int hot_y)
{
   fPic  = pic;
   fMask = mask;
   fHotx = hot_x;
   fHoty = hot_y;
}

////////////////////////////////////////////////////////////////////////////////
/// Start dragging.

Bool_t TGDNDManager::StartDrag(TGFrame *src, int x_root, int y_root,
                              Window_t grabWin)
{
   if (fDragging) return kTRUE;

   fLocalSource = src;

   if ((TGWindow *)fMain != src->GetMainFrame()) {
      fMain = (TGFrame *)src->GetMainFrame();
   }

   if (!gVirtualX->SetSelectionOwner(fMain->GetId(), fgDNDSelection)) {
      // hmmm... failed to acquire ownership of XdndSelection!
      return kFALSE;
   }

   if (grabWin == kNone) grabWin = fMain->GetId();

   gVirtualX->GrabPointer(grabWin, fGrabEventMask, kNone, fDNDNoDropCursor, kTRUE, kFALSE);

   fLocalTarget = 0;
   fDragging = kTRUE;
   fTarget = kNone;
   fTargetIsDNDAware = kFALSE;
   fStatusPending = kFALSE;
   if (fDropTimeout) delete fDropTimeout;
   fDropTimeout = 0;
   fDropAccepted = kFALSE;
   fAcceptedAction = kNone;
   fLocalAction = kNone;

   if (!fDragWin && fPic != kNone && fMask != kNone) {
      fDragWin = new TGDragWindow(gClient->GetDefaultRoot(), fPic, fMask);
      fDragWin->Move((x_root-fHotx)|1, (y_root-fHoty)|1);
      fDragWin->MapSubwindows();
      fDragWin->MapRaised();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop.

Bool_t TGDNDManager::Drop()
{
   if (!fDragging) return kFALSE;

   if (fTargetIsDNDAware) {
      if (fDropAccepted) {
         if (fStatusPending) {
            if (fDropTimeout) delete fDropTimeout;
            fDropTimeout = new TTimer(this, 5000);
         } else {
            SendDNDDrop(fTarget);
         }
      } else {
         SendDNDLeave(fTarget);
         fStatusPending = kFALSE;
      }
   }
   EndDrag();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// End dragging.

Bool_t TGDNDManager::EndDrag()
{
   if (!fDragging) return kFALSE;

   gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);

   if (fSource)
      SendDNDFinished(fSource);
   if (fLocalSource)
      fLocalSource->HandleDNDFinished();

   fDragging = kFALSE;
   if (fDragWin) {
      fDragWin->DeleteWindow();
      fDragWin = 0;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Process drag event.

Bool_t TGDNDManager::Drag(int x_root, int y_root, Atom_t action, Time_t timestamp)
{
   if (!fDragging) return kFALSE;

   Window_t newTarget = FindWindow(gVirtualX->GetDefaultRootWindow(),
                                   x_root, y_root, 15);

   if (newTarget == kNone) {
      Window_t t = GetRootProxy();
      if (t != kNone) newTarget = t;
   }

   if (fTarget != newTarget) {

      if (fTargetIsDNDAware) SendDNDLeave(fTarget);

      fTarget = newTarget;
      fTargetIsDNDAware = IsDNDAware(fTarget);
      fStatusPending = kFALSE;
      fDropAccepted = kFALSE;
      fAcceptedAction = kNone;

      if (fTargetIsDNDAware) SendDNDEnter(fTarget);

      if (fDragWin)
         gVirtualX->ChangeActivePointerGrab(fDragWin->GetId(), fGrabEventMask,
                                            fDNDNoDropCursor);
   }

   if (fTargetIsDNDAware && !fStatusPending) {
      SendDNDPosition(fTarget, x_root, y_root, action, timestamp);

      // this is to avoid sending XdndPosition messages over and over
      // if the target is not responding
      fStatusPending = kTRUE;
   }

   if (fDragWin) {
      fDragWin->RaiseWindow();
      fDragWin->Move((x_root-fHotx)|1, (y_root-fHoty)|1);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set root window proxy.

Bool_t TGDNDManager::SetRootProxy()
{
   Window_t mainw = fMain->GetId();
   int result = kFALSE;

   if (GetRootProxy() == kNone) {
      gVirtualX->ChangeProperties(gVirtualX->GetDefaultRootWindow(),
                                  fgDNDProxy, XA_WINDOW, 32,
                                  (unsigned char *) &mainw, 1);
      gVirtualX->ChangeProperties(mainw, fgDNDProxy, XA_WINDOW, 32,
                                  (unsigned char *) &mainw, 1);

      fProxyOurs = kTRUE;
      result = kTRUE;
   }
   // XSync(_dpy, kFALSE);
   gVirtualX->UpdateWindow(0);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove root window proxy.

Bool_t TGDNDManager::RemoveRootProxy()
{
   if (!fProxyOurs) return kFALSE;

   gVirtualX->DeleteProperty(fMain->GetId(), fgDNDProxy);
   gVirtualX->DeleteProperty(gVirtualX->GetDefaultRootWindow(), fgDNDProxy);
   // the following is to ensure that the properties
   // (specially the one on the root window) are deleted
   // in case the application is exiting...

   // XSync(_dpy, kFALSE);
   gVirtualX->UpdateWindow(0);

   fProxyOurs = kFALSE;

   return kTRUE;
}
