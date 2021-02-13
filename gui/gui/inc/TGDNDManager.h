// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   19/04/07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGDNDManager
#define ROOT_TGDNDManager

#include "TGFrame.h"

class TGMainFrame;
class TGDragWindow;
class TTimer;

//----------------------------------------------------------------------

class TGDragWindow : public TGFrame {

protected:
   static Cursor_t fgDefaultCursor; // Default Cursor

protected:
   virtual void DoRedraw();

   Window_t fInput;                 // Input Window
   Pixmap_t fPic, fMask;            // Pixmaps used as Window shape
   UInt_t   fPw, fPh;               // Hot point coordinates (x and y)

public:
   TGDragWindow(const TGWindow *p, Pixmap_t pic, Pixmap_t mask,
                UInt_t options = kChildFrame, Pixel_t back = GetWhitePixel());
   virtual ~TGDragWindow();

   virtual TGDimension GetDefaultSize() const { return TGDimension(fPw, fPh); }

   virtual void MapWindow();
   virtual void UnmapWindow();
   virtual void RaiseWindow();
   virtual void LowerWindow();
   virtual void MapRaised();

   virtual void Layout();

   Window_t GetInputId() const { return fInput; }
   Bool_t HasWindow(Window_t w) const { return (w == fId || w == fInput); }

   ClassDef(TGDragWindow, 0) // Window used for dragging
};

//----------------------------------------------------------------------

//_____________________________________________________________________________
//
// TDNDData
//
// Drag and drop data container.
//_____________________________________________________________________________

class TDNDData : public TObject {
private:
   TDNDData(const TDNDData&);            // Not implemented
   TDNDData& operator=(const TDNDData&); // Not implemented

public:
   TDNDData(Atom_t dt = kNone, void *d = 0, Int_t len = 0, Atom_t act = kNone) :
      fDataType(dt), fAction(act), fData(d), fDataLength(len) {}
   ~TDNDData() {}

   Atom_t    fDataType;       // Data type description
   Atom_t    fAction;         // Action description
   void     *fData;           // Actual data
   Int_t     fDataLength;     // Length of data

   ClassDef(TDNDData, 0) // Drag and drop specific data
};

//----------------------------------------------------------------------

class TGDNDManager : public TObject {

private:
   TGDNDManager(const TGDNDManager&);            // Not implemented
   TGDNDManager& operator=(const TGDNDManager&); // Not implemented

protected:
   TGFrame       *fMain;                         // pointer on TGMainFrame
   Atom_t         fVersion;                      // not really an Atom, but a long
   Atom_t        *fTypelist, *fDraggerTypes;     // lists of DND types
   Atom_t         fDropType;                     // drop type
   Atom_t         fAcceptedAction, fLocalAction; // accepted and local actions

   Bool_t         fDragging;                     // kTRUE while dragging
   Bool_t         fDropAccepted;                 // kTRUE if drop accepted
   Bool_t         fStatusPending;                // kTRUE if status is pending
   Bool_t         fUseVersion;                   // kTRUE if DND version is used
   Bool_t         fProxyOurs;                    // kTRUE if root proxy is ours
   Window_t       fSource, fTarget;              // source and target windows
   Bool_t         fTargetIsDNDAware;             // kTRUE if target is DND aware
   UInt_t         fGrabEventMask;                // pointer grab event mask
   TGFrame       *fLocalSource, *fLocalTarget;   // local source and target

   TTimer        *fDropTimeout;                  // drop timeout
   TGDragWindow  *fDragWin;                      // drag window

   Pixmap_t       fPic, fMask;                   // pixmap used for the drag window
   Int_t          fHotx, fHoty;                  // hot point coordinates
   Cursor_t       fDNDNoDropCursor;              // no drop cursor type

protected:
   static Atom_t  fgDNDAware, fgDNDSelection, fgDNDProxy;
   static Atom_t  fgDNDEnter, fgDNDLeave, fgDNDPosition, fgDNDStatus;
   static Atom_t  fgDNDDrop, fgDNDFinished;
   static Atom_t  fgDNDVersion;
   static Atom_t  fgDNDActionCopy, fgDNDActionMove, fgDNDActionLink;
   static Atom_t  fgDNDActionAsk, fgDNDActionPrivate;
   static Atom_t  fgDNDTypeList, fgDNDActionList, fgDNDActionDescrip;
   static Atom_t  fgXCDNDData;

   static Bool_t  fgInit;
   static Atom_t  fgXAWMState;

protected:
   void           InitAtoms();
   Window_t       GetRootProxy();
   Window_t       FindWindow(Window_t root, Int_t x, Int_t y, Int_t maxd);
   Bool_t         IsDNDAware(Window_t win, Atom_t *typelist = 0);
   Bool_t         IsTopLevel(Window_t win);

   void           SendDNDEnter(Window_t target);
   void           SendDNDLeave(Window_t target);
   void           SendDNDPosition(Window_t target, int x, int y,
                                  Atom_t action, Time_t timestamp);
   void           SendDNDStatus(Window_t target, Atom_t action);
   void           SendDNDDrop(Window_t target);
   void           SendDNDFinished(Window_t src);

   Bool_t         HandleDNDEnter(Window_t src, long vers, Atom_t dataTypes[3]);
   Bool_t         HandleDNDLeave(Window_t src);
   Bool_t         HandleDNDPosition(Window_t src, int x_root, int y_root, Atom_t action, Time_t timestamp);
   Bool_t         HandleDNDStatus(Window_t from, int accepted,
                                  Rectangle_t skip, Atom_t action);
   Bool_t         HandleDNDDrop(Window_t src, Time_t timestamp);
   Bool_t         HandleDNDFinished(Window_t target);

public:
   TGDNDManager(TGFrame *toplevel, Atom_t *typelist);
   virtual ~TGDNDManager();

   Bool_t         HandleClientMessage(Event_t *event);
   Bool_t         HandleSelectionRequest(Event_t *event);
   Bool_t         HandleSelection(Event_t *event);

   Bool_t         HandleTimer(TTimer *t);

  //--- called by widgets

   TGFrame       *GetMainFrame() const { return fMain; }
   void           SetMainFrame(TGFrame *main) { fMain = main; }
   void           SetDragPixmap(Pixmap_t pic, Pixmap_t mask, Int_t hot_x, Int_t hot_y);
   Bool_t         SetRootProxy();
   Bool_t         RemoveRootProxy();

   Bool_t         StartDrag(TGFrame *src, Int_t x_root, Int_t y_root,
                            Window_t grabWin = kNone);
   Bool_t         Drag(Int_t x_root, Int_t y_root, Atom_t action, Time_t timestamp);
   Bool_t         Drop();
   Bool_t         EndDrag();

   Bool_t         IsDragging() const { return fDragging; }
   Window_t       GetSource() const { return fSource; }
   Window_t       GetTarget() const { return fTarget; }
   Atom_t        *GetTypeList() const { return fTypelist; }

   static Atom_t  GetDNDAware();
   static Atom_t  GetDNDSelection();
   static Atom_t  GetDNDProxy();
   static Atom_t  GetDNDEnter();
   static Atom_t  GetDNDLeave();
   static Atom_t  GetDNDPosition();
   static Atom_t  GetDNDStatus();
   static Atom_t  GetDNDDrop();
   static Atom_t  GetDNDFinished();
   static Atom_t  GetDNDVersion();
   static Atom_t  GetDNDActionCopy();
   static Atom_t  GetDNDActionMove();
   static Atom_t  GetDNDActionLink();
   static Atom_t  GetDNDActionAsk();
   static Atom_t  GetDNDActionPrivate();
   static Atom_t  GetDNDTypeList();
   static Atom_t  GetDNDActionList();
   static Atom_t  GetDNDActionDescrip();
   static Atom_t  GetXCDNDData();

   ClassDef(TGDNDManager, 0) // The main Drag and Drop Manager
};

R__EXTERN TGDNDManager *gDNDManager; // global drag and drop manager

#endif  // ROOT_TGDNDManager

