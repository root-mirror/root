// @(#)root/guibuilder:$Id: d2f0a1966f9911570cafe9d356f2158a2773edd1 $
// Author: Valeriy Onuchin   12/09/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TRootGuiBuilder.h"
#include "TGuiBldDragManager.h"
#include "TGuiBldEditor.h"

#include "TGShutter.h"
#include "TGSplitter.h"
#include "TGLayout.h"
#include "TGResourcePool.h"
#include "TGButton.h"
#include "TROOT.h"
#include "TGDockableFrame.h"
#include "TGMdi.h"
#include "TGStatusBar.h"
#include "TG3DLine.h"
#include "TGLabel.h"
#include "TColor.h"
#include "TGToolBar.h"
#include "TGToolTip.h"
#include "KeySymbols.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TSystem.h"
#include "TRootHelpDialog.h"
#include "TGListTree.h"
#include "TImage.h"
#include "TGTextEdit.h"
#include "TGTab.h"
#include "TGListBox.h"
#include "TGComboBox.h"
#include "TGProgressBar.h"
#include "TVirtualX.h"


/** \class TRootGuiBuilder
    \ingroup guibuilder

### %ROOT GUI Builder principles

 With the GUI builder, we try to make the next step from WYSIWYG
 to embedded editing concept - WYSIWYE ("what you see is what you edit").
 The ROOT GUI Builder allows modifying real GUI objects.
 For example, one can edit the existing GUI application created by
 $ROOTSYS/tutorials/gui/guitest.C.
 GUI components can be added to a design area from a widget palette,
 or can be borrowed from another application.
 One can drag and and drop TCanvas's menu bar into the application.
 GUI objects can be resized and dragged, copied and pasted.
 ROOT GUI Builder allows changing the layout, snap to grid, change object's
 layout order via the GUI Builder toolbar, or by options in the right-click
 context menus.
 A final design can be immediatly tested and used, or saved as a C++ macro.
 For example, it's possible to rearrange buttons in control bar,
 add separators etc. and continue to use a new fancy control bar in the
 application.


 The following is a short description of the GUI Builder actions and key shortcuts:

  - Press Ctrl-Double-Click to start/stop edit mode
  - Press Double-Click to activate quick edit action (defined in root.mimes)
  - Warning: some shortcuts might not work if NumLock is enabled

### Selection, grabbing, dropping

   It is possible to select, drag any frame and drop it to any frame

  - Click left mouse button or Ctrl-Click to select an object to edit.
  - Press right mouse button to activate context menu
  - Multiple selection (grabbing):
     - draw lasso and press Return key
     - press Shift key and draw lasso
  - Dropping:
     - select frame and press Ctrl-Return key
  - Changing layout order:
     - select frame and use arrow keys to change layout order
  - Alignment:
     - draw lasso and press arrow keys (or Shift-Arrow key) to align frames

### Key shortcuts

  - Return    - grab selected frames
  - Ctrl-Return - drop frames
  - Del       - delete selected frame
  - Shift-Del - crop action
  - Ctrl-X    - cut action
  - Ctrl-C    - copy action
  - Ctrl-V    - paste frame into the last clicked position
  - Ctrl-L    - compact
  - Ctrl-B    - enable/disable layout
  - Ctrl-H    - switch horizontal-vertical layout
  - Ctrl-G    - switch on/off grid
  - Ctrl-S    - save action
  - Ctrl-O    - open and execute a ROOT macro file.  GUI components created
                after macro execution will be emebedded to currently edited
                design area.
  - Ctrl-N    - create new main frame

*/


const char gHelpBuilder[] = "\
               Start/Stop Edit Mode\n\
     ************************************************\n\
 o Select File menu / Edit\n\
 o Select Start Edit button on the toolbar\n\
 o Ctrl-Double-Click on the project frame\n\
 o Double-Click to activate quick edit action (defined in root.mimes)\n\
\n\
               Select, Grab, Drop\n\
     ************************************************\n\
  It is possible to select & drag any frame and drop it to another frame\n\
\n\
 o Press left mouse button Click or Ctrl-Click to select an object.\n\
 o Press right mouse button to activate context menu\n\
 o Multiple selection can be done in two ways (grabbing):\n\
      - draw lasso and press Return key\n\
      - press Shift key and draw lasso\n\
 o Dropping:\n\
      - select frame and press Ctrl-Return key\n\
 o Changing layout order of widgets:\n\
      - set broken layout mode via toolbar button or check button\n\
        \'Layout subframes\' in tab \'Layout\'\n\
      - select a widget and use arrow keys to change the layout order\n\
 o Alignment:\n\
      - remove the selection (if any) by using the space bar\n\
      - draw lasso and use the four toolbar buttons for widget alignment\n\
      - arrow keys align the frames too, if you prefer the keyboard\n\
\n\
                    Key shortcuts\n\
     ************************************************\n\
 o Return    - grab selected frames\n\
 o Ctrl-Return - drop frames\n\
 o Del       - delete selected frame\n\
 o Shift-Del - crop\n\
 o Ctrl-X    - cut\n\
 o Ctrl-C    - copy\n\
 o Ctrl-V    - paste frame into the last clicked position\n\
 o Ctrl-L    - compact frame\n\
 o Ctrl-B    - enable/disable layout\n\
 o Ctrl-H    - switch Horizontal-Vertical layout\n\
 o Ctrl-G    - switch ON/OFF grid\n\
 o Ctrl-S    - save\n\
 o Ctrl-O    - open and execute ROOT macro file\n\
 o Ctrl-N    - create new main frame\n\
 o Ctrl-Z    - undo last action (not implemented)\n\
 o Shift-Ctrl-Z - redo (not implemented)\n\
\n\
                    More information\n\
     ************************************************\n\
\n\
For more information, please see the GuiBuilder Howto page at:\n\
\n\
   http://root.cern.ch/root/HowtoGuiBuilder.html\n\
\n\
";

const char gHelpAboutBuilder[] = "\
                  ROOT Gui Builder\n\
\n\
************************************************************\n\
* Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.  *\n\
* All rights reserved.                                     *\n\
*                                                          *\n\
* For the licensing terms see $ROOTSYS/LICENSE.            *\n\
* For the list of contributors see $ROOTSYS/README/CREDITS.*\n\
************************************************************\n\
";

//----- Toolbar stuff...

static ToolBarData_t gToolBarData[] = {
   { "bld_edit.png",   "Start Edit (Ctrl-Dbl-Click)",   kFALSE, kEditableAct, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "bld_new.png",   "New (Ctrl-N)",   kFALSE, kNewAct, 0 },
   { "bld_open.png",   "Open (Ctrl-O)",   kFALSE, kOpenAct, 0 },
   { "bld_save.png",   "Save As (Ctrl-S)",   kFALSE, kSaveAct, 0 },
   { "",                 "",               kFALSE, -1, 0 },
//   { "bld_pointer.xpm",   "Selector (Ctrl-Click)",   kTRUE, kSelectAct, 0 },
//   { "bld_grab.xpm",   "Grab Selected Frames (Return)",   kTRUE, kGrabAct, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "bld_compact.png",   "Compact selected frame (Ctrl-L)",        kFALSE,  kCompactAct, 0 },
   { "bld_break.png",   "Disable/Enable layout (Ctrl-B)",        kFALSE,  kBreakLayoutAct, 0 },
   { "bld_hbox.png",  "Layout selected frame horizontally (Ctrl-H)",    kFALSE,  kLayoutHAct, 0 },
   { "bld_vbox.png",   "Layout selected frame vertically (Ctrl-H)",    kFALSE,  kLayoutVAct, 0 },
   { "bld_grid.png",   "On/Off grid (Ctrl+G)",     kFALSE,  kGridAct, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "bld_AlignTop.png",   "Align selected frames to the top line of lasso (Up  Arrow)",        kFALSE,  kUpAct, 0 },
   { "bld_AlignBtm.png",   "Align selected frames to the down line of lasso (Down Arrow)",        kFALSE,  kDownAct, 0 },
   { "bld_AlignLeft.png",   "Align selected frames to the left line of lasso (Left  Arrow)",        kFALSE,  kLeftAct, 0 },
   { "bld_AlignRight.png",   "Align selected frames to the right line of lasso (Right  Arrow)",        kFALSE,  kRightAct, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "bld_cut.png",   "Cut (Ctrl-X)",        kFALSE,  kCutAct, 0 },
   { "bld_copy.png",   "Copy (Ctrl-C)",        kFALSE,  kCopyAct, 0 },
   { "bld_paste.png",   "Paste frame into the last clicked position (Ctrl-V)",        kFALSE,  kPasteAct, 0 },
//   { "bld_paste_into.png",   "Paste with replacing of selected frame (Ctrl-R)",        kFALSE,  kReplaceAct, 0 },
   { "bld_delete.png",   "Delete (Del/Backspace)",        kFALSE,  kDeleteAct, 0 },
   { "bld_crop.png",   "Crop (Shift-Del)",        kFALSE,  kCropAct, 0 },
//   { "",                 "",               kFALSE, -1, 0 },
//   { "bld_undo.png",   "Undo (Ctrl-Z)",        kFALSE,  kUndoAct, 0 },
//   { "bld_redo.png",   "Redo (Shift-Ctrl-Z)",        kFALSE,  kRedoAct, 0 },
   { 0,                  0,                kFALSE, 0, 0 }
};


ClassImp(TRootGuiBuilder);


TGGC *TRootGuiBuilder::fgBgnd = 0;
TGGC *TRootGuiBuilder::fgBgndPopup = 0;
TGGC *TRootGuiBuilder::fgBgndPopupHlght = 0;


////////////////////////////////////////////////////////////////////////////////
//
// Here are few experimental GUI classes which give a nice&fancy appearence
// to GuiBuilder.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
class TGuiBldMenuTitle : public TGMenuTitle {

private:
   Pixel_t fBgndColor;

protected:
   void DoRedraw();

public:
   virtual ~TGuiBldMenuTitle() {}
   TGuiBldMenuTitle(const TGWindow *p, TGHotString *s, TGPopupMenu *menu) :
      TGMenuTitle(p, s, menu) {
         fEditDisabled = kEditDisable;
         fBgndColor = TRootGuiBuilder::GetBgnd();
         SetBackgroundColor(fBgndColor);
         AddInput(kEnterWindowMask | kLeaveWindowMask);
   }

   Bool_t HandleCrossing(Event_t *event);
};

////////////////////////////////////////////////////////////////////////////////
/// Handle  crossing events.

Bool_t TGuiBldMenuTitle::HandleCrossing(Event_t *event)
{
   if (event->fType == kEnterNotify) {
      fBgndColor = TRootGuiBuilder::GetPopupHlght();
   } else {
      fBgndColor = TRootGuiBuilder::GetBgnd();
   }
   DoRedraw();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw builder menu title.

void TGuiBldMenuTitle::DoRedraw()
{
   TGFrame::DoRedraw();

   int x, y, max_ascent, max_descent;
   x = y = 4;

   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   if (fState) {
      gVirtualX->SetForeground(fNormGC, GetDefaultSelectedBackground());
      gVirtualX->FillRectangle(fId,fNormGC, 0, 0, fWidth, fHeight);
      gVirtualX->SetForeground(fNormGC, GetForeground());
      fLabel->Draw(fId, fSelGC, x, y + max_ascent);
   } else {
      gVirtualX->SetForeground(fNormGC, fBgndColor);
      gVirtualX->FillRectangle(fId,fNormGC, 0, 0, fWidth, fHeight);
      gVirtualX->SetForeground(fNormGC, GetForeground());
      fLabel->Draw(fId, fNormGC, x, y + max_ascent);
   }
   if (fBgndColor == TRootGuiBuilder::GetPopupHlght()) {
      gVirtualX->DrawRectangle(fId, TGFrame::GetBlackGC()(),
                               0, 0, GetWidth()-1, GetHeight()-1);
   }
}


////////////////////////////////////////////////////////////////////////////////
class TGuiBldPopupMenu : public TGPopupMenu {

public:
   virtual ~TGuiBldPopupMenu() { }
   TGuiBldPopupMenu() :
      TGPopupMenu(gClient->GetDefaultRoot()) {
      fEditDisabled = kEditDisable;
      SetBackgroundColor(TRootGuiBuilder::GetPopupBgnd());
      fEntrySep = 8;
   }
   void DrawEntry(TGMenuEntry *entry);
};

////////////////////////////////////////////////////////////////////////////////
/// Draw popup menu entry.

void TGuiBldPopupMenu::DrawEntry(TGMenuEntry *entry)
{
   FontStruct_t  font;
   GCValues_t    gcval;

   if (entry->GetStatus() & kMenuHideMask)
      return;

   if (entry->GetStatus() & kMenuDefaultMask) {
      font = fHifontStruct;
      gcval.fMask = kGCFont;
      gcval.fFont = gVirtualX->GetFontHandle(font);
      gVirtualX->ChangeGC(fNormGC, &gcval);
      gVirtualX->ChangeGC(fSelGC, &gcval);
   } else {
      font = fFontStruct;
   }

   UInt_t tw = 0;
   UInt_t sep = fEntrySep;
   Int_t max_ascent, max_descent;
   gVirtualX->GetFontProperties(font, max_ascent, max_descent);
   if (entry->GetShortcut())
      tw = 7 + gVirtualX->TextWidth(fFontStruct, entry->GetShortcutText(),
                                    entry->GetShortcut()->Length());

   Int_t tx = entry->GetEx() + fXl;
   Int_t ty = entry->GetEy() + max_ascent + 2;
   UInt_t h = max_ascent + max_descent + sep;
   Int_t picposy = 0;
   if (entry->GetPic() != 0) {
      picposy = entry->GetEy() + h / 2;
      picposy -= entry->GetPic()->GetHeight() / 2;
   }

   switch (entry->GetType()) {
      case kMenuPopup:
      case kMenuLabel:
      case kMenuEntry:
         if ((entry->GetStatus() & kMenuActiveMask) &&
             entry->GetType() != kMenuLabel) {
            if (entry->GetStatus() & kMenuEnableMask) {
               gVirtualX->FillRectangle(fId,
                              TRootGuiBuilder::GetPopupHlghtGC()->GetGC(),
                              entry->GetEx()+1, entry->GetEy(),
                              fMenuWidth-6, h - 1);
               gVirtualX->DrawRectangle(fId,  TGFrame::GetBlackGC()(),
                                        entry->GetEx()+ 1, entry->GetEy()-1,
                                        fMenuWidth - entry->GetEx()- 6, h - 1);
            }

            if (entry->GetType() == kMenuPopup) {
               DrawTrianglePattern(fSelGC, fMenuWidth-10, entry->GetEy() + 3,
                                   fMenuWidth-6, entry->GetEy() + 11);
            }

            if (entry->GetStatus() & kMenuCheckedMask) {
               DrawCheckMark(fSelGC, 6, entry->GetEy()+sep, 14,
                             entry->GetEy()+11);
            }

            if (entry->GetStatus() & kMenuRadioMask) {
               DrawRCheckMark(fSelGC, 6, entry->GetEy()+sep, 14,
                              entry->GetEy()+11);
            }

            if (entry->GetPic() != 0) {
               entry->GetPic()->Draw(fId, fSelGC, 8, picposy);
            }

            entry->GetLabel()->Draw(fId,
                           (entry->GetStatus() & kMenuEnableMask) ? fSelGC :
                            GetShadowGC()(), tx, ty);
            if (entry->GetShortcut())
               entry->GetShortcut()->Draw(fId,
                           (entry->GetStatus() & kMenuEnableMask) ? fSelGC :
                           GetShadowGC()(), fMenuWidth - tw, ty);
         } else {
            if ( entry->GetType() != kMenuLabel) {
               gVirtualX->FillRectangle(fId,
                           TRootGuiBuilder::GetBgndGC()->GetGC(),
                           entry->GetEx()+1, entry->GetEy()-1, tx-4, h);

               gVirtualX->FillRectangle(fId,
                           TRootGuiBuilder::GetPopupBgndGC()->GetGC(),
                           tx-1, entry->GetEy()-1, fMenuWidth-tx-1, h);
            } else { // we need some special background for labels
               gVirtualX->FillRectangle(fId, TGFrame::GetBckgndGC()(),
                                       entry->GetEx()+1, entry->GetEy()-1,
                                       fMenuWidth - entry->GetEx()- 3, h);
            }

            if (entry->GetType() == kMenuPopup) {
               DrawTrianglePattern(fNormGC, fMenuWidth-10, entry->GetEy() + 3,
                                   fMenuWidth-6, entry->GetEy() + 11);
            }

            if (entry->GetStatus() & kMenuCheckedMask) {
               DrawCheckMark(fNormGC, 6, entry->GetEy()+sep, 14,
                             entry->GetEy()+11);
            }

            if (entry->GetStatus() & kMenuRadioMask) {
               DrawRCheckMark(fNormGC, 6, entry->GetEy()+sep, 14,
                              entry->GetEy()+11);
            }

            if (entry->GetPic() != 0) {
               entry->GetPic()->Draw(fId, fNormGC, 8, picposy);
            }

            if (entry->GetStatus() & kMenuEnableMask) {
               entry->GetLabel()->Draw(fId, fNormGC, tx, ty);
               if (entry->GetShortcut())
                  entry->GetShortcut()->Draw(fId, fNormGC, fMenuWidth - tw, ty);
            } else {
               entry->GetLabel()->Draw(fId, GetHilightGC()(), tx+1, ty+1);
               entry->GetLabel()->Draw(fId, GetShadowGC()(), tx, ty);
               if (entry->GetShortcut()) {
                  entry->GetShortcut()->Draw(fId, GetHilightGC()(),
                                             fMenuWidth - tw+1, ty+1);
                  entry->GetShortcut()->Draw(fId, GetShadowGC()(),
                                             fMenuWidth - tw, ty);
               }
            }
         }
         break;

      case kMenuSeparator:
         gVirtualX->FillRectangle(fId, TRootGuiBuilder::GetBgndGC()->GetGC(),
                                     entry->GetEx()+1, entry->GetEy()-1,
                                     tx-4, 4);
         gVirtualX->DrawLine(fId, TGFrame::GetBlackGC()(), tx+1,
                             entry->GetEy()+1, fMenuWidth-sep,
                             entry->GetEy()+1);
         break;
   }

   // restore font
   if (entry->GetStatus() & kMenuDefaultMask) {
      gcval.fFont = gVirtualX->GetFontHandle(fFontStruct);
      gVirtualX->ChangeGC(fNormGC, &gcval);
      gVirtualX->ChangeGC(fSelGC, &gcval);
   }
}

////////////////////////////////////////////////////////////////////////////////
class TGuiBldToolButton : public TGPictureButton {

private:
   Pixel_t fBgndColor;

protected:
   void  DoRedraw();

public:
   virtual ~TGuiBldToolButton() { }
   TGuiBldToolButton(const TGWindow *p, const TGPicture *pic, Int_t id = -1) :
         TGPictureButton(p, pic, id) {
      fBgndColor = TRootGuiBuilder::GetBgnd();
      ChangeOptions(GetOptions() & ~kRaisedFrame);
   }

   Bool_t IsDown() const { return (fOptions & kSunkenFrame); }
   void SetState(EButtonState state, Bool_t emit = kTRUE);
   Bool_t HandleCrossing(Event_t *event);
   void SetBackgroundColor(Pixel_t bgnd) { fBgndColor = bgnd; TGFrame::SetBackgroundColor(bgnd); }
};

////////////////////////////////////////////////////////////////////////////////
/// Redraw tool button.

void TGuiBldToolButton::DoRedraw()
{
   int x = (fWidth - fTWidth) >> 1;
   int y = (fHeight - fTHeight) >> 1;
   UInt_t w = GetWidth() - 1;
   UInt_t h = GetHeight()- 1;

   TGFrame::SetBackgroundColor(fBgndColor);

   TGFrame::DoRedraw();
   if (fState == kButtonDown || fState == kButtonEngaged) {
      ++x; ++y;
      w--; h--;
   }

   const TGPicture *pic = fPic;
   if (fState == kButtonDisabled) {
      if (!fPicD) CreateDisabledPicture();
      pic = fPicD ? fPicD : fPic;
   }
   if (fBgndColor == TRootGuiBuilder::GetPopupHlght()) {
      x--; y--;
      gVirtualX->DrawRectangle(fId, TGFrame::GetBlackGC()(), 0, 0, w, h);
   }
   pic->Draw(fId, fNormGC, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle crossing events.

Bool_t TGuiBldToolButton::HandleCrossing(Event_t *event)
{
   if (fTip) {
      if (event->fType == kEnterNotify) {
         fTip->Reset();
      } else {
         fTip->Hide();
      }
   }

   if ((event->fType == kEnterNotify) && (fState != kButtonDisabled)) {
      fBgndColor = TRootGuiBuilder::GetPopupHlght();
   } else {
      fBgndColor = TRootGuiBuilder::GetBgnd();
   }
   if (event->fType == kLeaveNotify) {
      fBgndColor = TRootGuiBuilder::GetBgnd();
      if (fState != kButtonDisabled && fState != kButtonEngaged)
         SetState(kButtonUp, kFALSE);
   }
   DoRedraw();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set state of tool bar button and emit a signal according
/// to passed arguments.

void TGuiBldToolButton::SetState(EButtonState state, Bool_t emit)
{
   Bool_t was = !IsDown();

   if (state != fState) {
      switch (state) {
         case kButtonEngaged:
         case kButtonDown:
            fOptions &= ~kRaisedFrame;
            fOptions |= kSunkenFrame;
            break;
         case kButtonDisabled:
         case kButtonUp:
            fOptions &= ~kRaisedFrame;
            fOptions &= ~kSunkenFrame;
            break;
      }
      fState = state;
      DoRedraw();
      if (emit) EmitSignals(was);
   }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Create GUI builder application.

TRootGuiBuilder::TRootGuiBuilder(const TGWindow *p) : TGuiBuilder(),
   TGMainFrame(p ? p : gClient->GetDefaultRoot(), 1, 1)
{
   SetCleanup(kDeepCleanup);
   gGuiBuilder  = this;
   fManager = 0;
   fEditor = 0;
   fActionButton = 0;
   fClosing = 0;

   if (gDragManager) {
      fManager = (TGuiBldDragManager *)gDragManager;
   } else {
      gDragManager = fManager = new TGuiBldDragManager();
   }
   fManager->SetBuilder(this);

   fMenuBar = new TGMdiMenuBar(this, 10, 10);
   AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
   InitMenu();

   AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsTop |
            kLHintsExpandX, 0,0,2,2));

   fToolDock = new TGDockableFrame(this);
   AddFrame(fToolDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
   fToolDock->SetWindowName("GuiBuilder ToolBar");

   fToolBar = new TGToolBar(fToolDock);
   fToolDock->AddFrame(fToolBar, new TGLayoutHints(kLHintsTop |
                       kLHintsExpandX));

   int spacing = 8;

   for (int i = 0; gToolBarData[i].fPixmap; i++) {
      if (strlen(gToolBarData[i].fPixmap) == 0) {
         spacing = 8;
         continue;
      }

      const TGPicture *pic = fClient->GetPicture(gToolBarData[i].fPixmap);
      TGuiBldToolButton *pb = new TGuiBldToolButton(fToolBar, pic,
                                                    gToolBarData[i].fId);
      pb->SetStyle(gClient->GetStyle());

      pb->SetToolTipText(gToolBarData[i].fTipText);

      TGToolTip *tip = pb->GetToolTip();
      tip->SetDelay(200);

      tip->Connect("Reset()", "TRootGuiBuilder", this, "UpdateStatusBar(=0)");
      tip->Connect("Hide()", "TRootGuiBuilder", this, "EraseStatusBar()");

      fToolBar->AddButton(this, pb, spacing);
      spacing = 0;

      if (gToolBarData[i].fId == kEditableAct) {
         fStartButton = pb;
         continue;
      }

      if ((gToolBarData[i].fId == kUndoAct) ||
          (gToolBarData[i].fId == kRedoAct)) {
         pb->SetState(kButtonDisabled);
      }
   }

   fToolBar->Connect("Clicked(Int_t)", "TGuiBldDragManager", fManager,
                     "HandleAction(Int_t)");

   AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsTop |
            kLHintsExpandX, 0,0,2,5));

   TGCompositeFrame *cf = new TGHorizontalFrame(this, 1, 1);
   AddFrame(cf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   //fShutterDock = new TGDockableFrame(cf);
   //cf->AddFrame(fShutterDock, new TGLayoutHints(kLHintsNormal ));
   //fShutterDock->SetWindowName("Widget Factory");
   //fShutterDock->EnableUndock(kTRUE);
   //fShutterDock->EnableHide(kTRUE);
   //fShutterDock->DockContainer();

   fShutter = new TGShutter(cf, kSunkenFrame);
   cf->AddFrame(fShutter, new TGLayoutHints(kLHintsNormal | kLHintsExpandY));
   fShutter->ChangeOptions(fShutter->GetOptions() | kFixedWidth);

   TGVSplitter *splitter = new TGVSplitter(cf);
   splitter->SetFrame(fShutter, kTRUE);
   cf->AddFrame(splitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   fMain = new TGMdiMainFrame(cf, fMenuBar, 1, 1);
   fMain->Connect("FrameClosed(Int_t)", "TRootGuiBuilder", this,
                  "HandleWindowClosed(Int_t)");

   TQObject::Connect("TGMdiFrame", "CloseWindow()", "TRootGuiBuilder", this,
                     "MaybeCloseWindow()");

   cf->AddFrame(fMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   fMain->GetContainer()->SetEditDisabled(kEditDisable);

   const TGPicture *pbg = fClient->GetPicture("bld_bg.png");
   if (pbg) {
      fMain->GetContainer()->SetBackgroundPixmap(pbg->GetPicture());
   }

   if (fManager) {
      fEditor = new TGuiBldEditor(cf);
      cf->AddFrame(fEditor, new TGLayoutHints(kLHintsNormal | kLHintsExpandY));
      fManager->SetPropertyEditor(fEditor);
      fEditor->SetEmbedded();
   }

   AddSection("Projects");
   AddSection("Buttons");
   AddSection("Containers");
   AddSection("Bars");
   AddSection("Input");
   AddSection("Complex Input");
   AddSection("Display");
   AddSection("Dialogs");

   // create an empty section
   AddSection("User's Macros");
   TGShutterItem *item = fShutter->GetItem("User's Macros");
   TGCompositeFrame *cont = (TGCompositeFrame *)item->GetContainer();
   cont->SetBackgroundColor(TColor::Number2Pixel(18));

   TGuiBldAction *act = new TGuiBldAction("TGMainFrame", "Empty Frame",
                                          kGuiBldProj);
   act->fAct = "empty";
   act->fPic = "bld_mainframe.xpm";
   AddAction(act, "Projects");

   act = new TGuiBldAction("TGMainFrame", "Horizontal Frame", kGuiBldProj);
   act->fAct = "horizontal";
   act->fPic = "bld_mainframe.xpm";
   AddAction(act, "Projects");

   act = new TGuiBldAction("TGMainFrame", "Vertical Frame", kGuiBldProj);
   act->fAct = "vertical";
   act->fPic = "bld_mainframe.xpm";
   AddAction(act, "Projects");

   // Standard
   act = new TGuiBldAction("TGTextButton", "Text Button", kGuiBldCtor);
   act->fAct = "new TGTextButton()";
   act->fPic = "bld_textbutton.xpm";
   AddAction(act, "Buttons");

   act = new TGuiBldAction("TGCheckButton", "Check Button", kGuiBldCtor);
   act->fAct = "new TGCheckButton()";
   act->fPic = "bld_checkbutton.xpm";
   AddAction(act, "Buttons");

   act = new TGuiBldAction("TGRadioButton", "Radio Button", kGuiBldCtor);
   act->fAct = "new TGRadioButton()";
   act->fPic = "bld_radiobutton.xpm";
   AddAction(act, "Buttons");

   act = new TGuiBldAction("TGPictureButton", "Picture Button", kGuiBldCtor);
   act->fAct = "new TGPictureButton()";
   act->fPic = "bld_image.xpm";
   AddAction(act, "Buttons");

   act = new TGuiBldAction("TGTextEntry", "Text Entry", kGuiBldCtor);
   act->fAct = "new TGTextEntry()";
   act->fPic = "bld_entry.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGTextEdit", "Text Edit", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildTextEdit()";
   act->fPic = "bld_text.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGNumberEntry", "Number Entry", kGuiBldCtor);
   act->fAct = "new TGNumberEntry()";
   act->fPic = "bld_numberentry.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGComboBox", "Combo Box", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildComboBox()";
   act->fPic = "bld_combobox.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGListBox", "List Box", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildListBox()";
   act->fPic = "bld_listbox.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGHSlider", "Horizontal Slider", kGuiBldCtor);
   act->fAct = "new TGHSlider()";
   act->fPic = "bld_hslider.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGVSlider", "Vertical Slider", kGuiBldCtor);
   act->fAct = "new TGVSlider()";
   act->fPic = "bld_vslider.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGHScrollBar", "HScrollbar", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildHScrollBar()";
   act->fPic = "bld_hscrollbar.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGVScrollBar", "VScrollbar", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildVScrollBar()";
   act->fPic = "bld_vscrollbar.xpm";
   AddAction(act, "Input");

   act = new TGuiBldAction("TGListTree", "List Tree", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildListTree()";
   act->fPic = "bld_listtree.xpm";
   AddAction(act, "Complex Input");

   act = new TGuiBldAction("TGLabel", "Text Label", kGuiBldCtor);
   act->fAct = "new TGLabel()";
   act->fPic = "bld_label.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TGIcon", "Icon", kGuiBldCtor);
   act->fAct = "new TGIcon()";
   act->fPic = "bld_image.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TGHorizontal3DLine", "Horizontal Line",
                           kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildH3DLine()";
   act->fPic = "bld_hseparator.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TGVertical3DLine", "Vertical Line", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildV3DLine()";
   act->fPic = "bld_vseparator.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TGStatusBar", "Status Bar", kGuiBldCtor);
   act->fAct = "new TGStatusBar()";
   act->fPic = "bld_statusbar.xpm";
   act->fHints = new TGLayoutHints(kLHintsBottom | kLHintsExpandX);
   AddAction(act, "Bars");

   act = new TGuiBldAction("TGHProgressBar", "HProgress Bar", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildHProgressBar()";
   act->fPic = "bld_hprogressbar.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TGVProgressBar", "VProgress Bar", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildVProgressBar()";
   act->fPic = "bld_vprogressbar.xpm";
   AddAction(act, "Display");

   act = new TGuiBldAction("TRootEmbeddedCanvas", "Embed Canvas", kGuiBldCtor);
   act->fAct = "new TRootEmbeddedCanvas()";
   act->fPic = "bld_embedcanvas.xpm";
   AddAction(act, "Display");

   // Containers
   act = new TGuiBldAction("TGHorizontalFrame", "Horizontal Frame",
                           kGuiBldCtor);
   act->fAct = "new TGHorizontalFrame(0,200,100)";
   act->fPic = "bld_hbox.xpm";
   AddAction(act, "Containers");

   act = new TGuiBldAction("TGVerticalFrame", "Vertical Frame", kGuiBldCtor);
   act->fAct = "new TGVerticalFrame(0,100,200)";
   act->fPic = "bld_vbox.xpm";
   AddAction(act, "Containers");

   act = new TGuiBldAction("TGGroupFrame", "Group Frame", kGuiBldCtor);
   act->fAct = "new TGGroupFrame()";
   act->fPic = "bld_groupframe.xpm";
   AddAction(act, "Containers");

   act = new TGuiBldAction("TGTab", "Tabbed Frame", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildTab()";
   act->fPic = "bld_tab.xpm";
   AddAction(act, "Containers");

   act = new TGuiBldAction("TGShutter", "Shutter", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildShutter()";
   act->fPic = "bld_shutter.png";
   AddAction(act, "Containers");


   act = new TGuiBldAction("TGCanvas", "Scrolled Canvas", kGuiBldCtor);
   act->fAct = "TRootGuiBuilder::BuildCanvas()";
   act->fPic = "bld_canvas.xpm";
   AddAction(act, "Containers");
/*
   act = new TGuiBldAction("TGVSplitter", "Horizontal Panes", kGuiBldFunc);
   act->fAct = "TRootGuiBuilder::VSplitter()";
   act->fPic = "bld_hpaned.xpm";
   AddAction(act, "Containers");

   act = new TGuiBldAction("TGHSplitter", "Vertical Panes", kGuiBldFunc);
   act->fAct = "TRootGuiBuilder::HSplitter()";
   act->fPic = "bld_vpaned.xpm";
   AddAction(act, "Containers");
*/
   act = new TGuiBldAction("TGColorSelect", "Color Selector", kGuiBldFunc);
   act->fAct = "new TGColorSelect()";
   act->fPic = "bld_colorselect.xpm";
   AddAction(act, "Dialogs");

   fShutter->Resize(140, fShutter->GetHeight());

   fStatusBar = new TGStatusBar(this, 40, 10);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
            0, 0, 3, 0));

   PropagateBgndColor(this, GetBgnd());
   SetEditDisabled(kEditDisable);   // disable editting to all subframes

   MapSubwindows();

   Int_t qq;
   UInt_t ww;
   UInt_t hh;
   gVirtualX->GetWindowSize(gVirtualX->GetDefaultRootWindow(), qq, qq, ww, hh);
   MoveResize(100, 100, ww - 200, hh - 200);
   SetWMPosition(100, 100);

   SetWindowName("ROOT GuiBuilder");
   SetIconName("ROOT GuiBuilder");
   fIconPic = SetIconPixmap("bld_rgb.xpm");
   SetClassHints("ROOT", "GuiBuilder");

   fSelected = 0;
   Update();

   fMenuFile->Connect("Activated(Int_t)", "TRootGuiBuilder", this,
                      "HandleMenu(Int_t)");
   fMenuWindow->Connect("Activated(Int_t)", "TRootGuiBuilder", this,
                        "HandleMenu(Int_t)");
   fMenuHelp->Connect("Activated(Int_t)", "TRootGuiBuilder", this,
                      "HandleMenu(Int_t)");

   // doesn't work properly on Windows...
   if (gVirtualX->InheritsFrom("TGX11"))
      BindKeys();
   UpdateStatusBar("Ready");
   MapRaised();

   fEditor->SwitchLayout();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TRootGuiBuilder::~TRootGuiBuilder()
{
   if (fIconPic) gClient->FreePicture(fIconPic);
   delete fMenuFile;
   delete fMenuWindow;
   delete fMenuHelp;
   gGuiBuilder = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Close GUI builder via window manager "Close" button.

void TRootGuiBuilder::CloseWindow()
{
   TGWindow *root = (TGWindow*)fClient->GetRoot();
   if (root) root->SetEditable(kFALSE);

   fEditor->Reset();

   if (fMain->GetNumberOfFrames() == 0) {
      fMenuFile->DisableEntry(kGUIBLD_FILE_CLOSE);
      fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
      fMenuFile->DisableEntry(kGUIBLD_FILE_START);
   } else {
      fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
      fMenuFile->EnableEntry(kGUIBLD_FILE_START);
      fMenuFile->EnableEntry(kGUIBLD_FILE_CLOSE);
   }

   TGButton *btn = fToolBar->GetButton(kGridAct);
   if (btn) {
      btn->SetState(!fClient->IsEditable() ? kButtonDisabled : kButtonUp);
   }
   fClosing = 1;
   fMain->CloseAll();
   if (fClosing == -1) {
      fClosing = 0;
      return;
   }
   SwitchToolbarButton();
   Hide();
}

////////////////////////////////////////////////////////////////////////////////
/// Find action by name

TGButton *TRootGuiBuilder::FindActionButton(const char *name, const char *sect)
{
   if (!name || !sect) return 0;

   TGShutterItem *item = fShutter->GetItem(sect);
   if (!item) return 0;

   TGCompositeFrame *cont = (TGCompositeFrame *)item->GetContainer();
   TGHorizontalFrame *hf;
   TGFrameElement *fe;

   TIter next(cont->GetList());
   TGLabel *lb;
   TGButton *btn;

   while ((fe = (TGFrameElement*)next())) {
      hf = (TGHorizontalFrame*)fe->fFrame;
      btn = (TGButton*)((TGFrameElement*)hf->GetList()->First())->fFrame;
      lb = (TGLabel*)((TGFrameElement*)hf->GetList()->Last())->fFrame;
      if (*(lb->GetText()) == name) {
         return (TGButton*)btn;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add new action to widget palette.

void TRootGuiBuilder::AddAction(TGuiBldAction *act, const char *sect)
{
   if (!act || !sect) return;

   TGShutterItem *item = fShutter->GetItem(sect);
   TGButton *btn = 0;

   if (!item) return;
   TGCompositeFrame *cont = (TGCompositeFrame *)item->GetContainer();
   cont->SetBackgroundColor(TColor::Number2Pixel(18));

   const TGPicture *pic = 0;
   if (!act->fPicture) {
      act->fPicture = fClient->GetPicture(act->fPic);
   }
   pic = act->fPicture;

   TGHorizontalFrame *hf = new TGHorizontalFrame(cont);

   if (pic) {
      btn = new TGPictureButton(hf, pic);
   } else {
      btn = new TGTextButton(hf, act->GetName());
   }

   btn->SetToolTipText(act->GetTitle(), 200);
   btn->SetUserData((void*)act);
   btn->Connect("Clicked()", "TRootGuiBuilder", this, "HandleButtons()");

   hf->AddFrame(btn, new TGLayoutHints(kLHintsTop | kLHintsCenterY,3,3,3,3));

   TGLabel *lb = new TGLabel(hf, act->fType != kGuiBldMacro ? act->GetTitle() :
                             act->GetName());
   lb->SetBackgroundColor(cont->GetBackground());
   hf->AddFrame(lb, new TGLayoutHints(kLHintsTop | kLHintsCenterY,3,3,3,3));
   hf->SetBackgroundColor(cont->GetBackground());

   // disable edit
   cont->SetEditDisabled(kEditDisable);
   hf->SetEditDisabled(kEditDisable);

   cont->AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 0));
   cont->MapSubwindows();
   cont->Resize();  // invoke Layout()
}

////////////////////////////////////////////////////////////////////////////////
/// Add new shutter item.

void TRootGuiBuilder::AddSection(const char *sect)
{
   static int id = 10000;
   TGShutterItem *item = new TGShutterItem(fShutter, new TGHotString(sect),
                                           id++);
   fShutter->AddItem(item);
   item->Connect("Selected()", "TRootGuiBuilder", this, "HandleMenu(=3)");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle buttons in the GUI builder's widget palette.

void TRootGuiBuilder::HandleButtons()
{
   TGFrame *parent;

   if (fActionButton) {
      parent = (TGFrame*)fActionButton->GetParent();
      parent->ChangeOptions(parent->GetOptions() & ~kSunkenFrame);
      fClient->NeedRedraw(parent, kTRUE);
   }

   if (!fClient->IsEditable()) {
      HandleMenu(kGUIBLD_FILE_START);
   }

   fActionButton = (TGButton *)gTQSender;
   TGuiBldAction *act  = (TGuiBldAction *)fActionButton->GetUserData();
   parent = (TGFrame*)fActionButton->GetParent();

   parent->ChangeOptions(parent->GetOptions() | kSunkenFrame);
   fClient->NeedRedraw(parent, kTRUE);

   if (act) {
      fAction = act;
      fManager->UngrabFrame();
      if (fAction->fType != kGuiBldCtor) ExecuteAction();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute an action.

TGFrame *TRootGuiBuilder::ExecuteAction()
{
   if (!fAction || fAction->fAct.IsNull()) return 0;

   TGFrame *ret = 0;

   if (!fClient->IsEditable() && (fAction->fType != kGuiBldMacro)) {
      TGMdiFrame *current = fMain->GetCurrent();
      if (current) current->SetEditable(kTRUE);
   }

   TString s = "";

   switch (fAction->fType) {
      case kGuiBldProj:
         s = fAction->fAct.Data();
         NewProject(s);
         fAction = 0;
         break;
      case kGuiBldMacro:
         {
         TGWindow *root = (TGWindow*)fClient->GetRoot();
         if (root) root->SetEditable(kFALSE);
         gROOT->Macro(fAction->fAct.Data());
         if (root) root->SetEditable(kTRUE);
         fAction = 0;
         break;
         }
      default:
         ret = (TGFrame *)gROOT->ProcessLineFast(fAction->fAct.Data());
         break;
   }

   Update();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Inititiate GUI Builder menus.

void TRootGuiBuilder::InitMenu()
{
   fMenuFile = new TGuiBldPopupMenu();
   fMenuFile->AddEntry(new TGHotString("&Edit (Ctrl+double-click)"),
                       kGUIBLD_FILE_START, 0,
                       fClient->GetPicture("bld_edit.png"));
   fMenuFile->AddEntry(new TGHotString("&Stop (Ctrl+double-click)"),
                       kGUIBLD_FILE_STOP, 0,
                       fClient->GetPicture("bld_stop.png"));
   fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
   fMenuFile->DisableEntry(kGUIBLD_FILE_START);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(new TGHotString("&New Project"), kGUIBLD_FILE_NEW,
                       0, fClient->GetPicture("bld_new.png"));
   fMenuFile->AddEntry(new TGHotString("&Open"), kGUIBLD_FILE_OPEN,
                       0, fClient->GetPicture("bld_open.png"));
   fMenuFile->AddEntry(new TGHotString("&Close"), kGUIBLD_FILE_CLOSE,
                        0, fClient->GetPicture("bld_delete.png"));
   fMenuFile->AddEntry(new TGHotString("&Save project as"), kGUIBLD_FILE_SAVE,
                       0, fClient->GetPicture("bld_save.png"));
   fMenuFile->DisableEntry(kGUIBLD_FILE_CLOSE);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(new TGHotString("E&xit"), kGUIBLD_FILE_EXIT,
                        0, fClient->GetPicture("bld_exit.png"));
/*
   fMenuEdit = new TGuiBldPopupMenu();
   fMenuEdit->AddSeparator();
   fMenuEdit->AddEntry(new TGHotString("&Preferences ..."), kGUIBLD_EDIT_PREF);
*/
   fMenuWindow = new TGuiBldPopupMenu();
   fMenuWindow->AddEntry(new TGHotString("Tile &Horizontally"),
                         kGUIBLD_WINDOW_HOR);
   fMenuWindow->AddEntry(new TGHotString("Tile &Vertically"),
                         kGUIBLD_WINDOW_VERT);
   fMenuWindow->AddEntry(new TGHotString("&Cascade"),
                         kGUIBLD_WINDOW_CASCADE);
   fMenuWindow->AddSeparator();
   //fMenuWindow->AddPopup(new TGHotString("&Windows"), fMain->GetWinListMenu());
   fMenuWindow->AddEntry(new TGHotString("&Arrange icons"),
                         kGUIBLD_WINDOW_ARRANGE);
   fMenuWindow->AddSeparator();
   fMenuWindow->AddEntry(new TGHotString("&Opaque resize"),
                         kGUIBLD_WINDOW_OPAQUE);
   fMenuWindow->CheckEntry(kGUIBLD_WINDOW_OPAQUE);

   fMenuHelp = new TGuiBldPopupMenu();
   fMenuHelp->AddEntry(new TGHotString("&Contents"), kGUIBLD_HELP_CONTENTS);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry(new TGHotString("&About"), kGUIBLD_HELP_ABOUT);
   //fMenuHelp->AddSeparator();
   //fMenuHelp->AddEntry(new TGHotString("&Send Bug Report"),kGUIBLD_HELP_BUG);

   TGMenuBar *bar = fMenuBar->GetMenuBar();

   TGuiBldMenuTitle *title;
   title = new TGuiBldMenuTitle(bar, new TGHotString("&File"), fMenuFile);
   bar->AddTitle(title, new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0));

   //title = new TGuiBldMenuTitle(bar, new TGHotString("&Edit"), fMenuEdit);
   //bar->AddTitle(title, new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0));

   title = new TGuiBldMenuTitle(bar, new TGHotString("&Windows"), fMenuWindow);
   bar->AddTitle(title, new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0));

   title = new TGuiBldMenuTitle(bar, new TGHotString("&Help"), fMenuHelp);
   bar->AddTitle(title, new TGLayoutHints(kLHintsTop | kLHintsRight,4,4,0,0));

   fMenuBar->SetEditDisabled(kEditDisable);
   PropagateBgndColor(fMenuBar, GetBgnd());
}

////////////////////////////////////////////////////////////////////////////////
/// Set selected frame.

void TRootGuiBuilder::ChangeSelected(TGFrame *f)
{
   fSelected = f;
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Enable toolbar buttons for alignment.

void TRootGuiBuilder::EnableLassoButtons(Bool_t on)
{
   TGButton *btn = 0;

   btn = fToolBar->GetButton(kUpAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kDownAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kRightAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kLeftAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kDeleteAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kGrabAct);
   if (btn) {
      btn->SetState(kButtonUp);
   }

   btn = fToolBar->GetButton(kCropAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable toolbar buttons according to the selected frame.

void TRootGuiBuilder::EnableSelectedButtons(Bool_t on)
{
   fSelected = fManager->GetSelected();

   TGButton *btn = 0;

   if (!fSelected) {
      btn = fToolBar->GetButton(kCompactAct);
      if (btn) btn->SetState(kButtonDisabled);

      btn = fToolBar->GetButton(kLayoutVAct);
      if (btn) btn->SetState(kButtonDisabled);

      btn = fToolBar->GetButton(kLayoutHAct);
      if (btn) btn->SetState(kButtonDisabled);

      btn = fToolBar->GetButton(kBreakLayoutAct);
      if (btn) btn->SetState(kButtonDisabled);
      return;
   }

   Bool_t comp = kFALSE;
   TGLayoutManager *lm = 0;
   Bool_t hor = kFALSE;
   Bool_t fixed = kFALSE;
   Bool_t enable = on;
   Bool_t compact_disable = kTRUE;

   if (fSelected->InheritsFrom(TGCompositeFrame::Class())) {
      lm = ((TGCompositeFrame*)fSelected)->GetLayoutManager();
      comp = kTRUE;
      hor = lm && lm->InheritsFrom(TGHorizontalLayout::Class());
      fixed = !fManager->CanChangeLayout(fSelected);
      compact_disable = !fManager->CanCompact(fSelected);
   } else {
      enable = kFALSE;
   }

   btn = fToolBar->GetButton(kCompactAct);
   if (btn) btn->SetState(enable && comp && !fixed && !compact_disable ?
                          kButtonUp : kButtonDisabled);

   btn = fToolBar->GetButton(kLayoutHAct);
   if (btn) {
      btn->SetState(enable && comp && !hor && !fixed ? kButtonUp :
                    kButtonDisabled);
   }

   btn = fToolBar->GetButton(kLayoutVAct);
   if (btn) {
      btn->SetState(enable && comp && hor && !fixed ? kButtonUp :
                    kButtonDisabled);
   }

   btn = fToolBar->GetButton(kBreakLayoutAct);
   if (btn) {
      btn->SetState(enable && comp && !fixed ? kButtonUp : kButtonDisabled);
   }
/*
   btn = fToolBar->GetButton(kGrabAct);
   if (btn) {
      btn->SetState(enable && comp ? kButtonDown : kButtonUp);
      TGToolTip *tt = btn->GetToolTip();
      tt->SetText(btn->IsDown() ? "Drop Frames (Ctrl-Return)" :
                                  "Grab Selected Frames (Return)");
   }
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable toolbar buttons according to the selected frame.

void TRootGuiBuilder::EnableEditButtons(Bool_t on)
{
   TGButton *btn = 0;

   Bool_t lasso = fManager->IsLassoDrawn() && on;

   btn = fToolBar->GetButton(kReplaceAct);
   if (btn) {
      btn->SetState(!on ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kCutAct);
   if (btn) {
      btn->SetState(!on || lasso ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kDropAct);
   if (btn) {
      btn->SetState(!on || lasso ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kCopyAct);
   if (btn) {
      btn->SetState(!on || lasso ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kPasteAct);
   if (btn) {
      btn->SetState(!on || !fManager->IsPasteFrameExist() ?
                    kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kCropAct);
   if (btn) {
      btn->SetState(!on && !lasso ? kButtonDisabled : kButtonUp);
   }

   btn = fToolBar->GetButton(kDeleteAct);
   if (btn) {
      btn->SetState(!on && !lasso ? kButtonDisabled : kButtonUp);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update gui builder.

void TRootGuiBuilder::Update()
{
   if (!fManager) {
      return;
   }

   EnableLassoButtons(fManager->IsLassoDrawn());
   fSelected = fManager->GetSelected();
   EnableSelectedButtons(fSelected);
   EnableEditButtons(fClient->IsEditable() && (fManager->IsLassoDrawn() ||
                     fManager->GetSelected() ||
                     fManager->IsPasteFrameExist()));

   if (fActionButton) {
      TGFrame *parent = (TGFrame*)fActionButton->GetParent();
      parent->ChangeOptions(parent->GetOptions() & ~kSunkenFrame);
      fClient->NeedRedraw(parent, kTRUE);
   }

   if (!fClient->IsEditable()) {
      UpdateStatusBar("");
      fMenuFile->EnableEntry(kGUIBLD_FILE_START);
      fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
      fEditable = 0;
      //fShutter->SetSelectedItem(fShutter->GetItem("Projects"));
   } else {
      fMenuFile->DisableEntry(kGUIBLD_FILE_START);
      fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
   }

   SwitchToolbarButton();
   fActionButton = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the status of the selected mode.

Bool_t TRootGuiBuilder::IsSelectMode() const
{
   TGButton *btn = 0;
   btn = fToolBar->GetButton(kSelectAct);

   if (!btn) return kFALSE;

   return btn->IsDown();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the status of the grab mode.

Bool_t TRootGuiBuilder::IsGrabButtonDown() const
{
   TGButton *btn = fToolBar->GetButton(kGrabAct);

   if (!btn) return kFALSE;

   return btn->IsDown();
}

////////////////////////////////////////////////////////////////////////////////
class TGuiBldSaveFrame : public TGMainFrame {

public:
   TGuiBldSaveFrame(const TGWindow *p, UInt_t w , UInt_t h) :
      TGMainFrame(p, w, h) {}
   void SetList(TList *li) { fList = li; }
};

static const char *gSaveMacroTypes[] = {
   "Macro files", "*.[C|c]*",
   "All files",   "*",
   0,             0
};

////////////////////////////////////////////////////////////////////////////////
/// Handle keys.

Bool_t TRootGuiBuilder::HandleKey(Event_t *event)
{
   if (event->fType == kGKeyPress) {
      UInt_t keysym;
      char str[2];
      gVirtualX->LookupString(event, str, sizeof(str), keysym);

      if (event->fState & kKeyControlMask) {
         if (str[0] == 19) {  // ctrl-s
            if (fMain->GetCurrent()) {
               return SaveProject(event);
            } else {
               return kFALSE; //TGMainFrame::HandleKey(event);
            }
         } else if (str[0] == 14) { //ctrl-n
            return NewProject();  //event not needed
         } else if (str[0] == 15) { // ctrl-o
            return OpenProject(event);
         }
      }
      fManager->HandleKey(event);
      return TGMainFrame::HandleKey(event);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new project.

Bool_t TRootGuiBuilder::NewProject(TString type)
{
   TGWindow *root = (TGWindow*)fClient->GetRoot();

   if (root) root->SetEditable(kFALSE);
   fEditable = new TGMdiFrame(fMain, 500, 400, kOwnBackground);
   fEditable->DontCallClose();
   fEditable->SetMdiHints(kMdiDefaultHints);
   fEditable->SetWindowName(fEditable->GetName());
   fEditable->SetEditDisabled(0);   // enable editting
   fEditable->MapRaised();
   fEditable->AddInput(kKeyPressMask | kButtonPressMask);
   fEditable->SetEditable(kTRUE);

   if (type == "horizontal") {
      TGHorizontalFrame *hor = new TGHorizontalFrame(fEditable, 100, 100);
      fEditable->AddFrame(hor, new TGLayoutHints( kLHintsExpandX |
                          kLHintsExpandY, 1, 1, 1, 1));
      hor->SetEditable(kTRUE);
      fClient->NeedRedraw(hor, kTRUE);
      fEditable->MapSubwindows();
      fEditable->MapWindow();
      fClient->NeedRedraw(fEditable, kTRUE);
      fEditable->SetLayoutBroken(kFALSE);
      fEditable->Layout();
   }
   else if (type == "vertical") {
      TGVerticalFrame *vert = new TGVerticalFrame(fEditable, 100, 100);
      fEditable->AddFrame(vert, new TGLayoutHints( kLHintsExpandX |
                          kLHintsExpandY,1,1,1,1));
      vert->SetEditable(kTRUE);
      fClient->NeedRedraw(vert, kTRUE);
      fEditable->MapSubwindows();
      fEditable->MapWindow();
      fClient->NeedRedraw(fEditable, kTRUE);
      fEditable->SetLayoutBroken(kFALSE);
      fEditable->Layout();

   } else {
      fEditable->SetLayoutBroken(kTRUE);
   }
   fManager->SetEditable(kTRUE);
   fMenuFile->EnableEntry(kGUIBLD_FILE_CLOSE);
   fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
   fEditable->SetCleanup(kDeepCleanup);


   SwitchToolbarButton();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open new gui builder project.

Bool_t TRootGuiBuilder::OpenProject(Event_t *event)
{

   TGButton *btn = fToolBar->GetButton(kOpenAct);
   if (btn) {
      btn->SetBackgroundColor(GetBgnd());
      fClient->NeedRedraw(btn, kTRUE);
   }

   TGFileInfo fi;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TString fname;

   fi.fFileTypes = gSaveMacroTypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;
   TGWindow *root = (TGWindow*)fClient->GetRoot();
   root->SetEditable(kFALSE);

   new TGFileDialog(fClient->GetDefaultRoot(), this, kFDOpen, &fi);

   if (!fi.fFilename) {
      root->SetEditable(kTRUE);
      return kFALSE;
   }

   dir    = fi.fIniDir;
   overwr = fi.fOverwrite;
   fname  = fi.fFilename;

   if (fname.EndsWith(".C", TString::kIgnoreCase) || fname.EndsWith(".cxx") ||
       fname.EndsWith(".cpp") || fname.EndsWith(".cc")) {
      NewProject();        // create new project
      gROOT->Macro(fname.Data()); // put content of the macro as child frame
   } else {
      Int_t retval;
      new TGMsgBox(fClient->GetDefaultRoot(), this, "Error...",
                   TString::Format("file (%s) must have source extension (.C, .c, .cxx, .cpp, .cc)",
                   fname.Data()), kMBIconExclamation, kMBRetry | kMBCancel,
                   &retval);

      if (retval == kMBRetry) {
         OpenProject(event);
      }
   }

   fMenuFile->EnableEntry(kGUIBLD_FILE_CLOSE);
   fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
   SwitchToolbarButton();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save selected project.

Bool_t TRootGuiBuilder::SaveProject(Event_t *event)
{
   TGButton *btn = fToolBar->GetButton(kSaveAct);
   if (btn) {
      btn->SetBackgroundColor(GetBgnd());
      fClient->NeedRedraw(btn, kTRUE);
   }

   TGMdiFrame *savfr = fMain->GetCurrent();
   if (!savfr) return kFALSE;

   static TImage *img = 0;

   if (!img) {
      img = TImage::Create();
   }
   img->FromWindow(savfr->GetParent()->GetId());

   TGWindow *root = (TGWindow*)fClient->GetRoot();
   TGFileInfo fi;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TString fname;
   root->SetEditable(kFALSE);

   fi.fFileTypes = gSaveMacroTypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;

   new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);

   if (!fi.fFilename) {
      root->SetEditable(kTRUE);
      SetEditable(kTRUE);
      return kFALSE;
   }

   dir = fi.fIniDir;
   overwr = fi.fOverwrite;
   fname = gSystem->UnixPathName(fi.fFilename);

   if (fname.EndsWith(".C", TString::kIgnoreCase) || fname.EndsWith(".cxx") ||
       fname.EndsWith(".cpp") || fname.EndsWith(".cc")) {
      TGuiBldSaveFrame *main = new TGuiBldSaveFrame(fClient->GetDefaultRoot(),
                                                    savfr->GetWidth(),
                                                    savfr->GetHeight());
      TList *list = main->GetList();
      TString name = savfr->GetName();
      savfr->SetName(main->GetName());
      main->SetList(savfr->GetList());
      main->SetLayoutBroken(savfr->IsLayoutBroken());
      main->SaveSource(fname.Data(), "keep_names");
      savfr->SetWindowName(fname.Data());
      main->SetList(list);

      main->SetMWMHints(kMWMDecorAll, kMWMFuncAll,
                        kMWMInputFullApplicationModal);
      main->SetWMSize(main->GetWidth(), main->GetHeight());
      main->SetWMSizeHints(main->GetDefaultWidth(), main->GetDefaultHeight(),
                           10000, 10000, 0, 0);
      main->SetWindowName(fname.Data());
      main->SetIconName(fname.Data());
      main->SetClassHints(fname.Data(), fname.Data());
      // some problems here under win32
      if (gVirtualX->InheritsFrom("TGX11")) main->SetIconPixmap("bld_rgb.xpm");

      savfr->SetName(name.Data());

      AddMacro(fname.Data(), img);
      delete main;
   } else {
      Int_t retval;
      new TGMsgBox(fClient->GetDefaultRoot(), this, "Error...",
                   TString::Format("file (%s) must have source extension (.C, .c, .cxx, .cpp, .cc)",
                   fname.Data()), kMBIconExclamation, kMBRetry | kMBCancel,
                   &retval);
      if (retval == kMBRetry) {
         SaveProject(event);
      }
      SwitchToolbarButton();
   }
   SwitchToolbarButton();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add macro to "User's Macros" section
/// Input img must be static - do not delete it.

void TRootGuiBuilder::AddMacro(const char *macro, TImage *img)
{
   if (!img || !img->GetWidth() || !img->GetHeight()) {
      return;
   }

   UInt_t w = 100;
   Float_t ratio = Float_t(w)/img->GetWidth();
   Float_t rh = img->GetHeight()*ratio;
   UInt_t h = UInt_t(rh);
   img->Scale(w, h);
   img->Merge(img, "overlay");

   static int i = 0;
   const TGPicture *pic = fClient->GetPicturePool()->GetPicture(
                                       TString::Format("%s;%d", macro, i++),
                                       img->GetPixmap(),img->GetMask());
   const char *name = gSystem->BaseName(macro);

   TGButton *btn = FindActionButton(name, "User's Macros");
   TGuiBldAction *act = 0;

   if (!btn) {
      act = new TGuiBldAction(name, macro, kGuiBldMacro);
      act->fAct = macro;
      act->fPic = macro;
      act->fPicture = pic;

      AddAction(act, "User's Macros");
   } else {
      act = (TGuiBldAction*)btn->GetUserData();
      act->fAct = macro;
      act->fPic = macro;
      act->fPicture = pic;

      if (btn->InheritsFrom(TGPictureButton::Class())) {
         btn->Resize(w, h);
         fClient->FreePicture(((TGPictureButton*)btn)->GetPicture());
         ((TGPictureButton*)btn)->SetPicture(pic);
      }
   }
   fClient->NeedRedraw(fShutter);
}

////////////////////////////////////////////////////////////////////////////////
/// Find the editable frame.

TGMdiFrame *TRootGuiBuilder::FindEditableMdiFrame(const TGWindow *win)
{
   const TGWindow *parent = win;
   TGMdiFrame *ret = 0;

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (parent->InheritsFrom(TGMdiFrame::Class())) {
         ret = (TGMdiFrame*)parent;
         return ret;
      }
      parent = parent->GetParent();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Switch image of toolbar edit button according to the current state.

void TRootGuiBuilder::SwitchToolbarButton()
{
   static const TGPicture *start = fClient->GetPicture("bld_edit.png");
   static const TGPicture *stop = fClient->GetPicture("bld_stop.png");

   if (fClient->IsEditable()) {
      fStartButton->SetEnabled(kTRUE);
      fStartButton->SetPicture(stop);
      fToolBar->SetId(fStartButton, kEndEditAct);
      fStartButton->SetToolTipText("Stop Edit (Ctrl-Dbl-Click)");
   } else {
      if (fMain->GetNumberOfFrames() < 1) {
         fStartButton->SetEnabled(kFALSE);
      } else {
         fStartButton->SetEnabled(kTRUE);
         fStartButton->SetPicture(start);
         fToolBar->SetId(fStartButton, kEditableAct);
         fStartButton->SetToolTipText("Start Edit (Ctrl-Dbl-Click)");
      }
   }

   fClient->NeedRedraw(fStartButton, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu items.

void TRootGuiBuilder::HandleMenu(Int_t id)
{
   TGWindow *root = (TGWindow*)fClient->GetRoot();
   TRootHelpDialog *hd;

   switch (id) {
      case kGUIBLD_FILE_START:
         if (fClient->IsEditable()) {
            break;
         }
         fEditable = fMain->GetCurrent();
         if (fEditable) {
            fEditable->SetEditable(kTRUE);
         } //else if (!fMain->GetCurrent()) {
            //NewProject();
         //}
         UpdateStatusBar("Start edit");
         fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
         fMenuFile->DisableEntry(kGUIBLD_FILE_START);
         SwitchToolbarButton();
         break;

      case kGUIBLD_FILE_STOP:
         if (!fClient->IsEditable()) {
            break;
         }
         fEditable = FindEditableMdiFrame(root);

         if (fEditable) {
            root->SetEditable(kFALSE);

            UpdateStatusBar("Stop edit");
            fMenuFile->EnableEntry(kGUIBLD_FILE_START);
            fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
            fEditable = 0;
            SwitchToolbarButton();
         }
         fEditor->Reset();
         break;

      case kGUIBLD_FILE_NEW:
         NewProject();
         SwitchToolbarButton();
         break;

      case kGUIBLD_FILE_CLOSE:
         fEditable = FindEditableMdiFrame(root);
         if (fEditable && (fEditable == fMain->GetCurrent())) {
            root->SetEditable(kFALSE);
         }
         fEditor->Reset();
         UpdateStatusBar("");
         fMain->Close(fMain->GetCurrent());

         if (fMain->GetNumberOfFrames() <= 1) {
            fMenuFile->DisableEntry(kGUIBLD_FILE_CLOSE);
            fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
            fMenuFile->DisableEntry(kGUIBLD_FILE_START);
         }

         if (fClient->IsEditable()) {
            fMenuFile->DisableEntry(kGUIBLD_FILE_START);
            fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
         } else {
            fMenuFile->EnableEntry(kGUIBLD_FILE_START);
            fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
         }
         SwitchToolbarButton();
         break;

      case kGUIBLD_FILE_EXIT:
         CloseWindow();
         break;

      case kGUIBLD_FILE_OPEN:
         OpenProject();
         SwitchToolbarButton();
         break;

      case kGUIBLD_FILE_SAVE:
         SaveProject();
         SwitchToolbarButton();
         break;

      case kGUIBLD_WINDOW_HOR:
         fMain->TileHorizontal();
         break;

      case kGUIBLD_WINDOW_VERT:
         fMain->TileVertical();
         break;

      case kGUIBLD_WINDOW_CASCADE:
         fMain->Cascade();
         break;

      case kGUIBLD_WINDOW_ARRANGE:
         fMain->ArrangeMinimized();
         break;

      case kGUIBLD_WINDOW_OPAQUE:
         if (fMenuWindow->IsEntryChecked(kGUIBLD_WINDOW_OPAQUE)) {
            fMenuWindow->UnCheckEntry(kGUIBLD_WINDOW_OPAQUE);
            fMain->SetResizeMode(kMdiNonOpaque);
         } else {
            fMenuWindow->CheckEntry(kGUIBLD_WINDOW_OPAQUE);
            fMain->SetResizeMode(kMdiOpaque);
         }
         break;
      case  kGUIBLD_HELP_CONTENTS:
         root->SetEditable(kFALSE);
         hd = new TRootHelpDialog(this, "Help on Gui Builder...", 600, 400);
         hd->SetText(gHelpBuilder);
         hd->SetEditDisabled();
         hd->Popup();
         root->SetEditable(kTRUE);
         break;

      case  kGUIBLD_HELP_ABOUT:
         root->SetEditable(kFALSE);
         hd = new TRootHelpDialog(this, "About Gui Builder...", 520, 160);
         hd->SetEditDisabled();
         hd->SetText(gHelpAboutBuilder);
         hd->Popup();
         root->SetEditable(kTRUE);
         break;

      default:
         fMain->SetCurrent(id);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handler before closing MDI frame.

void TRootGuiBuilder::MaybeCloseWindow()
{
   int retval;
   if (fClosing == -1)
      return;
   TGMdiFrame *mdiframe = (TGMdiFrame *)gTQSender;
   fManager->SetEditable(kFALSE);
   new TGMsgBox(gClient->GetDefaultRoot(), this,
                "Closing project", "Do you want to save the project before closing?",
                kMBIconExclamation, kMBYes | kMBNo | kMBCancel, &retval);

   fManager->SetEditable(kTRUE);
   if (retval == kMBYes) {
      SaveProject();
   }
   if (retval == kMBCancel) {
      fClosing = -1;
      if (!fClient->IsEditable())
         HandleMenu(kGUIBLD_FILE_START);
      return;
   }
   fEditor->RemoveFrame(mdiframe);
   mdiframe->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Handler for closed MDI frame.

void TRootGuiBuilder::HandleWindowClosed(Int_t )
{
   fEditable = 0;

   if (fClient->IsEditable()) {
      fManager->SetEditable(kFALSE);
      fMenuFile->DisableEntry(kGUIBLD_FILE_START);
      fMenuFile->EnableEntry(kGUIBLD_FILE_STOP);
   } else {
      fMenuFile->EnableEntry(kGUIBLD_FILE_START);
      fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
   }
   fEditor->Reset();
   UpdateStatusBar("");

   if (fMain->GetNumberOfFrames() == 0) {
      fMenuFile->DisableEntry(kGUIBLD_FILE_CLOSE);
      fMenuFile->DisableEntry(kGUIBLD_FILE_STOP);
      fMenuFile->DisableEntry(kGUIBLD_FILE_START);
      SwitchToolbarButton();
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update information shown on the status bar.

void TRootGuiBuilder::UpdateStatusBar(const char *txt)
{
   if (!fStatusBar) return;

   const char *text = 0;

   if (!txt) {
      TObject *o = (TObject *)gTQSender;

      if (o && o->InheritsFrom(TGToolTip::Class())) {
         TGToolTip *tip = (TGToolTip*)o;
         text = tip->GetText()->Data();
      }
   } else {
      text = txt;
   }
   fStatusBar->SetText(text);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear information shown in the status bar.

void TRootGuiBuilder::EraseStatusBar()
{
   if (!fStatusBar) return;

   fStatusBar->SetText("");
}

////////////////////////////////////////////////////////////////////////////////
/// Keyboard key binding.

void TRootGuiBuilder::BindKeys()
{
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_a),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_n),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_o),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Return),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Return),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Enter),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Enter),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_x),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_c),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_v),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_r),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_z),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_z),
                      kKeyControlMask | kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_b),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_l),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_g),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_h),
                      kKeyControlMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Delete),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Backspace),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Space),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Left),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Right),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Up),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Down),
                      0, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Left),
                      kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Right),
                      kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Up),
                      kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Down),
                      kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Delete),
                      kKeyShiftMask, kTRUE);

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Backspace),
                      kKeyShiftMask, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Create new vertical splitter (TGVSplitter).

TGFrame *TRootGuiBuilder::VSplitter()
{
   TGHorizontalFrame *ret = new TGHorizontalFrame();
   ret->SetCleanup(kDeepCleanup);
   TGVerticalFrame *v1 = new TGVerticalFrame(ret, 40, 10, kSunkenFrame |
                                             kFixedWidth);
   ret->AddFrame(v1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));
   //v1->SetEditDisabled(kEditDisableGrab);

   TGVSplitter *splitter = new TGVSplitter(ret);
   splitter->SetFrame(v1, kTRUE);
   ret->AddFrame(splitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));
   splitter->SetEditDisabled(kEditDisableBtnEnable);

   TGVerticalFrame *v2 = new TGVerticalFrame(ret, 10, 10, kSunkenFrame);
   v2->ChangeOptions(kSunkenFrame);
   ret->AddFrame(v2, new TGLayoutHints(kLHintsRight | kLHintsExpandX |
                 kLHintsExpandY));
   //v2->SetEditDisabled(kEditDisableGrab);
   ret->SetEditDisabled(kEditDisableLayout);

   ret->MapSubwindows();
   ret->SetLayoutBroken(kFALSE);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
///  Creates new horizontal splitter (TGHSplitter).

TGFrame *TRootGuiBuilder::HSplitter()
{
   TGVerticalFrame *ret = new TGVerticalFrame();
   ret->SetCleanup(kDeepCleanup);
   TGHorizontalFrame *v1 = new TGHorizontalFrame(ret, 10, 40, kSunkenFrame |
                                                 kFixedHeight);
   ret->AddFrame(v1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
   //v1->SetEditDisabled(kEditDisableGrab);

   TGHSplitter *splitter = new TGHSplitter(ret);
   splitter->SetFrame(v1, kTRUE);
   ret->AddFrame(splitter, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
   splitter->SetEditDisabled(kEditDisable);

   TGHorizontalFrame *v2 = new TGHorizontalFrame(ret, 10, 10);
   v2->ChangeOptions(kSunkenFrame);
   ret->AddFrame(v2, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                 kLHintsExpandY));
   //v2->SetEditDisabled(kEditDisableGrab);
   ret->SetEditDisabled(kEditDisableLayout);

   ret->MapSubwindows();
   ret->SetLayoutBroken(kFALSE);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Hide builder.

void TRootGuiBuilder::Hide()
{
   //fMain->CloseAll();
   UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default background color.

ULong_t TRootGuiBuilder::GetBgnd()
{
   return GetDefaultFrameBackground();

   static ULong_t gPixel = 0;

   if (gPixel) return gPixel;

   Float_t r, g, b;

   r = 232./255;
   g = 232./255;
   b = 222./255;

   gPixel = TColor::RGB2Pixel(r, g, b);
   return gPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Return background GC.

TGGC *TRootGuiBuilder::GetBgndGC()
{
   if (fgBgnd) return fgBgnd;

   fgBgnd = new TGGC(TGFrame::GetBckgndGC());

   Pixel_t back = GetBgnd();
   fgBgnd->SetBackground(back);
   fgBgnd->SetForeground(back);

   return fgBgnd;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a background color to frame and all its subframes.

void TRootGuiBuilder::PropagateBgndColor(TGFrame *frame, Pixel_t color)
{
   if (!frame) return;

   frame->SetBackgroundColor(color);
   if (!frame->InheritsFrom(TGCompositeFrame::Class())) return;

   TIter next(((TGCompositeFrame*)frame)->GetList());
   TGFrameElement *fe;

   while ((fe = (TGFrameElement*)next())) {
      if (fe->fFrame->GetBackground() == TGFrame::GetWhitePixel()) {
         continue;
      }
      PropagateBgndColor(fe->fFrame, color);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return background color for popup menus.

ULong_t TRootGuiBuilder::GetPopupBgnd()
{
   return GetDefaultFrameBackground();

   static ULong_t gPixel = 0;

   if (gPixel) return gPixel;

   Float_t r, g, b;

   r = 250./255;
   g = 250./255;
   b = 250./255;

   gPixel = TColor::RGB2Pixel(r, g, b);

   return gPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Return background GC for popup menus.

TGGC *TRootGuiBuilder::GetPopupBgndGC()
{
   if (fgBgndPopup) return fgBgndPopup;

   fgBgndPopup = new TGGC(TGFrame::GetBckgndGC());

   Pixel_t back = GetPopupBgnd();
   fgBgndPopup->SetBackground(back);
   fgBgndPopup->SetForeground(back);

   return fgBgndPopup;
}

////////////////////////////////////////////////////////////////////////////////
/// Return highlighted color for popup menu entry.

ULong_t TRootGuiBuilder::GetPopupHlght()
{
   return GetDefaultSelectedBackground();

   static ULong_t gPixel = 0;

   if (gPixel) return gPixel;

   Float_t r, g, b;

   r = 120./255;
   g = 120./255;
   b = 222./255;

   gPixel = TColor::RGB2Pixel(r, g, b);

   return gPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Return background GC for highlighted popup menu entry.

TGGC *TRootGuiBuilder::GetPopupHlghtGC()
{
   if (fgBgndPopupHlght) return fgBgndPopupHlght;

   fgBgndPopupHlght = new TGGC(TGFrame::GetHilightGC());

   Pixel_t back = GetPopupHlght();
   fgBgndPopupHlght->SetBackground(back);
   fgBgndPopupHlght->SetForeground(back);

   return fgBgndPopupHlght;
}

////////////////////////////////////////////////////////////////////////////////
/// Return style popup menu.

TGPopupMenu *TRootGuiBuilder::CreatePopup()
{
   return new TGuiBldPopupMenu();
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method used in guibuilding

TGFrame *TRootGuiBuilder::BuildListTree()
{
   TGCanvas *canvas = new TGCanvas(gClient->GetRoot(), 100, 100);
   TGListTree *lt = new TGListTree(canvas, kHorizontalFrame);
   lt->AddItem(0, "Entry 1");
   lt->AddItem(0, "Entry 2");
   lt->AddItem(0, "Entry 3");
   lt->AddItem(0, "Entry 4");
   lt->AddItem(0, "Entry 5");
   canvas->Resize(100, 60);
   canvas->MapSubwindows();

   return canvas;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method used in guibuilding to create TGCanvas widget

TGFrame *TRootGuiBuilder::BuildCanvas()
{
   TGCanvas *canvas = new TGCanvas(gClient->GetRoot(), 100, 100);
   TGCompositeFrame *cont = new TGCompositeFrame(canvas->GetViewPort(),
                                                 200, 200, kHorizontalFrame |
                                                 kOwnBackground);

   cont->SetCleanup(kDeepCleanup);
   cont->SetLayoutManager(new TGTileLayout(cont, 8));
   cont->AddFrame(new TGTextButton(cont, "Button1"));
   cont->AddFrame(new TGTextButton(cont, "Button2"));
   cont->AddFrame(new TGTextButton(cont, "Button3"));
   cont->AddFrame(new TGTextButton(cont, "Button4"));

   canvas->SetContainer(cont);
   return canvas;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method used in guibuilding to create TGShutter widget

TGFrame *TRootGuiBuilder::BuildShutter()
{
   TGShutterItem *item;
   TGCompositeFrame *container;
   const TGPicture  *buttonpic;
   TGPictureButton  *button;

   TGLayoutHints *l = new TGLayoutHints(kLHintsTop | kLHintsCenterX,5,5,5,0);
   TGShutter *shut = new TGShutter();

   item = shut->AddPage("Histograms");
   container = (TGCompositeFrame *)item->GetContainer();
   buttonpic = gClient->GetPicture("h1_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TH1");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("h2_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TH2");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("h3_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TH3");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("profile_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TProfile");
      container->AddFrame(button, l);
   }

   // new page
   item = shut->AddPage("Functions");
   container = (TGCompositeFrame *)item->GetContainer();
   buttonpic = gClient->GetPicture("f1_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TF1");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("f2_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TF2");
      container->AddFrame(button, l);
   }

   // new page
   item = shut->AddPage("Trees");
   container = (TGCompositeFrame *)item->GetContainer();
   buttonpic = gClient->GetPicture("ntuple_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TNtuple");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("tree_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TTree");
      container->AddFrame(button, l);
   }
   buttonpic = gClient->GetPicture("chain_s.xpm");

   if (buttonpic) {
      button = new TGPictureButton(container, buttonpic);
      button->SetToolTipText("TChain");
      container->AddFrame(button, l);
   }

   shut->MapSubwindows();
   return shut;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGTextEdit widget

TGFrame *TRootGuiBuilder::BuildTextEdit()
{
   TGTextEdit *te = new TGTextEdit();

   te->AddLine("all work and no play makes jack a pretty");
   te->AddLine("dull boy. all work and no play makes jack");
   te->AddLine("a pretty dull boy. all work and no play ");
   te->AddLine("makes jack a pretty dull boy. all work");
   te->AddLine("and no play makes jack a pretty dull boy.");

   te->MapSubwindows();
   te->Layout();
   te->Resize(100, 60);

   return te;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGTab widget

TGFrame *TRootGuiBuilder::BuildTab()
{
   TGTab *tab = new TGTab();

   tab->AddTab("Tab1");
   tab->AddTab("Tab2");
   tab->MapSubwindows();

   return tab;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGListBox widget

TGFrame *TRootGuiBuilder::BuildListBox()
{
   TGListBox *lb = new TGListBox();

   lb->AddEntry("Entry 1", 0);
   lb->AddEntry("Entry 2", 1);
   lb->AddEntry("Entry 3", 2);
   lb->AddEntry("Entry 4", 3);
   lb->AddEntry("Entry 5", 4);
   lb->AddEntry("Entry 6", 5);
   lb->AddEntry("Entry 7", 6);
   lb->MapSubwindows();

   lb->Resize(100,100);

   return lb;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGComboBox widget

TGFrame *TRootGuiBuilder::BuildComboBox()
{
   TGComboBox *cb = new TGComboBox();

   cb->AddEntry("Entry 1 ", 0);
   cb->AddEntry("Entry 2 ", 1);
   cb->AddEntry("Entry 3 ", 2);
   cb->AddEntry("Entry 4 ", 3);
   cb->AddEntry("Entry 5 ", 4);
   cb->AddEntry("Entry 6 ", 5);
   cb->AddEntry("Entry 7 ", 6);
   cb->MapSubwindows();

   FontStruct_t fs = TGTextLBEntry::GetDefaultFontStruct();
   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(fs, max_ascent, max_descent);

   cb->Resize(cb->GetListBox()->GetDefaultWidth(), max_ascent+max_descent+7);
   return cb;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGHorizontal3DLine widget.

TGFrame *TRootGuiBuilder::BuildH3DLine()
{
   TGHorizontal3DLine *l = new TGHorizontal3DLine(0, 100, 2);
   l->Resize(100, 2);

   return l;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGVertical3DLine widget.

TGFrame *TRootGuiBuilder::BuildV3DLine()
{
   TGVertical3DLine *l = new TGVertical3DLine();
   l->Resize(2, 100);

   return l;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGHScrollBar

TGFrame *TRootGuiBuilder::BuildHScrollBar()
{
   TGHScrollBar *b = new TGHScrollBar();

   b->Resize(100, b->GetDefaultHeight());
   b->SetRange(100, 20);
   b->MapSubwindows();

   return b;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGVScrollBar

TGFrame *TRootGuiBuilder::BuildVScrollBar()
{
   TGVScrollBar *b = new TGVScrollBar();

   b->Resize(b->GetDefaultWidth(), 100);
   b->MapSubwindows();
   b->SetRange(100, 20);

   return b;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGHProgressBar

TGFrame *TRootGuiBuilder::BuildHProgressBar()
{
   TGHProgressBar *b = new TGHProgressBar();

   b->Resize(100, b->GetDefaultHeight());
   b->SetPosition(25);
   b->Format("%.2f");
   b->SetFillType(TGProgressBar::kBlockFill);

   return b;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to create TGVProgressBar

TGFrame *TRootGuiBuilder::BuildVProgressBar()
{
   TGVProgressBar *b = new TGVProgressBar();

   b->Resize(b->GetDefaultWidth(), 100);
   b->SetPosition(25);
   b->SetFillType(TGProgressBar::kBlockFill);

   return b;
}


