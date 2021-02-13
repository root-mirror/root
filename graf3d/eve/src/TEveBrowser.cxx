// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveBrowser.h"

#include "TEveUtil.h"
#include "TEveElement.h"
#include "TEveManager.h"
#include "TEveSelection.h"
#include "TEveGedEditor.h"
#include "TEveWindow.h"
#include "TEveWindowManager.h"

#include "TGFileBrowser.h"
#include "TBrowser.h"

#include "TClass.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVirtualX.h"

#include "TClassMenuItem.h"

#include "TColor.h"

#include "TGCanvas.h"
#include "TGMenu.h"
#include "TGSplitter.h"
#include "TGTab.h"
#include <KeySymbols.h>
#include "TGLSAViewer.h"

#include "TGeoVolume.h"
#include "TGeoNode.h"

/** \class TEveListTreeItem
\ingroup TEve
Special list-tree-item for Eve.

Most state is picked directly from TEveElement, no need to store it
locally nor to manage its consistency.

Handles also selected/highlighted colors and, in the future, drag-n-drop.
*/

ClassImp(TEveListTreeItem);

////////////////////////////////////////////////////////////////////////////////
/// Warn about access to function members that should never be called.
/// TGListTree calls them in cases that are not used by Eve.

void TEveListTreeItem::NotSupported(const char* func) const
{
   Warning(Form("TEveListTreeItem::%s()", func), "not supported.");
}

////////////////////////////////////////////////////////////////////////////////
/// Return highlight color corresponding to current state of TEveElement.

Pixel_t TEveListTreeItem::GetActiveColor() const
{
   switch (fElement->GetSelectedLevel())
   {
      case 1: return TColor::Number2Pixel(kBlue - 2);
      case 2: return TColor::Number2Pixel(kBlue - 6);
      case 3: return TColor::Number2Pixel(kCyan - 2);
      case 4: return TColor::Number2Pixel(kCyan - 6);
   }
   return TGFrame::GetDefaultSelectedBackground();
}

////////////////////////////////////////////////////////////////////////////////
/// Item's check-box state has been toggled ... forward to element's
/// render-state.

void TEveListTreeItem::Toggle()
{
   fElement->SetRnrState(!IsChecked());
   fElement->ElementChanged(kTRUE, kTRUE);
}

/** \class TEveGListTreeEditorFrame
\ingroup TEve
Composite GUI frame for parallel display of a TGListTree and TEveGedEditor.
*/

ClassImp(TEveGListTreeEditorFrame);

TString TEveGListTreeEditorFrame::fgEditorClass("TEveGedEditor");

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGListTreeEditorFrame::TEveGListTreeEditorFrame(const TGWindow* p, Int_t width, Int_t height) :
   TGMainFrame (p ? p : gClient->GetRoot(), width, height),
   fFrame      (0),
   fLTFrame    (0),
   fListTree   (0),
   fSplitter   (0),
   fEditor     (0),
   fCtxMenu    (0),
   fSignalsConnected (kFALSE)
{
   SetCleanup(kNoCleanup);

   fFrame = new TGCompositeFrame(this, width, height, kVerticalFrame);

   // List-tree
   fLTFrame  = new TGCompositeFrame(fFrame, width, 3*height/7, kVerticalFrame);
   fLTCanvas = new TGCanvas(fLTFrame, 10, 10, kSunkenFrame | kDoubleBorder);
   fListTree = new TGListTree(fLTCanvas->GetViewPort(), 10, 10, kHorizontalFrame);
   fListTree->SetCanvas(fLTCanvas);
   fListTree->Associate(fFrame);
   fListTree->SetColorMode(TGListTree::EColorMarkupMode(TGListTree::kColorUnderline | TGListTree::kColorBox));
   fListTree->SetAutoCheckBoxPic(kFALSE);
   fListTree->SetUserControl(kTRUE);
   fLTCanvas->SetContainer(fListTree);
   fLTFrame->AddFrame(fLTCanvas, new TGLayoutHints
                      (kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1));
   fFrame  ->AddFrame(fLTFrame, new TGLayoutHints
                      (kLHintsNormal | kLHintsExpandX | kLHintsExpandY));

   // Splitter
   fSplitter = new TGHSplitter(fFrame);
   fFrame->AddFrame(fSplitter, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 1,1,2,2));

   // Editor
   fFrame->SetEditDisabled(kEditEnable);
   fFrame->SetEditable();
   fEditor = (TEveGedEditor*) gROOT->GetClass(fgEditorClass)->New();
   fEditor->SetGlobal(kFALSE);
   fEditor->ChangeOptions(fEditor->GetOptions() | kFixedHeight);
   fFrame->SetEditable(kEditDisable);
   fFrame->SetEditable(kFALSE);
   {
      TGFrameElement *el = 0;
      TIter next(fFrame->GetList());
      while ((el = (TGFrameElement *) next())) {
         if (el->fFrame == fEditor)
            if (el->fLayout) {
               el->fLayout->SetLayoutHints(kLHintsTop | kLHintsExpandX);
               el->fLayout->SetPadLeft(0); el->fLayout->SetPadRight(1);
               el->fLayout->SetPadTop(2);  el->fLayout->SetPadBottom(1);
               break;
            }
      }
   }
   fSplitter->SetFrame(fEditor, kFALSE);

   AddFrame(fFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));

   fCtxMenu = new TContextMenu("", "");

   Layout();
   MapSubwindows();
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveGListTreeEditorFrame::~TEveGListTreeEditorFrame()
{
   DisconnectSignals();

   delete fCtxMenu;

   // Should un-register editor, all items and list-tree from gEve ... eventually.

   delete fEditor;
   delete fSplitter;
   delete fListTree;
   delete fLTCanvas;
   delete fLTFrame;
   delete fFrame;
}

////////////////////////////////////////////////////////////////////////////////
/// Set GED editor class.

void TEveGListTreeEditorFrame::SetEditorClass(const char* edclass)
{
   fgEditorClass = edclass;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect list-tree signals.

void TEveGListTreeEditorFrame::ConnectSignals()
{
   fListTree->Connect("MouseOver(TGListTreeItem*, UInt_t)", "TEveGListTreeEditorFrame",
                      this, "ItemBelowMouse(TGListTreeItem*, UInt_t)");
   fListTree->Connect("Clicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)", "TEveGListTreeEditorFrame",
                      this, "ItemClicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)");
   fListTree->Connect("DoubleClicked(TGListTreeItem*, Int_t)", "TEveGListTreeEditorFrame",
                      this, "ItemDblClicked(TGListTreeItem*, Int_t)");
   fListTree->Connect("KeyPressed(TGListTreeItem*, ULong_t, ULong_t)", "TEveGListTreeEditorFrame",
                      this, "ItemKeyPress(TGListTreeItem*, UInt_t, UInt_t)");

   fSignalsConnected = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect list-tree signals.

void TEveGListTreeEditorFrame::DisconnectSignals()
{
   if (!fSignalsConnected) return;

   fListTree->Disconnect("MouseOver(TGListTreeItem*, UInt_t)",
                      this, "ItemBelowMouse(TGListTreeItem*, UInt_t)");
   fListTree->Disconnect("Clicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)",
                      this, "ItemClicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)");
   fListTree->Disconnect("DoubleClicked(TGListTreeItem*, Int_t)",
                      this, "ItemDblClicked(TGListTreeItem*, Int_t)");
   fListTree->Disconnect("KeyPressed(TGListTreeItem*, ULong_t, ULong_t)",
                      this, "ItemKeyPress(TGListTreeItem*, UInt_t, UInt_t)");

   fSignalsConnected = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Reconfigure to horizontal layout, list-tree and editor side by side.

void TEveGListTreeEditorFrame::ReconfToHorizontal()
{
   UnmapWindow();

   fFrame->ChangeOptions(kHorizontalFrame);
   fLTFrame->ChangeOptions(kHorizontalFrame);
   fListTree->ChangeOptions(kVerticalFrame);

   TGFrameElement *el = 0;
   TIter next(fFrame->GetList());
   while ((el = (TGFrameElement *) next()))
   {
      if (el->fFrame == fSplitter)
      {
         // This is needed so that splitter window gets destroyed on server.
         fSplitter->ReparentWindow(fClient->GetDefaultRoot());
         delete fSplitter;
         el->fFrame = fSplitter = new TGVSplitter(fFrame);
         el->fLayout->SetLayoutHints(kLHintsLeft | kLHintsExpandY);
         el->fLayout->SetPadLeft(2); el->fLayout->SetPadRight (2);
         el->fLayout->SetPadTop (1); el->fLayout->SetPadBottom(1);
      }
      else if (el->fFrame == fEditor)
      {
         fEditor->ChangeOptions(fEditor->GetOptions() & (~kFixedHeight));
         fEditor->ChangeOptions(fEditor->GetOptions() |   kFixedWidth);
         el->fLayout->SetLayoutHints(kLHintsLeft | kLHintsExpandY);
      }
   }

   fEditor->Resize(fEditor->GetWidth() / 2 - 1, fEditor->GetHeight());
   fSplitter->SetFrame(fEditor, kFALSE);

   Layout();
   MapSubwindows();
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Reconfigure to vertical layout, list-tree above the editor.

void TEveGListTreeEditorFrame::ReconfToVertical()
{
   UnmapWindow();

   fFrame->ChangeOptions(kVerticalFrame);
   fLTFrame->ChangeOptions(kVerticalFrame);
   fListTree->ChangeOptions(kHorizontalFrame);

   TGFrameElement *el = 0;
   TIter next(fFrame->GetList());
   while ((el = (TGFrameElement *) next()))
   {
      if (el->fFrame == fSplitter)
      {
         // This is needed so that splitter window gets destroyed on server.
         fSplitter->ReparentWindow(fClient->GetDefaultRoot());
         delete fSplitter;
         el->fFrame = fSplitter = new TGHSplitter(fFrame);
         el->fLayout->SetLayoutHints(kLHintsTop | kLHintsExpandX);
         el->fLayout->SetPadLeft(2); el->fLayout->SetPadRight (2);
         el->fLayout->SetPadTop (1); el->fLayout->SetPadBottom(1);
      }
      else if (el->fFrame == fEditor)
      {
         fEditor->ChangeOptions(fEditor->GetOptions() & (~kFixedWidth));
         fEditor->ChangeOptions(fEditor->GetOptions() |   kFixedHeight);
         el->fLayout->SetLayoutHints(kLHintsTop | kLHintsExpandX);
      }
   }

   fEditor->Resize(fEditor->GetWidth(), fEditor->GetHeight() / 2 - 1);
   fSplitter->SetFrame(fEditor, kFALSE);

   Layout();
   MapSubwindows();
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Different item is below mouse.

void TEveGListTreeEditorFrame::ItemBelowMouse(TGListTreeItem *entry, UInt_t /*mask*/)
{
   TEveElement* el = entry ? (TEveElement*) entry->GetUserData() : 0;
   gEve->GetHighlight()->UserPickedElement(el, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Item has been clicked, based on mouse button do:
///  - M1 - select, show in editor;
///  - M2 - paste (call gEve->ElementPaste();
///  - M3 - popup context menu.

void TEveGListTreeEditorFrame::ItemClicked(TGListTreeItem *item, Int_t btn, UInt_t mask, Int_t x, Int_t y)
{
   //printf("ItemClicked item %s List %d btn=%d, x=%d, y=%d\n",
   //  item->GetText(),fDisplayFrame->GetList()->GetEntries(), btn, x, y);

   static const TEveException eh("TEveGListTreeEditorFrame::ItemClicked ");

   TEveElement* el = (TEveElement*) item->GetUserData();
   if (el == 0) return;
   TObject* obj = el->GetObject(eh);

   switch (btn)
   {
      case 1:
         gEve->GetSelection()->UserPickedElement(el, mask & kKeyControlMask);
         break;

      case 2:
         if (gEve->ElementPaste(el))
            gEve->Redraw3D();
         break;

      case 3:
         // If control pressed, show menu for render-element itself.
         // event->fState & kKeyControlMask
         // ??? how do i get current event?
         // !!!!! Have this now ... fix.
         if (obj) fCtxMenu->Popup(x, y, obj);
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Item has been double-clicked, potentially expand the children.

void TEveGListTreeEditorFrame::ItemDblClicked(TGListTreeItem* item, Int_t btn)
{
   static const TEveException eh("TEveGListTreeEditorFrame::ItemDblClicked ");

   if (btn != 1) return;

   TEveElement* el = (TEveElement*) item->GetUserData();
   if (el == 0) return;

   el->ExpandIntoListTree(fListTree, item);

   TObject* obj = el->GetObject(eh);
   if (obj)
   {
      // Browse geonodes.
      if (obj->IsA()->InheritsFrom(TGeoNode::Class()))
      {
         TGeoNode* n = dynamic_cast<TGeoNode*>(obj);
         if (item->GetFirstChild() == 0 && n->GetNdaughters())
         {
            fListTree->DeleteChildren(item);
            for (Int_t i=0; i< n->GetNdaughters(); i++)
            {
               TString title;
               title.Form("%d : %s[%d]", i,
                          n->GetDaughter(i)->GetVolume()->GetName(),
                          n->GetDaughter(i)->GetNdaughters());

               TGListTreeItem* child = fListTree->AddItem(item, title.Data());
               child->SetUserData(n->GetDaughter(i));
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// A key has been pressed for an item.
///
/// Only <Delete>, <Enter> and <Return> keys are handled here,
/// otherwise the control is passed back to TGListTree.

void TEveGListTreeEditorFrame::ItemKeyPress(TGListTreeItem *entry, UInt_t keysym, UInt_t mask)
{
   static const TEveException eh("TEveGListTreeEditorFrame::ItemKeyPress ");

   entry = fListTree->GetCurrent();
   if (entry == 0) return;

   TEveElement* el = (TEveElement*) entry->GetUserData();

   fListTree->SetEventHandled(); // Reset back to false in default case.

   switch (keysym)
   {
      case kKey_Delete:
      {
         if (entry->GetParent())
         {
            if (el->GetDenyDestroy() > 0 && el->GetNItems() == 1)
               throw(eh + "DestroyDenied set for this item.");

            TEveElement* parent = (TEveElement*) entry->GetParent()->GetUserData();

            if (parent)
            {
               gEve->RemoveElement(el, parent);
               gEve->Redraw3D();
            }
         }
         else
         {
            if (el->GetDenyDestroy() > 0)
               throw(eh + "DestroyDenied set for this top-level item.");
            gEve->RemoveFromListTree(el, fListTree, entry);
            gEve->Redraw3D();
         }
         break;
      }

      case kKey_Enter:
      case kKey_Return:
      {
         gEve->GetSelection()->UserPickedElement(el, mask & kKeyControlMask);
         break;
      }

      default:
      {
         fListTree->SetEventHandled(kFALSE);
         break;
      }
   }
}


/** \class TEveBrowser
\ingroup TEve
Specialization of TRootBrowser for Eve.
*/

ClassImp(TEveBrowser);

////////////////////////////////////////////////////////////////////////////////
/// Add "Export to CINT" into context-menu for class cl.

void TEveBrowser::SetupCintExport(TClass* cl)
{
   TList* l = cl->GetMenuList();
   TClassMenuItem* n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction, cl,
                                          "Export to CINT", "ExportToCINT", this, "const char*,TObject*", 1);

   l->AddFirst(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate position of a widget for reparenting into parent.

void TEveBrowser::CalculateReparentXY(TGObject* parent, Int_t& x, Int_t& y)
{
   UInt_t   w, h;
   Window_t childdum;
   gVirtualX->GetWindowSize(parent->GetId(), x, y, w, h);
   gVirtualX->TranslateCoordinates(parent->GetId(),
                                   gClient->GetDefaultRoot()->GetId(),
                                   0, 0, x, y, childdum);
}

namespace
{
enum EEveMenu_e {
   kNewMainFrameSlot, kNewTabSlot,
   kNewViewer,  kNewScene,
   kNewBrowser, kNewCanvas, kNewCanvasExt, kNewTextEditor, kNewHtmlBrowser,
   kSel_PS_Ignore, kSel_PS_Element, kSel_PS_Projectable, kSel_PS_Compound,
   kSel_PS_PableCompound, kSel_PS_Master, kSel_PS_END,
   kHil_PS_Ignore, kHil_PS_Element, kHil_PS_Projectable, kHil_PS_Compound,
   kHil_PS_PableCompound, kHil_PS_Master, kHil_PS_END,
   kVerticalBrowser,
   kWinDecorNormal, kWinDecorHide, kWinDecorTitleBar, kWinDecorMiniBar
};

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveBrowser::TEveBrowser(UInt_t w, UInt_t h) :
   TRootBrowser(0, "Eve Main Window", w, h, "", kFALSE),
   fFileBrowser(0),
   fEvePopup   (0),
   fSelPopup   (0),
   fHilPopup   (0)
{
   // Construct Eve menu.

   fEvePopup = new TGPopupMenu(gClient->GetRoot());
   fEvePopup->AddEntry("New &MainFrame Slot", kNewMainFrameSlot);
   fEvePopup->AddEntry("New &Tab Slot",       kNewTabSlot);
   fEvePopup->AddSeparator();
   fEvePopup->AddEntry("New &Viewer",         kNewViewer);
   fEvePopup->AddEntry("New &Scene",          kNewScene);
   fEvePopup->AddSeparator();
   fEvePopup->AddEntry("New &Browser",        kNewBrowser);
   fEvePopup->AddEntry("New &Canvas",         kNewCanvas);
   fEvePopup->AddEntry("New Canvas Ext",      kNewCanvasExt);
   fEvePopup->AddEntry("New Text &Editor",    kNewTextEditor);
   // fEvePopup->AddEntry("New HTML Browser", kNewHtmlBrowser);
   fEvePopup->AddSeparator();

   {
      fSelPopup = new TGPopupMenu(gClient->GetRoot());
      fSelPopup->AddEntry("Ignore",      kSel_PS_Ignore);
      fSelPopup->AddEntry("Element",     kSel_PS_Element);
      fSelPopup->AddEntry("Projectable", kSel_PS_Projectable);
      fSelPopup->AddEntry("Compound",    kSel_PS_Compound);
      fSelPopup->AddEntry("Projectable and Compound",
                          kSel_PS_PableCompound);
      fSelPopup->AddEntry("Master",      kSel_PS_Master);
      fSelPopup->RCheckEntry(kSel_PS_Ignore + gEve->GetSelection()->GetPickToSelect(),
                             kSel_PS_Ignore, kSel_PS_END - 1);
      fEvePopup->AddPopup("Selection", fSelPopup);
   }
   {
      fHilPopup = new TGPopupMenu(gClient->GetRoot());
      fHilPopup->AddEntry("Ignore",      kHil_PS_Ignore);
      fHilPopup->AddEntry("Element",     kHil_PS_Element);
      fHilPopup->AddEntry("Projectable", kHil_PS_Projectable);
      fHilPopup->AddEntry("Compound",    kHil_PS_Compound);
      fHilPopup->AddEntry("Projectable and Compound",
                          kHil_PS_PableCompound);
      fHilPopup->AddEntry("Master",      kHil_PS_Master);
      fHilPopup->RCheckEntry(kHil_PS_Ignore + gEve->GetHighlight()->GetPickToSelect(),
                             kHil_PS_Ignore, kHil_PS_END - 1);
      fEvePopup->AddPopup("Highlight", fHilPopup);
   }

   fEvePopup->AddSeparator();
   fEvePopup->AddEntry("Vertical browser", kVerticalBrowser);
   fEvePopup->CheckEntry(kVerticalBrowser);
   {
      TGPopupMenu *wd = new TGPopupMenu(gClient->GetRoot());
      wd->AddEntry("Normal",     kWinDecorNormal);
      wd->AddEntry("Hide",       kWinDecorHide);
      wd->AddEntry("Title bars", kWinDecorTitleBar);
      wd->AddEntry("Mini bars",  kWinDecorMiniBar);
      fEvePopup->AddPopup("Window decorations", wd);
   }

   fEvePopup->Connect("Activated(Int_t)", "TEveBrowser",
                       this, "EveMenu(Int_t)");

   fMenuBar->AddPopup("&Eve", fEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));

   fPreMenuFrame->ChangeOptions(fPreMenuFrame->GetOptions() | kRaisedFrame);
   fTopMenuFrame->Layout();
   fTopMenuFrame->MapSubwindows();

   // Rename "Close Window" to "Close Eve"
   fMenuFile->GetEntry(kCloseWindow)->GetLabel()->SetString("Close Eve");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle events from Eve menu.

void TEveBrowser::EveMenu(Int_t id)
{
   switch (id)
   {
      case kNewMainFrameSlot: {
         TEveWindowSlot* ew_slot = TEveWindow::CreateWindowMainFrame(0);
         gEve->GetWindowManager()->SelectWindow(ew_slot);
         break;
      }
      case kNewTabSlot: {
         TEveWindowSlot* ew_slot = TEveWindow::CreateWindowInTab(GetTabRight(), 0);
         gEve->GetWindowManager()->SelectWindow(ew_slot);
         break;
      }
      case kNewViewer: {
         gEve->SpawnNewViewer("Viewer Pepe");
         break;
      }
      case kNewScene: {
         gEve->SpawnNewScene("Scena Mica");
         break;
      }
      case kNewBrowser: {
         gROOT->ProcessLineFast("new TBrowser");
         break;
      }
      case kNewCanvas: {
         StartEmbedding(1);
         gROOT->ProcessLineFast("new TCanvas");
         StopEmbedding();
         SetTabTitle("Canvas", 1);
         break;
      }
      case kNewCanvasExt: {
         gROOT->ProcessLineFast("new TCanvas");
         break;
      }
      case kNewTextEditor: {
         StartEmbedding(1);
         gROOT->ProcessLineFast(Form("new TGTextEditor((const char *)0, (const TGWindow *)0x%lx)", (ULong_t)gClient->GetRoot()));
         StopEmbedding();
         SetTabTitle("Editor", 1);
         break;
      }
      case kNewHtmlBrowser: {
         gSystem->Load("libGuiHtml");
         if (gSystem->Load("libRHtml") >= 0)
         {
            StartEmbedding(1);
            gROOT->ProcessLine(Form("new TGHtmlBrowser(\"http://root.cern.ch/root/html/ClassIndex.html\", \
                              (const TGWindow *)0x%lx)", (ULong_t)gClient->GetRoot()));
            StopEmbedding();
            SetTabTitle("HTML", 1);
         }
         break;
      }
      case kSel_PS_Ignore:
      case kSel_PS_Element:
      case kSel_PS_Projectable:
      case kSel_PS_Compound:
      case kSel_PS_PableCompound:
      case kSel_PS_Master: {
         gEve->GetSelection()->SetPickToSelect(id - kSel_PS_Ignore);
         fSelPopup->RCheckEntry(kSel_PS_Ignore + gEve->GetSelection()->GetPickToSelect(),
                                kSel_PS_Ignore, kSel_PS_END - 1);
         break;
      }
      case kHil_PS_Ignore:
      case kHil_PS_Element:
      case kHil_PS_Projectable:
      case kHil_PS_Compound:
      case kHil_PS_PableCompound:
      case kHil_PS_Master: {
         gEve->GetHighlight()->SetPickToSelect(id - kHil_PS_Ignore);
         fHilPopup->RCheckEntry(kHil_PS_Ignore + gEve->GetHighlight()->GetPickToSelect(),
                                kHil_PS_Ignore, kHil_PS_END - 1);
         break;
      }
      case kVerticalBrowser: {
         if (fEvePopup->IsEntryChecked(kVerticalBrowser)) {
            gEve->GetLTEFrame()->ReconfToHorizontal();
            fEvePopup->UnCheckEntry(kVerticalBrowser);
         } else {
            gEve->GetLTEFrame()->ReconfToVertical();
            fEvePopup->CheckEntry(kVerticalBrowser);
         }
         break;
      }
      case kWinDecorNormal: {
         gEve->GetWindowManager()->ShowNormalEveDecorations();
         break;
      }
      case kWinDecorHide: {
         gEve->GetWindowManager()->HideAllEveDecorations();
         break;
      }
      case kWinDecorTitleBar: {
         gEve->GetWindowManager()->SetShowTitleBars(kTRUE);
         break;
      }
      case kWinDecorMiniBar: {
         gEve->GetWindowManager()->SetShowTitleBars(kFALSE);
         break;
      }

      default: {
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize standard plugins.

void TEveBrowser::InitPlugins(Option_t *opt)
{
   TString o(opt);

   // File Browser plugin ... we have to process it here.
   if (o.Contains('F'))
   {
      StartEmbedding(0);
      TGFileBrowser *fb = MakeFileBrowser();
      fb->BrowseObj(gROOT);
      fb->Show();
      fFileBrowser = fb;
      StopEmbedding("Files");
      o.ReplaceAll("F", ".");
   }

   TRootBrowser::InitPlugins(o);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a file-browser. Caller should provide Start/StopEmbedding() calls
/// and populate the new browser.
///
/// If flag make_default is kTRUE, the default file-browser is set to the
/// newly created browser.

TGFileBrowser* TEveBrowser::MakeFileBrowser(Bool_t make_default)
{
   TBrowserImp    imp;
   TBrowser      *tb = new TBrowser("Pipi", "Strel", &imp);
   TGFileBrowser *fb = new TGFileBrowser(gClient->GetRoot(), tb, 200, 500);
   tb->SetBrowserImp((TBrowserImp *)this);
   fb->SetBrowser(tb);
   fb->SetNewBrowser(this);
   gROOT->GetListOfBrowsers()->Remove(tb);
   // This guy is never used and stays in list-of-cleanups after destruction.
   // So let's just delete it now.
   delete tb->GetContextMenu();

   if (make_default)
      fFileBrowser = fb;

   return fb;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the default file-browser.

TGFileBrowser* TEveBrowser::GetFileBrowser() const
{
  return fFileBrowser;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the default file browser.

void TEveBrowser::SetFileBrowser(TGFileBrowser* b)
{
  fFileBrowser = b;
}

////////////////////////////////////////////////////////////////////////////////
/// Override from TRootBrowser. We need to be more brutal as fBrowser is
/// not set in Eve case.

void TEveBrowser::ReallyDelete()
{
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TRootBrowser. Need to intercept closing of Eve tabs.

void TEveBrowser::CloseTab(Int_t id)
{
   // Check if this is an Eve window and destroy accordingly.
   TGCompositeFrame *pcf = fTabRight->GetTabContainer(id);
   if (pcf)
   {
      TGFrameElement *fe = (TGFrameElement *) pcf->GetList()->First();
      if (fe)
      {
         TEveCompositeFrame *ecf = dynamic_cast<TEveCompositeFrame*>(fe->fFrame);
         if (ecf)
         {
            ecf->GetEveWindow()->DestroyWindowAndSlot();
            return;
         }
      }
   }

   // Fallback to standard tab destruction
   TRootBrowser::CloseTab(id);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TGMainFrame. Calls TEveManager::Terminate().

void TEveBrowser::CloseWindow()
{
   TEveManager::Terminate();
}

////////////////////////////////////////////////////////////////////////////////
/// Hide the bottom tab (usually holding command-line widget).

void TEveBrowser::HideBottomTab()
{
   fV2->HideFrame(fHSplitter);
   fV2->HideFrame(fH2);
}

////////////////////////////////////////////////////////////////////////////////
/// TRootBrowser keeps (somewhat unnecessarily) counters for number ob tabs
/// on each position. Eve bastardizes the right tab so we have to fix the counters
/// when a new window is added ... it doesn't seem to be needed when it is removed.

void TEveBrowser::SanitizeTabCounts()
{
   fNbTab[TRootBrowser::kRight] = fTabRight->GetNumberOfTabs();
   fCrTab[TRootBrowser::kRight] = fTabRight->GetNumberOfTabs() - 1;
}
