// @(#)root/guibuilder:$Id$
// Author: Valeriy Onuchin   12/09/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGuiNameFrame                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TGuiBldNameFrame.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGuiBldEditor.h"
#include "TGLayout.h"
#include "TG3DLine.h"
#include "TColor.h"
#include "TROOT.h"
#include "TRootGuiBuilder.h"
#include "TGButton.h"
#include "TGFrame.h"
#include "TGMdiFrame.h"
#include "TGCanvas.h"
#include "TGListTree.h"
#include "TGuiBldDragManager.h"
#include "TGMsgBox.h"

ClassImp(TGuiBldNameFrame);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGuiBldNameFrame::TGuiBldNameFrame(const TGWindow *p, TGuiBldEditor *editor) :
                  TGCompositeFrame(p, 1, 1)
{
   fEditor = editor;
   fBuilder = (TRootGuiBuilder*)TRootGuiBuilder::Instance();
   fManager = fBuilder->GetManager();
   fEditDisabled = kEditDisable;
   SetCleanup(kDeepCleanup);
   TGFrame *frame = 0;
   TGFrame *fSelected = fEditor->GetSelected();
   if (fSelected) frame = fSelected;

   TGVerticalFrame *cf = new TGVerticalFrame(this, 180, 400);

   //list tree
   TGHorizontalFrame *f = new TGHorizontalFrame(cf);
   f->AddFrame(new TGLabel(f, "MDI Frame content"),
               new TGLayoutHints(kLHintsLeft, 0, 1, 0, 0));
   f->AddFrame(new TGHorizontal3DLine(f), new TGLayoutHints(kLHintsExpandX,
                                                            5, 5, 7, 7));
   cf->AddFrame(f, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   fCanvas = new TGCanvas(cf, 180, 110);
   fListTree = new TGListTree(fCanvas, 0);
   fCanvas->MapSubwindows();
   cf->AddFrame(fCanvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));

   //nameframe
   fTitleFrame = new TGHorizontalFrame(cf, 100, 30);
   fTitleFrame->AddFrame(new TGLabel(fTitleFrame, "Variable name"),
                         new TGLayoutHints(kLHintsLeft | kLHintsCenterY ,
                                           0, 1, 0, 0));
   fTitleFrame->AddFrame(new TGHorizontal3DLine(fTitleFrame),
                         new TGLayoutHints(kLHintsCenterY | kLHintsExpandX,
                                           1, 1, 1, 1));
   cf->AddFrame(fTitleFrame, new TGLayoutHints(kLHintsExpandX | kLHintsTop));

   TString name = "";
   if (frame) {
      name = frame->ClassName();
   }
   fLabel = new TGLabel(cf, name.Data());
   cf->AddFrame(fLabel, new TGLayoutHints(kLHintsCenterX, 1, 1, 0, 0));

   TGCompositeFrame *sub  = new TGHorizontalFrame(cf, 100, 30);
   fFrameName = new TGTextEntry(sub, frame ? frame->GetName() : "noname");
   fFrameName->SetAlignment(kTextLeft);
   fFrameName->Resize(120, fFrameName->GetHeight());
   sub->AddFrame(fFrameName, new TGLayoutHints(kLHintsTop | kLHintsCenterX,
                                               2, 2, 0, 0));
   fFrameName->SetEnabled(kTRUE);

   TGTextButton *btn = new TGTextButton(sub, "   Set Name   ");
   sub->AddFrame(btn, new TGLayoutHints(kLHintsTop));
   cf->AddFrame(sub, new TGLayoutHints(kLHintsTop | kLHintsCenterX,
                                       2, 2, 0, 0));

   AddFrame(cf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   btn->Connect("Clicked()", "TGuiBldNameFrame", this, "UpdateName()");
   btn->SetToolTipText("Set variable name");
   fListTree->Connect("Clicked(TGListTreeItem*, Int_t)", "TGuiBldNameFrame",
                      this, "SelectFrameByItem(TGListTreeItem*, Int_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw frame (just a prototype).

void TGuiBldNameFrame::DoRedraw()
{
   //TColor *col = gROOT->GetColor(29);
   //TRootGuiBuilder::SetBgndColor(fTitleFrame, col->GetPixel());
   TGCompositeFrame::DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset name frame.

void TGuiBldNameFrame::Reset()
{
   fFrameName->SetText("");
   fLabel->SetText("");
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Change selected frame.

void TGuiBldNameFrame::ChangeSelected(TGFrame *frame)
{
   fFrameName->Disconnect();

   if (!frame) {
      Reset();
      return;
   }

   TString name = frame->ClassName();

   fLabel->SetText(name.Data());
   fFrameName->SetText(frame->GetName());
   Resize();

   TGCompositeFrame *main = GetMdi(frame);

   if (main) {
      if (!fListTree->GetFirstItem())
         MapItems(main);
      else if ((fListTree->GetFirstItem()->GetUserData()) != main) {
         //different MDI
         //clear the list tree displayed
         while (fListTree->GetFirstItem()) {
            fListTree->DeleteItem(fListTree->GetFirstItem());
         }
         MapItems(main);
      }
      else // check if new items added or old ones reparented -> update tree
         CheckItems(main);
   }

   //highlight and open
   TGListTreeItem *item = 0;
   fListTree->OpenItem(fListTree->GetFirstItem()); //mdi
   item = fListTree->FindItemByObj(fListTree->GetFirstItem(), frame);
   if (item) {
      fListTree->HighlightItem(item);
      while (item->GetParent()) {
         item = item->GetParent();
         item->SetOpen(kTRUE);
      }
   }
   fClient->NeedRedraw(fListTree, kTRUE);
   fClient->NeedRedraw(fCanvas, kTRUE);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Set new name of frame, if it doesn't already exist in the same MDI frame.

void TGuiBldNameFrame::UpdateName()
{
   TGFrame *frame = fEditor->GetSelected();
   TString ch = fFrameName->GetText();

   if (!frame) {
      return;
   }

   if (FindItemByName(fListTree, ch, fListTree->GetFirstItem())) {
      fBuilder->UpdateStatusBar("Variable name already exists.");
      TGCompositeFrame *cf = (TGCompositeFrame*)frame->GetParent();
      int retval;
      fBuilder->GetManager()->SetEditable(kFALSE);
      new TGMsgBox(fClient->GetDefaultRoot(), fBuilder,
                   "Name conflict", "Variable name already exists.",
                   kMBIconExclamation, kMBOk, &retval);
      cf->SetEditable(kTRUE);
      // hack against selecting the message box itself
      fBuilder->GetManager()->SelectFrame(frame);
      frame->SetEditable(kTRUE);
   }
   else {
      fBuilder->UpdateStatusBar("Variable name changed.");
      frame->SetName(ch);
   }

   //clear the list tree displayed
   while (fListTree->GetFirstItem()) {
      fListTree->DeleteItem(fListTree->GetFirstItem());
   }

   TGCompositeFrame *main = GetMdi(frame);
   MapItems(main);

   fClient->NeedRedraw(fListTree, kTRUE);
   fClient->NeedRedraw(fFrameName);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Find the parent mdi frame

TGCompositeFrame *TGuiBldNameFrame::GetMdi(TGFrame *frame)
{
   TGFrame *p = frame;

   while (p && (p != fClient->GetDefaultRoot()) ) {
      if (p->InheritsFrom(TGMdiFrame::Class())) {
         return (TGCompositeFrame*)p;
      }
      else if (p->InheritsFrom(TGMainFrame::Class())) {
         return (TGCompositeFrame*)p;
      }
      p = (TGFrame*)p->GetParent();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Map all the frames and subframes in mdi frame to the list tree.

void TGuiBldNameFrame::MapItems(TGCompositeFrame *main)
{
   if (!main) {
     return;
   }

   TList *list = main->GetList(); //list of all elements in the frame
   TGFrameElement *el = 0;
   TIter next(list);

   while ((el = (TGFrameElement *) next())) {
      if (el->fFrame) {

         if (main->InheritsFrom(TGMdiFrame::Class()) ||
             main->InheritsFrom(TGMainFrame::Class())) {

            // first loop, we're in the main frame -> add items directly
            // to main frame folder of the tree list
            if (!fListTree->FindChildByData(0, main)) {
               // add main frame to root
               fListTree->AddItem(0, main->GetName(), main);
            }
             //add other items to mainframe
            fListTree->AddItem(fListTree->FindChildByData(0, main),
                               el->fFrame->GetName(), el->fFrame);

         } else { //means we're in recursion loop, browsing in subframe
            // result is the name of the tree folder to which we want to
            // place the element
            TGListTreeItem *result = 0;
            TGFrame *par = (TGFrame*)el->fFrame->GetParent();
            result = fListTree->FindItemByObj(fListTree->GetFirstItem(), par);
            if (result)
               fListTree->AddItem(result, el->fFrame->GetName(), el->fFrame);
         }

         if ( (el->fFrame->InheritsFrom(TGCompositeFrame::Class())) &&
              (!(el->fFrame->InheritsFrom(TGMdiFrame::Class()))) ) {
               //recursive call for composite subframes
            main = (TGCompositeFrame*)(el->fFrame);
            MapItems(main);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if items are in the list tree and at the same place.

Bool_t TGuiBldNameFrame::CheckItems(TGCompositeFrame *main)
{
   TList *list = main->GetList(); //list of all elements in the frame

   TGFrameElement *el = 0;
   TGListTreeItem *item = 0;
   TIter next(list);
   TGFrame *f = 0;
   TGListTreeItem *par = 0;

   while ((el = (TGFrameElement *) next())) {
      if (el && (el->fFrame)) {
         item = fListTree->FindItemByObj(fListTree->GetFirstItem(),
                                         el->fFrame);
         if (!item) {
            f = (TGFrame*)el->fFrame->GetParent();
            if (f) {
               par = fListTree->FindItemByObj(fListTree->GetFirstItem(), f);
               if (par)
                  fListTree->AddItem(par, el->fFrame->GetName(), el->fFrame);
            }
            //return kTRUE; //selected item not found = is newly created
         }
         else if (item->GetParent() && item->GetParent()->GetUserData() !=
                  el->fFrame->GetParent()) {
            f = (TGFrame*)el->fFrame->GetParent();
            if (f) {
               par = fListTree->FindItemByObj(fListTree->GetFirstItem(), f);
               if (par)
                  fListTree->Reparent(item, par);
            }
            //return kTRUE; //parent of the item changed
         }
         if (el->fFrame->InheritsFrom(TGCompositeFrame::Class())) {
            CheckItems((TGCompositeFrame*)el->fFrame);
         }
      }
   }
   return kFALSE; //treelist remains the same
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a frame.

void TGuiBldNameFrame::RemoveFrame(TGFrame *frame)
{
   TGListTreeItem *item;
   item = fListTree->FindItemByObj(fListTree->GetFirstItem(), frame);
   if (item) {
      fListTree->DeleteItem(item);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find item with GetText == name. Search tree downwards starting at item.

TGListTreeItem *TGuiBldNameFrame::FindItemByName(TGListTree *tree,
                                                 const char* name,
                                                 TGListTreeItem *item)
{
   TGListTreeItem *fitem;
   if (item && name) {
      if (!strcmp(item->GetText(), name)) { //if names are equal
         return item;
      }
      else {
         if (item->GetFirstChild()) {
            fitem = FindItemByName(tree, name, item->GetFirstChild());
            if (fitem) return fitem;
         }
         return FindItemByName(tree, name, item->GetNextSibling());
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// When list tree item is clicked, frame with that name is selected.

void TGuiBldNameFrame::SelectFrameByItem(TGListTreeItem* item, Int_t)
{
   TGFrame *frame = (TGFrame*)item->GetUserData();
   if (frame) {
      ((TGFrame*)frame->GetParent())->SetEditable(kTRUE);
      fManager->SelectFrame(frame);
      frame->SetEditable(kTRUE);
      fClient->NeedRedraw(frame);
   }
}


