// @(#)root/guibuilder:$Id$
// Author: Valeriy Onuchin, Lucie Flekova   12/09/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGuiBldEditor.h"
#include "TGuiBldHintsEditor.h"
#include "TGuiBldNameFrame.h"
#include "TGResourcePool.h"
#include "TGLabel.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"
#include "TG3DLine.h"
#include "TGuiBldGeometryFrame.h"
#include "TRootGuiBuilder.h"
#include "TGuiBldDragManager.h"
#include "TGFrame.h"

ClassImp(TGuiBldGeometryFrame);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGuiBldGeometryFrame::TGuiBldGeometryFrame(const TGWindow *p, TGuiBldEditor *ed)
   : TGVerticalFrame(p, 1, 1)
{
   fEditor = ed;
   fBuilder = (TRootGuiBuilder*)TRootGuiBuilder::Instance();
   fDragManager = (TGuiBldDragManager *)gDragManager;
   fSelected = fEditor->GetSelected();
   fEditDisabled = 1;
   SetCleanup(kDeepCleanup);

   TGGroupFrame *fGroupFrame = new TGGroupFrame(this, "Size");

   TGHorizontalFrame *hf = new TGHorizontalFrame(fGroupFrame);

   hf->AddFrame(new TGLabel(hf, " Width "), new TGLayoutHints(kLHintsLeft | 
                kLHintsCenterY, 2, 2, 2, 2));
   fNEWidth = new TGNumberEntry(hf, 0.0, 4, -1, (TGNumberFormat::EStyle)5);
   hf->AddFrame(fNEWidth, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                2, 2, 2, 2));

   hf->AddFrame(new TGLabel(hf, " Height "), new TGLayoutHints(kLHintsLeft | 
                kLHintsCenterY, 2, 2, 2, 2));
   fNEHeight = new TGNumberEntry(hf, 0.0, 4, -1, (TGNumberFormat::EStyle)5);
   hf->AddFrame(fNEHeight, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                2, 2, 2, 2));

   fGroupFrame->AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 
                         0, 0, 5, 0));

   AddFrame(fGroupFrame, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   
   fNEWidth->GetNumberEntry()->Connect("ReturnPressed()", "TGuiBldGeometryFrame",
                                       this, "ResizeSelected()");
   fNEWidth->Connect("ValueSet(Long_t)", "TGuiBldGeometryFrame", this,
                     "ResizeSelected()");
   fNEHeight->GetNumberEntry()->Connect("ReturnPressed()", "TGuiBldGeometryFrame",
                                        this, "ResizeSelected()");
   fNEHeight->Connect("ValueSet(Long_t)", "TGuiBldGeometryFrame", this,
                      "ResizeSelected()");

   if (!fSelected) {
      fNEWidth->SetNumber(0);
      fNEHeight->SetNumber(0);
   }
   else {
      fNEWidth->SetNumber(fSelected->GetWidth());
      fNEHeight->SetNumber(fSelected->GetHeight());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resize and redraw selected frame.

void TGuiBldGeometryFrame::ResizeSelected()
{
   if (!fEditor)
      return;

   fSelected = fEditor->GetSelected();

   if (!fSelected)
      return;

   Int_t w = fNEWidth->GetIntNumber();
   Int_t h = fNEHeight->GetIntNumber();

   if ((w > 0) && (h > 0)) {
      fSelected->MoveResize(fSelected->GetX(), fSelected->GetY(), w, h);
      fClient->NeedRedraw(fSelected, kTRUE);
      TGWindow *root = (TGWindow*)fClient->GetRoot();
      fClient->NeedRedraw(root, kTRUE);
      fDragManager->DrawGrabRectangles(fSelected);
      if (fBuilder) {
         fClient->NeedRedraw(fBuilder, kTRUE);
      }
   } else {
      fNEWidth->SetNumber(fSelected->GetWidth());
      fNEHeight->SetNumber(fSelected->GetHeight());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update number entries when new frame selected.

void TGuiBldGeometryFrame::ChangeSelected(TGFrame *frame)
{
   if (!frame) {
      fNEWidth->SetNumber(0);
      fNEHeight->SetNumber(0);
   } else {
      fNEWidth->SetNumber(frame->GetWidth());
      fNEHeight->SetNumber(frame->GetHeight());
   }
}


