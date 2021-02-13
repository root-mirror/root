// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Jun 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLClipSetEditor.h"
#include "TGLClip.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"

#include "TGLViewerEditor.h"

/** \class TGLClipSetSubEditor
\ingroup opengl
GUI sub-editor for TGLClipSet.
*/

ClassImp(TGLClipSetSubEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLClipSetSubEditor::TGLClipSetSubEditor(const TGWindow *p) :
   TGVerticalFrame(p),
   fM(0),
   fCurrentClip(TGLClip::kClipNone),
   fTypeButtons(0),
   fPlanePropFrame(0),
   fPlaneProp(),
   fBoxPropFrame(0),
   fBoxProp(),
   fClipInside(0),
   fAutoUpdate(0),
   fClipEdit(0),
   fClipShow(0),
   fApplyButton(0),
   fResetButton(0)
{
   fTypeButtons = new TGButtonGroup(this, "Clip Type", kChildFrame|kHorizontalFrame);
   new TGRadioButton(fTypeButtons, "None    ");
   new TGRadioButton(fTypeButtons, "Plane    ");
   new TGRadioButton(fTypeButtons, "Box");
   fTypeButtons->SetLayoutHints(new TGLayoutHints(kLHintsLeft|kLHintsBottom, 0, 0, 2, -10));
   AddFrame(fTypeButtons, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));
   // Clip inside / edit
   fClipInside = new TGCheckButton(this, "Clip away inside");
   AddFrame(fClipInside, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));
   fAutoUpdate = new TGCheckButton(this, "Auto update clip");
   AddFrame(fAutoUpdate, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));
   fClipEdit   = new TGCheckButton(this, "Edit In Viewer");
   AddFrame(fClipEdit, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));
   fClipShow   = new TGCheckButton(this, "Show In Viewer");
   AddFrame(fClipShow, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));

   // Plane properties
   fPlanePropFrame = new TGCompositeFrame(this);
   AddFrame(fPlanePropFrame, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));

   static const char * const planeStr[] = { "aX + ", "bY +", "cZ + ", "d = 0" };

   for (Int_t i = 0; i < 4; ++i)
      fPlaneProp[i] = TGLViewerEditor::MakeLabeledNEntry(fPlanePropFrame, planeStr[i] , 40);

   // Box properties
   fBoxPropFrame = new TGCompositeFrame(this);
   AddFrame(fBoxPropFrame, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));

   static const char * const boxStr[] = {"Center X", "Center Y", "Center Z", "Length X", "Length Y", "Length Z" };

   for (Int_t i = 0; i < 6; ++i)
      fBoxProp[i] = TGLViewerEditor::MakeLabeledNEntry(fBoxPropFrame, boxStr[i] , 60);

   // Apply button
   fApplyButton = new TGTextButton(this, "Apply");
   AddFrame(fApplyButton, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));

   fResetButton = new TGTextButton(this, "Reset");
   AddFrame(fResetButton, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 3, 2, 2));

   fTypeButtons->Connect("Clicked(Int_t)", "TGLClipSetSubEditor", this, "ClipTypeChanged(Int_t)");
   fClipInside->Connect("Clicked()", "TGLClipSetSubEditor", this, "UpdateViewerClip()");
   fAutoUpdate->Connect("Clicked()", "TGLClipSetSubEditor", this, "UpdateViewerClip()");
   fClipEdit->Connect("Clicked()", "TGLClipSetSubEditor", this, "UpdateViewerClip()");
   fClipShow->Connect("Clicked()", "TGLClipSetSubEditor", this, "UpdateViewerClip()");

   for (Int_t i = 0; i < 4; ++i)
      fPlaneProp[i]->Connect("ValueSet(Long_t)", "TGLClipSetSubEditor", this, "ClipValueChanged()");

   for (Int_t i = 0; i < 6; ++i)
      fBoxProp[i]->Connect("ValueSet(Long_t)", "TGLClipSetSubEditor", this, "ClipValueChanged()");

   fApplyButton->Connect("Pressed()", "TGLClipSetSubEditor", this, "UpdateViewerClip()");
   fResetButton->Connect("Pressed()", "TGLClipSetSubEditor", this, "ResetViewerClip()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TGLClipSetSubEditor::SetModel(TGLClipSet* m)
{
   fM = m;

   fTypeButtons->GetButton(fCurrentClip+1)->SetDown(kFALSE);
   fCurrentClip = fM->GetClipType();

   Double_t clip[6] = {0.};
   fM->GetClipState(fCurrentClip, clip);

   fApplyButton->SetState(kButtonDisabled);

   // Button ids run from 1
   if (TGButton *btn = fTypeButtons->GetButton(fCurrentClip+1)) {
      btn->SetDown();
      fTypeButtons->SetButton(fCurrentClip+1);
   }
   Bool_t active = (fCurrentClip != TGLClip::kClipNone);
   fClipInside->SetEnabled(active);
   fAutoUpdate->SetEnabled(active);
   fClipEdit  ->SetEnabled(active);
   fClipShow  ->SetEnabled(active);
   if (active) {
      fClipEdit->SetDown(fM->GetShowManip());
      fClipShow->SetDown(fM->GetShowClip());
      fClipInside->SetDown(fM->GetCurrentClip()->GetMode() == TGLClip::kInside);
      fAutoUpdate->SetDown(fM->GetAutoUpdate());

      if (fCurrentClip == TGLClip::kClipPlane) {
         HideFrame(fBoxPropFrame);
         ShowFrame(fPlanePropFrame);
         for (Int_t i = 0; i < 4; ++i)
            fPlaneProp[i]->SetNumber(clip[i]);
      } else if (fCurrentClip == TGLClip::kClipBox) {
         HideFrame(fPlanePropFrame);
         ShowFrame(fBoxPropFrame);
         for (Int_t i = 0; i < 6; ++i)
            fBoxProp[i]->SetNumber(clip[i]);
      }

      fResetButton->SetState(kButtonUp);
   } else {
      HideFrame(fPlanePropFrame);
      HideFrame(fBoxPropFrame);

      fResetButton->SetState(kButtonDisabled);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Changed signal.

void TGLClipSetSubEditor::Changed()
{
   Emit("Changed()");
}

////////////////////////////////////////////////////////////////////////////////
/// One of number entries was changed.

void TGLClipSetSubEditor::ClipValueChanged()
{
   fApplyButton->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Clip type radio button changed - update viewer.

void TGLClipSetSubEditor::ClipTypeChanged(Int_t id)
{
   switch (id)
   {
      case 2:  fCurrentClip = TGLClip::kClipPlane; break;
      case 3:  fCurrentClip = TGLClip::kClipBox;   break;
      default: fCurrentClip = TGLClip::kClipNone;  break;
   }
   fM->SetClipType(fCurrentClip);
   SetModel(fM);
   ((TGMainFrame*)GetMainFrame())->Layout();
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Change clipping volume.

void TGLClipSetSubEditor::UpdateViewerClip()
{
   Double_t data[6] = {0.};
   // Fetch GUI state for clip if 'type' into 'data' vector
   if (fCurrentClip == TGLClip::kClipPlane)
      for (Int_t i = 0; i < 4; ++i)
         data[i] = fPlaneProp[i]->GetNumber();
   else if (fCurrentClip == TGLClip::kClipBox)
      for (Int_t i = 0; i < 6; ++i)
         data[i] = fBoxProp[i]->GetNumber();

   fApplyButton->SetState(kButtonDisabled);

   fM->SetClipState(fCurrentClip, data);
   fM->SetShowManip(fClipEdit->IsDown());
   fM->SetShowClip (fClipShow->IsDown());
   if (fCurrentClip != TGLClip::kClipNone)
      fM->GetCurrentClip()->SetMode(fClipInside->IsDown() ? TGLClip::kInside : TGLClip::kOutside);
   fM->SetAutoUpdate(fAutoUpdate->IsDown());

   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset transformation of the current clip.

void TGLClipSetSubEditor::ResetViewerClip()
{
   if (fCurrentClip != TGLClip::kClipNone) {
      if (fM->GetClipType() == TGLClip::kClipPlane) {
         TGLPlane plane(0.0, -1.0, 0.0, 0.0);
         dynamic_cast<TGLClipPlane*>(fM->GetCurrentClip())->Set(plane);
      } else if (fM->GetClipType() == TGLClip::kClipBox) {
         fM->GetCurrentClip()->SetTransform(TGLMatrix());
      }
      Changed();
   }
}


/** \class TGLClipSetEditor
\ingroup opengl
GUI editor for TGLClipSet.
*/

ClassImp(TGLClipSetEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLClipSetEditor::TGLClipSetEditor(const TGWindow *p, Int_t width, Int_t height,
                                   UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM  (0),
   fSE (0)
{
   MakeTitle("TGLClipSet");

   fSE = new TGLClipSetSubEditor(this);
   AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
   fSE->Connect("Changed()", "TGLClipSetEditor", this, "Update()");
}


////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TGLClipSetEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TGLClipSet*>(obj);
   fSE->SetModel(fM);
}
