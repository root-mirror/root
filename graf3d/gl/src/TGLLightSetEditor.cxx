// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLLightSetEditor.h"
#include <TGLLightSet.h>

#include <TGButton.h>

/** \class TGLLightSetSubEditor
\ingroup opengl
Sub-editor for TGLLightSet.
*/

ClassImp(TGLLightSetSubEditor);

////////////////////////////////////////////////////////////////////////////////

TGLLightSetSubEditor::TGLLightSetSubEditor(const TGWindow *p) :
   TGVerticalFrame(p),
   fM             (0),

   fLightFrame    (0),
   fTopLight      (0),
   fRightLight    (0),
   fBottomLight   (0),
   fLeftLight     (0),
   fFrontLight    (0),
   fSpecularLight (0)
{
   // Constructor.

   fLightFrame = new TGGroupFrame(this, "Light sources:", kVerticalFrame);//, kLHintsTop | kLHintsCenterX);
   fLightFrame->SetTitlePos(TGGroupFrame::kLeft);
   AddFrame(fLightFrame, new TGLayoutHints(kLHintsTop| kLHintsExpandX, 1, 1, 1, 1));//-
   TGCompositeFrame* hf =0;

   hf = new TGHorizontalFrame(fLightFrame);
   fTopLight      = MakeLampButton("Top",      TGLLightSet::kLightTop, hf);
   fBottomLight   = MakeLampButton("Bottom",   TGLLightSet::kLightBottom, hf);
   fLightFrame->AddFrame(hf, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 0, 0, 2, 2));

   hf = new TGHorizontalFrame(fLightFrame);
   fLeftLight     = MakeLampButton("Left",     TGLLightSet::kLightLeft, hf);
   fRightLight    = MakeLampButton("Right",    TGLLightSet::kLightRight, hf);
   fLightFrame->AddFrame(hf, new TGLayoutHints(kLHintsTop|kLHintsExpandX , 0, 0, 0, 2));

   hf = new TGHorizontalFrame(fLightFrame);
   fFrontLight    = MakeLampButton("Front",    TGLLightSet::kLightFront, hf);
   fSpecularLight = MakeLampButton("Specular", TGLLightSet::kLightSpecular, hf);

   fLightFrame->AddFrame(hf, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 0, 0, 0, 2));
}

////////////////////////////////////////////////////////////////////////////////
/// Create a button for given lamp and set it up.

TGButton* TGLLightSetSubEditor::MakeLampButton(const char* name, Int_t wid,
                                               TGCompositeFrame* parent)
{
   TGButton* b = new TGCheckButton(parent, name, wid);
   parent->AddFrame(b, new TGLayoutHints(kLHintsNormal|kLHintsExpandX, -2, 0, 0, 2));
   b->Connect("Clicked()", "TGLLightSetSubEditor", this, "DoButton()");
   return b;
}

////////////////////////////////////////////////////////////////////////////////
/// New model was set, refresh data.

void TGLLightSetSubEditor::SetModel(TGLLightSet* m)
{
   fM = m;
   UInt_t als = fM->GetLightState();

   fTopLight   ->SetState((als & TGLLightSet::kLightTop)    ? kButtonDown : kButtonUp);
   fRightLight ->SetState((als & TGLLightSet::kLightRight)  ? kButtonDown : kButtonUp);
   fBottomLight->SetState((als & TGLLightSet::kLightBottom) ? kButtonDown : kButtonUp);
   fLeftLight  ->SetState((als & TGLLightSet::kLightLeft)   ? kButtonDown : kButtonUp);
   fFrontLight ->SetState((als & TGLLightSet::kLightFront)  ? kButtonDown : kButtonUp);

   fSpecularLight->SetState(fM->GetUseSpecular() ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Data in sub-editor has been changed, emit "Changed()" signal.

void TGLLightSetSubEditor::Changed()
{
   Emit("Changed()");
}

////////////////////////////////////////////////////////////////////////////////
/// Lights radio button was clicked.

void TGLLightSetSubEditor::DoButton()
{
   TGButton* b = (TGButton*) gTQSender;
   fM->SetLight(TGLLightSet::ELight(b->WidgetId()), b->IsOn());
   Changed();
}


//______________________________________________________________________________
// TGLLightSetEditor
//
// Editor for TGLLightSet.

ClassImp(TGLLightSetEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLLightSetEditor::TGLLightSetEditor(const TGWindow *p,
                                     Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM  (0),
   fSE (0)
{
   MakeTitle("TGLLightSet");

   fSE = new TGLLightSetSubEditor(this);
   AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
   fSE->Connect("Changed()", "TGLLightSetEditor", this, "Update()");
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLLightSetEditor::~TGLLightSetEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// SetModel ... forward to sub-editor.

void TGLLightSetEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TGLLightSet*>(obj);
   fSE->SetModel(fM);
}
