// @(#)root/eve:$Id$
// Author: Dmytro Kovalskyi, 28.2.2008

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveParamList.h"

#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"

#include "TEveGValuators.h"
#include "TGedEditor.h"

/** \class TEveParamList
\ingroup TEve
Collection of named parameters.
*/

ClassImp(TEveParamList);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveParamList::TEveParamList(const char* n, const char* t, Bool_t doColor) :
   TNamed(n, t),
   fColor(0)
{
   if (doColor) SetMainColorPtr(&fColor);
}

////////////////////////////////////////////////////////////////////////////////
/// Get config-struct for float parameter 'name'.

TEveParamList::FloatConfig_t TEveParamList::GetFloatParameter(const TString& name)
{
   static const TEveException eh("TEveParamList::GetFloatParameter ");

   for (FloatConfigVec_ci itr = fFloatParameters.begin(); itr != fFloatParameters.end(); ++itr)
      if (itr->fName.CompareTo(name)==0 ) return *itr;
   Error(eh, "parameter not found.");
   return FloatConfig_t();
}

////////////////////////////////////////////////////////////////////////////////
/// Get config-struct for int parameter 'name'.

TEveParamList::IntConfig_t TEveParamList::GetIntParameter(const TString& name)
{
   static const TEveException eh("TEveParamList::GetIntParameter ");

   for (IntConfigVec_ci itr = fIntParameters.begin(); itr != fIntParameters.end(); ++itr)
      if (itr->fName.CompareTo(name) == 0) return *itr;
   Error(eh, "parameter not found.");
   return IntConfig_t();
}

////////////////////////////////////////////////////////////////////////////////
/// Get value for bool parameter 'name'.

Bool_t TEveParamList::GetBoolParameter(const TString& name)
{
   static const TEveException eh("TEveParamList::GetBoolParameter ");

   for (BoolConfigVec_ci itr = fBoolParameters.begin(); itr != fBoolParameters.end(); ++itr)
      if ( itr->fName.CompareTo(name)==0 ) return itr->fValue;
   Error(eh, "parameter not found.");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit ParamChanged() signal.

void TEveParamList::ParamChanged(const char* name)
{
   Emit("ParamChanged(char*)", name);
}

/** \class TEveParamListEditor
\ingroup TEve
GUI editor for TEveParamList.

Slot methods from this object do not call Update, instead they call
their model's ParamChanged(const char* name) function which emits a
corresponding signal.

This helps in handling of parameter changes as they are probably
related to displayed objects in a more complicated way. Further,
the TGCheckButton::HandleButton() emits more signal after the
Clicked() signal and if model is reset in the editor, its contents
is removed. This results in a crash.
*/

ClassImp(TEveParamListEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveParamListEditor::TEveParamListEditor(const TGWindow *p, Int_t width, Int_t height,
                                         UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM          (0),
   fParamFrame (0)
{
   MakeTitle("TEveParamList");
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize widgets when a new object is selected.

void TEveParamListEditor::InitModel(TObject* obj)
{
   fM = dynamic_cast<TEveParamList*>(obj);

   if (fParamFrame) {
      fParamFrame->UnmapWindow();
      RemoveFrame(fParamFrame);
      fParamFrame->DestroyWindow();
      delete fParamFrame;
   }
   fParamFrame = new TGVerticalFrame(this);
   AddFrame(fParamFrame);

   // integers
   fIntParameters.clear();
   for (UInt_t i = 0; i < fM->fIntParameters.size(); ++i)
   {
      TGCompositeFrame* frame = new TGHorizontalFrame(fParamFrame);

      // number entry widget
      TGNumberEntry* widget = new TGNumberEntry
         (frame, fM->fIntParameters[i].fValue,
          5,                                 // number of digits
          i,                                 // widget ID
          TGNumberFormat::kNESInteger,       // style
          TGNumberFormat::kNEAAnyNumber,     // input value filter
          TGNumberFormat::kNELLimitMinMax,   // specify limits
          fM->fIntParameters[i].fMin,        // min value
          fM->fIntParameters[i].fMax);       // max value
      frame->AddFrame(widget, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 2,8,2,2));
      widget->Connect("ValueSet(Long_t)", "TEveParamListEditor", this, "DoIntUpdate()");
      fIntParameters.push_back(widget);

      // label
      frame->AddFrame(new TGLabel(frame,fM->fIntParameters[i].fName.Data()),
                      new TGLayoutHints(kLHintsLeft|kLHintsCenterY));

      fParamFrame->AddFrame(frame, new TGLayoutHints(kLHintsTop));
   }


   // floats
   fFloatParameters.clear();
   for (UInt_t i = 0; i < fM->fFloatParameters.size(); ++i)
   {
      TGCompositeFrame* frame = new TGHorizontalFrame(fParamFrame);

      // number entry widget
      TGNumberEntry* widget = new TGNumberEntry
         (frame, fM->fFloatParameters[i].fValue,
          5,                                // number of digits
          i,                                // widget ID
          TGNumberFormat::kNESReal,         // style
          TGNumberFormat::kNEAAnyNumber,    // input value filter
          TGNumberFormat::kNELLimitMinMax,  // specify limits
          fM->fFloatParameters[i].fMin,     // min value
          fM->fFloatParameters[i].fMax);    // max value
      frame->AddFrame(widget, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 2,8,2,2));
      widget->Connect("ValueSet(Long_t)", "TEveParamListEditor", this, "DoFloatUpdate()");
      fFloatParameters.push_back( widget );

      // label
      frame->AddFrame(new TGLabel(frame,fM->fFloatParameters[i].fName.Data()),
                      new TGLayoutHints(kLHintsLeft|kLHintsCenterY) );

      fParamFrame->AddFrame(frame, new TGLayoutHints(kLHintsTop));
   }

   // boolean
   fBoolParameters.clear();
   for (UInt_t i = 0; i < fM->fBoolParameters.size(); ++i)
   {
      TGCheckButton* widget = new TGCheckButton(fParamFrame,
                                                fM->fBoolParameters[i].fName.Data(),
                                                i);
      widget->Connect("Clicked()", "TEveParamListEditor", this, "DoBoolUpdate()");
      fBoolParameters.push_back(widget);

      fParamFrame->AddFrame(widget, new TGLayoutHints(kLHintsTop,2,0,1,1));
   }
   MapSubwindows();
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveParamListEditor::SetModel(TObject* obj)
{
   InitModel(obj);

   for (UInt_t i = 0; i < fIntParameters.size(); ++i)
      fIntParameters[i]->GetNumberEntry()->SetIntNumber(fM->fIntParameters[i].fValue);

   for (UInt_t i = 0; i < fFloatParameters.size(); ++i)
      fFloatParameters[i]->GetNumberEntry()->SetNumber(fM->fFloatParameters[i].fValue);

   for (UInt_t i = 0; i < fBoolParameters.size(); ++i)
      fBoolParameters[i]->SetState( fM->fBoolParameters[i].fValue ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for integer parameter update.

void TEveParamListEditor::DoIntUpdate()
{
   TGNumberEntry *widget = (TGNumberEntry*) gTQSender;
   Int_t id = widget->WidgetId();
   if (id < 0 || id >= (int) fM->fIntParameters.size()) return;
   fM->fIntParameters[id].fValue = widget->GetNumberEntry()->GetIntNumber();

   fM->ParamChanged(fM->fIntParameters[id].fName);
   gTQSender = (void*) widget;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for float parameter update.

void TEveParamListEditor::DoFloatUpdate()
{
   TGNumberEntry *widget = (TGNumberEntry*) gTQSender;
   Int_t id = widget->WidgetId();
   if (id < 0 || id >= (int) fM->fFloatParameters.size()) return;
   fM->fFloatParameters[id].fValue = widget->GetNumberEntry()->GetNumber();

   fM->ParamChanged(fM->fFloatParameters[id].fName);
   gTQSender = (void*) widget;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for bool parameter update.

void TEveParamListEditor::DoBoolUpdate()
{
   TGCheckButton *widget = (TGCheckButton*) gTQSender;
   Int_t id = widget->WidgetId();
   if (id < 0 || id >= (int) fM->fBoolParameters.size()) return;
   fM->fBoolParameters[id].fValue = widget->IsOn();

   fM->ParamChanged(fM->fBoolParameters[id].fName);
   gTQSender = (void*) widget;
}
