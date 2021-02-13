// @(#)root/gviz3d:$Id$
// Author: Tomasz Sosnicki   18/09/09

/************************************************************************
* Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/


#include "TStructNodeEditor.h"
#include "TStructNode.h"
#include "TStructNodeProperty.h"

#include <TGColorSelect.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TClass.h>

ClassImp(TStructNodeEditor);

//________________________________________________________________________
//////////////////////////////////////////////////////////////////////////
//
// TStructNodeEditor is an editor for changing node attributes such as
// maximum numbers of level or maximum number of objects diplayed if this
// node is our top node. We can also change color associated with a class
// or chagne color to default.
//
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Constructor of node attributes GUI.

TStructNodeEditor::TStructNodeEditor(TList* colors, const TGWindow *p, Int_t width, Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back), fColors(colors)
{
   MakeTitle("TStructNode");
   fInit = kFALSE;

   TGLayoutHints* expandX = new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 5,5,5,5);
   fNodeNameLabel = new TGLabel(this, "No node selected");
   this->AddFrame(fNodeNameLabel, expandX);

   fTypeName = new TGLabel(this);

   this->AddFrame(fTypeName, expandX);

   TGHorizontalFrame* maxObjectsFrame = new TGHorizontalFrame(this);
   TGLabel* fMaxObjectslabel = new TGLabel(maxObjectsFrame, "Max objects:");
   maxObjectsFrame->AddFrame(fMaxObjectslabel);

   fMaxObjectsNumberEntry = new TGNumberEntry(maxObjectsFrame, 0);
   fMaxObjectsNumberEntry->SetFormat(TGNumberEntry::kNESInteger);
   fMaxObjectsNumberEntry->SetLimits(TGNumberEntry::kNELLimitMin, 1);
   fMaxObjectsNumberEntry->SetState(kFALSE);
   fMaxObjectsNumberEntry->Connect("ValueSet(Long_t)", "TStructNodeEditor", this, "MaxObjectsValueSetSlot(Long_t)");
   maxObjectsFrame->AddFrame(fMaxObjectsNumberEntry);
   this->AddFrame(maxObjectsFrame, expandX);

   TGHorizontalFrame* maxLevelFrame = new TGHorizontalFrame(this);
   TGLabel* fMaxLevelsLabel = new TGLabel(maxLevelFrame, "Max levels:");
   maxLevelFrame->AddFrame(fMaxLevelsLabel);
   fMaxLevelsNumberEntry = new TGNumberEntry(maxLevelFrame, 0);
   fMaxLevelsNumberEntry->SetLimits(TGNumberEntry::kNELLimitMin, 1);
   fMaxLevelsNumberEntry->SetFormat(TGNumberEntry::kNESInteger);
   fMaxLevelsNumberEntry->SetState(kFALSE);
   fMaxLevelsNumberEntry->Connect("ValueSet(Long_t)", "TStructNodeEditor", this, "MaxLevelsValueSetSlot(Long_t)");
   maxLevelFrame->AddFrame(fMaxLevelsNumberEntry);
   this->AddFrame(maxLevelFrame, expandX);

   fNameEntry = new TGTextEntry(this, fName.Data());
   this->AddFrame(fNameEntry, expandX);
   fNameEntry->SetState(kFALSE);

   fColorSelect = new TGColorSelect(this);
   fColorSelect->Connect("ColorSelected(Pixel_t)", "TStructNodeEditor", this, "ColorSelectedSlot(Pixel_t)");
   this->AddFrame(fColorSelect, expandX);
   fColorSelect->SetEnabled(kFALSE);

   fAutoRefesh = new TGCheckButton(this, "Auto refesh");
   fAutoRefesh->SetOn();
   fAutoRefesh->Connect("Toggled(Bool_t)", "TStructNodeEditor", this, "AutoRefreshButtonSlot(Bool_t)");
   fAutoRefesh->SetEnabled(kFALSE);
   this->AddFrame(fAutoRefesh, expandX);

   fDefaultButton = new TGTextButton(this, "Default color");
   fDefaultButton->Connect("Clicked()", "TStructNodeEditor", this, "DefaultButtonSlot()");
   this->AddFrame(fDefaultButton, expandX);
   fDefaultButton->SetEnabled(kFALSE);


   fApplyButton = new TGTextButton(this, "Apply");
   fApplyButton->Connect("Clicked()", "TStructNodeEditor", this, "ApplyButtonSlot()");
   fApplyButton->SetEnabled(kFALSE);
   this->AddFrame(fApplyButton, expandX);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of node editor.

TStructNodeEditor::~TStructNodeEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// ApplyButton Slot. Activated when user press Apply button. Sets properties of a node

void TStructNodeEditor::ApplyButtonSlot()
{
   Bool_t needReset = false;

   if ((Int_t)(fNode->GetMaxLevel()) != fMaxLevelsNumberEntry->GetIntNumber()) {
      fNode->SetMaxLevel(fMaxLevelsNumberEntry->GetIntNumber());
      needReset = true;
   }

   if ((Int_t)(fNode->GetMaxObjects()) != fMaxObjectsNumberEntry->GetIntNumber()) {
      fNode->SetMaxObjects(fMaxObjectsNumberEntry->GetIntNumber());
      needReset = true;
   }

   if (fSelectedPropert) {
      fSelectedPropert->SetColor(fColorSelect->GetColor());
      fSelectedPropert->SetName(fNameEntry->GetText());
   }

   Update(needReset);
}

////////////////////////////////////////////////////////////////////////////////
/// Activated when user chage condition

void TStructNodeEditor::AutoRefreshButtonSlot(Bool_t on)
{
   if (on) {
      Update(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the fill area color.

void TStructNodeEditor::ColorSelectedSlot(Pixel_t color)
{
   if (fAvoidSignal) {
      return;
   }

   TStructNodeProperty* prop = FindNodeProperty(fNode);
   if (prop) {
      prop->SetColor(color);
   } else {
      // add property
      prop = new TStructNodeProperty(fNode->GetTypeName(), color);
      fColors->Add(prop);
      fColors->Sort();
      fSelectedPropert = prop;
      fNameEntry->SetText(fNode->GetTypeName());
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Defaulf button. Sets color of class to default

void TStructNodeEditor::DefaultButtonSlot()
{
   if (TStructNodeProperty* prop = FindNodeProperty(fNode)) {
      fColors->Remove(prop);
      fSelectedPropert = GetDefaultProperty();
      fNameEntry->SetText(fSelectedPropert->GetName());
      fColorSelect->SetColor(fSelectedPropert->GetPixel(), kFALSE);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Retruns property associated to the class of given node "node". If property isn't found
/// then returns NULL

TStructNodeProperty* TStructNodeEditor::FindNodeProperty(TStructNode* node)
{
   TIter it(fColors);
   TStructNodeProperty* prop;
   while ((prop = (TStructNodeProperty*) it() )) {
      TString propName(prop->GetName());
      if (propName.EndsWith("+")) {

         if (TClass* cl = TClass::GetClass(node->GetTypeName())) {
            propName.Remove(propName.Length()-1, 1);
            if (cl->InheritsFrom(propName.Data())) {
               return prop;
            }
         }
      } else {
         if (propName == TString(node->GetTypeName())) {
            return prop;
         }
      }
   }

   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns property with default color

TStructNodeProperty* TStructNodeEditor::GetDefaultProperty()
{
   return (TStructNodeProperty*)fColors->Last();
}

////////////////////////////////////////////////////////////////////////////////
/// Enables button and fields

void TStructNodeEditor::Init()
{
   fMaxObjectsNumberEntry->SetState(kTRUE);
   fMaxLevelsNumberEntry->SetState(kTRUE);
   fNameEntry->SetState(kTRUE);
   fColorSelect->SetEnabled(kTRUE);
   fDefaultButton->SetEnabled(kTRUE);
   fApplyButton->SetEnabled(kTRUE);
   fAutoRefesh->SetEnabled(kTRUE);
   fInit = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Emmited when user changes maximum number of levels

void TStructNodeEditor::MaxLevelsValueSetSlot(Long_t)
{
   fNode->SetMaxLevel(fMaxLevelsNumberEntry->GetIntNumber());

   if(fAutoRefesh->IsOn()) {
      Update(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emmited when user changes maximum number of objects

void TStructNodeEditor::MaxObjectsValueSetSlot(Long_t)
{
   fNode->SetMaxObjects(fMaxObjectsNumberEntry->GetIntNumber());

   if(fAutoRefesh->IsOn()) {
      Update(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used node attributes.

void TStructNodeEditor::SetModel(TObject* obj)
{
   fNode = dynamic_cast<TStructNode *>(obj);
   if (!fNode) return;

   // Add max level
   fMaxLevelsNumberEntry->SetIntNumber(fNode->GetMaxLevel());

   // Add max objects
   fMaxObjectsNumberEntry->SetIntNumber(fNode->GetMaxObjects());

   // Type label
   fTypeName->SetText(fNode->GetTypeName());

   // name label
   fNodeNameLabel->SetText(fNode->GetName());

   // Add color property
   fSelectedPropert = FindNodeProperty(fNode);
   if (!fSelectedPropert)
   {
      fSelectedPropert = GetDefaultProperty();
   }
   fNameEntry->SetText(fSelectedPropert->GetName());
   fColorSelect->SetColor(fSelectedPropert->GetPixel(), kFALSE);

   if (!fInit) {
      Init();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Signal emmited when color or other property like number of level is changed
/// without camera reset

void TStructNodeEditor::Update()
{
   Emit("Update(Bool_t)", false);
}

////////////////////////////////////////////////////////////////////////////////
/// Signal emmited when color or other property like number of level is changed.
/// If "resetCamera" is true, then current camera is reset.

void TStructNodeEditor::Update(Bool_t resetCamera)
{
   Emit("Update(Bool_t)", resetCamera);
}
