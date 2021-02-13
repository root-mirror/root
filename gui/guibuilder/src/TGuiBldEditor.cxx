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
// TGuiBldEditor - the property editor                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGuiBldEditor.h"
#include "TRootGuiBuilder.h"
#include "TGuiBldHintsEditor.h"
#include "TGuiBldNameFrame.h"
#include "TGuiBldGeometryFrame.h"
#include "TGResourcePool.h"
#include "TGTab.h"
#include "TGLabel.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"
#include "TG3DLine.h"
#include "TGColorSelect.h"
#include "TGuiBldDragManager.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGSplitter.h"
#include "TGTableLayout.h"

ClassImp(TGuiBldEditor);


////////////////////////////////////////////////////////////////////////////////
class TGuiBldBorderFrame : public TGVerticalFrame {

private:
   enum  EBldBorderFrameMode { 
      kBldBorderNone, kBldBorderSunken,
      kBldBorderPlain, kBldBorderRaised, 
      kBldBorderDouble };

private:
   TGuiBldEditor   *fEditor;
   TGFrame         *fSelected;
   TGButtonGroup   *fBtnGroup;
   TGColorSelect   *fBgndFrame;

public:
   TGuiBldBorderFrame(const TGWindow *p, TGuiBldEditor *editor);
   virtual ~TGuiBldBorderFrame() { }

   void  ChangeSelected(TGFrame*);
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGuiBldBorderFrame::TGuiBldBorderFrame(const TGWindow *p, TGuiBldEditor *editor) :
             TGVerticalFrame(p, 1, 1)
{
   fEditor = editor;
   fEditDisabled = 1;
   fBgndFrame = 0;

   SetCleanup(kDeepCleanup);

   fBtnGroup = new TGButtonGroup(this, "Border Mode");

   TGRadioButton *frame299 = new TGRadioButton(fBtnGroup," Sunken",kBldBorderSunken);
   frame299->SetToolTipText("Set a sunken border of the frame");
   TGRadioButton *frame302 = new TGRadioButton(fBtnGroup," None",kBldBorderPlain);
   frame302->SetToolTipText("Set no border of the frame");
   TGRadioButton *frame305 = new TGRadioButton(fBtnGroup," Raised",kBldBorderRaised);
   frame305->SetToolTipText("Set a raised border of the frame");
   frame305->SetState(kButtonDown);
   TGCheckButton *check = new TGCheckButton(fBtnGroup," Double",kBldBorderDouble);
   check->SetToolTipText("Set double border of the frame");
   //TQObject::Disconnect(check);

   fBtnGroup->SetRadioButtonExclusive(kTRUE);
   AddFrame(fBtnGroup, new TGLayoutHints(kLHintsCenterX | kLHintsTop));
   fBtnGroup->Connect("Pressed(Int_t)", "TGuiBldEditor", fEditor, "UpdateBorder(Int_t)");
   check->Connect("Pressed()", "TGuiBldEditor", fEditor, "UpdateBorder(=4)");
   check->Connect("Released()", "TGuiBldEditor", fEditor, "UpdateBorder(=5)");

   TGCompositeFrame *f = new TGGroupFrame(this, "Palette");
   TGHorizontalFrame *hf = new TGHorizontalFrame(f ,1, 1);
   fBgndFrame = new TGColorSelect(hf, 0, 1);
   fBgndFrame->SetEditDisabled();
   fBgndFrame->SetColor(GetDefaultFrameBackground());
   fBgndFrame->Connect("ColorSelected(Pixel_t)", "TGuiBldEditor", fEditor,
                       "UpdateBackground(Pixel_t)");
   hf->AddFrame(fBgndFrame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
   hf->AddFrame(new TGLabel(hf, "Backgrnd"), new TGLayoutHints(kLHintsTop | 
                kLHintsLeft, 2, 2, 2, 2));
   f->AddFrame(hf, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 2, 2, 2, 2));
   AddFrame(f, new TGLayoutHints(kLHintsCenterX | kLHintsTop));
}

////////////////////////////////////////////////////////////////////////////////
/// Perform actions when selected frame was changed.

void TGuiBldBorderFrame::ChangeSelected(TGFrame *frame)
{
   fSelected = frame;

   if (!frame) {
      return;
   }

   UInt_t opt = fSelected->GetOptions();

   fBtnGroup->SetButton(kBldBorderDouble, opt & kDoubleBorder);
   fBtnGroup->SetButton(kBldBorderSunken, opt & kSunkenFrame);
   fBtnGroup->SetButton(kBldBorderRaised, opt & kRaisedFrame);
   fBtnGroup->SetButton(kBldBorderPlain, !(opt & kRaisedFrame) && !(opt & kSunkenFrame));

   if (fBgndFrame) {
      TQObject::Disconnect(fBgndFrame);
      fBgndFrame->SetColor(fSelected->GetBackground());
      fBgndFrame->Connect("ColorSelected(Pixel_t)", "TGuiBldEditor", fEditor, "UpdateBackground(Pixel_t)");
   }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGuiBldEditor::TGuiBldEditor(const TGWindow *p) : TGVerticalFrame(p, 1, 1)
{
   TGHorizontalFrame *hf;
   TGVerticalFrame *vf;
   fSelected = 0;
   SetCleanup(kDeepCleanup);

   fNameFrame  = new TGuiBldNameFrame(this, this);
   AddFrame(fNameFrame,  new TGLayoutHints(kLHintsNormal | kLHintsExpandX,5,5,2,2));

   TGHSplitter *splitter = new TGHSplitter(this,100,5);
   AddFrame(splitter, new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,5,5));
   splitter->SetFrame(fNameFrame, kTRUE);

   //------------frame with layout switch
   hf = new TGHorizontalFrame(this);
   hf->AddFrame(new TGLabel(hf, "Composite Frame Layout"), 
                new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
   hf->AddFrame(new TGHorizontal3DLine(hf), new TGLayoutHints(kLHintsTop | 
                kLHintsExpandX, 2, 2, 2, 2));
   AddFrame(hf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
   
   vf = new TGVerticalFrame(this);
   fLayoutLabel = new TGLabel(vf, "Automatic Layout Disabled");
   vf->AddFrame(fLayoutLabel, new TGLayoutHints(kLHintsCenterX | kLHintsTop,
                2, 2, 2, 2));
   
   fLayoutButton = new TGTextButton(vf,"    Enable layout    ");
   fLayoutButton->SetEnabled(kFALSE);
   vf->AddFrame(fLayoutButton, new TGLayoutHints(kLHintsCenterX | kLHintsTop,
                2, 2, 2, 2));

   AddFrame(vf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

   AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsTop | 
            kLHintsExpandX, 2, 2, 2, 2));

   fLayoutButton->Connect("Clicked()", "TGuiBldEditor", this, "SwitchLayout()");
   fLayoutButton->SetToolTipText("If layout is on, all the frame \nelements get layouted automatically.");

   //-----------------------------

   fTab = new TGTab(this, 80, 40);
   AddFrame(fTab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
   fTablay = fTab->AddTab("Layout");
   TGCompositeFrame *tabcont = fTab->AddTab("Style");
   fLayoutId = 1; // 2nd tab
   fTab->Connect("Selected(Int_t)", "TGuiBldEditor", this, "TabSelected(Int_t)");

   fHintsFrame = new TGuiBldHintsEditor(fTablay, this);
   fTablay->AddFrame(fHintsFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
                     2, 2, 2, 2));

   fGeomFrame = new TGuiBldGeometryFrame(fTablay, this);
   fTablay->AddFrame(fGeomFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX,
                     2, 2, 2, 2));

   //----------------Position X,Y boxes---------------

   fPositionFrame = new TGGroupFrame(fTablay, "Position");

   hf = new TGHorizontalFrame(fPositionFrame);

   vf = new TGVerticalFrame(hf);
   vf->SetLayoutManager(new TGTableLayout(vf, 2, 2));

   vf->AddFrame(new TGLabel(vf, " X "), new TGTableLayoutHints(0, 1, 0, 1,
                kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));
   fXpos = new TGNumberEntry(vf, 0.0, 4, -1, (TGNumberFormat::EStyle)5);
   vf->AddFrame(fXpos, new TGTableLayoutHints(1, 2, 0, 1, kLHintsCenterY | 
                kLHintsLeft, 2, 2, 2, 2));

   vf->AddFrame(new TGLabel(vf, " Y "), new TGTableLayoutHints(0, 1, 1, 2,
                kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));
   fYpos = new TGNumberEntry(vf, 0.0, 4, -1, (TGNumberFormat::EStyle)5);
   vf->AddFrame(fYpos, new TGTableLayoutHints(1, 2, 1, 2, kLHintsCenterY | 
                kLHintsLeft, 2, 2, 2, 2));

   hf->AddFrame(vf, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));

   vf = new TGVerticalFrame(hf);
   vf->SetLayoutManager(new TGTableLayout(vf, 3, 3));

   TGTextButton *fTextButton6366 = new TGTextButton(vf, "^", -1, 
                TGButton::GetDefaultGC()(), 
                TGTextButton::GetDefaultFontStruct(),
                kRaisedFrame | kDoubleBorder | kFixedSize);
   fTextButton6366->Resize(20,20);
   vf->AddFrame(fTextButton6366, new TGTableLayoutHints(1, 2, 0, 1,
                kLHintsLeft | kLHintsTop, 1, 1, 1, 1));

   TGTextButton *fTextButton6367 = new TGTextButton(vf, "v", -1, 
                TGButton::GetDefaultGC()(), 
                TGTextButton::GetDefaultFontStruct(),
                kRaisedFrame | kDoubleBorder | kFixedSize);
   fTextButton6367->Resize(20,20);
   vf->AddFrame(fTextButton6367, new TGTableLayoutHints(1, 2, 2, 3,
                kLHintsLeft | kLHintsTop, 1, 1, 1, 1));

   TGTextButton *fTextButton6364 = new TGTextButton(vf, "<", -1, 
                TGButton::GetDefaultGC()(), 
                TGTextButton::GetDefaultFontStruct(),
                kRaisedFrame | kDoubleBorder | kFixedSize);
   fTextButton6364->Resize(20,20);
   vf->AddFrame(fTextButton6364, new TGTableLayoutHints(0, 1, 1, 2,
                kLHintsLeft | kLHintsTop, 1, 1, 1, 1));

   TGTextButton *fTextButton6365 = new TGTextButton(vf, ">", -1, 
                TGButton::GetDefaultGC()(), 
                TGTextButton::GetDefaultFontStruct(),
                kRaisedFrame | kDoubleBorder | kFixedSize);
   fTextButton6365->Resize(20,20);
   vf->AddFrame(fTextButton6365, new TGTableLayoutHints(2, 3, 1, 2,
                kLHintsLeft | kLHintsTop, 1, 1, 1, 1));

   hf->AddFrame(vf, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));

   fPositionFrame->AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   
   fTablay->AddFrame(fPositionFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   fXpos->Connect("ValueSet(Long_t)", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");
   fYpos->Connect("ValueSet(Long_t)", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");

   fTextButton6364->Connect("Clicked()", "TGNumberEntry", fXpos, "IncreaseNumber(TGNumberFormat::EStepSize=0,-1)");
   fTextButton6364->Connect("Clicked()", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");
   fTextButton6365->Connect("Clicked()", "TGNumberEntry", fXpos, "IncreaseNumber()");
   fTextButton6365->Connect("Clicked()", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");
   fTextButton6366->Connect("Clicked()", "TGNumberEntry", fYpos, "IncreaseNumber(TGNumberFormat::EStepSize=0,-1)");
   fTextButton6366->Connect("Clicked()", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");
   fTextButton6367->Connect("Clicked()", "TGNumberEntry", fYpos, "IncreaseNumber()");
   fTextButton6367->Connect("Clicked()", "TGuiBldHintsEditor", fHintsFrame, "SetPosition()");

   //----------------------------------------------------

   fBorderFrame = new TGuiBldBorderFrame(tabcont, this);
   tabcont->AddFrame(fBorderFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   MapSubwindows();
   SetWindowName("Frame Property Editor");
   SetEditDisabled(1);

   fEmbedded = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGuiBldEditor::~TGuiBldEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a frame.

void TGuiBldEditor::RemoveFrame(TGFrame *frame)
{
   fNameFrame->RemoveFrame(frame);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle  selected.

void TGuiBldEditor::TabSelected(Int_t id)
{
   if (id == fLayoutId) {
      //printf("%d\n", fSelected);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Hide editor.

void TGuiBldEditor::Hide()
{
   UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Change selected frame.

void TGuiBldEditor::ChangeSelected(TGFrame *frame)
{
   TGTabElement *tab = fTab->GetTabTab(fLayoutId);

   if (!frame) {
      fNameFrame->ChangeSelected(0);
      //fTab->SetTab(0);
      if (tab) {
         tab->SetEnabled(kFALSE);
         fClient->NeedRedraw(tab);
      }
      return;
   }

   fSelected = frame;
   TGWindow *parent = (TGWindow*)fSelected->GetParent();

   fNameFrame->ChangeSelected(fSelected);

   Bool_t enable_layout = kFALSE;
   enable_layout |= parent && !(parent->GetEditDisabled() & kEditDisableLayout);
   enable_layout |= !(fSelected->GetEditDisabled() & kEditDisableLayout);
   enable_layout |= parent && (parent->InheritsFrom(TGCompositeFrame::Class()) &&
                     !((TGCompositeFrame*)parent)->IsLayoutBroken());
   enable_layout |= (fSelected->InheritsFrom(TGCompositeFrame::Class()) &&
                     !((TGCompositeFrame*)fSelected)->IsLayoutBroken());

   if (enable_layout) {
      fHintsFrame->ChangeSelected(fSelected);

      if (tab) {
         tab->SetEnabled(kTRUE);
         fClient->NeedRedraw(tab);
      }
   } else {
      fHintsFrame->ChangeSelected(0);

      if (tab) {
         fTab->SetTab(0);
         tab->SetEnabled(kFALSE);
         fClient->NeedRedraw(tab);
      }
   }

   if ((frame->InheritsFrom(TGHorizontalFrame::Class())) ||
       (frame->InheritsFrom(TGVerticalFrame::Class())) ||
       (frame->InheritsFrom(TGGroupFrame::Class())) ) {

      fLayoutButton->SetEnabled(kTRUE);
      if (fSelected->IsLayoutBroken()) {
         fLayoutButton->SetText("    Enable layout    ");
         fLayoutLabel->SetText("Automatic layout disabled");
         if (fTablay) {
            fTablay->ShowFrame(fGeomFrame);
            fTablay->ShowFrame(fPositionFrame);
            fTablay->HideFrame(fHintsFrame);
         }
      } else {
         fLayoutButton->SetText("    Disable layout    ");
         fLayoutLabel->SetText("Automatic layout enabled");
         if (fTablay) {
            fTablay->HideFrame(fGeomFrame);
            fTablay->HideFrame(fPositionFrame);
            fTablay->ShowFrame(fHintsFrame);
         }
      }
   }
   else {
      fLayoutButton->SetEnabled(kFALSE);
      TGFrame *parentf = (TGFrame*)frame->GetParent();
      if (parentf->IsLayoutBroken()) {
         fLayoutButton->SetText("    Enable layout    ");
         fLayoutLabel->SetText("Automatic layout disabled");
         fTablay->ShowFrame(fGeomFrame);
         fTablay->ShowFrame(fPositionFrame);
         fTablay->HideFrame(fHintsFrame);
      } else {
         fLayoutButton->SetText("    Disable layout    ");
         fLayoutLabel->SetText("Automatic layout enabled");
         fTablay->HideFrame(fGeomFrame);
         fTablay->HideFrame(fPositionFrame);
         fTablay->ShowFrame(fHintsFrame);
      }
   }

   fYpos->SetIntNumber(frame->GetY());
   fXpos->SetIntNumber(frame->GetX());

   if (fBorderFrame) fBorderFrame->ChangeSelected(fSelected);
   if (fGeomFrame) fGeomFrame->ChangeSelected(fSelected);

   Emit("ChangeSelected(TGFrame*)", (long)fSelected);

   MapRaised();
}

////////////////////////////////////////////////////////////////////////////////
/// Update selected frame.

void TGuiBldEditor::UpdateSelected(TGFrame *frame)
{
   Emit("UpdateSelected(TGFrame*)", (long)frame);
}

////////////////////////////////////////////////////////////////////////////////
/// Update border of selected frame.

void TGuiBldEditor::UpdateBorder(Int_t b)
{
   if (!fSelected) return;

   UInt_t opt = fSelected->GetOptions();

   switch (b) {
      case 1:
         opt &= ~kRaisedFrame;
         opt |= kSunkenFrame;
         break;
      case 2:
         opt &= ~kSunkenFrame;
         opt &= ~kRaisedFrame;
         break;
      case 3:
         opt &= ~kSunkenFrame;
         opt |= kRaisedFrame;
         break;
      case 4:
         opt |= kDoubleBorder;
         break;
      case 5:
         opt &= ~kDoubleBorder;
         break;
      default:
         return;
   }
   fSelected->ChangeOptions(opt);
   fClient->NeedRedraw(fSelected, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Update background.

void TGuiBldEditor::UpdateBackground(Pixel_t col)
{
   if (!fSelected) return;

   fSelected->SetBackgroundColor(col);
   fClient->NeedRedraw(fSelected, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Update foreground.

void TGuiBldEditor::UpdateForeground(Pixel_t col)
{
   if (!fSelected) return;

   fSelected->SetForegroundColor(col);
   fClient->NeedRedraw(fSelected, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the editor.

void TGuiBldEditor::Reset()
{
   fSelected = 0;
   fNameFrame->Reset();
   TGTabElement *tab = fTab->GetTabTab(fLayoutId);
   if (tab) {
      fTab->SetTab(0);
      tab->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Popup dialog to set layout of editted frame off. If layout is on, all
/// the elements in the frame get layouted automatically.

void TGuiBldEditor::SwitchLayout()
{
   if (!fSelected) {
      fLayoutButton->SetText("    Enable layout    ");
      fLayoutButton->SetEnabled(kFALSE);
      fLayoutLabel->SetText("Automatic layout disabled");
      if (fTablay) {
         fTablay->ShowFrame(fGeomFrame);
         fTablay->ShowFrame(fPositionFrame);
         fTablay->HideFrame(fHintsFrame);
      }
      return;
   }

   TRootGuiBuilder *builder = (TRootGuiBuilder*)TRootGuiBuilder::Instance();
   TGFrame *frame = fSelected;
   TGCompositeFrame *cf = fNameFrame->GetMdi(frame);
   if (cf == 0)
      return;
   if (frame->IsLayoutBroken()) {
      Int_t retval;
      builder->GetManager()->SetEditable(kFALSE);
      new TGMsgBox(gClient->GetDefaultRoot(), builder, "Layout change",
                   "Enabling layout will automatically align and resize all the icons. \n Do you really want to layout them?",
                   kMBIconExclamation, kMBOk | kMBCancel, &retval);

      cf->SetEditable(kTRUE);
      // hack against selecting the message box itself
      builder->GetManager()->SelectFrame(frame);
      frame->SetEditable(kTRUE);

      if (retval == kMBOk) {
         frame->SetLayoutBroken(kFALSE);
         frame->Layout();
         fLayoutButton->SetText("    Disable layout    ");
         fLayoutLabel->SetText("Automatic layout enabled");
         if (fTablay) {
            fTablay->HideFrame(fGeomFrame);
            fTablay->HideFrame(fPositionFrame);
            fTablay->ShowFrame(fHintsFrame);
            fTablay->Resize(fHintsFrame->GetWidth(),fHintsFrame->GetHeight());
         }
      }
   } else {
      //set layout off - without dialog, because nothing "bad" can happen
      frame->SetLayoutBroken(kTRUE);
      fLayoutButton->SetText("    Enable layout    ");
      fLayoutLabel->SetText("Automatic layout disabled");
      if (fTablay) {
         fTablay->ShowFrame(fGeomFrame);
         fTablay->ShowFrame(fPositionFrame);
         fTablay->HideFrame(fHintsFrame);
      }
   }
   fClient->NeedRedraw(frame, kTRUE);
   if (fTablay) fClient->NeedRedraw(fTablay, kTRUE);
}



