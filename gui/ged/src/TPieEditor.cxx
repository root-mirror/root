// @(#)root/ged:$Id$
// Author: Guido Volpi 12/10/2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TPieEditor                                                          //
//                                                                      //
//  Implements GUI for pie-chart attributes.                            //
//                                                                      //
//  Title': set the title of the graph                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TPieEditor.gif">
*/
//End_Html

#include "TGedEditor.h"
#include "TGComboBox.h"
#include "TGButtonGroup.h"
#include "TPieEditor.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TPie.h"
#include "TVirtualPad.h"
#include "TGColorSelect.h"
#include "TColor.h"
#include "TBox.h"
#include "TPaveLabel.h"
#include "TVirtualX.h"
#include "snprintf.h"

ClassImp(TPieEditor);

enum EPieWid {
   kPie = 0,
   kPIE_HOR,
   kPIE_RAD,
   kPIE_TAN,
   kPIE_FILL,
   kPIE_OUTLINE,
   kPIE_TITLE,
   kPIE_3D,
   kPIE_3DANGLE,
   kPIE_3DTHICKNESS,
   kFONT_COLOR,
   kFONT_SIZE,
   kFONT_STYLE
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor of graph editor.

TPieEditor::TPieEditor(const TGWindow *p, Int_t width,
                         Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fPie = 0;
   // TextEntry to change the title
   MakeTitle("Pie Chart");

   fTitle = new TGTextEntry(this, new TGTextBuffer(50), kPIE_TITLE);
   fTitle->Resize(135, fTitle->GetDefaultHeight());
   fTitle->SetToolTipText("Enter the pie title string");
   // better take kLHintsLeft and Right - Right is not working at the moment
   AddFrame(fTitle, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   // Radio Buttons to change the draw options of the graph
   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   fgr = new TGButtonGroup(f2,3,1,3,5,"Label direction");
   fgr->SetRadioButtonExclusive(kTRUE);
   fLblDirH = new TGRadioButton(fgr,"Horizontal",kPIE_HOR);   // no draw option
   fLblDirH->SetToolTipText("Draw horizontal labels");
   fLblDirR = new TGRadioButton(fgr,"Radial",kPIE_RAD);  // option C
   fLblDirR->SetToolTipText("Draw labels radially");
   fLblDirT = new TGRadioButton(fgr,"Tangential",kPIE_TAN); // option L
   fLblDirT->SetToolTipText("Draw labels tangential to the piechart");

   fgr->SetLayoutHints(fShape1lh=new TGLayoutHints(kLHintsLeft, 0,3,0,0), fLblDirH);
   fgr->Show();
   fgr->ChangeOptions(kFitWidth|kChildFrame|kVerticalFrame);
   f2->AddFrame(fgr, new TGLayoutHints(kLHintsLeft, 4, 0, 0, 0));

   // CheckBox to activate/deactivate the drawing of the Marker
   fOutlineOnOff = new TGCheckButton(f2,"Outline",kPIE_OUTLINE);
   fOutlineOnOff->SetToolTipText("Draw a line to mark the pie");
   f2->AddFrame(fOutlineOnOff, new TGLayoutHints(kLHintsTop, 5, 1, 0, 3));
   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   // Exclusion zone parameters
   MakeTitle("3D options");
   TGCompositeFrame *f3 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f3, new TGLayoutHints(kLHintsTop, 1, 1, 5, 0));

   fIs3D = new TGCheckButton(f3,"3D",kPIE_3D);
   fIs3D->SetToolTipText("Draw a 3D charts");
   f3->AddFrame(fIs3D, new TGLayoutHints(kLHintsTop, 5, 1, 0, 0));

   f3DAngle = new TGNumberEntry(f3, 0, 2, kPIE_3DANGLE,  TGNumberEntry::kNESInteger, TGNumberEntry::kNEANonNegative,TGNumberFormat::kNELLimitMinMax, 0, 90);
   //f3DAngle->SetToolTipText("3D angle: 0-90")
   f3DAngle->Resize(50, 20);
   f3->AddFrame(f3DAngle, new TGLayoutHints(kLHintsLeft, 7, 1, 1, 1));
   f3DAngle->Associate(f3);

   f3DHeight = new TGNumberEntry(f3, 0, 3, kPIE_3DTHICKNESS,  TGNumberEntry::kNESReal, TGNumberEntry::kNEANonNegative);
   //f3DHeight->SetToolTipText("3D thick")
   f3DHeight->Resize(50, 20);
   f3->AddFrame(f3DHeight, new TGLayoutHints(kLHintsLeft, 7, 1, 1, 1));
   f3DHeight->Associate(f3);

   MakeTitle("Text");
   TGCompositeFrame *f4 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fColorSelect = new TGColorSelect(f4, 0, kFONT_COLOR);
   f4->AddFrame(fColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fColorSelect->Associate(this);
   fSizeCombo = BuildFontSizeComboBox(f4, kFONT_SIZE);
   f4->AddFrame(fSizeCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   fSizeCombo->Resize(91, 20);
   fSizeCombo->Associate(this);
   AddFrame(f4, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
   fTypeCombo = new TGFontTypeComboBox(this, kFONT_STYLE);
   fTypeCombo->Resize(137, 20);
   AddFrame(fTypeCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor of pie editor.

TPieEditor::~TPieEditor()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TPieEditor::ConnectSignals2Slots()
{
   fTitle->Connect("TextChanged(const char *)","TPieEditor",this,"DoTitle(const char *)");
   fgr->Connect("Clicked(Int_t)","TPieEditor",this,"DoShape()");
   fOutlineOnOff->Connect("Toggled(Bool_t)","TPieEditor",this,"DoMarkerOnOff(Bool_t)");
   f3DAngle->Connect("ValueSet(Long_t)", "TPieEditor", this, "DoChange3DAngle()");
   f3DHeight->Connect("ValueSet(Long_t)", "TPieEditor", this, "DoChange3DAngle()");
   fIs3D->Connect("Clicked()","TPieEditor",this,"DoGraphLineWidth()");

   // text attributes connection
   fTypeCombo->Connect("Selected(Int_t)","TPieEditor",this,"DoTextChange()");
   fSizeCombo->Connect("Selected(Int_t)","TPieEditor",this,"DoTextChange()");
   fColorSelect->Connect("ColorSelected(Pixel_t)","TPieEditor",this,"DoTextChange()");

   fInit = kFALSE;  // connect the slots to the signals only once
}


////////////////////////////////////////////////////////////////////////////////
/// Exclude TAttTextEditor from this interface.

void TPieEditor::ActivateBaseClassEditors(TClass* cl)
{
   TGedEditor *gedEditor = GetGedEditor();
   gedEditor->ExcludeClassEditor(TAttText::Class());
   TGedFrame::ActivateBaseClassEditors(cl);
}


////////////////////////////////////////////////////////////////////////////////
/// Pick up the used values of graph attributes.

void TPieEditor::SetModel(TObject* obj)
{
   fPie = (TPie *)obj;
   fAvoidSignal = kTRUE;

   // set the Title TextEntry
   const char *text = fPie->GetTitle();
   fTitle->SetText(text);

   TString soption = GetDrawOption();

   // For the label orientation there are 3 possibilities:
   //   0: horizontal
   //   1: radial
   //   2: tangent
   Int_t lblor(0);

   // Parse the options
   Int_t idx;
   // Paint the TPie in an existing canvas
   if ( (idx=soption.Index("same"))>=0 ) {
      // optionSame = kTRUE;
      soption.Remove(idx,4);
   }

   if ( (idx=soption.Index("nol"))>=0 ) {
      fOutlineOnOff->SetState(kButtonUp,kFALSE);
      soption.Remove(idx,3);
   }
   else {
      fOutlineOnOff->SetState(kButtonDown,kFALSE);
   }

   // check if is active the pseudo-3d
   if ( (idx=soption.Index("3d"))>=0 ) {
      fIs3D->SetState(kButtonDown, kFALSE);
      f3DAngle->SetNumber(fPie->GetAngle3D());
      f3DHeight->SetNumber(fPie->GetHeight());
      soption.Remove(idx,2);
   } else {
      fIs3D->SetState(kButtonUp, kFALSE);
   }

   // seek if have to draw the labels around the pie chart
   if ( (idx=soption.Index("t"))>=0 ) {
      lblor = 2;
      soption.Remove(idx,1);
   }

   // Seek if have to paint the labels along the radii
   if ( (idx=soption.Index("r"))>=0 ) {
      lblor = 1;
      soption.Remove(idx,1);
   }

   switch(lblor) {
   case 0:
      fLblDirH->SetState(kButtonDown,kTRUE);
      break;
   case 1:
      fLblDirR->SetState(kButtonDown,kTRUE);
      break;
   case 2:
      fLblDirT->SetState(kButtonDown,kTRUE);
      break;
   }

   // set text attributes
   fTypeCombo->Select(fPie->GetTextFont() / 10);

   Color_t c = fPie->GetTextColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fColorSelect->SetColor(p, kFALSE);

   Float_t s = fPie->GetTextSize();
   Float_t dy;

   if (obj->InheritsFrom(TPaveLabel::Class())) {
      TBox *pl = (TBox*)obj;
      dy = s * (pl->GetY2() - pl->GetY1());
   }
   else
      dy = s * (fGedEditor->GetPad()->GetY2() - fGedEditor->GetPad()->GetY1());

   Int_t size = fGedEditor->GetPad()->YtoPixel(0.0) - fGedEditor->GetPad()->YtoPixel(dy);
   if (size > 50) size = 50;
   if (size < 0)  size = 0;
   fSizeCombo->Select(size, kFALSE);

   if (fInit) ConnectSignals2Slots();
   fAvoidSignal = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Slot for setting the graph title.

void TPieEditor::DoTitle(const char *text)
{
   if (fAvoidSignal) return;
   fPie->SetTitle(text);
   Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the draw options.

void TPieEditor::DoShape()
{
   if (fAvoidSignal) return;

   TString opt = GetDrawOption();

   if (fLblDirH->GetState()==kButtonDown) {
      if (opt.Contains("t")) opt.Remove(opt.First("t"),1);
      if (opt.Contains("r")) opt.Remove(opt.First("r"),1);
   }
   else if (fLblDirR->GetState()==kButtonDown) {
      if (opt.Contains("t")) opt.Remove(opt.First("t"),1);
      if (!opt.Contains("r")) opt += "r";
   }
   else if (fLblDirT->GetState()==kButtonDown) {
      if (!opt.Contains("t")) opt += "t";
      if (opt.Contains("r")) opt.Remove(opt.First("r"),1);
   }

   SetDrawOption(opt);
   if (gPad && gPad->GetVirtCanvas()) gPad->GetVirtCanvas()->SetCursor(kPointer);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kPointer));
}


////////////////////////////////////////////////////////////////////////////////
/// Slot for setting markers as visible/invisible.

void TPieEditor::DoMarkerOnOff(Bool_t)
{
   if (fAvoidSignal) return;
   TString t = GetDrawOption();

   if (t.Contains("nol") && fOutlineOnOff->GetState() == kButtonDown) {
      t.Remove(t.First("nol"),3);
   }
   else if (!t.Contains("nol") && fOutlineOnOff->GetState() == kButtonUp) {
      t += "nol";
   }

   SetDrawOption(t);
}


////////////////////////////////////////////////////////////////////////////////
/// Slot for setting the 3D angle

void TPieEditor::DoChange3DAngle()
{
   if (fAvoidSignal) return;

   fPie->SetAngle3D(static_cast<Int_t>(f3DAngle->GetNumber()));
   fPie->SetHeight(f3DHeight->GetNumber());

   Update();

}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the graph line width.

void TPieEditor::DoGraphLineWidth()
{
   if (fAvoidSignal) return;

   TString opt = GetDrawOption();
   if (!opt.Contains("3d") && fIs3D->IsDown())
      opt += "3d";
   else if (opt.Contains("3d") && !fIs3D->IsDown())
      opt.Remove(opt.First("3d"),2);

   SetDrawOption(opt);

   Update();
}



////////////////////////////////////////////////////////////////////////////////
/// Change text.

void TPieEditor::DoTextChange()
{
   if (fAvoidSignal) return;

   // font color
   fPie->SetTextColor(TColor::GetColor(fColorSelect->GetColor()));

   // font type
   Int_t fontPrec = fPie->GetTextFont()%10;
   Int_t fontType = fTypeCombo->GetSelected();
   fPie->SetTextFont(fontType*10+fontPrec);

   // font size
   TVirtualPad* pad = fGedEditor->GetPad();

   Float_t val = TString(fSizeCombo->GetSelectedEntry()->GetTitle()).Atoi();

   Float_t dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY((Int_t)val);
   Float_t textSize;

   if (fGedEditor->GetModel()->InheritsFrom(TPaveLabel::Class())) {
      TBox *pl = (TBox*)fGedEditor->GetModel();
      textSize = dy/(pl->GetY2() - pl->GetY1());
   }
   else
      textSize = dy/(pad->GetY2() - pad->GetY1());

   fPie->SetTextSize(textSize);

   Update();

}


////////////////////////////////////////////////////////////////////////////////
/// Create text size combo box.

TGComboBox* TPieEditor::BuildFontSizeComboBox(TGFrame* parent, Int_t id)
{
   char a[100];
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("Default", 0);
   for (int i = 1; i <= 50; i++) {
      snprintf(a, 99, "%d", i);
      c->AddEntry(a, i);
   }

   return c;
}
