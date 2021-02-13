// @(#)root/treeviewer:$Id$
// Author: Bastien Dalla Piazza  02/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TParallelCoordEditor.h"
#include "TParallelCoord.h"
#include "TParallelCoordRange.h"
#include "TParallelCoordVar.h"

#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGComboBox.h"
#include "TGColorSelect.h"
#include "TColor.h"
#include "TG3DLine.h"
#include "TGSlider.h"
#include "TGDoubleSlider.h"
#include "TGedPatternSelect.h"
#include "TCanvas.h"

#include "TROOT.h"

ClassImp(TParallelCoordEditor);


/** \class TParallelCoordEditor

This is the TParallelCoord editor. It brings tools to explore datas
Using parallel coordinates. The main tools are:

  - Dots spacing : Set the dots spacing with whichone the lines
    must be drawn. This tool is useful to reduce the image
    cluttering.
  - The Selections section : Set the current edited selection and
    allows to apply it to the tree through a generated entry list.
  - The Entries section : Set how many events must be drawn.
    A weight cut can be defioned here (see TParallelCoord for a
    a description of the weight cut).
  - The Variables tab : To define the global settings to display
    the axes. It is also possible to add a variable from its
    expression or delete a selected one (also possible using right
    click on the pad.
*/

enum EParallelWid {
   kGlobalLineColor,
   kLineTypeBgroup,
   kLineTypePoly,
   kLineTypeCurves,
   kGlobalLineWidth,
   kDotsSpacing,
   kDotsSpacingField,
   kAlpha,
   kAlphaField,
   kSelectionSelect,
   kSelectLineColor,
   kSelectLineWidth,
   kActivateSelection,
   kDeleteSelection,
   kAddSelection,
   kAddSelectionEntry,
   kShowRanges,
   kPaintEntries,
   kEntriesToDraw,
   kFirstEntry,
   kNentries,
   kApplySelect,
   kUnApply,
   kDelayDrawing,
   kHideAllRanges,
   kVariables,
   kDeleteVar,
   kHistHeight,
   kHistWidth,
   kHistBinning,
   kRenameVar,
   kWeightCut,
   kHistColorSelect,
   kHistPatternSelect
};

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TParallelCoordEditor::TParallelCoordEditor(const TGWindow* /*p*/,
                                           Int_t/*width*/, Int_t /*height*/,
                                           UInt_t /*options*/, Pixel_t /*back*/)
{
   fPriority = 1;
   fDelay = kTRUE;

   // Line
   MakeTitle("Line");

   TGHorizontalFrame *f1 = new TGHorizontalFrame(this);
   fGlobalLineColor = new TGColorSelect(f1,0,kGlobalLineColor);
   f1->AddFrame(fGlobalLineColor,new TGLayoutHints(kLHintsLeft | kLHintsTop));
   fGlobalLineWidth = new TGLineWidthComboBox(f1, kGlobalLineWidth);
   fGlobalLineWidth->Resize(91, 20);
   f1->AddFrame(fGlobalLineWidth, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsTop));

   if (!TCanvas::SupportAlpha()) {

      AddFrame(new TGLabel(this,"Dots spacing"),
               new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

      TGHorizontalFrame *f2 = new TGHorizontalFrame(this);
      fDotsSpacing = new TGHSlider(f2,100,kSlider2|kScaleNo,kDotsSpacing);
      fDotsSpacing->SetRange(0,60);
      f2->AddFrame(fDotsSpacing,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
      fDotsSpacingField = new TGNumberEntryField(f2, kDotsSpacingField, 0,
                                                 TGNumberFormat::kNESInteger,
                                                 TGNumberFormat::kNEANonNegative);
      fDotsSpacingField->Resize(40,20);
      f2->AddFrame(fDotsSpacingField,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
      AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   }
   else {
      TGLabel *AlphaLabel = new TGLabel(this,"Opacity");
      AddFrame(AlphaLabel,
               new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
      TGHorizontalFrame *f2a = new TGHorizontalFrame(this);
      fAlpha = new TGHSlider(f2a,100,kSlider2|kScaleNo,kAlpha);
      fAlpha->SetRange(0,1000);
      f2a->AddFrame(fAlpha,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
      fAlphaField = new TGNumberEntryField(f2a, kAlphaField, 0,
                                           TGNumberFormat::kNESReal,
                                           TGNumberFormat::kNEANonNegative);
      fAlphaField->Resize(40,20);
      f2a->AddFrame(fAlphaField,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
      AddFrame(f2a, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   }

   fLineTypeBgroup = new TGButtonGroup(this,2,1,0,0, "Line type");
   fLineTypeBgroup->SetRadioButtonExclusive(kTRUE);
   fLineTypePoly = new TGRadioButton(fLineTypeBgroup,"Polyline", kLineTypePoly);
   fLineTypePoly->SetToolTipText("Draw the entries with a polyline");
   fLineTypeCurves = new TGRadioButton(fLineTypeBgroup,"Curves",
                                       kLineTypeCurves);
   fLineTypeCurves->SetToolTipText("Draw the entries with a curve");
   fLineTypeBgroup->ChangeOptions(kChildFrame|kVerticalFrame);
   AddFrame(fLineTypeBgroup, new TGLayoutHints(kLHintsCenterY | kLHintsLeft));

   // Selections
   MakeTitle("Selections");

   fHideAllRanges = new TGCheckButton(this,"Hide all ranges",kHideAllRanges);
   AddFrame(fHideAllRanges);

   fSelectionSelect = new TGComboBox(this,kSelectionSelect);
   fSelectionSelect->Resize(140,20);
   AddFrame(fSelectionSelect, new TGLayoutHints(kLHintsCenterY | kLHintsLeft,0,0,3,0));

   TGHorizontalFrame *f3 = new TGHorizontalFrame(this);
   fSelectLineColor = new TGColorSelect(f3,0,kSelectLineColor);
   f3->AddFrame(fSelectLineColor,new TGLayoutHints(kLHintsLeft | kLHintsTop));
   fSelectLineWidth = new TGLineWidthComboBox(f3, kSelectLineWidth);
   fSelectLineWidth->Resize(94, 20);
   f3->AddFrame(fSelectLineWidth, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,0,3,0));

   fActivateSelection = new TGCheckButton(this,"Activate",kActivateSelection);
   fActivateSelection->SetToolTipText("Activate the current selection");
   AddFrame(fActivateSelection);
   fShowRanges = new TGCheckButton(this,"Show ranges",kShowRanges);
   AddFrame(fShowRanges);

   TGHorizontalFrame *f5 = new TGHorizontalFrame(this);
   fAddSelectionField = new TGTextEntry(f5);
   fAddSelectionField->Resize(57,20);
   f5->AddFrame(fAddSelectionField, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   fAddSelection = new TGTextButton(f5,"Add");
   fAddSelection->SetToolTipText("Add a new selection (Right click on the axes to add a range).");
   f5->AddFrame(fAddSelection, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,3,0,0,0));
   fDeleteSelection = new TGTextButton(f5,"Delete",kDeleteSelection);
   fDeleteSelection->SetToolTipText("Delete the current selection");
   f5->AddFrame(fDeleteSelection, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,3,0,0,0));
   AddFrame(f5, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,3,0,0,0));

   TGHorizontalFrame *f7 = new TGHorizontalFrame(this);
   fApplySelect = new TGTextButton(f7,"Apply to tree",kApplySelect);
   fApplySelect->SetToolTipText("Generate an entry list for the current selection and apply it to the tree.");
   f7->AddFrame(fApplySelect);
   fUnApply = new TGTextButton(f7,"Reset tree",kUnApply);
   fUnApply->SetToolTipText("Reset the tree entry list");
   f7->AddFrame(fUnApply, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,10,0,0,0));
   AddFrame(f7, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,0,0,3,0));

   // Entries
   MakeTitle("Entries");

   fPaintEntries = new TGCheckButton(this,"Draw entries",kPaintEntries);
   AddFrame(fPaintEntries);
   fDelayDrawing = new TGCheckButton(this,"Delay Drawing", kDelayDrawing);
   AddFrame(fDelayDrawing);

   fEntriesToDraw = new TGDoubleHSlider(this,140,kDoubleScaleNo,kEntriesToDraw);
   AddFrame(fEntriesToDraw);

   TGHorizontalFrame *f6 = new TGHorizontalFrame(this);
   TGVerticalFrame *v1 = new TGVerticalFrame(f6);
   TGVerticalFrame *v2 = new TGVerticalFrame(f6);
   v1->AddFrame(new TGLabel(v1,"First entry:"));
   fFirstEntry = new TGNumberEntryField(v1, kFirstEntry, 0,
                                        TGNumberFormat::kNESInteger,
                                        TGNumberFormat::kNEANonNegative);
   fFirstEntry->Resize(68,20);
   v1->AddFrame(fFirstEntry);
   v2->AddFrame(new TGLabel(v2,"# of entries:"));
   fNentries = new TGNumberEntryField(v2, kFirstEntry, 0,
                                     TGNumberFormat::kNESInteger,
                                     TGNumberFormat::kNEANonNegative);
   fNentries->Resize(68,20);
   v2->AddFrame(fNentries);
   f6->AddFrame(v1);
   f6->AddFrame(v2, new TGLayoutHints(kLHintsLeft,4,0,0,0));
   AddFrame(f6);

   AddFrame(new TGLabel(this,"Weight cut:"));

   TGHorizontalFrame *f8 = new TGHorizontalFrame(this);
   fWeightCut = new TGHSlider(f8,100,kSlider2|kScaleNo,kDotsSpacing);
   fWeightCutField = new TGNumberEntryField(f8,kDotsSpacingField, 0,
                                            TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative);
   fWeightCutField->Resize(40,20);
   f8->AddFrame(fWeightCut);
   f8->AddFrame(fWeightCutField);
   AddFrame(f8);

   MakeVariablesTab();
}

////////////////////////////////////////////////////////////////////////////////
/// Make the "variable" tab.

void TParallelCoordEditor::MakeVariablesTab()
{
   fVarTab = CreateEditorTabSubFrame("Variables");
   // Variable

   TGHorizontalFrame *f9 = new TGHorizontalFrame(fVarTab);
   fAddVariable = new TGTextEntry(f9);
   fAddVariable->Resize(71,20);
   f9->AddFrame(fAddVariable, new TGLayoutHints(kLHintsCenterY));
   fButtonAddVar = new TGTextButton(f9,"Add");
   fButtonAddVar->SetToolTipText("Add a new variable from the tree (must be a valid expression).");
   f9->AddFrame(fButtonAddVar, new TGLayoutHints(kLHintsCenterY,4,0,0,0));
   fVarTab->AddFrame(f9);

   TGHorizontalFrame *f10 = new TGHorizontalFrame(fVarTab);
   fVariables = new TGComboBox(f10,kVariables);
   fVariables->Resize(105,20);
   f10->AddFrame(fVariables, new TGLayoutHints(kLHintsCenterY));
   fVarTab->AddFrame(f10,new TGLayoutHints(kLHintsLeft,0,0,2,0));

   TGHorizontalFrame *f12 = new TGHorizontalFrame(fVarTab);
   fDeleteVar = new TGTextButton(f12,"Delete",kDeleteVar);
   fDeleteVar->SetToolTipText("Delete the current selected variable");
   f12->AddFrame(fDeleteVar, new TGLayoutHints(kLHintsCenterY,1,0,0,0));
   fRenameVar = new TGTextButton(f12,"Rename",kRenameVar);
   fRenameVar->SetToolTipText("Rename the current selected variable");
   f12->AddFrame(fRenameVar, new TGLayoutHints(kLHintsCenterY,4,0,0,0));
   fVarTab->AddFrame(f12,new TGLayoutHints(kLHintsLeft,0,0,2,0));

   fVarTab->AddFrame(new TGLabel(fVarTab,"Axis histograms:"));

   TGHorizontalFrame *f11 = new TGHorizontalFrame(fVarTab);
   TGVerticalFrame *v3 = new TGVerticalFrame(f11);
   TGVerticalFrame *v4 = new TGVerticalFrame(f11);
   v3->AddFrame(new TGLabel(v3,"Binning:"));
   fHistBinning = new TGNumberEntryField(v3, kHistWidth, 0,
                                         TGNumberFormat::kNESInteger,
                                         TGNumberFormat::kNEANonNegative);
   fHistBinning->Resize(68,20);
   v3->AddFrame(fHistBinning);
   v4->AddFrame(new TGLabel(v4,"Width:"));
   fHistWidth = new TGNumberEntryField(v4, kHistWidth, 0,
                                       TGNumberFormat::kNESInteger,
                                       TGNumberFormat::kNEANonNegative);
   fHistWidth->Resize(68,20);
   v4->AddFrame(fHistWidth, new TGLayoutHints(kLHintsLeft,4,0,0,0));
   f11->AddFrame(v3);
   f11->AddFrame(v4);
   fVarTab->AddFrame(f11);

   fHistShowBoxes = new TGCheckButton(fVarTab,"Show box histograms");
   fVarTab->AddFrame(fHistShowBoxes);

   fVarTab->AddFrame(new TGLabel(fVarTab,"Bar histograms style:"));

   TGCompositeFrame *f13 = new TGCompositeFrame(fVarTab, 80, 20, kHorizontalFrame);
   fHistColorSelect = new TGColorSelect(f13, 0, kHistColorSelect);
   f13->AddFrame(fHistColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fHistColorSelect->Associate(this);
   fHistPatternSelect = new TGedPatternSelect(f13, 1, kHistPatternSelect);
   f13->AddFrame(fHistPatternSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fHistPatternSelect->Associate(this);
   fVarTab->AddFrame(f13, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TParallelCoordEditor::~TParallelCoordEditor()
{
   delete fLineTypePoly;
   delete fLineTypeCurves;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up the selection combo box.

void TParallelCoordEditor::CleanUpSelections()
{
   TList *list = fParallel->GetSelectList();
   fSelectionSelect->RemoveAll();
   Bool_t enable = list->GetSize() > 0;
   fSelectionSelect->SetEnabled(enable);
   fSelectLineColor->SetEnabled(enable);
   fSelectLineWidth->SetEnabled(enable);
   fActivateSelection->SetEnabled(enable);
   fShowRanges->SetEnabled(enable);
   fDeleteSelection->SetEnabled(enable);
   if (list->GetSize() > 0) {
      Int_t i = 0;
      TIter next(list);
      TParallelCoordSelect* sel;
      while ((sel = (TParallelCoordSelect*)next())) {
         fSelectionSelect->AddEntry(sel->GetTitle(),i);
         TGLBEntry *entry = fSelectionSelect->GetListBox()->GetEntry(i);
         if (entry)
            entry->SetBackgroundColor(TColor::Number2Pixel(sel->GetLineColor()));
         ++i;
      }
      sel = fParallel->GetCurrentSelection();
      if (sel) {
         fSelectionSelect->Select(list->IndexOf(sel),kFALSE);
         Color_t c;
         Pixel_t p;
         c = sel->GetLineColor();
         p = TColor::Number2Pixel(c);
         fSelectLineColor->SetColor(p);
         fSelectLineWidth->Select(sel->GetLineWidth());
         fActivateSelection->SetOn(sel->TestBit(TParallelCoordSelect::kActivated));
         fShowRanges->SetOn(sel->TestBit(TParallelCoordSelect::kShowRanges));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up the variables combo box.

void TParallelCoordEditor::CleanUpVariables()
{
   TList *list = fParallel->GetVarList();
   fVariables->RemoveAll();
   Bool_t enable = list->GetSize() > 0;
   fVariables->SetEnabled(enable);
   fDeleteVar->SetEnabled(enable);
   fHistShowBoxes->SetEnabled(enable);
   fHistWidth->SetEnabled(enable);
   fHistBinning->SetEnabled(enable);
   if (list->GetSize() > 0) {
      Int_t i = 0;
      TIter next(list);
      TParallelCoordVar* var;
      while ((var = (TParallelCoordVar*)next())) {
         fVariables->AddEntry(var->GetTitle(),i);
         ++i;
      }
      var = (TParallelCoordVar*)list->First();
      fVariables->Select(0,kFALSE);
      fHistShowBoxes->SetOn(var->TestBit(TParallelCoordVar::kShowBarHisto));
      fHistWidth->SetNumber(var->GetHistLineWidth());
      fHistBinning->SetNumber(var->GetHistBinning());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TParallelCoordEditor::ConnectSignals2Slots()
{
   fGlobalLineColor->Connect("ColorSelected(Pixel_t)","TParallelCoordEditor",
                             this, "DoGlobalLineColor(Pixel_t)");
   fGlobalLineWidth->Connect("Selected(Int_t)","TParallelCoordEditor",
                             this, "DoGlobalLineWidth(Int_t)");
   if (!TCanvas::SupportAlpha()) {
      fDotsSpacing->Connect("Released()","TParallelCoordEditor",
                           this, "DoDotsSpacing()");
      fDotsSpacing->Connect("PositionChanged(Int_t)","TParallelCoordEditor",
                           this, "DoLiveDotsSpacing(Int_t)");
      fDotsSpacingField->Connect("ReturnPressed()","TParallelCoordEditor",
                                 this, "DoDotsSpacingField()");
   }
   else {
      fAlpha->Connect("Released()","TParallelCoordEditor",
                           this, "DoAlpha()");
      fAlpha->Connect("PositionChanged(Int_t)","TParallelCoordEditor",
                           this, "DoLiveAlpha(Int_t)");
      fAlphaField->Connect("ReturnPressed()","TParallelCoordEditor",
                           this, "DoAlphaField()");
   }
   fLineTypeBgroup->Connect("Clicked(Int_t)", "TParallelCoordEditor",
                            this, "DoLineType()");
   fSelectionSelect->Connect("Selected(const char*)","TParallelCoordEditor",
                            this, "DoSelectionSelect(const char*)");
   fSelectLineColor->Connect("ColorSelected(Pixel_t)","TParallelCoordEditor",
                             this, "DoSelectLineColor(Pixel_t)");
   fSelectLineWidth->Connect("Selected(Int_t)","TParallelCoordEditor",
                             this, "DoSelectLineWidth(Int_t)");
   fActivateSelection->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                               this, "DoActivateSelection(Bool_t)");
   fShowRanges->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                        this, "DoShowRanges(Bool_t)");
   fDeleteSelection->Connect("Clicked()","TParallelCoordEditor",
                             this, "DoDeleteSelection()");
   fAddSelection->Connect("Clicked()","TParallelCoordEditor",
                             this, "DoAddSelection()");
   fPaintEntries->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                               this, "DoPaintEntries(Bool_t)");
   fEntriesToDraw->Connect("Released()","TParallelCoordEditor",
                          this, "DoEntriesToDraw()");
   fEntriesToDraw->Connect("PositionChanged()","TParallelCoordEditor",
                          this, "DoLiveEntriesToDraw()");
   fFirstEntry->Connect("ReturnPressed()","TParallelCoordEditor",
                        this, "DoFirstEntry()");
   fNentries->Connect("ReturnPressed()","TParallelCoordEditor",
                     this, "DoNentries()");
   fApplySelect->Connect("Clicked()","TParallelCoordEditor",
                         this, "DoApplySelect()");
   fUnApply->Connect("Clicked()","TParallelCoordEditor",
                     this, "DoUnApply()");
   fDelayDrawing->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                     this, "DoDelayDrawing(Bool_t)");
   fAddVariable->Connect("ReturnPressed()","TParallelCoordEditor",
                     this, "DoAddVariable()");
   fButtonAddVar->Connect("Clicked()","TParallelCoordEditor",
                     this, "DoAddVariable()");
   fHideAllRanges->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                           this, "DoHideAllRanges(Bool_t)");
   fVariables->Connect("Selected(const char*)","TParallelCoordEditor",
                      this, "DoVariableSelect(const char*)");
   fDeleteVar->Connect("Clicked()","TParallelCoordEditor",
                      this, "DoDeleteVar()");
   fHistWidth->Connect("ReturnPressed()","TParallelCoordEditor",
                       this, "DoHistWidth()");
   fHistBinning->Connect("ReturnPressed()","TParallelCoordEditor",
                         this, "DoHistBinning()");
   fWeightCut->Connect("Released()","TParallelCoordEditor",
                       this, "DoWeightCut()");
   fWeightCut->Connect("PositionChanged(Int_t)","TParallelCoordEditor",
                       this, "DoLiveWeightCut(Int_t)");
   fWeightCutField->Connect("ReturnPressed()","TParallelCoordEditor",
                            this, "DoWeightCut()");
   fHistColorSelect->Connect("ColorSelected(Pixel_t)", "TParallelCoordEditor",
                             this, "DoHistColorSelect(Pixel_t)");
   fHistPatternSelect->Connect("PatternSelected(Style_t)", "TParallelCoordEditor",
                               this, "DoHistPatternSelect(Style_t)");
   fHistShowBoxes->Connect("Toggled(Bool_t)","TParallelCoordEditor",
                               this, "DoHistShowBoxes(Bool_t)");

   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to activate or not a selection.

void TParallelCoordEditor::DoActivateSelection(Bool_t on)
{
   if (fAvoidSignal) return;

   TParallelCoordSelect* sel = fParallel->GetCurrentSelection();
   if (sel) {
      sel->SetActivated(on);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to add a selection.

void TParallelCoordEditor::DoAddSelection()
{
   TString title = fAddSelectionField->GetText();
   if (title == "") title = "Selection";
   TString titlebis = title;
   Bool_t found = kTRUE;
   Int_t i=1;
   while (found){
      if (fSelectionSelect->FindEntry(titlebis)) {
         titlebis = title;
         titlebis.Append(Form("(%d)",i));
      }
      else found = kFALSE;
      ++i;
   }

   fParallel->AddSelection(titlebis.Data());

   CleanUpSelections();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to add a variable.

void TParallelCoordEditor::DoAddVariable()
{
   if (fAvoidSignal) return;

   fParallel->AddVariable(fAddVariable->GetText());
   CleanUpVariables();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to apply a selection to the tree.

void TParallelCoordEditor::DoApplySelect()
{
   //FIXME I forgot to update the slider over the entries
   //      (nentries and firstentry might have changed after applying the selection)

   if (fAvoidSignal) return;

   fParallel->ApplySelectionToTree();
   Update();
   SetModel(fParallel);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to delay the drawing.

void TParallelCoordEditor::DoDelayDrawing(Bool_t on)
{
   if (fAvoidSignal) return;

   fDelay = on;
   fParallel->SetLiveRangesUpdate(!on);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to delete a selection.

void TParallelCoordEditor::DoDeleteSelection()
{
   if (fAvoidSignal) return;

   fParallel->DeleteSelection(fParallel->GetCurrentSelection());

   CleanUpSelections();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to delete a variable().

void TParallelCoordEditor::DoDeleteVar()
{
   if (fAvoidSignal) return;

   Bool_t hasDeleted = fParallel->RemoveVariable(((TGTextLBEntry*)fVariables->GetSelectedEntry())->GetTitle());
   CleanUpVariables();
   if (hasDeleted) Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the line dotspacing.

void TParallelCoordEditor::DoDotsSpacing()
{
   if (fAvoidSignal) return;

   fParallel->SetDotsSpacing(fDotsSpacing->GetPosition());
   fDotsSpacingField->SetNumber((Int_t)fDotsSpacing->GetPosition());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the line dotspacing from the entry field.

void TParallelCoordEditor::DoDotsSpacingField()
{
   if (fAvoidSignal) return;

   fParallel->SetDotsSpacing((Int_t)fDotsSpacingField->GetNumber());
   fDotsSpacing->SetPosition((Int_t)fDotsSpacingField->GetNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value from the entry field.

void TParallelCoordEditor::DoAlphaField()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fParallel->GetLineColor())) {
      color->SetAlpha((Float_t)fAlphaField->GetNumber());
      fAlpha->SetPosition((Int_t)fAlphaField->GetNumber()*1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value

void TParallelCoordEditor::DoAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fParallel->GetLineColor())) {
      color->SetAlpha((Float_t)fAlpha->GetPosition()/1000);
      fAlphaField->SetNumber((Float_t)fAlpha->GetPosition()/1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to select the entries to be drawn.

void TParallelCoordEditor::DoEntriesToDraw()
{
   if (fAvoidSignal) return;

   Long64_t nentries,firstentry;
   firstentry = fEntriesToDraw->GetMinPositionL();
   nentries = (Long64_t)(fEntriesToDraw->GetMaxPositionD() - fEntriesToDraw->GetMinPositionD() + 1);

   fParallel->SetCurrentFirst(firstentry);
   fParallel->SetCurrentN(nentries);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the first entry.

void TParallelCoordEditor::DoFirstEntry()
{
   if (fAvoidSignal) return;

   fParallel->SetCurrentFirst((Long64_t)fFirstEntry->GetNumber());
   fEntriesToDraw->SetPosition((Long64_t)fFirstEntry->GetNumber(),(Long64_t)fFirstEntry->GetNumber()+fParallel->GetCurrentN());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the global line color.

void TParallelCoordEditor::DoGlobalLineColor(Pixel_t a)
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fParallel->GetLineColor())) {
      color->SetAlpha(1);
      color = gROOT->GetColor(TColor::GetColor(a));
      if (color) {
         color->SetAlpha((Float_t)fAlphaField->GetNumber());
         fParallel->SetLineColor(color->GetNumber());
      }
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the global line width.

void TParallelCoordEditor::DoGlobalLineWidth(Int_t wid)
{
   if (fAvoidSignal) return;

   fParallel->SetLineWidth(wid);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to hide all the ranges.

void TParallelCoordEditor::DoHideAllRanges(Bool_t on)
{
   if (fAvoidSignal) return;

   TIter next(fParallel->GetSelectList());
   TParallelCoordSelect* sel;
   while((sel = (TParallelCoordSelect*)next())) sel->SetShowRanges(!on);
   fShowRanges->SetOn(!on);
   fShowRanges->SetEnabled(!on);
   fShowRanges->SetOn(!on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the axes histogram binning.

void TParallelCoordEditor::DoHistBinning()
{
   if (fAvoidSignal) return;

   fParallel->SetAxisHistogramBinning((Int_t)fHistBinning->GetNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the histograms color.

void TParallelCoordEditor::DoHistColorSelect(Pixel_t p)
{
   if (fAvoidSignal) return;

   Color_t col = TColor::GetColor(p);
   TIter next(fParallel->GetVarList());
   TParallelCoordVar *var = NULL;
   while ((var = (TParallelCoordVar*)next())) var->SetFillColor(col);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set histogram height.

void TParallelCoordEditor::DoHistShowBoxes(Bool_t s)
{
   if (fAvoidSignal) return;

   TIter next(fParallel->GetVarList());
   TParallelCoordVar* var;
   while ((var = (TParallelCoordVar*)next())) var->SetBit(TParallelCoordVar::kShowBarHisto,s);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the histograms fill style.

void TParallelCoordEditor::DoHistPatternSelect(Style_t sty)
{
   if (fAvoidSignal) return;

   TIter next(fParallel->GetVarList());
   TParallelCoordVar *var = NULL;
   while ((var = (TParallelCoordVar*)next())) var->SetFillStyle(sty);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set histogram width.

void TParallelCoordEditor::DoHistWidth()
{
   if (fAvoidSignal) return;

   fParallel->SetAxisHistogramLineWidth((Int_t)fHistWidth->GetNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the line type.

void TParallelCoordEditor::DoLineType()
{
   if (fAvoidSignal) return;

   if (fLineTypePoly->GetState() == kButtonDown) fParallel->SetCurveDisplay(kFALSE);
   else fParallel->SetCurveDisplay(kTRUE);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the dots spacing online.

void TParallelCoordEditor::DoLiveDotsSpacing(Int_t a)
{
   if (fAvoidSignal) return;
   fDotsSpacingField->SetNumber(a);
   fParallel->SetDotsSpacing(a);
   if (!fDelay) Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set alpha value online.

void TParallelCoordEditor::DoLiveAlpha(Int_t a)
{
   if (fAvoidSignal) return;
   fAlphaField->SetNumber((Float_t)a/1000);

   if (TColor *color = gROOT->GetColor(fParallel->GetLineColor())) color->SetAlpha((Float_t)a/1000);
   if (!fDelay) Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update the entries fields from the slider position.

void TParallelCoordEditor::DoLiveEntriesToDraw()
{
   if (fAvoidSignal) return;

   Long64_t nentries,firstentry;
   firstentry = fEntriesToDraw->GetMinPositionL();
   nentries = (Long64_t)(fEntriesToDraw->GetMaxPositionD() - fEntriesToDraw->GetMinPositionD() + 1);

   fFirstEntry->SetNumber(firstentry);
   fNentries->SetNumber(nentries);

   if (!fDelay) {
      fParallel->SetCurrentFirst(firstentry);
      fParallel->SetCurrentN(nentries);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update the wieght cut entry field from the slider position.

void TParallelCoordEditor::DoLiveWeightCut(Int_t n)
{
   if (fAvoidSignal) return;

   fWeightCutField->SetNumber(n);
   if (!fDelay) {
      fParallel->SetWeightCut(n);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the number of entries to display.

void TParallelCoordEditor::DoNentries()
{
   if (fAvoidSignal) return;

   fParallel->SetCurrentN((Long64_t)fNentries->GetNumber());
   fEntriesToDraw->SetPosition(fParallel->GetCurrentFirst(),fParallel->GetCurrentFirst()+fParallel->GetCurrentN());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to postpone the entries drawing.

void TParallelCoordEditor::DoPaintEntries(Bool_t on)
{
   if (fAvoidSignal) return;

   fParallel->SetBit(TParallelCoord::kPaintEntries,on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the line color of selection.

void TParallelCoordEditor::DoSelectLineColor(Pixel_t a)
{
   if (fAvoidSignal) return;

   TParallelCoordSelect* sel = fParallel->GetCurrentSelection();
   if (sel) sel->SetLineColor(TColor::GetColor(a));
   fSelectionSelect->GetSelectedEntry()->SetBackgroundColor(a);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the line width of selection.

void TParallelCoordEditor::DoSelectLineWidth(Int_t wid)
{
   if (fAvoidSignal) return;

   TParallelCoordSelect* sel = fParallel->GetCurrentSelection();
   if (sel) {
      sel->SetLineWidth(wid);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the selection beeing edited.

void TParallelCoordEditor::DoSelectionSelect(const char* title)
{
   if (fAvoidSignal) return;

   if (!fParallel->SetCurrentSelection(title)) return;

   Color_t c = fParallel->GetCurrentSelection()->GetLineColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fSelectLineColor->SetColor(p,kFALSE);

   fSelectLineWidth->Select(fParallel->GetCurrentSelection()->GetLineWidth(),kFALSE);

   fActivateSelection->SetOn(fParallel->GetCurrentSelection()->TestBit(TParallelCoordSelect::kActivated));
   fShowRanges->SetOn(fParallel->GetCurrentSelection()->TestBit(TParallelCoordSelect::kShowRanges));
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to show or not the ranges on the pad.

void TParallelCoordEditor::DoShowRanges(Bool_t s)
{
   if (fAvoidSignal) return;

   TParallelCoordSelect *select = fParallel->GetCurrentSelection();
   if (select) {
      select->SetShowRanges(s);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to reset the tree entry list to the original one.

void TParallelCoordEditor::DoUnApply()
{
   if (fAvoidSignal) return;

   fParallel->ResetTree();
   Update();
   SetModel(fParallel);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to select a variable.

void TParallelCoordEditor::DoVariableSelect(const char* /*var*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update the weight cut.

void TParallelCoordEditor::DoWeightCut()
{
   if (fAvoidSignal) return;

   Int_t n = (Int_t)fWeightCutField->GetNumber();
   fParallel->SetWeightCut(n);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used parallel coordinates plot attributes.

void TParallelCoordEditor::SetModel(TObject* obj)
{
   if (!obj) return;
   fParallel = dynamic_cast<TParallelCoord*>(obj);
   if (!fParallel) return;
   fAvoidSignal = kTRUE;

   Color_t c = fParallel->GetLineColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fGlobalLineColor->SetColor(p);

   fGlobalLineWidth->Select(fParallel->GetLineWidth());

   fPaintEntries->SetOn(fParallel->TestBit(TParallelCoord::kPaintEntries));

   if (!TCanvas::SupportAlpha()) {
      fDotsSpacing->SetPosition(fParallel->GetDotsSpacing());
      fDotsSpacingField->SetNumber(fParallel->GetDotsSpacing());
   }
   else {
      if (TColor *color = gROOT->GetColor(fParallel->GetLineColor())) {
         fAlpha->SetPosition((Int_t)color->GetAlpha()*1000);
         fAlphaField->SetNumber(color->GetAlpha());
      }
   }

   Bool_t cur = fParallel->GetCurveDisplay();
   if (cur) fLineTypeBgroup->SetButton(kLineTypeCurves,kTRUE);
   else     fLineTypeBgroup->SetButton(kLineTypePoly,kTRUE);

   if (fInit) fHideAllRanges->SetOn(kFALSE);

   CleanUpSelections();
   CleanUpVariables();

   if (fInit) fEntriesToDraw->SetRange(0LL,fParallel->GetNentries());
   fEntriesToDraw->SetPosition(fParallel->GetCurrentFirst(), fParallel->GetCurrentFirst()+fParallel->GetCurrentN());

   fFirstEntry->SetNumber(fParallel->GetCurrentFirst());
   fNentries->SetNumber(fParallel->GetCurrentN());

   fDelayDrawing->SetOn(fDelay);

   fWeightCut->SetRange(0,(Int_t)(fParallel->GetNentries()/10)); // Maybe search here for better boundaries.
   fWeightCut->SetPosition(fParallel->GetWeightCut());
   fWeightCutField->SetNumber(fParallel->GetWeightCut());

   fHistColorSelect->SetColor(TColor::Number2Pixel(((TParallelCoordVar*)fParallel->GetVarList()->Last())->GetFillColor()), kFALSE);
   fHistPatternSelect->SetPattern(((TParallelCoordVar*)fParallel->GetVarList()->Last())->GetFillStyle(),kFALSE);

   if (fInit) ConnectSignals2Slots();

   fAvoidSignal = kFALSE;
}
