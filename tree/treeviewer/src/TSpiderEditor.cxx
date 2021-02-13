// @(#)root/treeviewer:$Id$
// Author: Bastien Dalla Piazza  20/07/07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TSpiderEditor.h"
#include "TSpider.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGButtonGroup.h"
#include "TGLabel.h"
#include "TGPicture.h"
#include "TGTextEntry.h"
#include "TG3DLine.h"
#include "TGComboBox.h"
#include "TGColorSelect.h"
#include "TGedPatternSelect.h"
#include "TColor.h"

ClassImp(TSpiderEditor);

/** \class TSpiderEditor
The TSpider editor class.
Provides the graphical user interface to the spider plots.
*/

enum ESpiderWid {
   kAverage,
   kNx,
   kNy,
   kPolyLines,
   kSegment,
   kGotoEntry,
   kPicNext,
   kPicPrevious,
   kPicFollowing,
   kPicPreceding,
   kAddVar,
   kDeleteVar,
   kAvLineStyle,
   kAvLineColor,
   kAvLineWidth,
   kAvFillColor,
   kAvFillStyle
};

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TSpiderEditor::TSpiderEditor(const TGWindow* /*p*/, Int_t /*width*/, Int_t /*height*/, UInt_t /*options*/, Pixel_t /*back*/)
{
   fPriority = 1;
   MakeTitle("Spider");

   fBgroup = new TGButtonGroup(this,2,1,0,0, "Plot type");
   fBgroup->SetRadioButtonExclusive(kTRUE);
   fPolyLines = new TGRadioButton(fBgroup, "PolyLine", kPolyLines);
   fPolyLines->SetToolTipText("Set a polyline plot type");
   fSegment = new TGRadioButton(fBgroup, "Segment", kSegment);
   fSegment->SetToolTipText("Set a segment plot type");
   fBgroup->ChangeOptions(kFitWidth|kChildFrame);
   AddFrame(fBgroup, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 4, 1, 0, 0));

   TGHorizontalFrame *f1 = new TGHorizontalFrame(this);

   TGLabel *nxLabel = new TGLabel(f1,"Nx:");
   f1->AddFrame(nxLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 1, 1));
   fSetNx = new TGNumberEntryField(f1, kNx, 2,
                                      TGNumberFormat::kNESInteger,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELLimitMinMax, 0, 99);
   fSetNx->SetToolTipText("Set the X number of plots");
   fSetNx->Resize(30,20);
   f1->AddFrame(fSetNx, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

   TGLabel *nyLabel = new TGLabel(f1,"Ny:");
   f1->AddFrame(nyLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 1, 1));
   fSetNy = new TGNumberEntryField(f1, kNy, 2,
                                      TGNumberFormat::kNESInteger,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELLimitMinMax, 0, 99);
   fSetNy->SetToolTipText("Set the Y number of plots");
   fSetNy->Resize(30,20);
   f1->AddFrame(fSetNy, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft,1,1,1,1));

   MakeTitle("Average");

   fDisplayAverage = new TGCheckButton(this, "Average",kAverage);
   fDisplayAverage->SetToolTipText("Display average");
   AddFrame(fDisplayAverage, new TGLayoutHints(kLHintsTop, 3, 1, 1, 1));

   TGHorizontalFrame *f2 = new TGHorizontalFrame(this);

   fAvLineColorSelect = new TGColorSelect(f2, 0, kAvLineColor);
   f2->AddFrame(fAvLineColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

   fAvLineWidthCombo = new TGLineWidthComboBox(f2, kAvLineWidth);
   fAvLineWidthCombo->Resize(91, 20);
   f2->AddFrame(fAvLineWidthCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));

   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   fAvLineStyleCombo = new TGLineStyleComboBox(this, kAvLineStyle);
   fAvLineStyleCombo->Resize(137, 20);
   AddFrame(fAvLineStyleCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));

   TGHorizontalFrame *f2b = new TGHorizontalFrame(this);

   fAvFillColorSelect = new TGColorSelect(f2b, 0, kAvFillColor);
   f2b->AddFrame(fAvFillColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

   fAvFillPatternSelect = new TGedPatternSelect(f2b, 1, kAvFillStyle);
   f2b->AddFrame(fAvFillPatternSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

   AddFrame(f2b, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   MakeBrowse();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of the TSpiderEditor.

TSpiderEditor::~TSpiderEditor()
{
   delete fPolyLines;
   delete fSegment;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TSpiderEditor::ConnectSignals2Slots()
{
   fDisplayAverage->Connect("Toggled(Bool_t)", "TSpiderEditor", this, "DoDisplayAverage(Bool_t)");
   fSetNx->Connect("ReturnPressed()", "TSpiderEditor", this, "DoSetNx()");
   fSetNy->Connect("ReturnPressed()", "TSpiderEditor", this, "DoSetNy()");
   fBgroup->Connect("Clicked(Int_t)","TSpiderEditor",this,"DoSetPlotType()");
   fGotoEntry->Connect("ReturnPressed()", "TSpiderEditor", this, "DoGotoEntry()");
   fGotoNext->Connect("Clicked()","TSpiderEditor",this,"DoGotoNext()");
   fGotoPrevious->Connect("Clicked()","TSpiderEditor",this,"DoGotoPrevious()");
   fGotoFollowing->Connect("Clicked()","TSpiderEditor",this,"DoGotoFollowing()");
   fGotoPreceding->Connect("Clicked()","TSpiderEditor",this,"DoGotoPreceding()");
   fAddVar->Connect("ReturnPressed()","TSpiderEditor",this,"DoAddVar()");
   fDeleteVar->Connect("ReturnPressed()","TSpiderEditor",this,"DoDeleteVar()");
   fAvLineStyleCombo->Connect("Selected(Int_t)", "TSpiderEditor", this, "DoAvLineStyle(Int_t)");
   fAvLineWidthCombo->Connect("Selected(Int_t)", "TSpiderEditor", this, "DoAvLineWidth(Int_t)");
   fAvLineColorSelect->Connect("ColorSelected(Pixel_t)", "TSpiderEditor", this, "DoAvLineColor(Pixel_t)");
   fAvFillColorSelect->Connect("ColorSelected(Pixel_t)", "TSpiderEditor", this, "DoAvFillColor(Pixel_t)");
   fAvFillPatternSelect->Connect("PatternSelected(Style_t)", "TSpiderEditor", this, "DoAvFillPattern(Style_t)");

   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make the Browse tab.

void TSpiderEditor::MakeBrowse()
{
   fBrowse = CreateEditorTabSubFrame("Browse");

   TGHorizontalFrame *title1 = new TGHorizontalFrame(fBrowse);
   title1->AddFrame(new TGLabel(title1, "Entries"),
                    new TGLayoutHints(kLHintsLeft, 3, 1, 0, 0));
   title1->AddFrame(new TGHorizontal3DLine(title1),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   fBrowse->AddFrame(title1, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));

   TGHorizontalFrame *f3 = new TGHorizontalFrame(fBrowse);

   TGLabel *gotoEntryLabel = new TGLabel(f3,"Go to:");
   f3->AddFrame(gotoEntryLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 18, 1, 5));
   fGotoEntry = new TGNumberEntryField(f3, kGotoEntry, 0,
                                      TGNumberFormat::kNESInteger,
                                      TGNumberFormat::kNEANonNegative);
   fGotoEntry->SetToolTipText("Jump to a specified entry");
   fGotoEntry->Resize(60,20);
   f3->AddFrame(fGotoEntry, new TGLayoutHints(kLHintsRight, 1, 1, 1, 1));

   fBrowse->AddFrame(f3,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   TGHorizontalFrame *f2 = new TGHorizontalFrame(fBrowse);

   fPicPrevious = gClient->GetPicture("first_t.xpm");
   fGotoPrevious = new TGPictureButton(f2,fPicPrevious,kPicPrevious);
   fGotoPrevious->SetToolTipText("Jump to the last entries");
   f2->AddFrame(fGotoPrevious,new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 1, 1, 1, 1));

   fPicPreceding = gClient->GetPicture("previous_t.xpm");
   fGotoPreceding = new TGPictureButton(f2,fPicPreceding,kPicPreceding);
   fGotoPreceding->SetToolTipText("Jump to the last entries");
   f2->AddFrame(fGotoPreceding,new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 1, 1, 1, 1));

   fPicFollowing = gClient->GetPicture("next_t.xpm");
   fGotoFollowing = new TGPictureButton(f2,fPicFollowing,kPicFollowing);
   fGotoFollowing->SetToolTipText("Jump to the last entries");
   f2->AddFrame(fGotoFollowing,new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 1, 1, 1, 1));

   fPicNext = gClient->GetPicture("last_t.xpm");
   fGotoNext = new TGPictureButton(f2,fPicNext,kPicNext);
   fGotoNext->SetToolTipText("Jump to the next entries");
   f2->AddFrame(fGotoNext,new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 1, 1, 1, 1));

   fBrowse->AddFrame(f2,new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));

   TGHorizontalFrame *title2 = new TGHorizontalFrame(fBrowse);

   title2->AddFrame(new TGLabel(title2, "Variables"),
                    new TGLayoutHints(kLHintsLeft, 3, 1, 0, 0));
   title2->AddFrame(new TGHorizontal3DLine(title2),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   fBrowse->AddFrame(title2, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));

   TGHorizontalFrame *f4 = new TGHorizontalFrame(fBrowse);

   TGVerticalFrame *v1 = new TGVerticalFrame(f4);

   TGLabel *addVar = new TGLabel(v1,"Add:");
   v1->AddFrame(addVar,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   TGLabel *deleteVar = new TGLabel(v1,"Delete:");
   v1->AddFrame(deleteVar,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   f4->AddFrame(v1, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,0,0));

   TGVerticalFrame *v2 = new TGVerticalFrame(f4);

   fAddVar = new TGTextEntry(v2, new TGTextBuffer(50), kAddVar);
   fAddVar->Resize(60,20);
   fAddVar->SetToolTipText("Add a variable");
   v2->AddFrame(fAddVar, new TGLayoutHints(kLHintsRight | kLHintsTop));

   fDeleteVar = new TGTextEntry(v2, new TGTextBuffer(50), kAddVar);
   fDeleteVar->Resize(60,20);
   fDeleteVar->SetToolTipText("Delete a variable");
   v2->AddFrame(fDeleteVar, new TGLayoutHints(kLHintsRight | kLHintsTop));

   f4->AddFrame(v2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,0,0));

   fBrowse->AddFrame(f4, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used spider attributes.

void TSpiderEditor::SetModel(TObject* obj)
{
   if (!obj) return;
   fSpider = dynamic_cast<TSpider*>(obj);
   if (!fSpider) return;
   fAvoidSignal = kTRUE;

   Bool_t av = fSpider->GetDisplayAverage();
   if(av) fDisplayAverage->SetState(kButtonDown);
   else fDisplayAverage->SetState(kButtonUp);

   fSetNx->SetNumber(fSpider->GetNx());
   fSetNy->SetNumber(fSpider->GetNy());

   Bool_t seg = fSpider->GetSegmentDisplay();

   if(seg) fBgroup->SetButton(kSegment,kTRUE);
   else fBgroup->SetButton(kPolyLines,kTRUE);

   fGotoEntry->SetNumber(fSpider->GetCurrentEntry());

   fAddVar->SetText("");
   fDeleteVar->SetText("");

   fAvLineStyleCombo->Select(fSpider->GetAverageLineStyle());
   fAvLineWidthCombo->Select(fSpider->GetAverageLineWidth());
   Color_t c = fSpider->GetAverageLineColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fAvLineColorSelect->SetColor(p);
   c = fSpider->GetAverageFillColor();
   p = TColor::Number2Pixel(c);
   fAvFillColorSelect->SetColor(p,kFALSE);
   fAvFillPatternSelect->SetPattern(fSpider->GetAverageFillStyle(),kFALSE);

   if (fInit) ConnectSignals2Slots();

   fAvoidSignal = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to add a variable.

void TSpiderEditor::DoAddVar()
{
   if (fAvoidSignal) return;

   const char * var = fAddVar->GetText();
   fSpider->AddVariable(var);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the average LineStyle.

void TSpiderEditor::DoAvLineStyle(Int_t a)
{
   if (fAvoidSignal) return;

   fSpider->SetAverageLineStyle(a);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the average LineWidth.

void TSpiderEditor::DoAvLineWidth(Int_t a)
{
   if (fAvoidSignal) return;

   fSpider->SetAverageLineWidth(a);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the average LineColor.

void TSpiderEditor::DoAvLineColor(Pixel_t a)
{
   if (fAvoidSignal) return;

   fSpider->SetAverageLineColor(TColor::GetColor(a));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the average Fill Color.

void TSpiderEditor::DoAvFillColor(Pixel_t a)
{
   if (fAvoidSignal) return;

   fSpider->SetAverageFillColor(TColor::GetColor(a));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the average FillStyle.

void TSpiderEditor::DoAvFillPattern(Style_t a)
{
   if (fAvoidSignal) return;

   fSpider->SetAverageFillStyle(a);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to delete a variable.

void TSpiderEditor::DoDeleteVar()
{
   if (fAvoidSignal) return;

   const char * var = fDeleteVar->GetText();
   fSpider->DeleteVariable(var);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot Connected to the average display.

void TSpiderEditor::DoDisplayAverage(Bool_t av)
{
   if (fAvoidSignal) return;

   fSpider->SetDisplayAverage(av);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to select an entry by number.

void  TSpiderEditor::DoGotoEntry()
{
   if (fAvoidSignal) return;
   Long64_t ev = (Long64_t)fGotoEntry->GetNumber();
   fSpider->GotoEntry(ev);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to Go to next entries.

void  TSpiderEditor::DoGotoNext()
{
   if (fAvoidSignal) return;
   fSpider->GotoNext();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to go to previous entries.

void  TSpiderEditor::DoGotoPrevious()
{
   if (fAvoidSignal) return;
   fSpider->GotoPrevious();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to go to next entry.

void  TSpiderEditor::DoGotoFollowing()
{
   if (fAvoidSignal) return;
   fSpider->GotoFollowing();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to go to last entry.

void  TSpiderEditor::DoGotoPreceding()
{
   if (fAvoidSignal) return;
   fSpider->GotoPreceding();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the nx setting.

void  TSpiderEditor::DoSetNx()
{
   if (fAvoidSignal) return;
   UInt_t nx = (UInt_t)fSetNx->GetNumber();
   fSpider->SetNx(nx);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the nx setting.

void  TSpiderEditor::DoSetNy()
{
   if (fAvoidSignal) return;
   UInt_t ny = (UInt_t)fSetNy->GetNumber();
   fSpider->SetNy(ny);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the plot type.

void  TSpiderEditor::DoSetPlotType()
{
   if(fSegment->GetState() == kButtonDown) fSpider->SetSegmentDisplay(kTRUE);
   else fSpider->SetSegmentDisplay(kFALSE);
   Update();
}
