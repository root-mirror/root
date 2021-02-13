// @(#)root/ged:$Id$
// Author: Ilka  Antcheva 24/04/06

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TLineEditor                                                        //
//                                                                      //
//  Implements GUI for editing line attributes: shape, size, angle.    //                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TLineEditor.gif">
*/
//End_Html


#include "TLineEditor.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TLine.h"

ClassImp(TLineEditor);

enum ELineWid {
   kLine_STAX,
   kLine_STAY,
   kLine_ENDX,
   kLine_ENDY,
   kLine_VERTICAL,
   kLine_HORIZONTAL
};


////////////////////////////////////////////////////////////////////////////////
/// Constructor of line GUI.

TLineEditor::TLineEditor(const TGWindow *p, Int_t width,
                           Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fLine = 0;

   MakeTitle("Points");

   TGCompositeFrame *f3 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f3, new TGLayoutHints(kLHintsTop, 1, 1, 3, 0));

   TGCompositeFrame *f3a = new TGCompositeFrame(f3, 80, 20);
   f3->AddFrame(f3a, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   TGLabel *fStartPointXLabel = new TGLabel(f3a, "Start X:");
   f3a->AddFrame(fStartPointXLabel, new TGLayoutHints(kLHintsNormal, 8, 0, 5, 5));

   TGLabel *fStartPointYLabel = new TGLabel(f3a, "Y:");
   f3a->AddFrame(fStartPointYLabel, new TGLayoutHints(kLHintsNormal, 37, 0, 5, 5));

   TGLabel *fEndPointXLabel = new TGLabel(f3a, "End X:");
   f3a->AddFrame(fEndPointXLabel, new TGLayoutHints(kLHintsNormal, 10, 0, 5, 5));

   TGLabel *fEndPointYLabel = new TGLabel(f3a, "Y:");
   f3a->AddFrame(fEndPointYLabel, new TGLayoutHints(kLHintsNormal, 37, 0, 5, 5));

   TGCompositeFrame *f3b = new TGCompositeFrame(f3, 80, 20, kFixedWidth);
   f3->AddFrame(f3b, new TGLayoutHints(kLHintsNormal, 8, 0, 0, 0));

   fStartPointX = new TGNumberEntry(f3b, 0.0, 8, kLine_STAX,
                                      TGNumberFormat::kNESRealThree,
                                      TGNumberFormat::kNEAAnyNumber,
                                      TGNumberFormat::kNELNoLimits);
   fStartPointX->GetNumberEntry()->SetToolTipText("Set start point X coordinate of Line.");
   f3b->AddFrame(fStartPointX, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

   fStartPointY = new TGNumberEntry(f3b, 0.0, 8, kLine_STAY,
                                      TGNumberFormat::kNESRealThree,
                                      TGNumberFormat::kNEAAnyNumber,
                                      TGNumberFormat::kNELNoLimits);
   fStartPointY->GetNumberEntry()->SetToolTipText("Set start point Y coordinate of Line.");
   f3b->AddFrame(fStartPointY, new TGLayoutHints(kLHintsExpandX, 1, 1, 3, 1));

   fEndPointX = new TGNumberEntry(f3b, 0.0, 8, kLine_ENDX,
                                    TGNumberFormat::kNESRealThree,
                                    TGNumberFormat::kNEAAnyNumber,
                                    TGNumberFormat::kNELNoLimits);
   fEndPointX->GetNumberEntry()->SetToolTipText("Set end point X xoordinate of Line.");
   f3b->AddFrame(fEndPointX, new TGLayoutHints(kLHintsExpandX, 1, 1, 3, 1));

   fEndPointY = new TGNumberEntry(f3b, 0.0, 8, kLine_ENDY,
                                    TGNumberFormat::kNESRealThree,
                                    TGNumberFormat::kNEAAnyNumber,
                                    TGNumberFormat::kNELNoLimits);
   fEndPointY->GetNumberEntry()->SetToolTipText("Set end point Y coordinate of Line.");
   f3b->AddFrame(fEndPointY, new TGLayoutHints(kLHintsExpandX, 1, 1, 3, 1));

   fVertical = new TGCheckButton(this,"Vertical",kLine_VERTICAL);
   fVertical->SetToolTipText("Set vertical");
   AddFrame(fVertical, new TGLayoutHints(kLHintsTop, 8, 1, 5, 0));

   fHorizontal = new TGCheckButton(this,"Horizontal",kLine_HORIZONTAL);
   fHorizontal->SetToolTipText("Set horizontal");
   AddFrame(fHorizontal, new TGLayoutHints(kLHintsTop, 8, 1, 3, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of line editor.

TLineEditor::~TLineEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TLineEditor::ConnectSignals2Slots()
{
   fStartPointX->Connect("ValueSet(Long_t)", "TLineEditor", this, "DoStartPoint()");
   (fStartPointX->GetNumberEntry())->Connect("ReturnPressed()", "TLineEditor", this, "DoStartPoint()");
   fStartPointY->Connect("ValueSet(Long_t)", "TLineEditor", this, "DoStartPoint()");
   (fStartPointY->GetNumberEntry())->Connect("ReturnPressed()", "TLineEditor", this, "DoStartPoint()");
   fEndPointX->Connect("ValueSet(Long_t)", "TLineEditor", this, "DoEndPoint()");
   (fEndPointX->GetNumberEntry())->Connect("ReturnPressed()", "TLineEditor", this, "DoEndPoint()");
   fEndPointY->Connect("ValueSet(Long_t)", "TLineEditor", this, "DoEndPoint()");
   (fEndPointY->GetNumberEntry())->Connect("ReturnPressed()", "TLineEditor", this, "DoEndPoint()");
   fVertical->Connect("Clicked()","TLineEditor",this,"DoLineVertical()");
   fHorizontal->Connect("Clicked()","TLineEditor",this,"DoLineHorizontal()");

   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used line attributes.

void TLineEditor::SetModel(TObject* obj)
{
   fLine = (TLine *)obj;
   fAvoidSignal = kTRUE;

   Float_t val = fLine->GetX1();
   fStartPointX->SetNumber(val);

   val = fLine->GetX2();
   fEndPointX->SetNumber(val);

   val = fLine->GetY1();
   fStartPointY->SetNumber(val);

   val = fLine->GetY2();
   fEndPointY->SetNumber(val);

   if (fLine->IsHorizontal()) fHorizontal->SetState(kButtonDown, kFALSE);
   else fHorizontal->SetState(kButtonUp, kFALSE);

   if (fLine->IsVertical()) fVertical->SetState(kButtonDown, kFALSE);
   else fVertical->SetState(kButtonUp, kFALSE);

   if (fInit) ConnectSignals2Slots();

   fAvoidSignal = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the line start point.

void TLineEditor::DoStartPoint()
{
   if (fAvoidSignal) return;
   fLine->SetX1((Double_t)fStartPointX->GetNumber());
   fLine->SetY1((Double_t)fStartPointY->GetNumber());
   fLine->Paint(fLine->GetDrawOption());
   Update();
}
////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the line EndPoint.

void TLineEditor::DoEndPoint()
{
   if (fAvoidSignal) return;
   fLine->SetX2((Double_t)fEndPointX->GetNumber());
   fLine->SetY2((Double_t)fEndPointY->GetNumber());
   fLine->Paint(fLine->GetDrawOption());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot so set the line vertical

void TLineEditor::DoLineVertical()
{
   if (fAvoidSignal) return;
   if (fVertical->GetState() == kButtonDown) {
      fLine->SetVertical();
      fHorizontal->SetState(kButtonUp, kFALSE);
   } else {
      fLine->SetVertical(kFALSE);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot so set the line horizontal

void TLineEditor::DoLineHorizontal()
{
   if (fAvoidSignal) return;
   if (fHorizontal->GetState() == kButtonDown) {
      fLine->SetHorizontal();
      fVertical->SetState(kButtonUp, kFALSE);
   } else {
      fLine->SetHorizontal(kFALSE);
   }
   Update();
}
