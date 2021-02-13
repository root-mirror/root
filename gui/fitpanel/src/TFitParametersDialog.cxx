// @(#)root/fitpanel:$Id$
// Author: Ilka Antcheva, Lorenzo Moneta 03/10/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TFitParametersDialog                                                //
//                                                                      //
//  Create a dialog for fit function parameter settings.                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TFitParametersDialog.h"
#include "TF1.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGLayout.h"
#include "TGTextEntry.h"
#include "TGMsgBox.h"
#include "TGNumberEntry.h"
#include "TGTripleSlider.h"
#include "TVirtualPad.h"
#include "TMath.h"

#include <limits>


enum EParametersDialogWid {
   kNAME,
   kFIX = 10,
   kBND = 20,
   kVAL = 30,
   kMIN = 40,
   kMAX = 50,
   kSLD = 60,
   kSTP = 70,
   kERR = 80,
   kUPDATE = 8888,
   kRESET,
   kAPPLY,
   kOK,
   kCANCEL
};

ClassImp(TFitParametersDialog);

////////////////////////////////////////////////////////////////////////////////

TFitParametersDialog::TFitParametersDialog(const TGWindow *p,
                                           const TGWindow *main,
                                           TF1 *func,
                                           TVirtualPad *pad,
                                           Int_t *ret_code) :
   TGTransientFrame(p, main, 10, 10, kVerticalFrame),
   fFunc           (func),
   fFpad           (pad),
   fHasChanges     (kFALSE),
   fImmediateDraw  (kTRUE),
   fRetCode        (ret_code)

{
   // Create a dialog for fit function parameters' settings.

   SetCleanup(kDeepCleanup);

   fFunc->GetRange(fRangexmin, fRangexmax);
   fNP = fFunc->GetNpar();
   fPmin = new Double_t[fNP];
   fPmax = new Double_t[fNP];
   fPval = new Double_t[fNP];
   fPerr = new Double_t[fNP];
   fPstp = new Double_t[fNP];

   for (Int_t i = 0; i < fNP; i++) {
      fFunc->GetParLimits(i, fPmin[i], fPmax[i]);
      fPval[i] = fFunc->GetParameter(i);
      fPerr[i] = fFunc->GetParError(i);
      if (fPerr[i] > 1E-16)
         fPstp[i] = TMath::Power(10, TMath::Floor(TMath::Log10(fPerr[i])));
      else { 
         if (TMath::Abs(fPval[i]) > 1.)
            // if error is zero use as step approx 10% of current value
            fPstp[i] = TMath::Power(10, TMath::Floor(TMath::Log10(fPval[i])) - 1);
         else
            fPstp[i] = 0.1;
      }
   }
   fParNam = new TGTextEntry*[fNP];
   fParFix = new TGCheckButton*[fNP];
   fParBnd = new TGCheckButton*[fNP];
   fParVal = new TGNumberEntry*[fNP];
   fParMin = new TGNumberEntryField*[fNP];
   fParMax = new TGNumberEntryField*[fNP];
   fParSld = new TGTripleHSlider*[fNP];
   fParStp = new TGNumberEntry*[fNP];
   fParErr = new TGNumberEntryField*[fNP];

   memset(fParNam, 0, sizeof(TGTextEntry*)*fNP);
   memset(fParFix, 0, sizeof(TGCheckButton*)*fNP);
   memset(fParBnd, 0, sizeof(TGCheckButton*)*fNP);
   memset(fParVal, 0, sizeof(TGNumberEntry*)*fNP);
   memset(fParMin, 0, sizeof(TGNumberEntryField*)*fNP);
   memset(fParMax, 0, sizeof(TGNumberEntryField*)*fNP);
   memset(fParSld, 0, sizeof(TGTripleHSlider*)*fNP);
   memset(fParStp, 0, sizeof(TGNumberEntry*)*fNP);
   memset(fParErr, 0, sizeof(TGNumberEntryField*)*fNP);

   TGCompositeFrame *f1 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f1, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   // column 'Name'
   fContNam = new TGCompositeFrame(f1, 80, 20, kVerticalFrame | kFixedWidth);
   fContNam->AddFrame(new TGLabel(fContNam,"Name"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParNam[i] = new TGTextEntry(fContNam, new TGTextBuffer(80), kNAME+i);
      fParNam[i]->SetText(Form("%s", fFunc->GetParName(i)));
      fParNam[i]->SetEnabled(kFALSE);
      fContNam->AddFrame(fParNam[i],
                         new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
   }
   f1->AddFrame(fContNam, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Fix'
   fContFix = new TGCompositeFrame(f1, 20, 20, kVerticalFrame | kFixedWidth);
   fContFix->AddFrame(new TGLabel(fContFix,"Fix"),
                      new TGLayoutHints(kLHintsTop, 2, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParFix[i] = new TGCheckButton(fContFix, "", kFIX*fNP+i);
      fParFix[i]->SetToolTipText(Form("Set %s to fixed", fFunc->GetParName(i)));
      fContFix->AddFrame(fParFix[i], new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                                       5, 5, 10, 7));
      if ((fPmin[i] == fPmax[i]) && (fPmin[i] || fPmax[i]))
         fParFix[i]->SetState(kButtonDown);
      else
         fParFix[i]->SetState(kButtonUp);
      fParFix[i]->Connect("Toggled(Bool_t)", "TFitParametersDialog", this, "DoParFix(Bool_t)");
   }
   f1->AddFrame(fContFix, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Bound'
   fContBnd = new TGCompositeFrame(f1, 40, 20, kVerticalFrame | kFixedWidth);
   fContBnd->AddFrame(new TGLabel(fContBnd,"Bound"),
                      new TGLayoutHints(kLHintsTop, 2, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParBnd[i] = new TGCheckButton(fContBnd, "", kBND*fNP+i);
      fParBnd[i]->SetToolTipText(Form("Set bound to %s", fFunc->GetParName(i)));
      fContBnd->AddFrame(fParBnd[i], new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                                       15, 5, 10, 7));
      fParBnd[i]->Connect("Toggled(Bool_t)", "TFitParametersDialog", this, "DoParBound(Bool_t)");
      if ( ((fPmin[i] != fPmax[i]) && (fPmin[i] || fPmax[i])) || (fParMin[i] < fParMax[i]) )
         fParBnd[i]->SetState(kButtonDown, kFALSE);
      else
         fParBnd[i]->SetState(kButtonUp, kFALSE);
   }
   f1->AddFrame(fContBnd, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Value'
   fContVal = new TGCompositeFrame(f1, 100, 20, kVerticalFrame | kFixedWidth);
   fContVal->AddFrame(new TGLabel(fContVal,"Value"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParVal[i] = new TGNumberEntry(fContVal, 1.2E-12, 15, kVAL*fNP+i,
                                     TGNumberFormat::kNESReal);
      fParVal[i]->SetNumber(fPval[i]);
      fParVal[i]->SetFormat(TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber); //tbs
      fContVal->AddFrame(fParVal[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
      (fParVal[i]->GetNumberEntry())->SetToolTipText(Form("%s", fFunc->GetParName(i)));
      (fParVal[i]->GetNumberEntry())->Connect("ReturnPressed()", "TFitParametersDialog",
                                              this, "DoParValue()");
      fParVal[i]->Connect("ValueSet(Long_t)", "TFitParametersDialog", this, "DoParValue()");
      (fParVal[i]->GetNumberEntry())->Connect("TabPressed()", "TFitParametersDialog", this, "HandleTab()");
      (fParVal[i]->GetNumberEntry())->Connect("ShiftTabPressed()", "TFitParametersDialog", this, "HandleShiftTab()");
      fTextEntries.Add(fParVal[i]->GetNumberEntry());
   }
   f1->AddFrame(fContVal, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Min'
   fContMin = new TGCompositeFrame(f1, 100, 20, kVerticalFrame | kFixedWidth);
   fContMin->AddFrame(new TGLabel(fContMin,"Min"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParMin[i] = new TGNumberEntryField(fContMin, kMIN*fNP+i, 0.0,
                                          TGNumberFormat::kNESReal,
                                          TGNumberFormat::kNEAAnyNumber);
      ((TGTextEntry*)fParMin[i])->SetToolTipText(Form("Lower limit of %s",
                                                 fFunc->GetParName(i)));
      fContMin->AddFrame(fParMin[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
      fParMin[i]->SetNumber(fPmin[i]);
      fParMin[i]->Connect("ReturnPressed()", "TFitParametersDialog", this,
                          "DoParMinLimit()");
      fParMin[i]->Connect("TabPressed()", "TFitParametersDialog", this, "HandleTab()");
      fParMin[i]->Connect("ShiftTabPressed()", "TFitParametersDialog", this, "HandleShiftTab()");
      fTextEntries.Add(fParMin[i]);
   }
   f1->AddFrame(fContMin, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Set Range'
   fContSld = new TGCompositeFrame(f1, 120, 20, kVerticalFrame | kFixedWidth);
   fContSld->AddFrame(new TGLabel(fContSld,"Set Range"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParSld[i] = new TGTripleHSlider(fContSld, 100, kDoubleScaleBoth, kSLD*fNP+i,
                                       kHorizontalFrame, GetDefaultFrameBackground(),
                                       kFALSE, kFALSE, kFALSE, kFALSE);
      fContSld->AddFrame(fParSld[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 5));
      fParSld[i]->SetConstrained(kTRUE);
   }
   f1->AddFrame(fContSld,  new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Max'
   fContMax = new TGCompositeFrame(f1, 100, 20, kVerticalFrame | kFixedWidth);
   fContMax->AddFrame(new TGLabel(fContMax,"Max"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParMax[i] = new TGNumberEntryField(fContMax, kMAX*fNP+i, 0.0,
                                          TGNumberFormat::kNESReal,
                                          TGNumberFormat::kNEAAnyNumber);
      ((TGTextEntry*)fParMax[i])->SetToolTipText(Form("Upper limit of %s",
                                                 fFunc->GetParName(i)));
      fContMax->AddFrame(fParMax[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
      fParMax[i]->SetNumber(fPmax[i]);
      fParMax[i]->Connect("ReturnPressed()", "TFitParametersDialog", this, "DoParMaxLimit()");
      fParMax[i]->Connect("TabPressed()", "TFitParametersDialog", this, "HandleTab()");
      fParMax[i]->Connect("ShiftTabPressed()", "TFitParametersDialog", this, "HandleShiftTab()");
      fTextEntries.Add(fParMax[i]);
   }
   f1->AddFrame(fContMax, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Step'
   fContStp = new TGCompositeFrame(f1, 100, 20, kVerticalFrame | kFixedWidth);
   fContStp->AddFrame(new TGLabel(fContStp,"Step"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParStp[i] = new TGNumberEntry(fContStp, 1.2E-12, 15, kSTP*fNP+i,
                                     TGNumberFormat::kNESReal);
      fParStp[i]->SetNumber(fPstp[i]);
      fParStp[i]->SetFormat(TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber); //tbs
      fContStp->AddFrame(fParStp[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
      (fParStp[i]->GetNumberEntry())->SetToolTipText(Form("%s", fFunc->GetParName(i)));
      (fParStp[i]->GetNumberEntry())->Connect("ReturnPressed()", "TFitParametersDialog",
                                              this, "DoParStep()");
      fParStp[i]->Connect("ValueSet(Long_t)", "TFitParametersDialog", this, "DoParStep()");
      (fParStp[i]->GetNumberEntry())->Connect("TabPressed()", "TFitParametersDialog", this, "HandleTab()");
      (fParStp[i]->GetNumberEntry())->Connect("ShiftTabPressed()", "TFitParametersDialog", this, "HandleShiftTab()");
      fTextEntries.Add(fParStp[i]->GetNumberEntry());
   }
   f1->AddFrame(fContStp, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   // column 'Error'
   fContErr = new TGCompositeFrame(f1, 80, 20, kVerticalFrame | kFixedWidth);
   fContErr->AddFrame(new TGLabel(fContErr,"Errors"),
                      new TGLayoutHints(kLHintsTop, 5, 0, 0, 0));
   for (Int_t i = 0; i < fNP; i++ ) {
      fParErr[i] = new TGNumberEntryField(fContErr, kERR*fNP+i, 0.0,
                                          TGNumberFormat::kNESReal,
                                          TGNumberFormat::kNEAAnyNumber);
      ((TGTextEntry*)fParErr[i])->SetToolTipText(Form("Error of %s",
                                                 fFunc->GetParName(i)));
      fContErr->AddFrame(fParErr[i], new TGLayoutHints(kLHintsExpandX, 2, 2, 7, 5));
      fParErr[i]->SetEnabled(kFALSE);
      if (fPerr[i])
         fParErr[i]->SetNumber(fPerr[i]);
      else
         ((TGTextEntry *)fParErr[i])->SetText("-");
   }
   f1->AddFrame(fContErr, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

   TGCompositeFrame *f2 = new TGCompositeFrame(this, 270, 20, kHorizontalFrame);
   AddFrame(f2, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));

   fUpdate = new TGCheckButton(f2, "&Immediate preview", kUPDATE);
   fUpdate->SetToolTipText("Immediate function redrawing");
   fUpdate->SetState(kButtonDown);
   f2->AddFrame(fUpdate, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 5, 5));
   fUpdate->Connect("Toggled(Bool_t)", "TFitParametersDialog", this, "HandleButtons(Bool_t)");

   TGCompositeFrame *f3 = new TGCompositeFrame(f2, 270, 20, kHorizontalFrame | kFixedWidth);
   f2->AddFrame(f3, new TGLayoutHints(kLHintsRight));

   fReset = new TGTextButton(f3, "&Reset", kRESET);
   f3->AddFrame(fReset, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,5,5));
   fReset->SetToolTipText("Reset the parameter settings");
   fReset->SetState(kButtonDisabled);
   fReset->Connect("Clicked()", "TFitParametersDialog", this, "DoReset()");

   fApply = new TGTextButton(f3, "&Apply", kAPPLY);
   f3->AddFrame(fApply, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,5,5));
   fApply->SetState(kButtonDisabled);
   fApply->Connect("Clicked()", "TFitParametersDialog", this, "DoApply()");
   fApply->SetToolTipText("Apply parameter settings and redraw the function");

   fOK = new TGTextButton(f3, "&OK", kOK);
   f3->AddFrame(fOK, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,5,5));
   fOK->SetToolTipText("Apply parameter settings, redraw function and close this dialog");
   fOK->Connect("Clicked()", "TFitParametersDialog", this, "DoOK()");

   fCancel = new TGTextButton(f3, "&Cancel", kCANCEL);
   f3->AddFrame(fCancel, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2,2,5,5));
   fCancel->SetToolTipText("Close this dialog with no parameter changes");
   fCancel->Connect("Clicked()", "TFitParametersDialog", this, "DoCancel()");
   *fRetCode = kFPDNoChange; // default setting

   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
   CenterOnParent(kFALSE, kBottomLeft);
   SetWindowName(Form("Set Parameters of %s", fFunc->GetTitle()));

   for (Int_t i = 0; i < fNP; i++ ) {
      if (fParFix[i]->GetState() == kButtonDown) {
         fParVal[i]->SetState(kFALSE);
         fParMin[i]->SetEnabled(kFALSE);
         fParMax[i]->SetEnabled(kFALSE);
         fParSld[i]->UnmapWindow();
      } else {
         if (fPmin[i]*fPmax[i] == 0 && fPmin[i] >= fPmax[i]) { //init
            // round again the values on the percent level
            Double_t u = TMath::Power(10, TMath::Floor(TMath::Log10(TMath::Abs(fPval[i])) )-2 );
            Double_t roundVal = int(fPval[i]/ u) * u; 
            // set min at +/- 100 step size 
            fParMin[i]->SetNumber( roundVal - 100* fPstp[i]); 
            fParMax[i]->SetNumber( roundVal + 100* fPstp[i]); 
            // if (!fPval[i]) {
            //    fParMin[i]->SetNumber(-10);
            //    fParMax[i]->SetNumber(10);
            // } else {
            //    fParMin[i]->SetNumber(-3*TMath::Abs(fPval[i]));
            //    fParMax[i]->SetNumber(3*TMath::Abs(fPval[i]));
            // }
         }
         fParSld[i]->SetRange(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
         fParSld[i]->SetPosition(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
         fParSld[i]->SetPointerPosition(fParVal[i]->GetNumber());
         fParSld[i]->Connect("PointerPositionChanged()", "TFitParametersDialog",
                             this, "DoSlider()");
         fParSld[i]->Connect("PositionChanged()", "TFitParametersDialog",
                             this, "DoSlider()");
      }
   }

   gClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFitParametersDialog::~TFitParametersDialog()
{
   DisconnectSlots();
   fTextEntries.Clear();
   Cleanup();
   delete [] fPval;
   delete [] fPmin;
   delete [] fPmax;
   delete [] fPerr;
   delete [] fPstp;

   delete [] fParNam;
   delete [] fParFix;
   delete [] fParBnd;
   delete [] fParVal;
   delete [] fParMin;
   delete [] fParMax;
   delete [] fParSld;
   delete [] fParStp;
   delete [] fParErr;
}

////////////////////////////////////////////////////////////////////////////////
/// Close parameters' dialog.

void TFitParametersDialog::CloseWindow()
{
   if (fHasChanges) {
      Int_t ret;
      const char *txt;
      txt = "Do you want to apply last parameters' setting?";
      new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                   "Parameters Have Been Changed", txt, kMBIconExclamation,
                   kMBYes | kMBNo | kMBCancel, &ret);
      if (ret == kMBYes) {
         SetParameters();
         *fRetCode = kFPDChange;
      }
      else if (ret == kMBNo)
         DoReset();
      else return;
   }

   DisconnectSlots();
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the Cancel button.

void TFitParametersDialog::DoCancel()
{
   if (fHasChanges)
      DoReset();
   // for (Int_t i = 0; i < fNP; i++ ) {
   //    if (fParBnd[i]->GetState() == kButtonDown)
   //       *fRetCode = kFPDBounded;
   // }
   CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the Bound check button.

void TFitParametersDialog::DoParBound(Bool_t on)
{
   TGButton *bt = (TGButton *) gTQSender;
   Int_t id = bt->WidgetId();
   fHasChanges = kTRUE;

   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kBND*fNP+i) {
         if (on) {
            if (fParMin[i]->GetNumber() >= fParMax[i]->GetNumber()) {
               Int_t ret;
               const char *txt;
               txt = "'Min' value cannot be bigger or equal to 'Max' - set the limits first!";
               new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                            "Parameter Limits", txt, kMBIconExclamation,kMBOk,&ret);

               fParBnd[i]->SetState(kButtonUp, kFALSE);
               return;
            }
            if ((fParVal[i]->GetNumber() < fParMin[i]->GetNumber()) ||
                (fParVal[i]->GetNumber() > fParMax[i]->GetNumber())) {
               Double_t v = (fParMax[i]->GetNumber()+fParMin[i]->GetNumber())/2.;
               fParVal[i]->SetNumber(v);
               fFunc->SetParameter(i, v);
               fClient->NeedRedraw(fParVal[i]);
            }
            fParVal[i]->SetLimits(TGNumberFormat::kNELLimitMinMax,
                                  fParMin[i]->GetNumber(),
                                  fParMax[i]->GetNumber());
            fClient->NeedRedraw(fParVal[i]);
            fFunc->SetParLimits(i, fParMin[i]->GetNumber(),
                                   fParMax[i]->GetNumber());
         } else {
            fParVal[i]->SetLimits(TGNumberFormat::kNELNoLimits);
            fFunc->ReleaseParameter(i);
            fFunc->GetParLimits(i, fPmin[i], fPmax[i]);
            fPval[i] = fFunc->GetParameter(i);
            if (fPmin[i]*fPmax[i] == 0 && fPmin[i] >= fPmax[i]) { //init
               if (!fPval[i]) {
                  fParMin[i]->SetNumber(-10);
                  fParMax[i]->SetNumber(10);
               } else {
                  fParMin[i]->SetNumber(-10*TMath::Abs(fPval[i]));
                  fParMax[i]->SetNumber(10*TMath::Abs(fPval[i]));
               }
            }
            fParSld[i]->SetRange(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPosition(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPointerPosition(fPval[i]);
         }
      }
   }
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to parameter step setting.

void TFitParametersDialog::DoParStep()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the Fix check button.

void TFitParametersDialog::DoParFix(Bool_t on)
{
   fReset->SetState(kButtonUp);

   TGButton *bt = (TGButton *) gTQSender;
   Int_t id = bt->WidgetId();
   fHasChanges = kTRUE;

   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kFIX*fNP+i) {
         if (on) {
            // no bound available
            fParBnd[i]->Disconnect("Toggled(Bool_t)");
            fParBnd[i]->SetEnabled(kFALSE);
            fParBnd[i]->SetToolTipText(Form("DISABLED - %s is fixed", fFunc->GetParName(i)));
            if (fParVal[i]->GetNumber() != 0) {
               fParMin[i]->SetNumber(fParVal[i]->GetNumber());
               fParMin[i]->SetEnabled(kFALSE);
               fParMax[i]->SetNumber(fParVal[i]->GetNumber());
               fParMax[i]->SetEnabled(kFALSE);
            } else {
               fParMin[i]->SetNumber(1.);
               fParMin[i]->SetEnabled(kFALSE);
               fParMax[i]->SetNumber(1.);
               fParMax[i]->SetEnabled(kFALSE);
            }
            fParVal[i]->SetState(kFALSE);
            fParStp[i]->SetState(kFALSE);
            fParSld[i]->Disconnect("PointerPositionChanged()");
            fParSld[i]->Disconnect("PositionChanged()");
            fParSld[i]->UnmapWindow();
            fFunc->FixParameter(i, fParVal[i]->GetNumber());
         } else if (!fParMin[i]->IsEnabled()) {
            if (fPmin[i] != fPmax[i]) {
               if (fPmin[i])
                  fParMin[i]->SetNumber(fPmin[i]);
               else if (fPerr[i])
                  fParMin[i]->SetNumber(fPval[i]-3*fPerr[i]);
               else if (fPval[i])
                  fParMin[i]->SetNumber(fPval[i]-0.1*fPval[i]);
               else
                  fParMin[i]->SetNumber(1.0);
               if (fPmax[i])
                  fParMax[i]->SetNumber(fPmax[i]);
               else if (fPerr[i])
                  fParMax[i]->SetNumber(fPval[i]+3*fPerr[i]);
               else if (fPval[i])
                  fParMax[i]->SetNumber(fPval[i]+0.1*fPval[i]);
               else
                  fParMax[i]->SetNumber(1.0);
            } else if (fPval[i]) {
               fParMin[i]->SetNumber(fPval[i]-0.1*fPval[i]);
               fParMax[i]->SetNumber(fPval[i]+0.1*fPval[i]);
            } else {
               fParMin[i]->SetNumber(1.0);
               fParMax[i]->SetNumber(1.0);
            }
            if (fParMax[i]->GetNumber() < fParMin[i]->GetNumber()){
               Double_t temp;
               temp = fParMax[i]->GetNumber();
               fParMax[i]->SetNumber(fParMin[i]->GetNumber());
               fParMin[i]->SetNumber(temp);
            }
            fParBnd[i]->SetEnabled(kTRUE);
            fParBnd[i]->Connect("Toggled(Bool_t)",  "TFitParametersDialog",
                                this, "DoParBound(Bool_t)");
            fParBnd[i]->SetState(kButtonUp);
            fParMax[i]->SetEnabled(kTRUE);
            fParMin[i]->SetEnabled(kTRUE);
            fParSld[i]->MapWindow();
            fParVal[i]->SetState(kTRUE);
            fParStp[i]->SetState(kTRUE);
            fParSld[i]->SetRange(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPosition(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPointerPosition(fPval[i]);
            fParSld[i]->Connect("PointerPositionChanged()", "TFitParametersDialog",
                                this, "DoSlider()");
            fParSld[i]->Connect("PositionChanged()", "TFitParametersDialog",
                                this, "DoSlider()");
            fFunc->SetParLimits(i, fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
         }
      }
   }
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the parameter values inside the function

void TFitParametersDialog::SetParameters()
{
   fFunc->SetRange(fRangexmin, fRangexmax);
   for (Int_t i = 0; i < fNP; i++ ) {
      // first make sure the current value is up to date
      fParVal[i]->GetNumberEntry()->ReturnPressed();
      if (fParFix[i]->GetState() == kButtonDown) {
         fFunc->SetParameter(i, fParVal[i]->GetNumber());
         fFunc->FixParameter(i, fParVal[i]->GetNumber());
         *fRetCode = kFPDBounded;
      } else {
         if (fParBnd[i]->GetState() == kButtonDown) {
            fFunc->SetParameter(i, fParVal[i]->GetNumber());
            fFunc->SetParLimits(i, fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            *fRetCode = kFPDBounded;
         } else {
            fFunc->ReleaseParameter(i);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the OK button.

void TFitParametersDialog::DoOK()
{
   if (fHasChanges) 
      DrawFunction();

   SetParameters();

   // we want here to confirm the parameters settings so
   // it is like having changed them
   *fRetCode = kFPDChange;

   
   CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the Preview button.

void TFitParametersDialog::DoApply()
{
   DrawFunction();
   fApply->SetState(kButtonDisabled);
   if (fReset->GetState() == kButtonDisabled)
      fReset->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the Reset button.

void TFitParametersDialog::DoReset()
{
   fHasChanges = kTRUE;
   Int_t k = fNP;
   for (Int_t i = 0; i < fNP; i++) {
      if (fParVal[i]->GetNumber() == fPval[i])
         k--;
      else
         break;
   }

   if (!k) {
      if (fReset->GetState() == kButtonUp)
         fReset->SetState(kButtonDisabled);
      fHasChanges = kFALSE;
      return;
   }
   for (Int_t i = 0; i < fNP; i++) {
      fFunc->SetParameter(i, fPval[i]);
      fFunc->SetParLimits(i, fPmin[i], fPmax[i]);
      fFunc->SetParError(i, fPerr[i]);

      if (fPmin[i])
         fParMin[i]->SetNumber(fPmin[i]);
      else if (fPerr[i])
         fParMin[i]->SetNumber(fPval[i]-3*fPerr[i]);
      else if (fPval[i])
         fParMin[i]->SetNumber(-3*TMath::Abs(fPval[i]));
      else
         fParMin[i]->SetNumber(1.0);

      if (fPmax[i])
         fParMax[i]->SetNumber(fPmax[i]);
      else if (fPerr[i])
         fParMax[i]->SetNumber(fPval[i]+3*fPerr[i]);
      else if (fPval[i])
         fParMax[i]->SetNumber(3*TMath::Abs(fPval[i]));
      else
         fParMax[i]->SetNumber(1.0);
      if (fParMax[i]->GetNumber() < fParMin[i]->GetNumber()){
         Double_t temp;
         temp = fParMax[i]->GetNumber();
         fParMax[i]->SetNumber(fParMin[i]->GetNumber());
         fParMin[i]->SetNumber(temp);
      }
      if (fParMin[i]->GetNumber() == fParMax[i]->GetNumber()) {
         fParVal[i]->SetState(kFALSE);
         fParMin[i]->SetEnabled(kFALSE);
         fParMax[i]->SetEnabled(kFALSE);
         fParStp[i]->SetState(kFALSE);
         fParSld[i]->Disconnect("PointerPositionChanged()");
         fParSld[i]->Disconnect("PositionChanged()");
         fParSld[i]->UnmapWindow();
         fParBnd[i]->Disconnect("Toggled(Bool_t)");
         fParBnd[i]->SetEnabled(kFALSE);
         fFunc->FixParameter(i, fParVal[i]->GetNumber());
         fParFix[i]->SetState(kButtonDown);
      } else {
         fParFix[i]->SetState(kButtonUp);
         if (!fParMax[i]->IsEnabled()) {
            fParMax[i]->SetEnabled(kTRUE);
            fParMin[i]->SetEnabled(kTRUE);
            fParVal[i]->SetState(kTRUE);
            fParStp[i]->SetState(kTRUE);
            fParSld[i]->SetRange(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPosition(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
            fParSld[i]->SetPointerPosition(fPval[i]);
            fParSld[i]->MapWindow();
            fParSld[i]->Connect("PointerPositionChanged()", "TFitParametersDialog",
                                this, "DoSlider()");
            fParSld[i]->Connect("PositionChanged()", "TFitParametersDialog",
                                this, "DoSlider()");
            fParBnd[i]->SetEnabled(kTRUE);
            fParBnd[i]->Connect("Toggled(Bool_t)", "TFitParametersDialog",
                                this, "DoParBound(Bool_t)");
         }
      }
      fParVal[i]->SetNumber(fPval[i]);

      fParSld[i]->SetRange(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
      fParSld[i]->SetPosition(fParMin[i]->GetNumber(), fParMax[i]->GetNumber());
      fParSld[i]->SetPointerPosition(fPval[i]);
   }

   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
   fHasChanges = kFALSE;
   *fRetCode = kFPDNoChange;
   fReset->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the parameters' value settings.

void TFitParametersDialog::DoSlider()
{
   TGTripleHSlider *sl = (TGTripleHSlider *) gTQSender;
   Int_t id = sl->WidgetId();

   fHasChanges = kTRUE;
   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kSLD*fNP+i) {
         fFunc->SetParameter(i,fParSld[i]->GetPointerPosition());
         fFunc->SetParLimits(i,fParSld[i]->GetMinPosition(),
                               fParSld[i]->GetMaxPosition());
         fParMin[i]->SetNumber(fParSld[i]->GetMinPosition());
         fParMax[i]->SetNumber(fParSld[i]->GetMaxPosition());
         fParVal[i]->SetNumber(fParSld[i]->GetPointerPosition());
      }
   }
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
   if (fReset->GetState() == kButtonDisabled)
      fReset->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the parameter value settings.

void TFitParametersDialog::DoParValue()
{
   TGNumberEntry *ne = (TGNumberEntry *) gTQSender;
   Int_t id = ne->WidgetId();

   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kVAL*fNP+i)  {
         if (fParVal[i]->GetNumber() < fParMin[i]->GetNumber()) {
            Double_t extraIncrement = (fParMax[i]->GetNumber() - fParMin[i]->GetNumber()) / 4;
            fParMin[i]->SetNumber(fParVal[i]->GetNumber() - extraIncrement );
            fClient->NeedRedraw(fParMin[i]);
            fParSld[i]->SetRange(fParMin[i]->GetNumber(),
                                 fParMax[i]->GetNumber());
            fParSld[i]->SetPosition(fParMin[i]->GetNumber(),
                                    fParMax[i]->GetNumber());
         }
         if (fParVal[i]->GetNumber() > fParMax[i]->GetNumber()) {
            Double_t extraIncrement = (fParMax[i]->GetNumber() - fParMin[i]->GetNumber()) / 4;
            fParMax[i]->SetNumber(fParVal[i]->GetNumber()  + extraIncrement );
            fClient->NeedRedraw(fParMax[i]);
            fParSld[i]->SetRange(fParMin[i]->GetNumber(),
                                 fParMax[i]->GetNumber());
            fParSld[i]->SetPosition(fParMin[i]->GetNumber(),
                                    fParMax[i]->GetNumber());
         }
         fParSld[i]->SetPointerPosition(fParVal[i]->GetNumber());
         fClient->NeedRedraw(fParSld[i]);
         fFunc->SetParameter(i,fParSld[i]->GetPointerPosition());
         if (fParBnd[i]->GetState() == kButtonDown)
            fFunc->SetParLimits(i,fParSld[i]->GetMinPosition(),
                                  fParSld[i]->GetMaxPosition());
         else
            fFunc->ReleaseParameter(i);
      }
   }
   fHasChanges = kTRUE;
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
   if (fReset->GetState() == kButtonDisabled)
      fReset->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the minumum parameter limit settings.

void TFitParametersDialog::DoParMinLimit()
{
   TGNumberEntryField *ne = (TGNumberEntryField *) gTQSender;
   Int_t id = ne->WidgetId();

   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kMIN*fNP+i) {
         if ((fParMin[i]->GetNumber() >= fParMax[i]->GetNumber()) &&
             (fParBnd[i]->GetState() == kButtonDown)) {
            Int_t ret;
            const char *txt;
            txt = "'Min' cannot be bigger then 'Max' if this parameter is bounded.";
            new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                         "Parameter Limits", txt, kMBIconExclamation, kMBOk, &ret);
            fParMin[i]->SetNumber(fParVal[i]->GetNumber()-fParStp[i]->GetNumber());
            return;
         }
         if (fParBnd[i]->GetState() == kButtonDown) {
            Double_t val = (fParMax[i]->GetNumber()+fParMin[i]->GetNumber())/2.;
            fParVal[i]->SetNumber(val);
            fParVal[i]->SetLimitValues(fParMin[i]->GetNumber(),
                                       fParMax[i]->GetNumber());
         }
         fParSld[i]->SetRange(fParMin[i]->GetNumber(),
                              fParMax[i]->GetNumber());
         fParSld[i]->SetPosition(fParMin[i]->GetNumber(),
                                 fParMax[i]->GetNumber());
         fParSld[i]->SetPointerPosition(fParVal[i]->GetNumber());
         fClient->NeedRedraw(fParSld[i]);
      }
   }
   fHasChanges = kTRUE;
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
   if (fReset->GetState() == kButtonDisabled)
      fReset->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot related to the maximum parameter limit settings.

void TFitParametersDialog::DoParMaxLimit()
{
   TGNumberEntryField *ne = (TGNumberEntryField *) gTQSender;
   Int_t id = ne->WidgetId();

   for (Int_t i = 0; i < fNP; i++ ) {
      if (id == kMAX*fNP+i) {
         if ((fParMin[i]->GetNumber() >= fParMax[i]->GetNumber()) &&
             (fParBnd[i]->GetState() == kButtonDown)) {
            Int_t ret;
            const char *txt;
            txt = "'Min' cannot be bigger then 'Max' if this parameter is bounded.";
            new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                         "Parameter Limits", txt, kMBIconExclamation, kMBOk, &ret);
            fParMax[i]->SetNumber(fParVal[i]->GetNumber()+fParStp[i]->GetNumber());
            return;
         }
         if (fParBnd[i]->GetState() == kButtonDown) {
            Double_t val = (fParMax[i]->GetNumber()+(fParMin[i]->GetNumber()))/2.;
            fParVal[i]->SetNumber(val);
            fParVal[i]->SetLimitValues(fParMin[i]->GetNumber(),
                                       fParMax[i]->GetNumber());
         }
         fParSld[i]->SetRange(fParMin[i]->GetNumber(),
                              fParMax[i]->GetNumber());
         fParSld[i]->SetPosition(fParMin[i]->GetNumber(),
                                 fParMax[i]->GetNumber());
         fParSld[i]->SetPointerPosition(fParVal[i]->GetNumber());
         fClient->NeedRedraw(fParSld[i]);
      }
   }
   fHasChanges = kTRUE;
   if (fUpdate->GetState() == kButtonDown)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges)
      fApply->SetState(kButtonUp);
   if (fReset->GetState() == kButtonDisabled)
      fReset->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw function graphics.

void TFitParametersDialog::DrawFunction()
{
   if ( !fFpad ) return;
   TVirtualPad *save = 0;
   save = gPad;
   gPad = fFpad;
   gPad->cd();

   Style_t st = fFunc->GetLineStyle();
   fFunc->SetLineStyle(2);

   TString opt = fFunc->GetDrawOption();
   opt.ToUpper();
   if (!opt.Contains("SAME"))
      opt += "SAME";
   //fFunc->SetRange(fRXmin, fRXmax);
   fFunc->Draw(opt);
   gPad->Modified();
   gPad->Update();
   fHasChanges = kFALSE;

   fFunc->SetLineStyle(st);
   if (save) gPad = save;
   *fRetCode = kFPDBounded;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle the button dependent states in this dialog.

void TFitParametersDialog::HandleButtons(Bool_t update)
{
   if (update && fHasChanges)
      DrawFunction();
   else if ((fApply->GetState() == kButtonDisabled) && fHasChanges) {
      fApply->SetState(kButtonUp);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect signals from slot methods.

void TFitParametersDialog::DisconnectSlots()
{
   for (Int_t i = 0; i < fNP; i++ ) {
      fParFix[i]->Disconnect("Toggled(Bool_t)");
      fParBnd[i]->Disconnect("Toggled(Bool_t)");
      fParVal[i]->Disconnect("ValueSet(Long_t)");
      fParMin[i]->Disconnect("ReturnPressed()");
      fParMax[i]->Disconnect("ReturnPressed()");
      fParSld[i]->Disconnect("PointerPositionChanged()");
      fParSld[i]->Disconnect("PositionChanged()");
      fParStp[i]->Disconnect("ValueSet(Long_t)");
      fParVal[i]->Disconnect("TabPressed(Long_t)");
      fParVal[i]->Disconnect("ShiftTabPressed(Long_t)");
      fParMin[i]->Disconnect("TabPressed(Long_t)");
      fParMin[i]->Disconnect("ShiftTabPressed(Long_t)");
      fParMax[i]->Disconnect("TabPressed(Long_t)");
      fParMax[i]->Disconnect("ShiftTabPressed(Long_t)");
      fParStp[i]->Disconnect("TabPressed(Long_t)");
      fParStp[i]->Disconnect("ShiftTabPressed(Long_t)");
   }
   fUpdate->Disconnect("Toggled(Bool_t)");
   fReset->Disconnect("Clicked()");
   fApply->Disconnect("Clicked()");
   fOK->Disconnect("Clicked()");
   fCancel->Disconnect("Clicked()");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Shift+Tab key event (set focus to the previous number entry field)

void TFitParametersDialog::HandleShiftTab()
{
   TGNumberEntryField *next, *sender = (TGNumberEntryField *)gTQSender;
   next = (TGNumberEntryField *)fTextEntries.Before((TObject *)sender);
   if (next == 0)
      next = (TGNumberEntryField *)fTextEntries.Last();
   if (next) {
      next->SetFocus();
      next->Home();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Tab key event (set focus to the next number entry field)

void TFitParametersDialog::HandleTab()
{
   TGNumberEntryField *next, *sender = (TGNumberEntryField *)gTQSender;
   next = (TGNumberEntryField *)fTextEntries.After((TObject *)sender);
   if (next == 0)
      next = (TGNumberEntryField *)fTextEntries.First();
   if (next) {
      next->SetFocus();
      next->Home();
   }
}



