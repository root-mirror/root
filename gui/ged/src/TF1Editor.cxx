// @(#)root/ged:$Id$
// Author: Ilka Antcheva 21/03/06

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TF1Editor                                                           //
//                                                                      //
//  GUI for TF1 attributes and parameters.                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TF1Editor.h"
#include "TGedEditor.h"
#include "TH1.h"
#include "TF1.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"
#include "TString.h"
#include "TGNumberEntry.h"
#include "TG3DLine.h"
#include "TFunctionParametersDialog.h"
#include "TVirtualPad.h"


ClassImp(TF1Editor);

enum ETF1Wid {
   kTF1_TIT,  kTF1_NPX,
   kTF1_XSLD, kTF1_XMIN, kTF1_XMAX,
   kTF1_PAR,  kTF1_DRW
};

TF1Editor::TF1Editor(const TGWindow *p, Int_t width, Int_t height,
                     UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   // Constructor of TF1 editor.

   MakeTitle("Function");

   fTitle = new TGTextEntry(this, new TGTextBuffer(50), kTF1_TIT);
   fTitle->Resize(137, fTitle->GetDefaultHeight());
   fTitle->SetEnabled(kFALSE);
   fTitle->SetToolTipText(Form("Function expression or predefined name"));
   AddFrame(fTitle, new TGLayoutHints(kLHintsLeft,3, 2, 2, 3));

   TGCompositeFrame *f3a = new TGCompositeFrame(this, 137, 20, kHorizontalFrame);
   AddFrame(f3a, new TGLayoutHints(kLHintsTop, 0, 1, 3, 0));
   fDrawMode = new TGCheckButton(f3a, "Update", kTF1_DRW);
   fDrawMode->SetToolTipText("Immediate function redrawing");
   f3a->AddFrame(fDrawMode, new TGLayoutHints(kLHintsLeft | kLHintsBottom, 3, 1, 1, 0));
   fParLabel = new TGLabel(f3a, "");
   f3a->AddFrame(fParLabel, new TGLayoutHints(kLHintsRight | kLHintsBottom, 25, 2, 1, 0));

   TGCompositeFrame *f3 = new TGCompositeFrame(this, 137, 20, kHorizontalFrame | kFixedWidth);
   fSetPars = new TGTextButton(f3, "Set Parameters...", kTF1_PAR);
   f3->AddFrame(fSetPars, new TGLayoutHints(kLHintsRight | kLHintsTop | kLHintsExpandX,
                                            0, 1, 5, 0));
   fSetPars->SetToolTipText("Open a dialog for parameter(s) settings");
   AddFrame(f3, new TGLayoutHints(kLHintsTop | kLHintsLeft, 3, 2, 2, 3));

   MakeTitle("X-Range");

   TGCompositeFrame *f4 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   TGCompositeFrame *f4a = new TGCompositeFrame(f4, 66, 20, kVerticalFrame | kFixedWidth);
   TGLabel *fNpxLabel = new TGLabel(f4a, "Points: ");
   f4a->AddFrame(fNpxLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 5, 1));
   f4->AddFrame(f4a, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 1, 0, 0));

   TGCompositeFrame *f4b = new TGCompositeFrame(f4, 40, 20, kVerticalFrame);
   fNXpoints = new TGNumberEntry(f4b, 100, 7, kTF1_NPX,
                                 TGNumberFormat::kNESInteger,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax,4,100000);
   fNXpoints->GetNumberEntry()->SetToolTipText("Points along x-axis (4-100 000)");
   f4b->AddFrame(fNXpoints, new TGLayoutHints(kLHintsLeft, 0, 0, 1, 0));
   f4->AddFrame(f4b, new TGLayoutHints(kLHintsTop | kLHintsRight, 0, 1, 0, 0));
   AddFrame(f4, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 1, 0, 0));

   TGCompositeFrame *f5 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fSliderX = new TGDoubleHSlider(f5, 1, 2);
   fSliderX->Resize(137,20);
   f5->AddFrame(fSliderX, new TGLayoutHints(kLHintsLeft));
   AddFrame(f5, new TGLayoutHints(kLHintsTop, 3, 7, 4, 1));

   TGCompositeFrame *f6 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fSldMinX = new TGNumberEntryField(f6, kTF1_XMIN, 0.0,
                                     TGNumberFormat::kNESRealFour,
                                     TGNumberFormat::kNEAAnyNumber);
   ((TGTextEntry*)fSldMinX)->SetToolTipText("Lower bound along x-axis");
   fSldMinX->Resize(65,20);
   fSldMinX->SetState(kFALSE);
   f6->AddFrame(fSldMinX, new TGLayoutHints(kLHintsLeft));
   fSldMaxX = new TGNumberEntryField(f6, kTF1_XMAX, 0.0,
                                     TGNumberFormat::kNESRealFour,
                                     TGNumberFormat::kNEAAnyNumber);
   ((TGTextEntry*)fSldMaxX)->SetToolTipText("Upper bound along x-axis");
   fSldMaxX->SetState(kFALSE);
   fSldMaxX->Resize(65,20);
   f6->AddFrame(fSldMaxX, new TGLayoutHints(kLHintsLeft, 4, 0, 0, 0));
   AddFrame(f6, new TGLayoutHints(kLHintsTop, 3, 3, 5, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of TF1 editor.

TF1Editor::~TF1Editor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TF1Editor::ConnectSignals2Slots()
{
   fNXpoints->Connect("ValueSet(Long_t)", "TF1Editor", this, "DoXPoints()");
   (fNXpoints->GetNumberEntry())->Connect("ReturnPressed()", "TF1Editor",
                                          this, "DoXPoints()");
   fSetPars->Connect("Clicked()", "TF1Editor", this, "DoParameterSettings()");
   fSliderX->Connect("Pressed()","TF1Editor", this,"DoSliderXPressed()");
   fSliderX->Connect("Released()","TF1Editor", this,"DoSliderXReleased()");
   fSliderX->Connect("PositionChanged()","TF1Editor", this,"DoSliderXMoved()");

   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the function parameters and options.

void TF1Editor::SetModel(TObject* obj)
{
   if (obj == 0 || !obj->InheritsFrom(TF1::Class())) {
      return;
   }

   fF1 = (TF1*)obj;
   fAvoidSignal = kTRUE;

   const char *text = fF1->GetTitle();
   fTitle->SetText(text);

   fNP = fF1->GetNpar();
   fParLabel->SetText(Form("Npar: %d", fNP));
   fClient->NeedRedraw(fParLabel);

   fNXpoints->SetNumber(fF1->GetNpx());

   if (!fNP)
      fSetPars->SetState(kButtonDisabled, kFALSE);
   else
      fSetPars->SetState(kButtonUp, kFALSE);

   TAxis *x = fF1->GetHistogram()->GetXaxis();
   Int_t nx = x->GetNbins();
   Int_t nxbinmin = x->GetFirst();
   Int_t nxbinmax = x->GetLast();
   fSliderX->SetRange(1,nx);
   fSliderX->SetPosition((Double_t)nxbinmin,(Double_t)nxbinmax);
   fSldMinX->SetNumber(x->GetBinLowEdge(nxbinmin));
   fSldMaxX->SetNumber(x->GetBinUpEdge(nxbinmax));

   if (fInit) ConnectSignals2Slots();
   fAvoidSignal = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the function parameter(s) settings.

void TF1Editor::DoParameterSettings()
{
   TGMainFrame *main =  (TGMainFrame *)GetMainFrame();
   Double_t rmin = fSldMinX->GetNumber();
   Double_t rmax = fSldMaxX->GetNumber();
   new TFunctionParametersDialog(gClient->GetDefaultRoot(), main,
                                 fF1, fGedEditor->GetPad(), rmin, rmax);

}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the number of points setting.

void TF1Editor::DoXPoints()
{
   if (fAvoidSignal) return;
   Double_t rmin, rmax;
   fF1->GetRange(rmin, rmax);
   fF1->SetRange(fSldMinX->GetNumber(), fSldMaxX->GetNumber());
   fF1->SetNpx((Int_t)fNXpoints->GetNumber());
   fF1->GetHistogram()->GetXaxis()->Set((Int_t)fNXpoints->GetNumber(),
                                          fSldMinX->GetNumber(),
                                          fSldMaxX->GetNumber());
   Update();
   fF1->SetRange(rmin, rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-Slider range for function redrawing.

void TF1Editor::DoSliderXMoved()
{
   if (fAvoidSignal) return;

   TVirtualPad *save = 0;
   save = gPad;
   gPad = fGedEditor->GetPad();
   fGedEditor->GetPad()->cd();

   fF1->SetNpx((Int_t)fNXpoints->GetNumber());
   TAxis *x = fF1->GetHistogram()->GetXaxis();

   if (fDrawMode->GetState() == kButtonDown) {
      TString opt = fF1->GetDrawOption();
      opt.ToUpper();
      if (!opt.Contains("SAME"))
         opt += "SAME";
      fF1->Draw(opt);

      x->SetRange((Int_t)((fSliderX->GetMinPosition())+0.5),
                  (Int_t)((fSliderX->GetMaxPosition())+0.5));
      fSldMinX->SetNumber(x->GetBinLowEdge(x->GetFirst()));
      fSldMaxX->SetNumber(x->GetBinUpEdge(x->GetLast()));
      fClient->NeedRedraw(fSliderX,kTRUE);
      fClient->NeedRedraw(fSldMinX,kTRUE);
      fClient->NeedRedraw(fSldMaxX,kTRUE);
      Update();

   } else {
      x->SetRange((Int_t)((fSliderX->GetMinPosition())+0.5),
                  (Int_t)((fSliderX->GetMaxPosition())+0.5));
      fSldMinX->SetNumber(x->GetBinLowEdge(x->GetFirst()));
      fSldMaxX->SetNumber(x->GetBinUpEdge(x->GetLast()));
      fClient->NeedRedraw(fSliderX,kTRUE);
      fClient->NeedRedraw(fSldMinX,kTRUE);
      fClient->NeedRedraw(fSldMaxX,kTRUE);
   }
   if(save) gPad = save;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-Slider.

void TF1Editor::DoSliderXPressed()
{
   if (fAvoidSignal || (fDrawMode->GetState() == kButtonDown)) return;

   TVirtualPad *save = 0;
   save = gPad;
   gPad = fGedEditor->GetPad();
   fGedEditor->GetPad()->cd();

   fF1->SetNpx((Int_t)fNXpoints->GetNumber());
   TAxis *x = fF1->GetHistogram()->GetXaxis();
   TString opt = fF1->GetDrawOption();
   opt.ToUpper();
   if (!opt.Contains("SAME"))
      opt += "SAME";
   fF1->Draw(opt);

   x->SetRange((Int_t)((fSliderX->GetMinPosition())+0.5),
               (Int_t)((fSliderX->GetMaxPosition())+0.5));
   fSldMinX->SetNumber(x->GetBinLowEdge(x->GetFirst()));
   fSldMaxX->SetNumber(x->GetBinUpEdge(x->GetLast()));
   fClient->NeedRedraw(fSliderX,kTRUE);
   fClient->NeedRedraw(fSldMinX,kTRUE);
   fClient->NeedRedraw(fSldMaxX,kTRUE);
   Update();

   if(save) gPad = save;

}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-Slider.

void TF1Editor::DoSliderXReleased()
{
   if (fAvoidSignal || (fDrawMode->GetState() == kButtonDown)) return;

   TVirtualPad *save = 0;
   save = gPad;
   gPad = fGedEditor->GetPad();
   fGedEditor->GetPad()->cd();

   fF1->SetNpx((Int_t)fNXpoints->GetNumber());
   TAxis *x = fF1->GetHistogram()->GetXaxis();
   TString opt = fF1->GetDrawOption();
   opt.ToUpper();
   if (!opt.Contains("SAME"))
      opt += "SAME";
   fF1->Draw(opt);

   x->SetRange((Int_t)((fSliderX->GetMinPosition())+0.5),
               (Int_t)((fSliderX->GetMaxPosition())+0.5));
   fSldMinX->SetNumber(x->GetBinLowEdge(x->GetFirst()));
   fSldMaxX->SetNumber(x->GetBinUpEdge(x->GetLast()));
   fClient->NeedRedraw(fSliderX,kTRUE);
   fClient->NeedRedraw(fSldMinX,kTRUE);
   fClient->NeedRedraw(fSldMaxX,kTRUE);
   Update();

   if(save) gPad = save;
}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to min/max settings of the slider range.

void TF1Editor::DoXRange()
{
   if (fAvoidSignal) return;
   TAxis *x = fF1->GetHistogram()->GetXaxis();
   Int_t nx = x->GetNbins();
   Double_t width = x->GetBinWidth(1);
   Double_t lowLimit = x->GetBinLowEdge(1);
   Double_t upLimit = x->GetBinUpEdge(nx);
   if ((fSldMinX->GetNumber()+width/2) < (lowLimit))
      fSldMinX->SetNumber(lowLimit);
   if ((fSldMaxX->GetNumber()-width/2) > (upLimit))
      fSldMaxX->SetNumber(upLimit);
   x->SetRangeUser(fSldMinX->GetNumber()+width/2,
                   fSldMaxX->GetNumber()-width/2);
   Int_t nxbinmin = x->GetFirst();
   Int_t nxbinmax = x->GetLast();
   fSliderX->SetPosition((Double_t)(nxbinmin),(Double_t)(nxbinmax));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Exclude TAttFillEditor from this interface.

void TF1Editor::ActivateBaseClassEditors(TClass* cl)
{
   fGedEditor->ExcludeClassEditor(TAttFill::Class());
   TGedFrame::ActivateBaseClassEditors(cl);
}

