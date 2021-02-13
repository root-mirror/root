// @(#)root/sessionviewer:$Id$
// Author: Fons Rademakers   21/03/03

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProofProgressDialog                                                 //
//                                                                      //
// This class provides a query progress bar.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TProofProgressDialog.h"
#include "TProofProgressLog.h"
#include "TProofProgressMemoryPlot.h"
#include "TEnv.h"
#include "TError.h"
#include "TGLabel.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGProgressBar.h"
#include "TGSpeedo.h"
#include "TProof.h"
#include "TSlave.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TH1F.h"
#include "THLimitsFinder.h"
#include "TVirtualX.h"


#ifdef PPD_SRV_NEWER
#undef PPD_SRV_NEWER
#endif
#define PPD_SRV_NEWER(v) (fProof && fProof->GetRemoteProtocol() > v)

Bool_t TProofProgressDialog::fgKeepDefault = kTRUE;
Bool_t TProofProgressDialog::fgLogQueryDefault = kFALSE;
TString TProofProgressDialog::fgTextQueryDefault = "last";

//static const Int_t gSVNMemPlot = 25090;

ClassImp(TProofProgressDialog);

////////////////////////////////////////////////////////////////////////////////
/// Create PROOF processing progress dialog.

TProofProgressDialog::TProofProgressDialog(TProof *proof, const char *selector,
                                           Int_t files, Long64_t first,
                                           Long64_t entries) : fDialog(0),
   fBar(0), fClose(0), fStop(0), fAbort(0), fAsyn(0), fLog(0), fRatePlot(0),
   fMemPlot(0), fKeepToggle(0), fLogQueryToggle(0), fTextQuery(0), fEntry(0),
   fTitleLab(0), fFilesEvents(0), fTimeLab(0), fProcessed(0), fEstim(0),
   fTotal(0), fRate(0), fInit(0), fSelector(0), fSpeedo(0), fSmoothSpeedo(0)
{
   fProof         = proof;
   fFiles         = files;
   fFirst         = first;
   fEntries       = entries;
   fPrevProcessed = 0;
   fPrevTotal     = 0;
   fLogWindow     = 0;
   fMemWindow     = 0;
   fStatus        = kRunning;
   fKeep          = fgKeepDefault;
   fLogQuery      = fgLogQueryDefault;
   fRatePoints    = 0;
   fRateGraph     = 0;
   fMBRtGraph     = 0;
   fActWGraph     = 0;
   fTotSGraph     = 0;
   fEffSGraph     = 0;
   fProcTime      = 0.;
   fInitTime      = 0.;
   fAvgRate       = 0.;
   fAvgMBRate     = 0.;
   fRightInfo     = 0;
   fSpeedoEnabled = kFALSE;
   fSpeedo        = 0;
   fUpdtSpeedo    = 0;
   fSmoothSpeedo  = 0;

   // Make sure we are attached to a good instance
   if (!proof || !(proof->IsValid())) {
      Error("TProofProgressDialog", "proof instance is invalid (%p, %s): protocol error?",
                                    proof, (proof && !(proof->IsValid())) ? "invalid" : "undef");
      return;
   }

   // Have to save this information here, in case gProof is dead when
   // the logs are requested
   fSessionUrl = (proof && proof->GetManager()) ? proof->GetManager()->GetUrl() : "";

   if (PPD_SRV_NEWER(25)) {
      fRatePoints = new TNtuple("RateNtuple","Rate progress info","tm:evr:mbr:act:tos:efs");
   } else if (PPD_SRV_NEWER(11)) {
      fRatePoints = new TNtuple("RateNtuple","Rate progress info","tm:evr:mbr");
   }

   fDialog = new TGTransientFrame(0, 0, 10, 10);
   fDialog->Connect("CloseWindow()", "TProofProgressDialog", this, "DoClose()");
   fDialog->DontCallClose();
   fDialog->SetCleanup(kDeepCleanup);

//=======================================================================================

   TGHorizontalFrame *hf4 = new TGHorizontalFrame(fDialog, 100, 100);

   TGVerticalFrame *vf4 = new TGVerticalFrame(hf4, 100, 100);

   // Title label
   TString buf;
   buf.Form("Executing on PROOF cluster \"%s\" with %d parallel workers:",
            fProof ? fProof->GetMaster() : "<dummy>",
            fProof ? fProof->GetParallel() : 0);
   fTitleLab = new TGLabel(vf4, buf);
   fTitleLab->SetTextJustify(kTextTop | kTextLeft);
   vf4->AddFrame(fTitleLab, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));
   buf.Form("Selector: %s", selector);
   fSelector = new TGLabel(vf4, buf);
   fSelector->SetTextJustify(kTextTop | kTextLeft);
   vf4->AddFrame(fSelector, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));
   buf.Form("%d files, number of events %lld, starting event %lld",
            fFiles, fEntries, fFirst);
   fFilesEvents = new TGLabel(vf4, buf);
   fFilesEvents->SetTextJustify(kTextTop | kTextLeft);
   vf4->AddFrame(fFilesEvents, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));

   // Progress bar
   fBar = new TGHProgressBar(vf4, TGProgressBar::kFancy, 200);
   fBar->SetBarColor("green");
   fBar->Percent(kTRUE);
   fBar->ShowPos(kTRUE);
   vf4->AddFrame(fBar, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                     kLHintsExpandX, 10, 10, 5, 5));

   // Status labels
   if (PPD_SRV_NEWER(11)) {
      TGHorizontalFrame *hf0 = new TGHorizontalFrame(vf4, 0, 0);
      TGCompositeFrame *cf0 = new TGCompositeFrame(hf0, 110, 0, kFixedWidth);
      cf0->AddFrame(new TGLabel(cf0, "Initialization time:"));
      hf0->AddFrame(cf0);
      fInit = new TGLabel(hf0, "- secs");
      fInit->SetTextJustify(kTextTop | kTextLeft);
      hf0->AddFrame(fInit, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 0, 0));
      vf4->AddFrame(hf0, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));
   }

   TGHorizontalFrame *hf1 = new TGHorizontalFrame(vf4, 0, 0);
   TGCompositeFrame *cf1 = new TGCompositeFrame(hf1, 110, 0, kFixedWidth);
   fTimeLab = new TGLabel(cf1, "Estimated time left:");
   fTimeLab->SetTextJustify(kTextTop | kTextLeft);
   cf1->AddFrame(fTimeLab, new TGLayoutHints(kLHintsLeft));
   hf1->AddFrame(cf1);
   fEstim = new TGLabel(hf1, "- sec");
   fEstim->SetTextJustify(kTextTop | kTextLeft);
   hf1->AddFrame(fEstim, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 0, 0));
   vf4->AddFrame(hf1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));

   hf1 = new TGHorizontalFrame(vf4, 0, 0);
   cf1 = new TGCompositeFrame(hf1, 110, 0, kFixedWidth);
   fProcessed = new TGLabel(cf1, "Processing status:");
   fProcessed->SetTextJustify(kTextTop | kTextLeft);
   cf1->AddFrame(fProcessed, new TGLayoutHints(kLHintsLeft));
   hf1->AddFrame(cf1);
   fTotal= new TGLabel(hf1, "- / - events");
   fTotal->SetTextJustify(kTextTop | kTextLeft);
   hf1->AddFrame(fTotal, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 0, 0));

   vf4->AddFrame(hf1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));

   TGHorizontalFrame *hf2 = new TGHorizontalFrame(vf4, 0, 0);
   TGCompositeFrame *cf2 = new TGCompositeFrame(hf2, 110, 0, kFixedWidth);
   cf2->AddFrame(new TGLabel(cf2, "Processing rate:"));
   hf2->AddFrame(cf2);
   fRate = new TGLabel(hf2, "- events/sec \n");
   fRate->SetTextJustify(kTextTop | kTextLeft);
   hf2->AddFrame(fRate, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 0, 0));
   vf4->AddFrame(hf2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 10, 10, 5, 0));

   // Keep toggle button
   fKeepToggle = new TGCheckButton(vf4,
                    new TGHotString("Close dialog when processing is complete"));
   if (!fKeep) fKeepToggle->SetState(kButtonDown);
   fKeepToggle->Connect("Toggled(Bool_t)",
                        "TProofProgressDialog", this, "DoKeep(Bool_t)");
   vf4->AddFrame(fKeepToggle, new TGLayoutHints(kLHintsBottom, 10, 10, 10, 5));

   hf4->AddFrame(vf4, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));

   TGVerticalFrame *vf51 = new TGVerticalFrame(hf4, 20, 20);

   Int_t enablespeedo = gEnv->GetValue("Proof.EnableSpeedo", 0);
   if (enablespeedo) fSpeedoEnabled = kTRUE;

   fSpeedo = new TGSpeedo(vf51, 0.0, 1.0, "", "  Ev/s");
   if (fSpeedoEnabled) {
      fSpeedo->Connect("OdoClicked()", "TProofProgressDialog", this, "ToggleOdometerInfos()");
      fSpeedo->Connect("LedClicked()", "TProofProgressDialog", this, "ToggleThreshold()");
   }
   vf51->AddFrame(fSpeedo);
   fSpeedo->SetDisplayText("Init Time", "[ms]");
   fSpeedo->EnablePeakMark();
   fSpeedo->SetThresholds(0.0, 25.0, 50.0);
   fSpeedo->SetThresholdColors(TGSpeedo::kRed, TGSpeedo::kOrange, TGSpeedo::kGreen);
   fSpeedo->SetOdoValue(0);
   fSpeedo->EnableMeanMark();

   fSmoothSpeedo = new TGCheckButton(vf51, new TGHotString("Smooth speedometer update"));
   if (fSpeedoEnabled) {
      fSmoothSpeedo->SetState(kButtonDown);
      fSmoothSpeedo->SetToolTipText("Control smoothness in refreshing the speedo");
   } else {
      fSmoothSpeedo->SetToolTipText("Speedo refreshing is disabled");
      fSmoothSpeedo->SetState(kButtonDisabled);
   }
   vf51->AddFrame(fSmoothSpeedo, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 0));

   hf4->AddFrame(vf51, new TGLayoutHints(kLHintsBottom, 5, 5, 5, 5));

   fDialog->AddFrame(hf4, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

//==========================================================================================

   // Stop, cancel and close buttons
   TGHorizontalFrame *hf3 = new TGHorizontalFrame(fDialog, 60, 20);

   UInt_t  nb1 = 0, width1 = 0, height1 = 0;

   fAsyn = new TGTextButton(hf3, "&Run in background");
   if (fProof->GetRemoteProtocol() >= 22 && fProof->IsSync()) {
      fAsyn->SetToolTipText("Continue running in the background (asynchronous mode), releasing the ROOT prompt");
   } else {
      fAsyn->SetToolTipText("Switch to asynchronous mode disabled: functionality not supported by the server");
      fAsyn->SetState(kButtonDisabled);
   }
   fAsyn->Connect("Clicked()", "TProofProgressDialog", this, "DoAsyn()");
   hf3->AddFrame(fAsyn, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));
   height1 = TMath::Max(height1, fAsyn->GetDefaultHeight());
   width1  = TMath::Max(width1, fAsyn->GetDefaultWidth()); ++nb1;

   fStop = new TGTextButton(hf3, "&Stop");
   fStop->SetToolTipText("Stop processing, Terminate() will be executed");
   fStop->Connect("Clicked()", "TProofProgressDialog", this, "DoStop()");
   hf3->AddFrame(fStop, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));
   height1 = TMath::Max(height1, fStop->GetDefaultHeight());
   width1  = TMath::Max(width1, fStop->GetDefaultWidth()); ++nb1;

   fAbort = new TGTextButton(hf3, "&Cancel");
   fAbort->SetToolTipText("Cancel processing, Terminate() will NOT be executed");
   fAbort->Connect("Clicked()", "TProofProgressDialog", this, "DoAbort()");
   hf3->AddFrame(fAbort, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));
   height1 = TMath::Max(height1, fAbort->GetDefaultHeight());
   width1  = TMath::Max(width1, fAbort->GetDefaultWidth()); ++nb1;

   fClose = new TGTextButton(hf3, "&Close");
   fClose->SetToolTipText("Close this dialog");
   fClose->SetState(kButtonDisabled);
   fClose->Connect("Clicked()", "TProofProgressDialog", this, "DoClose()");
   hf3->AddFrame(fClose, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));
   height1 = TMath::Max(height1, fClose->GetDefaultHeight());
   width1  = TMath::Max(width1, fClose->GetDefaultWidth()); ++nb1;

   fDialog->AddFrame(hf3, new TGLayoutHints(kLHintsBottom | kLHintsCenterX | kLHintsExpandX, 5, 5, 5, 5));

   TGHorizontalFrame *hf5 = new TGHorizontalFrame(fDialog, 60, 20);

   fLog = new TGTextButton(hf5, "&Show Logs");
   fLog->SetToolTipText("Show query log messages");
   fLog->Connect("Clicked()", "TProofProgressDialog", this, "DoLog()");
   hf5->AddFrame(fLog, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));

   if (PPD_SRV_NEWER(11)) {
      fRatePlot = new TGTextButton(hf5, "&Performance plot");
      fRatePlot->SetToolTipText("Show rates, chunck sizes, cluster activities ... vs time");
      fRatePlot->SetState(kButtonDisabled);
      fRatePlot->Connect("Clicked()", "TProofProgressDialog", this, "DoPlotRateGraph()");
      hf5->AddFrame(fRatePlot, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));
   }

   fMemPlot = new TGTextButton(hf5, "&Memory Plot");
   fMemPlot->Connect("Clicked()", "TProofProgressDialog", this, "DoMemoryPlot()");
   fMemPlot->SetToolTipText("Show memory consumption vs entry / merging phase");
   hf5->AddFrame(fMemPlot, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));

   fUpdtSpeedo = new TGTextButton(hf5, "&Enable speedometer");
   fUpdtSpeedo->Connect("Clicked()", "TProofProgressDialog", this, "DoEnableSpeedo()");
   if (fSpeedoEnabled) {
      fUpdtSpeedo->ChangeText("&Disable speedometer");
      fUpdtSpeedo->SetToolTipText("Disable speedometer");
   } else {
      fUpdtSpeedo->ChangeText("&Enable speedometer");
      fUpdtSpeedo->SetToolTipText("Enable speedometer (may have an impact on performance)");
   }
   hf5->AddFrame(fUpdtSpeedo, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 7, 7, 0, 0));

   fDialog->AddFrame(hf5, new TGLayoutHints(kLHintsBottom | kLHintsCenterX | kLHintsExpandX, 5, 5, 5, 5));

   // Only enable if master supports it
   if (!PPD_SRV_NEWER(18)) {
      fMemPlot->SetState(kButtonDisabled);
      TString tip = TString::Format("Not supported by the master: required protocol 19 > %d",
                                    (fProof ? fProof->GetRemoteProtocol() : -1));
      fMemPlot->SetToolTipText(tip.Data());
   } else {
      fMemPlot->SetToolTipText("Show memory consumption");
   }

   // Connect slot to proof progress signal
   if (fProof) {
      fProof->Connect("Progress(Long64_t,Long64_t)", "TProofProgressDialog",
                      this, "Progress(Long64_t,Long64_t)");
      fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                      "TProofProgressDialog", this,
                      "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                      "TProofProgressDialog", this,
                      "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fProof->Connect("StopProcess(Bool_t)", "TProofProgressDialog", this,
                      "IndicateStop(Bool_t)");
      fProof->Connect("ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)",
                      "TProofProgressDialog", this,
                      "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)");
      fProof->Connect("CloseProgressDialog()", "TProofProgressDialog", this, "DoClose()");
      fProof->Connect("DisableGoAsyn()", "TProofProgressDialog", this, "DisableAsyn()");
   }

   // Set dialog title
   if (fProof) {
      if (strlen(fProof->GetUser()) > 0)
         fDialog->SetWindowName(Form("PROOF Query Progress: %s@%s",
                                     fProof->GetUser(), fProof->GetMaster()));
      else
         fDialog->SetWindowName(Form("PROOF Query Progress: %s", fProof->GetMaster()));
   } else
      fDialog->SetWindowName("PROOF Query Progress: <dummy>");

   // Map all widgets and calculate size of dialog
   fDialog->MapSubwindows();

   fDialog->Resize(fDialog->GetDefaultSize());

   const TGWindow *main = gClient->GetRoot();
   // Position relative to the parent window (which is the root window)
   Window_t wdum;
   int      ax, ay;
   Int_t    mw = ((TGFrame *) main)->GetWidth();
   Int_t    mh = ((TGFrame *) main)->GetHeight();
   Int_t    width  = fDialog->GetDefaultWidth();
   Int_t    height = fDialog->GetDefaultHeight();

   gVirtualX->TranslateCoordinates(main->GetId(), main->GetId(),
                          (mw - width), (mh - height) >> 1, ax, ay, wdum);

   // Make the message box non-resizable
   fDialog->SetWMSize(width, height);
   fDialog->SetWMSizeHints(width, height, width, height, 0, 0);

   fDialog->SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                                       kMWMDecorMinimize | kMWMDecorMenu,
                        kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                                       kMWMFuncMinimize,
                        kMWMInputModeless);

   fDialog->Move(ax-10, ay - mh/4);
   fDialog->SetWMPosition(ax-10, ay - mh/4);
   // Popup dialog and wait till user replies
   fDialog->MapWindow();

   fStartTime = gSystem->Now();
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle information displayed in Analog Meter

void TProofProgressDialog::ToggleOdometerInfos()
{
   if (fRightInfo < 1)
      fRightInfo++;
   else
      fRightInfo = 0;
   if (fRightInfo == 0) {
      fSpeedo->SetDisplayText("Init Time", "[ms]");
      fSpeedo->SetOdoValue((Int_t)(fInitTime * 1000.0));
   }
   else if (fRightInfo == 1) {
      fSpeedo->SetDisplayText("Proc Time", "[ms]");
      fSpeedo->SetOdoValue((Int_t)(fProcTime * 1000.0));
   }
}

////////////////////////////////////////////////////////////////////////////////

void TProofProgressDialog::ToggleThreshold()
{
   if (fSpeedo->IsThresholdActive()) {
      fSpeedo->DisableThreshold();
      fSpeedo->Glow(TGSpeedo::kNoglow);
   }
   else
      fSpeedo->EnableThreshold();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset dialog box preparing for new query

void TProofProgressDialog::ResetProgressDialog(const char *selec,
                                               Int_t files, Long64_t first,
                                               Long64_t entries)
{
   TString buf;

   // Update title
   buf.Form("Executing on PROOF cluster \"%s\" with %d parallel workers:",
            fProof ? fProof->GetMaster() : "<dummy>",
            fProof ? fProof->GetParallel() : 0);
   fTitleLab->SetText(buf);

   // Reset members
   fFiles         = files;
   fFirst         = first;
   fEntries       = entries;
   fPrevProcessed = 0;
   fPrevTotal     = 0;
   fStatus        = kRunning;

   // Update selector name
   buf.Form("Selector: %s", selec);
   fSelector->SetText(buf);

   // Reset 'estim' and 'processed' text
   fTimeLab->SetText("Estimated time left:");
   fProcessed->SetText("Processing status:");

   // Update numbers
   buf.Form("%d files, number of events %lld, starting event %lld",
            fFiles, fEntries, fFirst);
   fFilesEvents->SetText(buf);

   // Reset progress bar
   fBar->SetBarColor("green");
   fBar->Reset();

   // Reset speedo
   fSpeedo->SetMinMaxScale(0.0, 1.0);
   fSpeedo->SetMeanValue(0.0);
   fSpeedo->ResetPeakVal();

   // Reset buttons
   fStop->SetState(kButtonUp);
   fAbort->SetState(kButtonUp);
   fClose->SetState(kButtonDisabled);
   if (fProof && fProof->IsSync() && fProof->GetRemoteProtocol() >= 22) {
      fAsyn->SetState(kButtonUp);
   } else {
      fAsyn->SetState(kButtonDisabled);
   }

   // Reconnect the slots
   if (fProof) {
      fProof->Connect("Progress(Long64_t,Long64_t)", "TProofProgressDialog",
                      this, "Progress(Long64_t,Long64_t)");
      fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                      "TProofProgressDialog", this,
                      "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                      "TProofProgressDialog", this,
                      "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fProof->Connect("StopProcess(Bool_t)", "TProofProgressDialog", this,
                      "IndicateStop(Bool_t)");
      fProof->Connect("DisableGoAsyn()", "TProofProgressDialog", this, "DisableAsyn()");
   }

   // Reset start time
   fStartTime = gSystem->Now();

   // Clear the list of performances points
   if (PPD_SRV_NEWER(11))
      fRatePoints->Reset();
   SafeDelete(fRateGraph);
   SafeDelete(fMBRtGraph);
   SafeDelete(fActWGraph);
   SafeDelete(fTotSGraph);
   SafeDelete(fEffSGraph);
   fAvgRate = 0.;
   fAvgMBRate = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Update progress bar and status labels.
/// Use "processed == total" or "processed < 0" to indicate end of processing.

void TProofProgressDialog::Progress(Long64_t total, Long64_t processed)
{
   Long_t tt;
   UInt_t hh=0, mm=0, ss=0;
   TString buf;
   TString stm;
   static const char *cproc[] = { "running", "done",
                                  "STOPPED", "ABORTED", "***EVENTS SKIPPED***"};

   // Update title
   buf.Form("Executing on PROOF cluster \"%s\" with %d parallel workers:",
            fProof ? fProof->GetMaster() : "<dummy>",
            fProof ? fProof->GetParallel() : 0);
   fTitleLab->SetText(buf);

   if (total < 0)
      total = fPrevTotal;
   else
      fPrevTotal = total;

   // Nothing to update
   if (fPrevProcessed == processed)
      return;

   // Number of processed events
   Long64_t evproc = (processed >= 0) ? processed : fPrevProcessed;

   if (fEntries != total) {
      fEntries = total;
      buf.Form("%d files, number of events %lld, starting event %lld",
               fFiles, fEntries, fFirst);
      fFilesEvents->SetText(buf);
   }

   // Update position
   Float_t pos = Float_t(Double_t(evproc * 100)/Double_t(total));
   fBar->SetPosition(pos);

   // get current time
   fEndTime = gSystem->Now();
   TTime tdiff = fEndTime - fStartTime;
   Float_t eta = 0;
   if (evproc > 0)
      eta = ((Float_t)((Long64_t)tdiff)*total/Float_t(evproc) - Long64_t(tdiff))/1000.;

   if (processed >= 0 && processed >= total) {
      tt = (Long_t)Long64_t(tdiff)/1000;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm.Form("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm.Form("%d min %d sec", mm, ss);
      else
         stm.Form("%d sec", ss);
      fProcessed->SetText("Processed:");
      buf.Form("%lld events in %s\n", total, stm.Data());
      fTotal->SetText(buf);

      fEstim->SetText("0 sec");

      if (fProof) {
         fProof->Disconnect("Progress(Long64_t,Long64_t)", this,
                            "Progress(Long64_t,Long64_t)");
         fProof->Disconnect("StopProcess(Bool_t)", this,
                            "IndicateStop(Bool_t)");
         fProof->Disconnect("DisableGoAsyn()", this, "DisableAsyn()");
      }

      // Set button state
      fAsyn->SetState(kButtonDisabled);
      fStop->SetState(kButtonDisabled);
      fAbort->SetState(kButtonDisabled);
      fClose->SetState(kButtonUp);
      if (!fKeep) DoClose();

      // Set the status to done
      fStatus = kDone;

   } else {
      // A negative value for process indicates that we are finished,
      // no matter whether the processing was complete
      Bool_t incomplete = (processed < 0 &&
                          (fPrevProcessed < total || fPrevProcessed == 0))
                        ? kTRUE : kFALSE;
      if (incomplete) {
         fStatus = kIncomplete;
         // We use a different color to highlight incompletion
         fBar->SetBarColor("magenta");
      }
      tt = (Long_t)eta;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm.Form("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm.Form("%d min %d sec", mm, ss);
      else
         stm.Form("%d sec", ss);

      fEstim->SetText(stm.Data());
      buf.Form("%lld / %lld events", evproc, total);
      if (fStatus > kDone) {
         buf += TString::Format(" - %s", cproc[fStatus]);
      }
      fTotal->SetText(buf);

      buf.Form("%.1f events/sec\n", Float_t(evproc)/Long64_t(tdiff)*1000.);
      fRate->SetText(buf);

      if (processed < 0) {
         // And we disable the buttons
         fAsyn->SetState(kButtonDisabled);
         fStop->SetState(kButtonDisabled);
         fAbort->SetState(kButtonDisabled);
         fClose->SetState(kButtonUp);

         // Set the status to done
         fStatus = kDone;
      }
   }
   fPrevProcessed = evproc;
}

////////////////////////////////////////////////////////////////////////////////
/// Update progress bar and status labels.
/// Use "processed == total" or "processed < 0" to indicate end of processing.

void TProofProgressDialog::Progress(Long64_t total, Long64_t processed,
                                    Long64_t bytesread,
                                    Float_t initTime, Float_t procTime,
                                    Float_t evtrti, Float_t mbrti,
                                    Int_t actw, Int_t tses, Float_t eses)
{
   Double_t BinLow, BinHigh;
   Int_t nbins;
   Long_t tt;
   UInt_t hh=0, mm=0, ss=0;
   TString buf;
   TString stm;
   static const char *cproc[] = { "running", "done",
                                  "STOPPED", "ABORTED", "***EVENTS SKIPPED***"};

   // Update title
   buf.Form("Executing on PROOF cluster \"%s\" with %d parallel workers:",
            fProof ? fProof->GetMaster() : "<dummy>",
            fProof ? fProof->GetParallel() : 0);
   fTitleLab->SetText(buf);

   if (gDebug > 1)
      Info("Progress","t: %lld, p: %lld, itm: %f, ptm: %f", total, processed, initTime, procTime);

   if (initTime >= 0.) {
      // Set init time
      fInitTime = initTime;
      buf.Form("%.1f secs", initTime);
      fInit->SetText(buf);
      if (fSpeedoEnabled && fRightInfo == 0)
         fSpeedo->SetOdoValue((Int_t)(fInitTime * 1000.0));
   }

   Bool_t over = kFALSE;
   if (total < 0) {
      total = fPrevTotal;
      over = kTRUE;
   } else {
      fPrevTotal = total;
   }

   // Show proc time by default when switching from init to proc
   if (processed > 0 && fPrevProcessed <= 0)
      while (fRightInfo != 1)
         ToggleOdometerInfos();

   // Nothing to update
   if (fPrevProcessed == processed)
      return;

   // Number of processed events
   Long64_t evproc = (processed >= 0) ? processed : fPrevProcessed;
   Float_t mbsproc = bytesread / TMath::Power(2.,20.);

   if (fEntries != total) {
      fEntries = total;
      buf.Form("%d files, number of events %lld, starting event %lld",
               fFiles, fEntries, fFirst);
      fFilesEvents->SetText(buf);
   }

   // Update position
   Float_t pos = Float_t(Double_t(evproc * 100)/Double_t(total));
   fBar->SetPosition(pos);

   Float_t eta = 0;
   if (evproc > 0 && procTime > 0.)
      eta = (Float_t) (total - evproc) / (Double_t)evproc * procTime;

   // Update average rates
   if (procTime > 0.) {
      fProcTime = procTime;
      fAvgRate = Float_t(evproc) / procTime;
      fAvgMBRate = mbsproc / procTime;
   }

   if (fSpeedoEnabled) {
      if (fRightInfo == 0)
         fSpeedo->SetOdoValue((Int_t)(fInitTime * 1000.0));
      else if (fRightInfo == 1)
         fSpeedo->SetOdoValue((Int_t)(fProcTime * 1000.0));
   }

   if (over || (processed >= 0 && processed >= total)) {

      // A negative value for process indicates that we are finished,
      // no matter whether the processing was complete
      Bool_t incomplete = (processed < 0 &&
                          (fPrevProcessed < total || fPrevProcessed == 0))
                        ? kTRUE : kFALSE;
      TString st = "";
      if (incomplete) {
         fStatus = kIncomplete;
         // We use a different color to highlight incompletion
         fBar->SetBarColor("magenta");
         st = TString::Format(" %s", cproc[fStatus]);
      }

      tt = (Long_t)fProcTime;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm.Form("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm.Form("%d min %d sec", mm, ss);
      else
         stm.Form("%d sec", ss);
      fProcessed->SetText("Processed:");
      TString sf("MB");
      Float_t xb = fAvgMBRate*fProcTime;
      xb = AdjustBytes(xb, sf);
      buf.Form("%lld events (%.2f %s)\n",
               std::max(fPrevProcessed, processed), xb, sf.Data());
      fTotal->SetText(buf);
      buf.Form("%s %s\n", stm.Data(), st.Data());
      fTimeLab->SetText("Processing time:");
      fEstim->SetText(buf);
      buf.Form("%.1f evts/sec (%.1f MB/sec)\n", fAvgRate, fAvgMBRate);
      fRate->SetText(buf);
      // Fill rate graph
      Bool_t useAvg = gEnv->GetValue("Proof.RatePlotUseAvg", 0);
      if (useAvg) {
         if (fAvgRate > 0.) {
            fRatePoints->Fill(procTime, fAvgRate, fAvgMBRate);
            fRatePlot->SetState(kButtonUp);
         }
      } else {
         if (evtrti > 0.) {
            fRatePoints->Fill(procTime, evtrti, mbrti, (Float_t)actw, (Float_t)tses, eses);
            fRatePlot->SetState(kButtonUp);
         }
      }

      if (fProof) {
         fProof->Disconnect("Progress(Long64_t,Long64_t)", this,
                            "Progress(Long64_t,Long64_t)");
         fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                            this,
                            "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
         fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                            this,
                            "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
         fProof->Disconnect("StopProcess(Bool_t)", this, "IndicateStop(Bool_t)");
         fProof->Disconnect("DisableGoAsyn()", this, "DisableAsyn()");
      }

      // Set button state
      fAsyn->SetState(kButtonDisabled);
      fStop->SetState(kButtonDisabled);
      fAbort->SetState(kButtonDisabled);
      fClose->SetState(kButtonUp);

      if (fSmoothSpeedo->GetState() == kButtonDown)
         fSpeedo->SetScaleValue(0.0, 0);
      else
         fSpeedo->SetScaleValue(0.0);
      fSpeedo->Glow(TGSpeedo::kNoglow);

      if (!fKeep) DoClose();

      // Set the status to done
      fStatus = kDone;

   } else {
      // A negative value for process indicates that we are finished,
      // no matter whether the processing was complete
      Bool_t incomplete = (processed < 0 &&
                          (fPrevProcessed < total || fPrevProcessed == 0))
                        ? kTRUE : kFALSE;
      if (incomplete) {
         fStatus = kIncomplete;
         // We use a different color to highlight incompletion
         fBar->SetBarColor("magenta");
      }
      tt = (Long_t)eta;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm.Form("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm.Form("%d min %d sec", mm, ss);
      else
         stm.Form("%d sec", ss);

      fEstim->SetText(stm.Data());
      TString sf("MB");
      Float_t xb = AdjustBytes(mbsproc, sf);
      buf.Form("%lld / %lld events - %.2f %s", evproc, total, xb, sf.Data());
      if (fStatus > kDone) {
         buf += TString::Format(" - %s", cproc[fStatus]);
      }
      fTotal->SetText(buf);

      // Post
      if (evtrti > 0.) {
         buf.Form("%.1f evts/sec \navg: %.1f evts/sec (%.1f MB/sec)",
                  evtrti, fAvgRate, fAvgMBRate);
         fRatePoints->Fill(procTime, evtrti, mbrti, (Float_t)actw, (Float_t)tses, eses);
         fRatePlot->SetState(kButtonUp);
         if (fSpeedoEnabled) {
            if (evtrti > fSpeedo->GetScaleMax()) {
               nbins = 4;
               BinLow = fSpeedo->GetScaleMin();
               BinHigh = 1.5 * evtrti;
               THLimitsFinder::OptimizeLimits(4, nbins, BinLow, BinHigh, kFALSE);
               fSpeedo->SetMinMaxScale(fSpeedo->GetScaleMin(), BinHigh);
            }
            if (fSmoothSpeedo->GetState() == kButtonDown)
               fSpeedo->SetScaleValue(evtrti, 0);
            else
               fSpeedo->SetScaleValue(evtrti);
            fSpeedo->SetMeanValue(fAvgRate);
         }
      } else {
         buf.Form("avg: %.1f evts/sec (%.1f MB/sec)", fAvgRate, fAvgMBRate);
      }
      fRate->SetText(buf);

      if (processed < 0) {
         // And we disable the buttons
         fAsyn->SetState(kButtonDisabled);
         fStop->SetState(kButtonDisabled);
         fAbort->SetState(kButtonDisabled);
         fClose->SetState(kButtonUp);

         if (fSpeedoEnabled) {
            if (fSmoothSpeedo->GetState() == kButtonDown)
               fSpeedo->SetScaleValue(0.0, 0);
            else
               fSpeedo->SetScaleValue(0.0);
            fSpeedo->Glow(TGSpeedo::kNoglow);
         }

         // Set the status to done
         fStatus = kDone;
      }
   }
   fPrevProcessed = evproc;
}

////////////////////////////////////////////////////////////////////////////////
/// Transform MBs to GBs ot TBs and get the correct suffix

Float_t TProofProgressDialog::AdjustBytes(Float_t mbs, TString &sf)
{
   Float_t xb = mbs;
   sf = "MB";
   if (xb > 1024.) {
      xb = xb / 1024.;
      sf = "GB";
   }
   if (xb > 1024.) {
      xb = xb / 1024.;
      sf = "TB";
   }
   // Done
   return xb;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup dialog.

TProofProgressDialog::~TProofProgressDialog()
{
   if (fProof) {
      fProof->Disconnect("Progress(Long64_t,Long64_t)", this,
                         "Progress(Long64_t,Long64_t)");
      fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                         this,
                         "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                         this,
                         "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fProof->Disconnect("StopProcess(Bool_t)", this, "IndicateStop(Bool_t)");
      fProof->Disconnect("DisableGoAsyn()", this, "DisableAsyn()");
      fProof->Disconnect("ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)",
                         this,
                         "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)");
      fProof->Disconnect("CloseProgressDialog()", this, "CloseProgressDialog()");
      fProof->ResetProgressDialogStatus();
      // We are called after a TProofDetach: we delete the instance
      if (!fProof->IsValid())
         SafeDelete(fProof);
   }
   if (fLogWindow)
      delete fLogWindow;
   if (fMemWindow)
      delete fMemWindow;
   fDialog->Cleanup();
   delete fDialog;
}

////////////////////////////////////////////////////////////////////////////////
/// Called when dialog is closed.

void TProofProgressDialog::CloseWindow()
{
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Disable the asyn switch when an external request for going asynchronous is issued

void TProofProgressDialog::DisableAsyn()
{
   fProof->Disconnect("DisableGoAsyn()", this, "DisableAsyn()");
   fAsyn->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Indicate that Cancel or Stop was clicked.

void TProofProgressDialog::IndicateStop(Bool_t aborted)
{
   if (aborted == kTRUE)
      fBar->SetBarColor("red");
   else
      fBar->SetBarColor("yellow");

   if (fProof) {
      fProof->Disconnect("Progress(Long64_t,Long64_t)", this,
                         "Progress(Long64_t,Long64_t)");
      fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                         this,
                         "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                         this,
                         "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fProof->Disconnect("StopProcess(Bool_t)", this, "IndicateStop(Bool_t)");
      fProof->Disconnect("DisableGoAsyn()", this, "DisableAsyn()");
      // These buttons are meaningless at this point
      fAsyn->SetState(kButtonDisabled);
      fStop->SetState(kButtonDisabled);
      fAbort->SetState(kButtonDisabled);
   }

   fClose->SetState(kButtonUp);
   if (!fKeep)
      DoClose();
}

////////////////////////////////////////////////////////////////////////////////
/// Load/append a log msg in the log frame, if open

void TProofProgressDialog::LogMessage(const char *msg, Bool_t all)
{
   if (fLogWindow) {
      if (all) {
         // load buffer
         fLogWindow->LoadBuffer(msg);
      } else {
         // append
         fLogWindow->AddBuffer(msg);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close dialog.

void TProofProgressDialog::DoClose()
{
   fClose->SetState(kButtonDisabled);
   TTimer::SingleShot(50, "TProofProgressDialog", this, "CloseWindow()");
}

////////////////////////////////////////////////////////////////////////////////
/// Ask proof session for logs

void TProofProgressDialog::DoLog()
{
   if (fProof) {
      if (!fLogWindow) {
         fLogWindow = new TProofProgressLog(this);
         if (!fLogWindow->TestBit(TObject::kInvalidObject))
            fLogWindow->DoLog();
      } else {
         // Clear window
         if (!fLogWindow->TestBit(TObject::kInvalidObject)) {
            fLogWindow->Clear();
            fLogWindow->DoLog();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keep toggle button.

void TProofProgressDialog::DoKeep(Bool_t)
{
   fKeep = !fKeep;

   // Last choice will be the default for the future
   fgKeepDefault = fKeep;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle log-current-query-only toggle button.

void TProofProgressDialog::DoSetLogQuery(Bool_t)
{
   fLogQuery = !fLogQuery;
   fEntry->SetEnabled(fLogQuery);
   if (fLogQuery)
      fEntry->SetToolTipText("Enter the query number ('last' for the last query)",50);
   else
      fEntry->SetToolTipText(0);

   // Last choice will be the default for the future
   fgLogQueryDefault = fLogQuery;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Stop button.

void TProofProgressDialog::DoStop()
{
   // Do not wait for ever, but al least 10 seconds
   Long_t timeout = gEnv->GetValue("Proof.ShutdownTimeout", 60) / 2;
   timeout = (timeout > 10) ? timeout : 10;
   fProof->StopProcess(kFALSE, timeout);
   fStatus = kStopped;

   // Set buttons states
   fAsyn->SetState(kButtonDisabled);
   fStop->SetState(kButtonDisabled);
   fAbort->SetState(kButtonDisabled);
   fClose->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Cancel button.

void TProofProgressDialog::DoAbort()
{
   fProof->StopProcess(kTRUE);
   fStatus = kAborted;

   // Set buttons states
   fAsyn->SetState(kButtonDisabled);
   fStop->SetState(kButtonDisabled);
   fAbort->SetState(kButtonDisabled);
   fClose->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Asyn button.

void TProofProgressDialog::DoAsyn()
{
   fProof->GoAsynchronous();

   // Set buttons states
   fAsyn->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Plot Rate Graph.

void TProofProgressDialog::DoPlotRateGraph()
{
   // We must have some point to plot
   if (!fRatePoints || fRatePoints->GetEntries() <= 0) {
      Info("DoPlotRateGraph","list is empty!");
      return;
   }

   // Fill the graphs
   Int_t np = (Int_t)fRatePoints->GetEntries();
   Double_t eymx = -1., bymx = -1., wymx = -1., tymx=-1., symx = -1.;
   SafeDelete(fRateGraph);
   SafeDelete(fMBRtGraph);
   SafeDelete(fActWGraph);
   SafeDelete(fTotSGraph);
   SafeDelete(fEffSGraph);
   fRateGraph = new TGraph(np);
   fMBRtGraph = new TGraph(np);
   if (PPD_SRV_NEWER(25)) {
      fActWGraph = new TGraph(np);
      fTotSGraph = new TGraph(np);
      fEffSGraph = new TGraph(np);
   }
   Float_t *nar = fRatePoints->GetArgs();
   Int_t ii = 0;
   for ( ; ii < np; ++ii) {
      fRatePoints->GetEntry(ii);
      if (!(nar[1] > 0.)) continue;
      // Evts/s
      fRateGraph->SetPoint(ii, (Double_t) nar[0], (Double_t) nar[1]);
      eymx = (nar[1] > eymx) ? nar[1] : eymx;
      // MBs/s
      fMBRtGraph->SetPoint(ii, (Double_t) nar[0], (Double_t) nar[2]);
      bymx = (nar[2] > bymx) ? nar[2] : bymx;
      // Active workers
      if (PPD_SRV_NEWER(25)) {
         fActWGraph->SetPoint(ii, (Double_t) nar[0], (Double_t) nar[3]);
         wymx = (nar[3] > wymx) ? nar[3] : wymx;
      }
      // Sessions info
      if (PPD_SRV_NEWER(25)) {
         fTotSGraph->SetPoint(ii, (Double_t) nar[0], (Double_t) nar[4]);
         tymx = (nar[4] > tymx) ? nar[4] : tymx;
         fEffSGraph->SetPoint(ii, (Double_t) nar[0], (Double_t) nar[5]);
         symx = (nar[5] > symx) ? nar[5] : symx;
      }
   }

   // Pad numbering
   Int_t npads = 4;
   Int_t kEvrt = 1;
   Int_t kMBrt = 2;
   Int_t kActW = 3;
   Int_t kSess = 4;
   if (bymx <= 0.) {
      SafeDelete(fMBRtGraph);
      npads--;
      kActW--;
      kSess--;
   }
   if (wymx <= 0.) {
      SafeDelete(fActWGraph);
      npads--;
      kSess--;
   }
   // Plot only if more than one active session during the query
   if (tymx <= 1.) {
      SafeDelete(fTotSGraph);
      SafeDelete(fEffSGraph);
      npads--;
      kSess--;
   }
   if (tymx <= 0.) SafeDelete(fTotSGraph);
   if (symx <= 0.) SafeDelete(fEffSGraph);

   // Create a canvas
   Int_t jsz = 200*npads;
   TCanvas *c1 = new TCanvas("c1","Rate vs Time",200,10,700,jsz);
   c1->SetFillColor(0);
   c1->SetGrid();
   c1->SetBorderMode(0);
   c1->SetFrameBorderMode(0);

   // Padding
   c1->Divide(1, npads);

   // Event Rate plot
   TPad *cpad = (TPad *) c1->GetPad(kEvrt);
   if (cpad) {
      cpad->cd();
      cpad->SetFillColor(0);
      cpad->SetBorderMode(20);
      cpad->SetFrameBorderMode(0);
   }
   fRateGraph->SetMinimum(0.);
   fRateGraph->SetMaximum(eymx*1.1);
   fRateGraph->SetLineColor(50);
   fRateGraph->SetLineWidth(2);
   fRateGraph->SetMarkerColor(38);
   fRateGraph->SetMarkerStyle(25);
   fRateGraph->SetMarkerSize(0.8);
   fRateGraph->SetTitle("Processing rate (evts/sec)");
   fRateGraph->GetXaxis()->SetTitle("elapsed time (sec)");
   fRateGraph->Draw("ALP");

   // Line with average
   TLine *line = new TLine(fRateGraph->GetXaxis()->GetXmin(),fAvgRate,
                           fRateGraph->GetXaxis()->GetXmax(),fAvgRate);
   line->SetLineColor(8);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();

   // Label
   Double_t xax0 = fRateGraph->GetXaxis()->GetXmin();
   Double_t xax1 = fRateGraph->GetXaxis()->GetXmax();
   Double_t yax0 = 0.;
   Double_t yax1 = eymx*1.1;
   Double_t x0 = xax0 + 0.05 * (xax1 - xax0);
   Double_t x1 = xax0 + 0.60 * (xax1 - xax0);
   Double_t y0 = yax0 + 0.10 * (yax1 - yax0);
   Double_t y1 = yax0 + 0.20 * (yax1 - yax0);
   TPaveText *pt = new TPaveText(x0, y0, x1, y1, "br");
   pt->SetFillColor(0);
   pt->AddText(Form("Global average: %.2f evts/sec", fAvgRate));
   pt->Draw();

   // MB Rate plot
   if (fMBRtGraph) {
      cpad = (TPad *) c1->GetPad(kMBrt);
      if (cpad) {
         cpad->cd();
         cpad->SetFillColor(0);
         cpad->SetBorderMode(0);
         cpad->SetFrameBorderMode(0);
      }
      fMBRtGraph->SetFillColor(38);
      TH1F *graph2 = new TH1F("graph2","Average read chunck size (MBs/request)",100,
                               fRateGraph->GetXaxis()->GetXmin(),fRateGraph->GetXaxis()->GetXmax());
      graph2->SetMinimum(0);
      graph2->SetMaximum(1.1*bymx);
      graph2->SetDirectory(0);
      graph2->SetStats(0);
      graph2->GetXaxis()->SetTitle("elapsed time (sec)");
      fMBRtGraph->SetHistogram(graph2);
      fMBRtGraph->Draw("AB");
   }

   // MB Rate plot
   if (fActWGraph) {
      cpad = (TPad *) c1->GetPad(kActW);
      if (cpad) {
         cpad->cd();
         cpad->SetFillColor(0);
         cpad->SetBorderMode(0);
         cpad->SetFrameBorderMode(0);
      }
      fActWGraph->SetMinimum(0.);
      fActWGraph->SetMaximum(wymx*1.1);
      fActWGraph->SetLineColor(50);
      fActWGraph->SetLineWidth(2);
      fActWGraph->SetMarkerColor(38);
      fActWGraph->SetMarkerStyle(25);
      fActWGraph->SetMarkerSize(0.8);
      fActWGraph->SetTitle("Active workers");
      fActWGraph->GetXaxis()->SetTitle("elapsed time (sec)");
      fActWGraph->Draw("ALP");
   }

   // MB Rate plot
   if (fTotSGraph) {
      cpad = (TPad *) c1->GetPad(kSess);
      if (cpad) {
         cpad->cd();
         cpad->SetFillColor(0);
         cpad->SetBorderMode(0);
         cpad->SetFrameBorderMode(0);
      }
      fTotSGraph->SetMinimum(0.);
      fTotSGraph->SetMaximum(tymx*1.1);
      fTotSGraph->SetLineColor(50);
      fTotSGraph->SetLineWidth(2);
      fTotSGraph->SetMarkerColor(38);
      fTotSGraph->SetMarkerStyle(25);
      fTotSGraph->SetMarkerSize(0.8);
      fTotSGraph->SetTitle("Active, Effective sessions");
      fTotSGraph->GetXaxis()->SetTitle("elapsed time (sec)");
      fTotSGraph->Draw("ALP");

      // Effective sessions
      if (fEffSGraph) {
         fEffSGraph->SetMinimum(0.);
         fEffSGraph->SetMaximum(tymx*1.1);
         fEffSGraph->SetLineColor(38);
         fEffSGraph->SetLineWidth(2);
         fEffSGraph->SetMarkerColor(50);
         fEffSGraph->SetMarkerStyle(21);
         fEffSGraph->SetMarkerSize(0.6);
         fEffSGraph->Draw("SLP");
      }
   }

   c1->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Do a memory plot

void TProofProgressDialog::DoMemoryPlot()
{
   if (!fMemWindow) {
      fMemWindow = new TProofProgressMemoryPlot(this, 500, 300);
      fMemWindow->DoPlot();
   } else {
      // Clear window
      fMemWindow->Clear();
      fMemWindow->DoPlot();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable speedometer

void TProofProgressDialog::DoEnableSpeedo()
{
   if (!fSpeedoEnabled) {
      // Enable and connect
      fSpeedoEnabled = kTRUE;
      fSpeedo->Connect("OdoClicked()", "TProofProgressDialog", this, "ToggleOdometerInfos()");
      fSpeedo->Connect("LedClicked()", "TProofProgressDialog", this, "ToggleThreshold()");
      fUpdtSpeedo->ChangeText("&Disable speedometer");
      fUpdtSpeedo->SetToolTipText("Disable speedometer");
      fSmoothSpeedo->SetState(kButtonDown);
      fSmoothSpeedo->SetToolTipText("Control smoothness in refreshing the speedo");
   } else {
      // Disable and disconnect
      fSpeedoEnabled = kFALSE;
      // Reset speedo
      fSpeedo->SetScaleValue(0);
      fUpdtSpeedo->ChangeText("&Enable speedometer");
      fUpdtSpeedo->SetToolTipText("Enable speedometer (may have an impact on performance)");
      fSmoothSpeedo->SetToolTipText("Speedo refreshing is disabled");
      fSmoothSpeedo->SetState(kButtonDisabled);
   }
}
