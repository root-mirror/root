// @(#)root/proofplayer:$Id$
// Author: G. Ganis May 2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TStatsFeedback
\ingroup proofkernel

Utility class to display PROOF stats feedback histos during queries

*/

#include "TStatsFeedback.h"

#include "TError.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TProof.h"
#include "TProofDebug.h"
#include "TROOT.h"
#include "TSeqCollection.h"
#include "TStyle.h"
#include "TVirtualPad.h"

ClassImp(TStatsFeedback);


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TStatsFeedback::TStatsFeedback(TProof *proof)
{
   if (proof == 0) proof = gProof;

   TProof *p = dynamic_cast<TProof*>(proof);
   if (p == 0) {
      Error("TStatsFeedback", "no valid proof session found");
      SetBit(TObject::kInvalidObject);
      return;
   }
   fProof = p;
   fName = fProof->GetSessionTag();

   if (!(proof->Connect("Feedback(TList*)", "TStatsFeedback",
                  this, "Feedback(TList*)"))) {
      Error("TStatsFeedback", "Connect() failed");
      SetBit(TObject::kInvalidObject);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TStatsFeedback::~TStatsFeedback()
{
   // Required since we overload TObject::Hash.
   ROOT::CallRecursiveRemoveIfNeeded(*this);

   fProof->Disconnect("Feedback(TList*)", this, "Feedback(TList*");
}

////////////////////////////////////////////////////////////////////////////////
/// Display feedback

void TStatsFeedback::Feedback(TList *objs)
{
   TSeqCollection *canvases = gROOT->GetListOfCanvases();

   PDB(kFeedback,1) Info("Feedback", "%d Objects", objs->GetSize());

   // Attach to the histograms we want to plot
   TH1D *hevt = 0, *hpck = 0;
   TH1I *hass = 0;
   TIter next(objs);
   TObject *o = 0;
   while((o = next())) {
      if (!strcmp(o->GetName(), "PROOF_EventsHist")) {
         hevt = dynamic_cast<TH1D *>(o);
      } else if (!strcmp(o->GetName(), "PROOF_PacketsHist")) {
         hpck = dynamic_cast<TH1D *>(o);
      } else if (!strcmp(o->GetName(), "PROOF_ProcPcktHist")) {
         hass = dynamic_cast<TH1I *>(o);
      }
      if (hevt && hpck && hass) break;
   }
   if (!hevt && !hpck && !hass) {
      Warning("Feedback", "none of the requested histograms has been found!");
      return;
   }

   // Number of histograms
   Int_t nh = 3;
   if (!hass) nh = 2;
   // Create or attach to canvas
   TString cvnm = TString::Format("Stats: %s", fProof->GetSessionTag());
   TVirtualPad *cv = 0;
   if (gROOT->GetListOfCanvases())
      cv = (TVirtualPad *) canvases->FindObject(cvnm.Data());
   if (cv && nh == 3 && !cv->GetPad(3)) SafeDelete(cv);
   if (!cv) {
      Int_t h = (nh == 3) ? 600 : 400;
      TString cvcmd = TString::Format("new TCanvas(\"%s\", \"Feedback Stats\",10,300,600,%d)",
                                      cvnm.Data(), h);
      if (!(cv = (TVirtualPad *) gROOT->ProcessLine(cvcmd))) {
         Warning("Feedback", "could not create canvas!");
         return;
      }
      PDB(kFeedback,2) Info("Feedback", "created canvas %s", cvnm.Data());
      // Create pads
      cv->Divide(1, nh);
   } else {
      cv->cd();
      PDB(kFeedback,2) Info("Feedback", "using canvas %s", cvnm.Data());
   }
   TVirtualPad *pd1 = (TVirtualPad *) cv->GetPad(1);
   TVirtualPad *pd2 = (TVirtualPad *) cv->GetPad(2);
   TVirtualPad *pd3 = (nh == 3) ? (TVirtualPad *) cv->GetPad(3) : 0;

   UInt_t optstat = gStyle->GetOptStat();
   gStyle->SetOptStat(11);
   // Plot
   if (hevt) {
      if (pd1) pd1->cd();
      hevt->SetFillColor(kGreen);
      hevt->DrawCopy();
   }
   if (hpck) {
      if (pd2) pd2->cd();
      hpck->SetFillColor(kAzure-5);
      hpck->DrawCopy();
   }
   if (hass) {
      if (pd3) pd3->cd();
      hass->SetFillColor(kGray);
      hass->SetMaximum(2);
      hass->DrawCopy();
   }

   cv->cd();
   cv->Update();
   gStyle->SetOptStat(optstat);
}
