// @(#)root/minuit2:$Id$
// Author:  L. Moneta 2012

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2012 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "TMinuit2TraceObject.h"
#include "TH1.h"
#include "TVirtualPad.h"
#include "TList.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MinimumState.h"

ClassImp(TMinuit2TraceObject);

TMinuit2TraceObject::TMinuit2TraceObject(int parNumber)
   : ROOT::Minuit2::MnTraceObject(parNumber), TNamed("Minuit2TraceObject", "ROOT Trace Object for Minuit2"),
     fIterOffset(0), fHistoFval(0), fHistoEdm(0), fHistoParList(0), fOldPad(0), fMinuitPad(0)
{
}

TMinuit2TraceObject::~TMinuit2TraceObject()
{
   // rest previous pad but do not delete histograms
   if (fOldPad && gPad && fOldPad != gPad)
      gPad = fOldPad;
   int niter = -1;
   if (fHistoFval) {
      niter = int(fHistoFval->GetEntries() + 0.5);
      fHistoFval->GetXaxis()->SetRange(1, niter);
   }
   if (fHistoEdm)
      fHistoEdm->GetXaxis()->SetRange(1, niter);
   if (fHistoParList) {
      for (int i = 0; i < fHistoParList->GetSize(); ++i) {
         TH1 *h1 = (TH1 *)fHistoParList->At(i);
         if (h1)
            h1->GetXaxis()->SetRange(1, niter);
      }
   }
}

void TMinuit2TraceObject::Init(const ROOT::Minuit2::MnUserParameterState &state)
{

   ROOT::Minuit2::MnTraceObject::Init(state);

   fIterOffset = 0;

   // build debug histogram
   if (fHistoFval)
      delete fHistoFval;
   if (fHistoEdm)
      delete fHistoEdm;
   if (fHistoParList) {
      fHistoParList->Delete();
      delete fHistoParList;
   }
   if (fMinuitPad)
      delete fMinuitPad;

   fHistoFval = new TH1D("minuit2_hist_fval", "Function Value/iteration", 2, 0, 1);
   fHistoEdm = new TH1D("minuit2_hist_edm", "Edm/iteration", 2, 0, 1);
   fHistoFval->SetCanExtend(TH1::kAllAxes);
   fHistoEdm->SetCanExtend(TH1::kAllAxes);

   // create histos for all parameters
   fHistoParList = new TList();
   for (unsigned int ipar = 0; ipar < state.Params().size(); ++ipar) {
      if (state.Parameter(ipar).IsFixed() || state.Parameter(ipar).IsConst())
         continue;
      TH1D *h1 = new TH1D(TString::Format("minuit2_hist_par%d", ipar),
                          TString::Format("Value of %s/iteration", state.Name(ipar)), 2, 0, 1);
      h1->SetCanExtend(TH1::kAllAxes);
      fHistoParList->Add(h1);
   }

   if (gPad)
      fOldPad = gPad;

   // fMinuitPad = new TCanvas("c1_minuit2","TMinuit2 Progress",2);
   // fMinuitPad->Divide(1,3);
   // fMinuitPad->cd(1); fHistoFval->Draw();
   // fMinuitPad->cd(2); fHistoEdm->Draw();
   // fMinuitPad->cd(3); fHistoPar->Draw();
   fHistoFval->Draw("hist");
   fMinuitPad = gPad;
}

void TMinuit2TraceObject::operator()(int iter, const ROOT::Minuit2::MinimumState &state)
{
   // action for each iteration: fill histograms
   // if iteration number is < 0 add at the end of current histograms
   // if offset is > 0 start filling from end of previous histogram

   int lastIter = int(fHistoFval->GetEntries() + 0.5);
   if (iter < 0)
      iter = lastIter;
   else {
      if (iter == 0 && lastIter > 0)
         fIterOffset = lastIter;

      iter += fIterOffset;
   }

   ROOT::Minuit2::MnTraceObject::operator()(iter, state);

   fHistoFval->SetBinContent(iter + 1, state.Fval());
   fHistoEdm->SetBinContent(iter + 1, state.Edm());

   for (unsigned int ipar = 0; ipar < state.Vec().size(); ++ipar) {
      double eval = UserState().Trafo().Int2ext(ipar, state.Vec()(ipar));
      TH1 *histoPar = (TH1 *)fHistoParList->At(ipar);
      histoPar->SetBinContent(iter + 1, eval);
   }

   if (fMinuitPad) {
      if (ParNumber() == -2)
         fHistoEdm->Draw();
      else if (ParNumber() >= 0 && ParNumber() < fHistoParList->GetSize()) {
         fHistoParList->At(ParNumber())->Draw();
      } else
         fHistoFval->Draw();
   }
}
