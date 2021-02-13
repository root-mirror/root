// @(#)root/proof:$Id$
// Author: Sangsu Ryu 22/06/2010

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProofBenchRunCPU
\ingroup proofbench

CPU-intensive PROOF benchmark test generates events and fill 1, 2, or 3-D histograms.
No I/O activity is involved.  

*/

#include "RConfigure.h"

#include "TProofBenchRunCPU.h"
#include "TProofNodes.h"
#include "TProofPerfAnalysis.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TProof.h"
#include "TString.h"
#include "Riostream.h"
#include "TMap.h"
#include "TEnv.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TPerfStats.h"
#include "TQueryResult.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"

ClassImp(TProofBenchRunCPU);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TProofBenchRunCPU::TProofBenchRunCPU(TPBHistType *histtype, Int_t nhists,
                                     TDirectory* dirproofbench, TProof* proof,
                                     TProofNodes* nodes, Long64_t nevents, Int_t ntries,
                                     Int_t start, Int_t stop, Int_t step, Int_t draw,
                                     Int_t debug)
                  : TProofBenchRun(proof, kPROOF_BenchSelCPUDef),
                    fHistType(histtype), fNHists(nhists),
                    fNEvents(nevents), fNTries(ntries), fStart(start), fStop(stop),
                    fStep(step), fDraw(draw), fDebug(debug), fDirProofBench(dirproofbench),
                    fNodes(nodes), fListPerfPlots(0),
                    fCanvas(0), fProfile_perfstat_event(0), fHist_perfstat_event(0),
                    fProfile_perfstat_evtmax(0), fNorm_perfstat_evtmax(0),
                    fProfile_queryresult_event(0), fNorm_queryresult_event(0), fProfile_cpu_eff(0),
                    fProfLegend(0), fNormLegend(0), fName(0)
{
   if (TestBit(kInvalidObject)) {
      Error("TProofBenchRunCPU", "problems validating PROOF session or enabling selector PAR");
      return;
   }

   fName = "CPU";

   if (!fNodes) fNodes = new TProofNodes(fProof);

   if (stop == -1) fStop = fNodes->GetNWorkersCluster();

   fListPerfPlots = new TList;

   gEnv->SetValue("Proof.StatsTrace",1);
   gStyle->SetOptStat(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofBenchRunCPU::~TProofBenchRunCPU()
{
   fProof=0;
   fDirProofBench=0;
   SafeDelete(fListPerfPlots);
   SafeDelete(fCanvas);
   SafeDelete(fNodes);
   SafeDelete(fProfLegend);
   SafeDelete(fNormLegend);
}

////////////////////////////////////////////////////////////////////////////////
/// Build histograms, profiles and graphs needed for this run

void TProofBenchRunCPU::BuildHistos(Int_t start, Int_t stop, Int_t step, Bool_t nx)
{
   TObject *o = 0;
   Int_t quotient = (stop - start) / step;
   Int_t ndiv = quotient + 1;
   Double_t ns_min = start - step/2.;
   Double_t ns_max = quotient*step + start + step/2.;

   fProfLegend = new TLegend(0.1, 0.8, 0.3, 0.9);
   fNormLegend = new TLegend(0.7, 0.8, 0.9, 0.9);

   TString axtitle("Active Workers"), namelab(GetName()), sellab(GetSelName());
   if (nx) {
      axtitle = "Active Workers/Node";
      namelab.Form("x_%s", GetName());
   }
   if (fSelName == kPROOF_BenchSelCPUDef)
      sellab.Form("%s_%s", GetSelName(), GetNameStem().Data());

   TString name, title;

   // Book perfstat profile (max evts)
   name.Form("Prof_%s_PS_MaxEvts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s PerfStat Event - %s", namelab.Data(), sellab.Data());
   fProfile_perfstat_evtmax = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_perfstat_evtmax->SetDirectory(fDirProofBench);
   fProfile_perfstat_evtmax->GetYaxis()->SetTitle("Events/sec");
   fProfile_perfstat_evtmax->GetXaxis()->SetTitle(axtitle);
   fProfile_perfstat_evtmax->SetMarkerStyle(23);
   fProfile_perfstat_evtmax->SetMarkerColor(2);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_perfstat_evtmax);
   fProfLegend->AddEntry(fProfile_perfstat_evtmax, "Maximum");

   // Book perfstat profile
   name.Form("Prof_%s_PS_Evts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s PerfStat Event - %s", namelab.Data(), sellab.Data());
   fProfile_perfstat_event = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_perfstat_event->SetDirectory(fDirProofBench);
   fProfile_perfstat_event->GetYaxis()->SetTitle("Events/sec");
   fProfile_perfstat_event->GetXaxis()->SetTitle(axtitle);
   fProfile_perfstat_event->SetMarkerStyle(21);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_perfstat_event);
   fProfLegend->AddEntry(fProfile_perfstat_event, "Average");

   // Book perfstat histogram
   name.Form("Hist_%s_PS_Evts_%s", namelab.Data(), sellab.Data());
   title.Form("Histogram %s PerfStat Event - %s", namelab.Data(), sellab.Data());
   fHist_perfstat_event = new TH2D(name, title, ndiv, ns_min, ns_max, 100, 0, 0);
   fHist_perfstat_event->SetDirectory(fDirProofBench);
   fHist_perfstat_event->GetYaxis()->SetTitle("Events/sec");
   fHist_perfstat_event->GetXaxis()->SetTitle(axtitle);
   fHist_perfstat_event->SetMarkerStyle(7);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fHist_perfstat_event);

   // Book normalized perfstat profile (max evts)
   name.Form("Norm_%s_PS_MaxEvts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s Normalized PerfStat Event - %s", namelab.Data(), sellab.Data());
   fNorm_perfstat_evtmax = new TProfile(name, title, ndiv, ns_min, ns_max);
   fNorm_perfstat_evtmax->SetDirectory(fDirProofBench);
   fNorm_perfstat_evtmax->GetYaxis()->SetTitle("Events/sec");
   fNorm_perfstat_evtmax->GetXaxis()->SetTitle(axtitle);
   fNorm_perfstat_evtmax->SetMarkerStyle(23);
   fNorm_perfstat_evtmax->SetMarkerColor(2);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fNorm_perfstat_evtmax);
   fNormLegend->AddEntry(fNorm_perfstat_evtmax, "Maximum");

   // Book queryresult profile
   name.Form("Prof_%s_QR_Evts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s QueryResult Event - %s", namelab.Data(), sellab.Data());
   fProfile_queryresult_event = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_queryresult_event->SetDirectory(fDirProofBench);
   fProfile_queryresult_event->GetYaxis()->SetTitle("Events/sec");
   fProfile_queryresult_event->GetXaxis()->SetTitle(axtitle);
   fProfile_queryresult_event->SetMarkerStyle(22);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_queryresult_event);

   // Book normalized queryresult profile
   name.Form("Norm_%s_QR_Evts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s Normalized QueryResult Event - %s", namelab.Data(), sellab.Data());
   fNorm_queryresult_event = new TProfile(name, title, ndiv, ns_min, ns_max);
   fNorm_queryresult_event->SetDirectory(fDirProofBench);
   fNorm_queryresult_event->GetYaxis()->SetTitle("Events/sec");
   fNorm_queryresult_event->GetXaxis()->SetTitle(axtitle);
   fNorm_queryresult_event->SetMarkerStyle(22);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fNorm_queryresult_event);
   fNormLegend->AddEntry(fNorm_queryresult_event, "Average");

   // Book CPU efficiency profile
   name.Form("Prof_%s_CPU_eff_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s CPU efficiency - %s", namelab.Data(), sellab.Data());
   fProfile_cpu_eff = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_cpu_eff->SetDirectory(fDirProofBench);
   fProfile_cpu_eff->GetYaxis()->SetTitle("Efficiency");
   fProfile_cpu_eff->GetXaxis()->SetTitle(axtitle);
   fProfile_cpu_eff->SetMarkerStyle(22);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_cpu_eff);
}

////////////////////////////////////////////////////////////////////////////////
/// Run benchmark
/// Input parameters
///   nevents:   Number of events to run per file. When it is -1, use data member fNEvents.
///   start: Start scan with 'start' workers. When it is -1, use data member fStart.
///          When 0, the same number of workers are activated on all nodes.
///   stop: Stop scan at 'stop' workers. When it is -1 , use data member fStop.
///   step: Scan every 'step' workers. When it is -1, use data member fStep.
///   ntries: Number of repetitions.  When it is -1, use data member fNTries.
///   debug: debug switch. When it is -1, use data member fDebug.
///   draw: draw switch. When it is -1, use data member fDraw.
/// Returns
///    Nothing

void TProofBenchRunCPU::Run(Long64_t nevents, Int_t start, Int_t stop,
                            Int_t step, Int_t ntries, Int_t debug, Int_t draw)
{
   if (!fProof){
      Error("Run", "Proof not set");
      return;
   }

   nevents = (nevents == -1) ? fNEvents : nevents;
   start = (start == -1) ? fStart : start;
   stop = (stop == -1) ? fStop : stop;
   step = (step == -1) ? fStep : step;
   ntries = (ntries == -1) ? fNTries : ntries;
   debug = (debug == -1) ? fDebug : debug;
   draw = (draw == -1) ? fDraw : draw;

   Bool_t nx = kFALSE;
   if (step == -2){
      nx = kTRUE;
      start = fStart;
      step = 1;
   }

   if (nx){
      Int_t minnworkersanode = fNodes->GetMinWrksPerNode();
      if (stop > minnworkersanode) stop = minnworkersanode;
   }

   // Load the selector, if needed
   if (!TClass::GetClass(fSelName)) {
      // Is it the default selector?
      if (fSelName == kPROOF_BenchSelCPUDef) {
         // Load the parfile
         TString par = TString::Format("%s/%s%s.par", TROOT::GetEtcDir().Data(), kPROOF_BenchParDir, kPROOF_BenchCPUSelPar);
         Info("Run", "Uploading '%s' ...", par.Data());
         if (fProof->UploadPackage(par) != 0) {
            Error("Run", "problems uploading '%s' - cannot continue", par.Data());
            return;
         }
         Info("Run", "Enabling '%s' ...", kPROOF_BenchCPUSelPar);
         if (fProof->EnablePackage(kPROOF_BenchCPUSelPar) != 0) {
            Error("Run", "problems enabling '%s' - cannot continue", kPROOF_BenchCPUSelPar);
            return;
         }
      } else {
         if (fParList.IsNull()) {
            Error("Run", "you should load the class '%s' before running the benchmark", fSelName.Data());
            return;
         } else {
            TString par;
            Int_t from = 0;
            while (fParList.Tokenize(par, from, ",")) {
               Info("Run", "Uploading '%s' ...", par.Data());
               if (fProof->UploadPackage(par) != 0) {
                  Error("Run", "problems uploading '%s' - cannot continue", par.Data());
                  return;
               }
               Info("Run", "Enabling '%s' ...", par.Data());
               if (fProof->EnablePackage(par) != 0) {
                  Error("Run", "problems enabling '%s' - cannot continue", par.Data());
                  return;
               }
            }
         }
      }
      // Check
      if (!TClass::GetClass(fSelName)) {
         Error("Run", "failed to load '%s'", fSelName.Data());
         return;
      }
   }

   // Build histograms, profiles and graphs needed for this run
   BuildHistos(start, stop, step, nx);

   // Get pad
   if (!fCanvas) fCanvas = new TCanvas("Canvas");
   // Cleanup up the canvas
   fCanvas->Clear();

   // Divide the canvas as many as the number of profiles in the list
   fCanvas->Divide(2,1);

   TString perfstats_name = "PROOF_PerfStats";

   SetParameters();

   if (nx){
      Info("Run", "Running CPU-bound tests; %d ~ %d active worker(s)/node,"
                  " every %d worker(s)/node.", start, stop, step);
   } else {
      Info("Run", "Running CPU-bound tests; %d ~ %d active worker(s),"
                  " every %d worker(s).", start, stop, step);
   }

   Int_t npad = 1; //pad number

   Int_t nnodes = fNodes->GetNNodes(); // Number of machines
   Int_t ncores = fNodes->GetNCores(); // Number of cores

   Double_t ymi = -1., ymx = -1., emx = -1.;
   for (Int_t nactive = start; nactive <= stop; nactive += step) {

      // For CPU effectiveness (ok for lite; should do it properly for standard clusters)
      Int_t ncoren = (nactive < ncores) ? nactive : ncores;

      // Actvate the wanted workers
      Int_t nw = -1;
      if (nx) {
         TString workers;
         workers.Form("%dx", nactive);
         nw = fNodes->ActivateWorkers(workers);
      } else {
         nw = fNodes->ActivateWorkers(nactive);
      }
      if (nw < 0){
         Error("Run", "could not activate the requested number of"
                      " workers/node on the cluster; skipping the test point"
                      " (%d workers/node)", nactive);
         continue;
      }

      for (Int_t j = 0; j < ntries; j++) {

         if (nx){
            Info("Run", "Running CPU-bound tests with %d active worker(s)/node;"
                        " trial %d/%d", nactive, j + 1, ntries);
         } else {
            Info("Run", "Running CPU-bound tests with %d active worker(s);"
                        " trial %d/%d", nactive, j + 1, ntries);
         }

         Int_t nevents_all=0;
         if (nx){
            nevents_all=nevents*nactive*nnodes;
         } else {
            nevents_all=nevents*nactive;
         }

         // Process
         fProof->Process(fSelName, nevents_all, fSelOption);

         TList *l = fProof->GetOutputList();

         // Save perfstats
         TTree *t = 0;
         if (l) t = dynamic_cast<TTree*>(l->FindObject(perfstats_name.Data()));
         if (t) {

            //FillPerfStatPerfPlots(t, profile_perfstat_event, nactive);
            FillPerfStatPerfPlots(t, nactive);

            TProofPerfAnalysis pfa(t);
            Double_t pf_eventrate = pfa.GetEvtRateAvgMax();
//            if (pf_eventrate > emx) emx = pf_eventrate;
            fProfile_perfstat_evtmax->Fill(nactive, pf_eventrate);
            fCanvas->cd(npad);
            fProfile_perfstat_evtmax->SetMaximum(emx*1.6);
            fProfile_perfstat_evtmax->SetMinimum(0.);
            fProfile_perfstat_evtmax->Draw("L");
            fProfLegend->Draw();
            gPad->Update();
            // The normalised histos
            // Use the first bin to set the Y range for the histo
            Double_t nert = nx ? pf_eventrate/nactive/nnodes : pf_eventrate/nactive;
            fNorm_perfstat_evtmax->Fill(nactive, nert);
            Double_t y1 = fNorm_perfstat_evtmax->GetBinContent(1);
            Double_t e1 = fNorm_perfstat_evtmax->GetBinError(1);
            Double_t dy = 5 * e1;
            if (dy / y1 < 0.2) dy = y1 * 0.2;
            if (dy > y1) dy = y1*.999999;
            if (ymi < 0.) ymi = y1 - dy;
            if (fNorm_perfstat_evtmax->GetBinContent(nactive) < ymi)
               ymi = fNorm_perfstat_evtmax->GetBinContent(nactive) / 2.;
            if (ymx < 0.) ymx = y1 + dy;
            if (fNorm_perfstat_evtmax->GetBinContent(nactive) > ymx)
               ymx = fNorm_perfstat_evtmax->GetBinContent(nactive) * 1.5;
            fNorm_perfstat_evtmax->SetMaximum(ymx);
            fNorm_perfstat_evtmax->SetMinimum(ymi);
            fCanvas->cd(npad + 1);
            fNorm_perfstat_evtmax->Draw("L");
            fNormLegend->Draw();
            gPad->Update();

            // Build up new name
            TString newname = TString::Format("%s_%s_%dwrks%dthtry", t->GetName(), GetName(), nactive, j);
            t->SetName(newname);

            if (debug && fDirProofBench->IsWritable()){
               TDirectory *curdir = gDirectory;
               TString dirn = nx ? "RunCPUx" : "RunCPU";
               if (!fDirProofBench->GetDirectory(dirn))
                  fDirProofBench->mkdir(dirn, "RunCPU results");
               if (fDirProofBench->cd(dirn)) {
                  t->SetDirectory(fDirProofBench);
                  t->Write();
                  l->Remove(t);
               } else {
                  Warning("Run", "cannot cd to subdirectory '%s' to store the results!", dirn.Data());
               }
               curdir->cd();
            }

         } else {
            if (l)
               Warning("Run", "%s: tree not found", perfstats_name.Data());
            else
               Error("Run", "PROOF output list is empty!");
         }

         // Performance measures from TQueryResult

         const char *drawopt = t ? "LSAME" : "L";
         TQueryResult *queryresult = fProof->GetQueryResult();
         if (queryresult) {
            queryresult->Print("F");
            TDatime qr_start = queryresult->GetStartTime();
            TDatime qr_end = queryresult->GetEndTime();
            Float_t qr_proc = queryresult->GetProcTime();

            Long64_t qr_entries = queryresult->GetEntries();

            // Calculate event rate
            Double_t qr_eventrate = qr_entries / Double_t(qr_proc);
            if (qr_eventrate > emx) emx = qr_eventrate;

            // Calculate and fill CPU efficiency
            Float_t qr_cpu_eff = -1.;
            if (qr_proc > 0.) {
               qr_cpu_eff = queryresult->GetUsedCPU() / ncoren / qr_proc ;
               fProfile_cpu_eff->Fill(nactive, qr_cpu_eff);
               Printf("cpu_eff: %f", qr_cpu_eff);
            }

            // Fill and draw
            fProfile_queryresult_event->Fill(nactive, qr_eventrate);
            fCanvas->cd(npad);
            fProfile_queryresult_event->Draw(drawopt);
            fProfLegend->Draw();
            gPad->Update();
            // The normalised histo
            Double_t nert = nx ? qr_eventrate/nactive/nnodes : qr_eventrate/nactive;
            fNorm_queryresult_event->Fill(nactive, nert);
            // Use the first bin to set the Y range for the histo
            Double_t y1 = fNorm_queryresult_event->GetBinContent(1);
            Double_t e1 = fNorm_queryresult_event->GetBinError(1);
            Double_t dy = 5 * e1;
            if (dy / y1 < 0.2) dy = y1 * 0.2;
            if (dy > y1) dy = y1*.999999;
            if (ymi < 0.) ymi = y1 - dy;
            if (fNorm_queryresult_event->GetBinContent(nactive) < ymi)
               ymi = fNorm_queryresult_event->GetBinContent(nactive) / 2.;
            if (ymx < 0.) ymx = y1 + dy;
            if (fNorm_queryresult_event->GetBinContent(nactive) > ymx)
               ymx = fNorm_queryresult_event->GetBinContent(nactive) * 1.5;
            fNorm_queryresult_event->SetMaximum(ymx);
//            fNorm_queryresult_event->SetMinimum(ymi);
            fNorm_queryresult_event->SetMinimum(0.);
            fCanvas->cd(npad+1);
            fNorm_queryresult_event->Draw(drawopt);
            fNormLegend->Draw();
         } else {
            Warning("Run", "TQueryResult not found!");
         }
         gPad->Update();

      } // for iterations
   } // for number of workers

   // Make the result persistent
   fCanvas->cd(npad);
   fProfile_queryresult_event->SetMaximum(1.6*emx);
   fProfile_queryresult_event->DrawCopy("L");
   fProfile_perfstat_evtmax->DrawCopy("LSAME");
   fProfLegend->Draw();
   fCanvas->cd(npad + 1);
   fNorm_queryresult_event->DrawCopy("L");
   fNorm_perfstat_evtmax->DrawCopy("LSAME");
   fNormLegend->Draw();
   gPad->Update();

   //save performance profiles to file
   if (fDirProofBench && fDirProofBench->IsWritable()){
      TDirectory *curdir = gDirectory;
      TString dirn = nx ? "RunCPUx" : "RunCPU";
      if (!fDirProofBench->GetDirectory(dirn))
         fDirProofBench->mkdir(dirn, "RunCPU results");
      if (fDirProofBench->cd(dirn)) {
         fListPerfPlots->Write(0, kOverwrite);
         fListPerfPlots->SetOwner(kFALSE);
         fListPerfPlots->Clear();
      } else {
         Warning("Run", "cannot cd to subdirectory '%s' to store the results!", dirn.Data());
      }
      curdir->cd();
   }
}

////////////////////////////////////////////////////////////////////////////////

void TProofBenchRunCPU::FillPerfStatPerfPlots(TTree* t, Int_t nactive)
{
   // Fill performance profiles using tree 't'(PROOF_PerfStats).
   // Input parameters
   //    t: Proof output tree (PROOF_PerfStat) containing performance statistics.
   //    profile: Profile to be filled up with information from tree 't'.
   //    nactive: Number of active workers processed the query.
   // Return
   //    Nothing

   // find perfstat profile
   if (!fProfile_perfstat_event){
      Error("FillPerfStatPerfPlots", "no perfstat profile found");
      return;
   }

   // find perfstat histogram
   if (!fHist_perfstat_event){
      Error("FillPerfStatPerfPlots", "no perfstat histogram found");
      return;
   }

   // extract timing information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   t->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = t->GetEntries();

   Double_t event_rate_packet = 0;

   for (Long64_t k=0; k<entries; k++) {

      t->GetEntry(k);

      // Skip information from workers
      if (pe.fEvtNode.Contains(".")) continue;

      if (pe.fType == TVirtualPerfStats::kPacket){
         if (pe.fProcTime != 0.0) {
            event_rate_packet = pe.fEventsProcessed / pe.fProcTime;
            fHist_perfstat_event->Fill(Double_t(nactive), event_rate_packet);
         }
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Show settings

void TProofBenchRunCPU::Print(Option_t* option) const
{
   Printf("+++ TProofBenchRunCPU +++++++++++++++++++++++++++++++++++++++++");
   Printf("Name      = %s", fName.Data());
   if (fProof) fProof->Print(option);
   Printf("fHistType = k%s", GetNameStem().Data());
   Printf("fNHists   = %d", fNHists);
   Printf("fNEvents  = %lld", fNEvents);
   Printf("fNTries   = %d", fNTries);
   Printf("fStart    = %d", fStart);
   Printf("fStop     = %d", fStop);
   Printf("fStep     = %d", fStep);
   Printf("fDraw     = %d", fDraw);
   Printf("fDebug    = %d", fDebug);
   if (fDirProofBench)
      Printf("fDirProofBench = %s", fDirProofBench->GetPath());
   if (fNodes) fNodes->Print(option);
   if (fListPerfPlots) fListPerfPlots->Print(option);
   if (fCanvas)
      Printf("Performance Canvas: Name = %s Title = %s",
              fCanvas->GetName(), fCanvas->GetTitle());
   Printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw Performance plots

void TProofBenchRunCPU::DrawPerfPlots()
{
   // Get canvas
   if (!fCanvas) fCanvas = new TCanvas("Canvas");

   fCanvas->Clear();

   // Divide the canvas as many as the number of profiles in the list
   Int_t nprofiles = fListPerfPlots->GetSize();
   if (nprofiles <= 2){
      fCanvas->Divide(1,nprofiles);
   } else {
      Int_t nside = (Int_t)TMath::Sqrt((Float_t)nprofiles);
      nside = (nside*nside<nprofiles)?nside+1:nside;
      fCanvas->Divide(nside,nside);
   }

   Int_t npad=1;
   TIter nxt(fListPerfPlots);
   TProfile* profile=0;
   while ((profile=(TProfile*)(nxt()))){
      fCanvas->cd(npad++);
      profile->Draw();
      gPad->Update();
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Set histogram type

void TProofBenchRunCPU::SetHistType(TPBHistType *histtype)
{
   fHistType = histtype;
   fName.Form("%sCPU", GetNameStem().Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Get name for this run

TString TProofBenchRunCPU::GetNameStem() const
{
   TString namestem("+++undef+++");
   if (fHistType) {
      switch (fHistType->GetType()) {
         case TPBHistType::kHist1D:
            namestem = "Hist1D";
            break;
         case TPBHistType::kHist2D:
            namestem = "Hist2D";
            break;
         case TPBHistType::kHist3D:
            namestem = "Hist3D";
            break;
         case TPBHistType::kHistAll:
            namestem = "HistAll";
            break;
         default:
            break;
      }
   }
   return namestem;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameters

Int_t TProofBenchRunCPU::SetParameters()
{
   if (!fProof) {
      Error("SetParameters", "proof not set; Doing nothing");
      return 1;
   }

   if (!fHistType) fHistType = new TPBHistType(TPBHistType::kHist1D);
   fProof->AddInput(fHistType);
   fProof->SetParameter("PROOF_BenchmarkNHists", fNHists);
   fProof->SetParameter("PROOF_BenchmarkDraw", Int_t(fDraw));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete parameters set for this run

Int_t TProofBenchRunCPU::DeleteParameters()
{
   if (!fProof){
      Error("DeleteParameters", "proof not set; Doing nothing");
      return 1;
   }
   if (fProof->GetInputList()) {
      TObject *type = fProof->GetInputList()->FindObject("PROOF_Benchmark_HistType");
      if (type) fProof->GetInputList()->Remove(type);
   }
   fProof->DeleteParameters("PROOF_BenchmarkNHists");
   fProof->DeleteParameters("PROOF_BenchmarkDraw");
   return 0;
}

