// @(#)root/proof:$Id$
// Author: Sangsu Ryu 22/06/2010

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProofBenchRunDataRead
\ingroup proofbench

I/O-intensive PROOF benchmark test reads in event files distributed
on the cluster. Number of events processed per second and size of
events processed per second are plotted against number of active
workers. Performance rate for unit packets and performance rate
for query are plotted.

*/

#include "RConfigure.h"

#include "TProofBenchRunDataRead.h"
#include "TProofBenchDataSet.h"
#include "TProofPerfAnalysis.h"
#include "TProofNodes.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TProof.h"
#include "TString.h"
#include "Riostream.h"
#include "TMap.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TPerfStats.h"
#include "THashList.h"
#include "TSortedList.h"
#include "TVirtualPad.h"
#include "TEnv.h"
#include "TLeaf.h"
#include "TQueryResult.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"

ClassImp(TProofBenchRunDataRead);

////////////////////////////////////////////////////////////////////////////////

TProofBenchRunDataRead::TProofBenchRunDataRead(TProofBenchDataSet *pbds, TPBReadType *readtype,
                                               TDirectory* dirproofbench, TProof* proof,
                                               TProofNodes* nodes, Long64_t nevents, Int_t ntries,
                                               Int_t start, Int_t stop, Int_t step, Int_t debug)
                       : TProofBenchRun(proof, kPROOF_BenchSelDataDef), fProof(proof),
                         fReadType(readtype), fDS(pbds),
                         fNEvents(nevents), fNTries(ntries), fStart(start), fStop(stop), fStep(step),
                         fDebug(debug), fFilesPerWrk(2), fReleaseCache(kTRUE),
                         fDirProofBench(dirproofbench), fNodes(nodes),
                         fListPerfPlots(0), fProfile_perfstat_event(0), fHist_perfstat_event(0),
                         fProfile_perfstat_evtmax(0), fNorm_perfstat_evtmax(0),
                         fProfile_queryresult_event(0), fNorm_queryresult_event(0),
                         fProfile_perfstat_IO(0), fHist_perfstat_IO(0),
                         fProfile_perfstat_IOmax(0), fNorm_perfstat_IOmax(0),
                         fProfile_queryresult_IO(0), fNorm_queryresult_IO(0), fProfile_cpu_eff(0),
                         fProfLegend_evt(0), fNormLegend_evt(0), fProfLegend_mb(0), fNormLegend_mb(0),
                         fCPerfProfiles(0), fName(0)
{
   // Default constructor

   if (!fProof) fProof = gProof;
   if (!fDS) fDS = new TProofBenchDataSet(fProof);

   // Set name
   fName = "DataRead";

   if (!fNodes) fNodes = new TProofNodes(fProof);
   fNodes->GetMapOfActiveNodes()->Print();

   if (stop == -1) fStop = fNodes->GetNWorkersCluster();

   fListPerfPlots = new TList;

   gEnv->SetValue("Proof.StatsTrace",1);
   gStyle->SetOptStat(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofBenchRunDataRead::~TProofBenchRunDataRead()
{
   fProof=0;
   fDirProofBench=0;
   SafeDelete(fListPerfPlots);
   if (fCPerfProfiles) delete fCPerfProfiles;
   SafeDelete(fProfLegend_evt);
   SafeDelete(fProfLegend_mb);
   SafeDelete(fNormLegend_evt);
   SafeDelete(fNormLegend_mb);
}

////////////////////////////////////////////////////////////////////////////////
/// Run benchmark
/// Input parameters
///    dset:    Dataset on which to run
///    start: Start scan with 'start' workers.
///    stop: Stop scan at 'stop workers.
///    step: Scan every 'step' workers.
///    ntries: Number of tries. When it is -1, data member fNTries is used.
///    debug: debug switch.
///    Int_t: Ignored
/// Returns
///    Nothing

void TProofBenchRunDataRead::Run(const char *dset, Int_t start, Int_t stop,
                                 Int_t step, Int_t ntries, Int_t debug, Int_t)
{
   if (!fProof){
      Error("Run", "Proof not set");
      return;
   }
   if (!dset || (dset && strlen(dset) <= 0)){
      Error("Run", "dataset name not set");
      return;
   }
   // Check if the dataset exists
   if (!fProof->ExistsDataSet(dset)) {
      Error("Run", "no such data set found; %s", dset);
      return;
   }

   start = (start == -1) ? fStart : start;
   stop = (stop == -1) ? fStop : stop;
   step = (step == -1) ? fStep : step;
   ntries = (ntries == -1) ? fNTries : ntries;
   debug = (debug == -1) ? fDebug : debug;

   Int_t fDebug_sav = fDebug;
   fDebug = debug;

   Bool_t nx = kFALSE;
   if (step == -2){
      nx = kTRUE;
      step = 1;
   }

   if (nx){
      Int_t minnworkersanode = fNodes->GetMinWrksPerNode();
      if (stop > minnworkersanode) stop = minnworkersanode;
   }

   // Load the selector, if needed
   if (!TClass::GetClass(fSelName) || !fDS->IsProof(fProof)) {
      // Is it the default selector?
      if (fSelName == kPROOF_BenchSelDataDef) {
         // Load the parfile
         TString par = TString::Format("%s/%s%s.par", TROOT::GetEtcDir().Data(), kPROOF_BenchParDir, kPROOF_BenchDataSelPar);
         Info("Run", "Uploading '%s' ...", par.Data());
         if (fProof->UploadPackage(par) != 0) {
            Error("Run", "problems uploading '%s' - cannot continue", par.Data());
            return;
         }
         Info("Run", "Enabling '%s' ...", kPROOF_BenchDataSelPar);
         if (fProof->EnablePackage(kPROOF_BenchDataSelPar) != 0) {
            Error("Run", "problems enabling '%s' - cannot continue", kPROOF_BenchDataSelPar);
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

   TString dsname(dset);
   TString dsbasename = gSystem->BaseName(dset);

   // Get pad
   if (!fCPerfProfiles){
      TString canvasname = TString::Format("Performance Profiles %s", GetName());
      fCPerfProfiles = new TCanvas(canvasname.Data(), canvasname.Data());
   }

   // Cleanup up the canvas
   fCPerfProfiles->Clear();

   fCPerfProfiles->Divide(2,2);

   Info("Run", "Running IO-bound tests on dataset '%s'; %d ~ %d active worker(s),"
               " every %d worker(s).", dset, start, stop, step);

   Int_t npad = 1; //pad number

   Int_t nnodes = fNodes->GetNNodes(); // Number of machines
   Int_t ncores = fNodes->GetNCores(); // Number of cores

   Bool_t drawpf = kFALSE;
   Double_t ymi = -1., ymx = -1., emx =- 1, ymiio = -1., ymxio = -1., mbmx = -1.;
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

      // Prepare the dataset for this run. possibly a subsample of
      // the total one
      TFileCollection *fc = GetDataSet(dsname, nactive, nx);
      if (!fc) {
         Error("Run", "could not retrieve dataset '%s'", dsname.Data());
         continue;
      }
      fc->Print("F");
      TString dsn = TString::Format("%s_%d_%d", dsbasename.Data(), nactive, (Int_t)nx);
      fProof->RegisterDataSet(dsn, fc, "OT");
      fProof->ShowDataSet(dsn, "F");

      for (Int_t j=0; j<ntries; j++) {

         if (nx){
            Info("Run", "Running IO-bound tests with %d active worker(s)/node;"
                        " trial %d/%d", nactive, j + 1, ntries);
         } else {
            Info("Run", "Running IO-bound tests with %d active worker(s);"
                        " trial %d/%d", nactive, j + 1, ntries);
         }

         // Cleanup run
         const char *dsnr = (fDS->IsProof(fProof)) ? dsn.Data() : dsname.Data();
         if (fReleaseCache) fDS->ReleaseCache(dsnr);

         DeleteParameters();
         SetParameters();

         Info("Run", "Processing data set %s with"
                     " %d active worker(s).", dsn.Data(), nactive);

         TTime starttime = gSystem->Now();
         fProof->Process(dsn, fSelName, fSelOption);

         DeleteParameters();

         TTime endtime = gSystem->Now();

         TList *l = fProof->GetOutputList();

         //save perfstats
         TString perfstats_name = "PROOF_PerfStats";
         TTree *t = 0;
         if (l) t = dynamic_cast<TTree*>(l->FindObject(perfstats_name.Data()));
         if (t) {
            drawpf = kTRUE;
            TTree* tnew=(TTree*)t->Clone("tnew");

            FillPerfStatProfiles(tnew, nactive);

            TProofPerfAnalysis pfa(tnew);
            Double_t pf_eventrate = pfa.GetEvtRateAvgMax();
            Double_t pf_IOrate = pfa.GetMBRateAvgMax();
            fProfile_perfstat_evtmax->Fill(nactive, pf_eventrate);
            fCPerfProfiles->cd(npad);
            fProfile_perfstat_evtmax->SetMaximum(1.6*emx);
            fProfile_perfstat_evtmax->SetMinimum(0.);
            fProfile_perfstat_evtmax->Draw();
            fProfLegend_evt->Draw();
            gPad->Update();
            fProfile_perfstat_IOmax->Fill(nactive, pf_IOrate);
            fCPerfProfiles->cd(npad + 2);
            fProfile_perfstat_IOmax->SetMaximum(1.6*mbmx);
            fProfile_perfstat_IOmax->SetMinimum(0.);
            fProfile_perfstat_IOmax->Draw();
            fProfLegend_mb->Draw();
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
            fCPerfProfiles->cd(npad + 1);
            fNorm_perfstat_evtmax->Draw();
            gPad->Update();
            //
            Double_t niort = nx ? pf_IOrate/nactive/nnodes : pf_IOrate/nactive;
            fNorm_perfstat_IOmax->Fill(nactive, niort);
            y1 = fNorm_perfstat_IOmax->GetBinContent(1);
            e1 = fNorm_perfstat_IOmax->GetBinError(1);
            dy = 5 * e1;
            if (dy / y1 < 0.2) dy = y1 * 0.2;
            if (dy > y1) dy = y1*.999999;
            if (ymiio < 0.) ymiio = y1 - dy;
            if (fNorm_perfstat_IOmax->GetBinContent(nactive) < ymiio)
               ymiio = fNorm_perfstat_IOmax->GetBinContent(nactive) / 2.;
            if (ymxio < 0.) ymxio = y1 + dy;
            if (fNorm_perfstat_IOmax->GetBinContent(nactive) > ymxio)
               ymxio = fNorm_perfstat_IOmax->GetBinContent(nactive) * 1.5;
            fNorm_perfstat_IOmax->SetMaximum(ymxio);
            fNorm_perfstat_IOmax->SetMinimum(ymiio);
            fCPerfProfiles->cd(npad + 3);
            fNorm_perfstat_IOmax->Draw();
            gPad->Update();

            //change the name
            TString newname = TString::Format("%s_%s_%dwrks%dthtry", t->GetName(), GetName(), nactive, j);
            tnew->SetName(newname);

            if (debug && fDirProofBench->IsWritable()){
               TDirectory *curdir = gDirectory;
               TString dirn = nx ? "RunDataReadx" : "RunDataRead";
               if (!fDirProofBench->GetDirectory(dirn))
                  fDirProofBench->mkdir(dirn, "RunDataRead results");
               if (fDirProofBench->cd(dirn)) {
                  tnew->SetDirectory(fDirProofBench);
                  tnew->Write();
                  l->Remove(tnew);
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

         //
         const char *drawopt = t ? "SAME" : "";
         // Performance measures from TQueryResult
         TQueryResult *queryresult = fProof->GetQueryResult();
         if (queryresult) {
            TDatime qr_start = queryresult->GetStartTime();
            TDatime qr_end = queryresult->GetEndTime();
            Float_t qr_proc = queryresult->GetProcTime();
            Long64_t qr_bytes = queryresult->GetBytes();

            Long64_t qr_entries = queryresult->GetEntries();

            // Calculate and fill CPU efficiency
            Float_t qr_cpu_eff = -1.;
            if (qr_proc > 0.) {
               qr_cpu_eff = queryresult->GetUsedCPU() / ncoren / qr_proc ;
               fProfile_cpu_eff->Fill(nactive, qr_cpu_eff);
               Printf("cpu_eff: %f", qr_cpu_eff);
            }

            // Calculate event rate, fill and draw
            Double_t qr_eventrate=0;

            qr_eventrate = qr_entries / Double_t(qr_proc);
            if (qr_eventrate > emx) emx = qr_eventrate;

            fProfile_queryresult_event->Fill(nactive, qr_eventrate);
            fCPerfProfiles->cd(npad);
            fProfile_queryresult_event->SetMinimum(0.);
            fProfile_queryresult_event->Draw(drawopt);
            fProfLegend_evt->Draw();
            gPad->Update();

            // Calculate IO rate, fill and draw
            Double_t qr_IOrate = 0;

            const Double_t Dmegabytes = 1024*1024;

            qr_IOrate = qr_bytes / Dmegabytes / Double_t(qr_proc);
            if (qr_IOrate > mbmx) mbmx = qr_IOrate;

            fProfile_queryresult_IO->Fill(nactive, qr_IOrate);
            fCPerfProfiles->cd(npad + 2);
            fProfile_queryresult_IO->SetMinimum(0.);
            fProfile_queryresult_IO->Draw(drawopt);
            fProfLegend_mb->Draw();
            gPad->Update();

            // The normalised histos
            // Use the first bin to set the Y range for the histo
            Double_t nert = nx ? qr_eventrate/nactive/nnodes : qr_eventrate/nactive;
            fNorm_queryresult_event->Fill(nactive, nert);
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
//            fNorm_queryresult_event->SetMaximum(ymx);
            fNorm_queryresult_event->SetMinimum(ymi);
            fCPerfProfiles->cd(npad + 1);
            fNorm_queryresult_event->Draw(drawopt);
            fNormLegend_evt->Draw();
            gPad->Update();
            //
            Double_t niort = nx ? qr_IOrate/nactive/nnodes : qr_IOrate/nactive;
            fNorm_queryresult_IO->Fill(nactive, niort);
            y1 = fNorm_queryresult_IO->GetBinContent(1);
            e1 = fNorm_queryresult_IO->GetBinError(1);
            dy = 5 * e1;
            if (dy / y1 < 0.2) dy = y1 * 0.2;
            if (dy > y1) dy = y1*.999999;
            if (ymiio < 0.) ymiio = y1 - dy;
            if (fNorm_queryresult_IO->GetBinContent(nactive) < ymiio)
               ymiio = fNorm_queryresult_IO->GetBinContent(nactive) / 2.;
            if (ymxio < 0.) ymxio = y1 + dy;
            if (fNorm_queryresult_IO->GetBinContent(nactive) > ymxio)
               ymxio = fNorm_queryresult_IO->GetBinContent(nactive) * 1.5;
//            fNorm_queryresult_IO->SetMaximum(ymxio);
            fNorm_queryresult_IO->SetMinimum(ymiio);
            fCPerfProfiles->cd(npad + 3);
            fNorm_queryresult_IO->Draw(drawopt);
            fNormLegend_mb->Draw();
            gPad->Update();
         }
         fCPerfProfiles->cd(0);
      }
      // Remove temporary dataset
      fProof->RemoveDataSet(dsn);
      SafeDelete(fc);
   }

   // Make the result persistent
   fCPerfProfiles->cd(npad);
   fProfile_queryresult_event->SetMaximum(1.6*emx);
   fProfile_queryresult_event->DrawCopy();
   if (drawpf) fProfile_perfstat_evtmax->DrawCopy("SAME");
   fProfLegend_evt->Draw();
   fCPerfProfiles->cd(npad + 2);
   fProfile_queryresult_IO->SetMaximum(1.6*mbmx);
   fProfile_queryresult_IO->DrawCopy();
   if (drawpf) fProfile_perfstat_IOmax->DrawCopy("SAME");
   fProfLegend_mb->Draw();
   fCPerfProfiles->cd(npad + 1);
   fNorm_queryresult_event->DrawCopy();
   if (drawpf) fNorm_perfstat_evtmax->DrawCopy("SAME");
   fNormLegend_evt->Draw();
   fCPerfProfiles->cd(npad + 3);
   fNorm_queryresult_IO->DrawCopy();
   if (drawpf) fNorm_perfstat_IOmax->DrawCopy("SAME");
   fProfLegend_mb->Draw();
   gPad->Update();

   //save performance profiles to file
   if (fDirProofBench->IsWritable()){
      TDirectory *curdir = gDirectory;
      TString dirn = nx ? "RunDataReadx" : "RunDataRead";
      if (!fDirProofBench->GetDirectory(dirn))
         fDirProofBench->mkdir(dirn, "RunDataRead results");
      if (fDirProofBench->cd(dirn)) {
         fListPerfPlots->Write(0, kOverwrite);
         fListPerfPlots->SetOwner(kFALSE);
         fListPerfPlots->Clear();
      } else {
         Warning("Run", "cannot cd to subdirectory '%s' to store the results!", dirn.Data());
      }
      curdir->cd();
   }
   // Restore member data
   fDebug = fDebug_sav;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a subsample of dsname suited to run with 'nact' and option 'nx'.

TFileCollection *TProofBenchRunDataRead::GetDataSet(const char *dset,
                                                    Int_t nact, Bool_t nx)
{
   TFileCollection *fcsub = 0;

   // Dataset must exists
   if (!fProof || (fProof && !fProof->ExistsDataSet(dset))) {
      Error("GetDataSet", "dataset '%s' does not exist", dset);
      return fcsub;
   }

   // Get the full collection
   TFileCollection *fcref = fProof->GetDataSet(dset);
   if (!fcref) {
      Error("GetDataSet", "dataset '%s' could not be retrieved", dset);
      return fcsub;
   }
   // Is it remote ?
   Bool_t remote = (fcref->TestBit(TFileCollection::kRemoteCollection)) ? kTRUE : kFALSE;

   // Separate info per server
   TMap *mpref = fcref->GetFilesPerServer(fProof->GetMaster(), kTRUE);
   if (!mpref) {
      SafeDelete(fcref);
      Error("GetDataSet", "problems classifying info on per-server base");
      return fcsub;
   }
   mpref->Print();

   // Get Active node information
   TMap *mpnodes = fNodes->GetMapOfActiveNodes();
   if (!mpnodes) {
      SafeDelete(fcref);
      SafeDelete(mpref);
      Error("GetDataSet", "problems getting map of active nodes");
      return fcsub;
   }
   mpnodes->Print();

   // Number of files: fFilesPerWrk per active worker
   Int_t nf = fNodes->GetNActives() * fFilesPerWrk;
   Printf(" number of files needed (ideally): %d (%d per worker)", nf, fFilesPerWrk);

   // The output dataset
   fcsub = new TFileCollection(TString::Format("%s_%d_%d", fcref->GetName(), nact, nx),
                                                           fcref->GetTitle());

   // Order reference sub-collections
   TIter nxnd(mpnodes);
   TObject *key = 0;
   TFileInfo *fi = 0;
   TFileCollection *xfc = 0;
   TList *lswrks = 0;
   while ((key = nxnd())) {
      TIter nxsrv(mpref);
      TObject *ksrv = 0;
      while ((ksrv = nxsrv())) {
         TUrl urlsrv(ksrv->GetName());
         if (TString(urlsrv.GetHostFQDN()).IsNull())
            urlsrv.SetHost(TUrl(gProof->GetMaster()).GetHostFQDN());
         if (remote ||
             !strcmp(urlsrv.GetHostFQDN(), TUrl(key->GetName()).GetHostFQDN())) {
            if ((xfc = dynamic_cast<TFileCollection *>(mpref->GetValue(ksrv)))) {
               if ((lswrks = dynamic_cast<TList *>(mpnodes->GetValue(key)))) {
                  Int_t nfnd = fFilesPerWrk * lswrks->GetSize();
                  while (nfnd-- && xfc->GetList()->GetSize() > 0) {
                     if ((fi = (TFileInfo *) xfc->GetList()->First())) {
                        xfc->GetList()->Remove(fi);
                        fcsub->Add(fi);
                     }
                  }
               } else {
                  Warning("GetDataSet", "could not attach to worker list for node '%s'",
                                        key->GetName());
               }
            } else {
               Warning("GetDataSet", "could not attach to file collection for server '%s'",
                                     ksrv->GetName());
            }
         }
      }
   }

   // Update counters
   fcsub->Update();
   fcsub->Print();

   // Make sure that the tree name is the one of the original dataset
   if (fcref) {
      TString dflt(fcref->GetDefaultTreeName());
      if (!dflt.IsNull()) fcsub->SetDefaultTreeName(dflt);
   }

   // Cleanup
   SafeDelete(fcref);
   SafeDelete(mpref);
   // Done
   return fcsub;
}

////////////////////////////////////////////////////////////////////////////////

void TProofBenchRunDataRead::FillPerfStatProfiles(TTree *t, Int_t nactive)
{
   // Fill performance profiles using tree 't'(PROOF_PerfStats).
   // Input parameters
   //    t: Proof output tree (PROOF_PerfStat) containing performance
   //       statistics.
   //    nactive: Number of active workers processed the query.
   // Return
   //    Nothing

   // extract timing information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   t->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = t->GetEntries();

   const Double_t Dmegabytes = 1024.*1024.;
   Double_t event_rate_packet = 0;
   Double_t IO_rate_packet = 0;

   for (Long64_t k=0; k<entries; k++) {
      t->GetEntry(k);

      // Skip information from workers
      if (pe.fEvtNode.Contains(".")) continue;

      if (pe.fType==TVirtualPerfStats::kPacket){
         if (pe.fProcTime != 0.0) {
            event_rate_packet = pe.fEventsProcessed / pe.fProcTime;
            fHist_perfstat_event->Fill(Double_t(nactive), event_rate_packet);
            IO_rate_packet = pe.fBytesRead / Dmegabytes / pe.fProcTime;
            fHist_perfstat_IO->Fill(Double_t(nactive), IO_rate_packet);
         }
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the content of this object

void TProofBenchRunDataRead::Print(Option_t* option) const
{
   Printf("Name         = %s", fName.Data());
   if (fProof) fProof->Print(option);
   Printf("fReadType    = %s%s", "k", GetNameStem().Data());
   Printf("fNEvents     = %lld", fNEvents);
   Printf("fNTries      = %d", fNTries);
   Printf("fStart       = %d", fStart);
   Printf("fStop        = %d", fStop);
   Printf("fStep        = %d", fStep);
   Printf("fDebug       = %d", fDebug);
   if (fDirProofBench)
      Printf("fDirProofBench = %s", fDirProofBench->GetPath());
   if (fNodes) fNodes->Print(option);
   if (fListPerfPlots) fListPerfPlots->Print(option);

   if (fCPerfProfiles)
      Printf("Performance Profiles Canvas: Name = %s Title = %s",
              fCPerfProfiles->GetName(), fCPerfProfiles->GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Get canvas

void TProofBenchRunDataRead::DrawPerfProfiles()
{
   if (!fCPerfProfiles){
      TString canvasname = TString::Format("Performance Profiles %s", GetName());
      fCPerfProfiles = new TCanvas(canvasname.Data(), canvasname.Data());
   }

   fCPerfProfiles->Clear();

   // Divide the canvas as many as the number of profiles in the list
   Int_t nprofiles = fListPerfPlots->GetSize();
   if (nprofiles <= 2){
      fCPerfProfiles->Divide(nprofiles);
   } else {
      Int_t nside = (Int_t)TMath::Sqrt((Float_t)nprofiles);
      nside = (nside*nside < nprofiles) ? nside + 1 : nside;
      fCPerfProfiles->Divide(nside,nside);
   }

   Int_t npad=1;
   TIter nxt(fListPerfPlots);
   TProfile* profile=0;
   while ((profile=(TProfile*)(nxt()))){
      fCPerfProfiles->cd(npad++);
      profile->Draw();
      gPad->Update();
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Get name for this run

TString TProofBenchRunDataRead::GetNameStem()const
{
   TString namestem("+++undef+++");
   if (fReadType) {
      switch (fReadType->GetType()) {
         case TPBReadType::kReadFull:
            namestem="Full";
            break;
         case TPBReadType::kReadOpt:
            namestem="Opt";
            break;
         case TPBReadType::kReadNo:
            namestem="No";
            break;
         default:
            break;
      }
   }
   return namestem;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameters

Int_t TProofBenchRunDataRead::SetParameters()
{
   if (!fProof){
      Error("SetParameters", "Proof not set; Doing nothing");
      return 1;
   }
   if (!fReadType) fReadType = new TPBReadType(TPBReadType::kReadOpt);
   fProof->AddInput(fReadType);
   fProof->SetParameter("PROOF_BenchmarkDebug", Int_t(fDebug));
   // For Mac Os X only: do not OS cache the files read
   fProof->SetParameter("PROOF_DontCacheFiles", Int_t(1));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete parameters set for this run

Int_t TProofBenchRunDataRead::DeleteParameters()
{
   if (!fProof){
      Error("DeleteParameters", "Proof not set; Doing nothing");
      return 1;
   }
   if (fProof->GetInputList()) {
      TObject *type = fProof->GetInputList()->FindObject("PROOF_Benchmark_ReadType");
      if (type) fProof->GetInputList()->Remove(type);
   }
   fProof->DeleteParameters("PROOF_BenchmarkDebug");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Build histograms, profiles and graphs needed for this run

void TProofBenchRunDataRead::BuildHistos(Int_t start, Int_t stop, Int_t step, Bool_t nx)
{
   TObject *o = 0;
   Int_t quotient = (stop - start) / step;
   Int_t ndiv = quotient + 1;
   Double_t ns_min = start - step/2.;
   Double_t ns_max = quotient*step + start + step/2.;

   fProfLegend_evt = new TLegend(0.1, 0.8, 0.3, 0.9);
   fNormLegend_evt = new TLegend(0.7, 0.8, 0.9, 0.9);
   fProfLegend_mb = new TLegend(0.1, 0.8, 0.3, 0.9);
   fNormLegend_mb = new TLegend(0.7, 0.8, 0.9, 0.9);

   TString axtitle("Active Workers"), namelab(GetName()), sellab(GetSelName());
   if (nx) {
      axtitle = "Active Workers/Node";
      namelab.Form("x_%s", GetName());
   }
   if (fSelName == kPROOF_BenchSelDataDef)
      sellab.Form("%s_%s", GetSelName(), GetNameStem().Data());

   TString name, title;

   // Book perfstat profile (max evts)
   name.Form("Prof_%s_PS_MaxEvts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s PerfStat Event- %s", namelab.Data(), sellab.Data());
   fProfile_perfstat_evtmax = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_perfstat_evtmax->SetDirectory(fDirProofBench);
   fProfile_perfstat_evtmax->GetYaxis()->SetTitle("Events/sec");
   fProfile_perfstat_evtmax->GetXaxis()->SetTitle(axtitle);
   fProfile_perfstat_evtmax->SetMarkerStyle(23);
   fProfile_perfstat_evtmax->SetMarkerStyle(2);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_perfstat_evtmax);
   fProfLegend_evt->AddEntry(fProfile_perfstat_evtmax, "Maximum");

   // Book perfstat profile (evts)
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

   // Book perfstat histogram (evts)
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
   name.Form("Norm_%s_PF_MaxEvts_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s Normalized QueryResult Event - %s", namelab.Data(), sellab.Data());
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
   fNormLegend_evt->AddEntry(fNorm_perfstat_evtmax, "Maximum");

   // Book queryresult profile (evts)
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
   fProfLegend_evt->AddEntry(fProfile_queryresult_event, "Average");

   // Book normalized queryresult profile (evts)
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
   fNormLegend_evt->AddEntry(fNorm_queryresult_event, "Average");

   // Book perfstat profile (mbs)
   name.Form("Prof_%s_PS_IO_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s PerfStat I/O %s", namelab.Data(), sellab.Data());
   fProfile_perfstat_IO = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_perfstat_IO->SetDirectory(fDirProofBench);
   fProfile_perfstat_IO->GetYaxis()->SetTitle("MB/sec");
   fProfile_perfstat_IO->GetXaxis()->SetTitle(axtitle);
   fProfile_perfstat_IO->SetMarkerStyle(21);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_perfstat_IO);

   // Book perfstat histogram (mbs)
   name.Form("Hist_%s_PS_IO_%s", namelab.Data(), sellab.Data());
   title.Form("Histogram %s PerfStat I/O - %s", namelab.Data(), sellab.Data());
   fHist_perfstat_IO = new TH2D(name, title, ndiv, ns_min, ns_max, 100, 0, 0);
   fHist_perfstat_IO->SetDirectory(fDirProofBench);
   fHist_perfstat_IO->GetYaxis()->SetTitle("MB/sec");
   fHist_perfstat_IO->GetXaxis()->SetTitle(axtitle);
   fHist_perfstat_IO->SetMarkerStyle(7);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fHist_perfstat_IO);

   // Book perfstat profile (max mbs)
   name.Form("Prof_%s_PS_MaxIO_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s PerfStat I/O - %s", namelab.Data(), sellab.Data());
   fProfile_perfstat_IOmax = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_perfstat_IOmax->SetDirectory(fDirProofBench);
   fProfile_perfstat_IOmax->GetYaxis()->SetTitle("MB/sec");
   fProfile_perfstat_IOmax->GetXaxis()->SetTitle(axtitle);
   fProfile_perfstat_IOmax->SetMarkerStyle(21);
   fProfile_perfstat_IOmax->SetMarkerColor(2);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_perfstat_IOmax);
   fProfLegend_mb->AddEntry(fProfile_perfstat_IOmax, "Maximum");

   // Book normalized perfstat profile (max mbs)
   name.Form("Norm_%s_PS_MaxIO_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s Normalized PerfStat I/O - %s", namelab.Data(), sellab.Data());
   fNorm_perfstat_IOmax = new TProfile(name, title, ndiv, ns_min, ns_max);
   fNorm_perfstat_IOmax->SetDirectory(fDirProofBench);
   fNorm_perfstat_IOmax->GetYaxis()->SetTitle("MB/sec");
   fNorm_perfstat_IOmax->GetXaxis()->SetTitle(axtitle);
   fNorm_perfstat_IOmax->SetMarkerStyle(23);
   fNorm_perfstat_IOmax->SetMarkerColor(2);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fNorm_perfstat_IOmax);
   fNormLegend_mb->AddEntry(fNorm_perfstat_IOmax, "Maximum");

   // Book queryresult profile (mbs)
   name.Form("Prof_%s_QR_IO_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s QueryResult I/O - %s", namelab.Data(), sellab.Data());
   fProfile_queryresult_IO = new TProfile(name, title, ndiv, ns_min, ns_max);
   fProfile_queryresult_IO->SetDirectory(fDirProofBench);
   fProfile_queryresult_IO->GetYaxis()->SetTitle("MB/sec");
   fProfile_queryresult_IO->GetXaxis()->SetTitle(axtitle);
   fProfile_queryresult_IO->SetMarkerStyle(22);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fProfile_queryresult_IO);
   fProfLegend_mb->AddEntry(fProfile_queryresult_IO, "Average");

   // Book normalized queryresult profile (mbs)
   name.Form("Norm_%s_QR_IO_%s", namelab.Data(), sellab.Data());
   title.Form("Profile %s Normalized QueryResult I/O - %s", namelab.Data(), sellab.Data());
   fNorm_queryresult_IO = new TProfile(name, title, ndiv, ns_min, ns_max);
   fNorm_queryresult_IO->SetDirectory(fDirProofBench);
   fNorm_queryresult_IO->GetYaxis()->SetTitle("MB/sec");
   fNorm_queryresult_IO->GetXaxis()->SetTitle(axtitle);
   fNorm_queryresult_IO->SetMarkerStyle(22);
   if ((o = fListPerfPlots->FindObject(name))) {
      fListPerfPlots->Remove(o);
      delete o;
   }
   fListPerfPlots->Add(fNorm_queryresult_IO);
   fNormLegend_mb->AddEntry(fNorm_queryresult_IO, "Average");

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
