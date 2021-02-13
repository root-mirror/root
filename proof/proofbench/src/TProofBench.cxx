// @(#)root/proof:$Id$
// Author: G.Ganis, S.Ryu Feb 2011

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**
  \defgroup proofbench PROOF benchmark utilities
  \ingroup proof

  Set of utilities to benchmark a PROOF facility.
  See also https://root.cern.ch/proof-benchmark-framework-tproofbench .

*/

/** \class TProofBench
\ingroup proofbench

  Steering class for PROOF benchmarks

*/

#include "RConfigure.h"

#include "TProofBench.h"
#include "Getline.h"
#include "TProofBenchRunCPU.h"
#include "TProofBenchRunDataRead.h"
#include "TProofBenchDataSet.h"
#include "TProofNodes.h"
#include "TClass.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "THashList.h"
#include "TKey.h"
#include "TMap.h"
#include "TObjString.h"
#include "TProof.h"
#include "TROOT.h"
#include "TSortedList.h"
#include "TTimeStamp.h"
#include "TUrl.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
#ifdef WIN32
#include <io.h>
#endif

ClassImp(TProofBench);

// Functions for fitting

TF1 *TProofBench::fgFp1 = 0;
TF1 *TProofBench::fgFp1n = 0;
TF1 *TProofBench::fgFp2 = 0;
TF1 *TProofBench::fgFp2n = 0;
TF1 *TProofBench::fgFp3 = 0;
TF1 *TProofBench::fgFp3n = 0;
TF1 *TProofBench::fgFio = 0;
TF1 *TProofBench::fgFioV = 0;
static Int_t gFioVn0 = -1;             // Number of real cores for fgFioV
static Int_t gFioVn1 = -1;             // Number of real+hyper cores for fgFioV

TList *TProofBench::fgGraphs = new TList;

////////////////////////////////////////////////////////////////////////////////
/// Simple polynomial 1st degree

Double_t funp1(Double_t *xx, Double_t *par)
{
   Double_t res = par[0] + par[1] * xx[0];
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Simple polynomial 2nd degree

Double_t funp2(Double_t *xx, Double_t *par)
{
   Double_t res = par[0] + par[1] * xx[0] + par[2] * xx[0] * xx[0];
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Normalized 1st degree

Double_t funp1n(Double_t *xx, Double_t *par)
{
   Double_t res = par[0] / xx[0] + par[1];
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Normalized 2nd degree

Double_t funp2n(Double_t *xx, Double_t *par)
{
   Double_t res = par[0] / xx[0] + par[1] + par[2] * xx[0];
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// I/O saturated rate function

Double_t funio(Double_t *xx, Double_t *par)
{
   Double_t sat = par[0] / par[1] * (xx[0] * par[1] / par[2] - 1.);
   if (xx[0] < par[2] / par[1]) sat = 0.;
   Double_t res = par[0] * xx[0] / (1. + sat);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// I/O saturated rate function with varying Rcpu

Double_t funiov(Double_t *xx, Double_t *par)
{
   // par[0] = rio
   // par[1] = b1
   // par[2] = b2
   // par[3] = nc
   // par[4] = ri

   Double_t rio = par[0] / par[3] * xx[0];
   if (xx[0] > par[3]) rio = par[0];

   Double_t rcpu = par[1] * xx[0];
   if (xx[0] > gFioVn0) rcpu = par[1]*gFioVn0 + par[2]*(xx[0] - gFioVn0);
   if (xx[0] > gFioVn1) rcpu = par[1]*gFioVn0 + par[2]*(gFioVn1 - gFioVn0);

   Double_t res = 1. / (1./par[4] + 1./rio + 1./rcpu);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Function with varying Rcpu

Double_t funcpuv(Double_t *xx, Double_t *par)
{
   // par[0] = offset
   // par[1] = rate contribution from real cores
   // par[2] = rate contribution from hyper cores

   Double_t n = (xx[0] - par[0]);
   Double_t rcpu = par[1] * n;
   if (xx[0] > gFioVn0) rcpu = par[1]*gFioVn0 + par[2]*(n - gFioVn0);
   if (xx[0] > gFioVn1) rcpu = par[1]*gFioVn0 + par[2]*(gFioVn1 - gFioVn0);

   return rcpu;
}

////////////////////////////////////////////////////////////////////////////////
/// Function with varying Rcpu normalized

Double_t funcpuvn(Double_t *xx, Double_t *par)
{
   // par[0] = offset
   // par[1] = rate contribution from real cores
   // par[2] = rate contribution from hyper cores

   Double_t n = (xx[0] - par[0]);
   Double_t rcpu = par[1] * n;
   if (xx[0] > gFioVn0) rcpu = par[1]*gFioVn0 + par[2]*(n - gFioVn0);
   if (xx[0] > gFioVn1) rcpu = par[1]*gFioVn0 + par[2]*(gFioVn1 - gFioVn0);

   return rcpu / xx[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor: check PROOF and load selectors PAR

TProofBench::TProofBench(const char *url, const char *outfile, const char *proofopt)
            : fUnlinkOutfile(kFALSE), fProofDS(0), fOutFile(0),
              fNtries(4), fHistType(0), fNHist(16), fReadType(0),
              fDataSet("BenchDataSet"), fNFilesWrk(2), fReleaseCache(kTRUE),
              fDataGenSel(kPROOF_BenchSelDataGenDef),
              fRunCPU(0), fRunDS(0), fDS(0), fDebug(kFALSE), fDescription(0)
{
   SetBit(kInvalidObject);
   if (!url) {
      Error("TProofBench", "specifying a PROOF master url is mandatory - cannot continue");
      return;
   }
   if (!(fProof = TProof::Open(url, proofopt)) || (fProof && !fProof->IsValid())) {
      Error("TProofBench", "could not open a valid PROOF session - cannot continue");
      return;
   }
   // Get the size of the cluster
   fNumWrkMax = fProof->GetParallel();
   if (fProof->UseDynamicStartup() && TProof::GetEnvVars()) {
      // It must be passed as PROOF option 'workers=N' and recorded in the envs vars
      TNamed *n = (TNamed *) TProof::GetEnvVars()->FindObject("PROOF_NWORKERS");
      if (!n) {
         Error("TProofBench", "dynamic mode: you must specify the max number of workers");
         fProof->Close();
         SafeDelete(fProof);
         return;
      }
      TString sn(n->GetTitle());
      if (sn.IsDigit()) fNumWrkMax = sn.Atoi();
      if (!sn.IsDigit()) {
         Error("TProofBench", "dynamic mode: wrong specification of the max number of"
                              " workers ('%s')", n->GetTitle());
         fProof->Close();
         SafeDelete(fProof);
         return;
      }
   }
   if (fNumWrkMax <= 0) {
      Error("TProofBench", "wrong max number of workers ('%d')", fNumWrkMax);
      fProof->Close();
      SafeDelete(fProof);
      return;
   }
   // By default we use the same instance for dataset actions
   fProofDS = fProof;
   // The object is now valid
   ResetBit(kInvalidObject);
   // Identifying string
   TUrl u(url);
   TString host(TString::Format("PROOF at %s", u.GetHost()));
   if (!strcmp(u.GetProtocol(), "lite")) host.Form("PROOF-Lite on %s", gSystem->HostName());
   fDescription = new TNamed("PB_description",
                             TString::Format("%s, %d workers", host.Data(), fNumWrkMax).Data());
   Printf(" Run description: %s", fDescription->GetTitle());
   // Set output file
   if (SetOutFile(outfile, kFALSE) != 0)
      Warning("TProofBench", "problems opening '%s' - ignoring: use SetOutFile to try"
                                   " again or with another file", outfile);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TProofBench::~TProofBench()
{
   CloseOutFile();
   if (fUnlinkOutfile) gSystem->Unlink(fOutFileName);
   SafeDelete(fReadType);
   SafeDelete(fRunCPU);
   SafeDelete(fRunDS);
   SafeDelete(fDescription);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the otuput file
/// Return 0 on success, -1 on error

Int_t TProofBench::OpenOutFile(Bool_t wrt, Bool_t verbose)
{
   // Remove any bad file
   if (fOutFile && fOutFile->IsZombie()) SafeDelete(fOutFile);

   Int_t rc = 0;
   if (!fOutFile && fOutFileName.Length() > 0) {
      const char *mode = 0;
      if (wrt)
         mode = gSystem->AccessPathName(fOutFileName) ? "RECREATE" : "UPDATE";
      else
         mode = "READ";
      if (!(fOutFile = TFile::Open(fOutFileName, mode)) || (fOutFile && fOutFile->IsZombie())) {
         if (verbose)
            Warning("OpenOutFile", "problems opening '%s' - ignoring: use SetOutFile to try"
                                   " again or with another file", fOutFileName.Data());
         rc = -1;
      }
      if (fOutFile) {
         gROOT->GetListOfFiles()->Remove(fOutFile);
         if (!strcmp(mode, "RECREATE")) {
            // Save the description string
            fOutFile->cd();
            fDescription->Write();
         }
      }
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the output file
/// Return 0 on success, -1 on error

Int_t TProofBench::SetOutFile(const char *outfile, Bool_t verbose)
{
   Int_t rc = 0;
   // Close existing file, if any
   if (fOutFile) {
      if (!fOutFile->IsZombie()) fOutFile->Close();
      SafeDelete(fOutFile);
   }

   fOutFileName = outfile;
   if (fOutFileName == "<default>") {
      // Default output file: proofbench-<master>-<DayMonthYear-hhmm>.root
      TDatime dat;
      const char *lite = (fProof->IsLite()) ? "-lite" : "";
      fOutFileName.Form("proofbench-%s%s-%dw-%d-%.2d%.2d.root",
                        fProof->GetMaster(), lite, fNumWrkMax,
                        dat.GetDate(), dat.GetHour(), dat.GetMinute());
      Info("SetOutFile", "using default output file: '%s'", fOutFileName.Data());
      fUnlinkOutfile = kTRUE;
   }
   if (!fOutFileName.IsNull()) {
      if ((rc = OpenOutFile(kTRUE, kFALSE)) != 0 && verbose)
         Warning("SetOutFile", "problems opening '%s' - ignoring: use SetOutFile to try"
                               " again or with another file", outfile);
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Close output file

void TProofBench::CloseOutFile()
{
   if (SetOutFile(0) != 0)
      Warning("CloseOutFile", "problems closing '%s'", fOutFileName.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the CPU run
/// Return 0 on success, -1 on error

Int_t TProofBench::RunCPU(Long64_t nevents, Int_t start, Int_t stop, Int_t step)
{
   // Open the file for the results
   if (OpenOutFile(kTRUE) != 0) {
      Error("RunCPU", "problems opening '%s' to save the result", fOutFileName.Data());
      return -1;
   }
   fUnlinkOutfile = kFALSE;

   SafeDelete(fRunCPU);
   TPBHistType *htype = new TPBHistType(TPBHistType::kHist1D); // Owned by the input list
   fRunCPU = new TProofBenchRunCPU(htype, fNHist, fOutFile);
   if (!fCPUSel.IsNull()) fRunCPU->SetSelName(fCPUSel);
   if (!fSelOption.IsNull()) fRunDS->SetSelOption(fSelOption);
   if (!fCPUPar.IsNull()) fRunCPU->SetParList(fCPUPar);
   fRunCPU->Run(nevents, start, stop, step, fNtries, fDebug, -1);

   // Close the file
   if (SetOutFile(0) != 0)
      Warning("RunCPU", "problems closing '%s'", fOutFileName.Data());

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the CPU run scanning over the number of workers per node
/// Return 0 on success, -1 on error

Int_t TProofBench::RunCPUx(Long64_t nevents, Int_t start, Int_t stop)
{
   // Open the file for the results
   if (OpenOutFile(kTRUE) != 0) {
      Error("RunCPUx", "problems opening '%s' to save the result", fOutFileName.Data());
      return -1;
   }
   fUnlinkOutfile = kFALSE;

   SafeDelete(fRunCPU);
   TPBHistType *htype = new TPBHistType(TPBHistType::kHist1D); // Owned by the input list
   fRunCPU = new TProofBenchRunCPU(htype, fNHist, fOutFile);
   if (!fCPUSel.IsNull()) fRunCPU->SetSelName(fCPUSel);
   if (!fSelOption.IsNull()) fRunDS->SetSelOption(fSelOption);
   if (!fCPUPar.IsNull()) fRunCPU->SetParList(fCPUPar);
   fRunCPU->Run(nevents, start, stop, -2, fNtries, fDebug, -1);

   // Close the file
   if (SetOutFile(0) != 0)
      Warning("RunCPUx", "problems closing '%s'", fOutFileName.Data());

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the CPU speedup plot.
///  opt = 'typewhat', e.g. 'std:max:'
///    type = 'std:'      draw standard evt/s plot
///           'stdx:'     draw standard evt/s plot, 1 worker per node
///           'norm:'     draw normalized plot
///           'normx:'    draw normalized plot, 1 worker per node
///    what = 'max:'      draw max rate
///           'avg:'      draw average rate
///           'all:'      draw max and average rate on same plot (default)
///  dofit =  0           no fit
///           1           fit with the relevant '1st degree related' function
///           2           fit with the relevant '2nd degree related' function
///           3           fit with varying rcpu function
///     n0 = for dofit == 3, number of real cores
///     n1 = for dofit == 3, number of total cores (real + hyperthreaded)
///

void TProofBench::DrawCPU(const char *outfile, const char *opt, Bool_t verbose,
                          Int_t dofit, Int_t n0, Int_t n1)
{
   // Get the TProfile an create the graphs
   TFile *fout = TFile::Open(outfile, "READ");
   if (!fout || (fout && fout->IsZombie())) {
      ::Error("DrawCPU", "could not open file '%s' ...", outfile);
      return;
   }

   // Get description
   TString description("<not available>");
   TNamed *nmdesc = (TNamed *) fout->Get("PB_description");
   if (nmdesc) description = nmdesc->GetTitle();

   // Parse option
   TString oo(opt);
   Bool_t isNorm = (oo.Contains("norm")) ? kTRUE : kFALSE;
   Bool_t isX = (oo.Contains("stdx:") || oo.Contains("normx:")) ? kTRUE : kFALSE;
   Bool_t doAvg = kTRUE, doMax = kTRUE;
   if (oo.Contains("avg:")) doMax = kFALSE;
   if (oo.Contains("max:")) doAvg = kFALSE;

   const char *dirn = (isX) ? "RunCPUx" : "RunCPU";
   TDirectory *d = (TDirectory *) fout->Get(dirn);
   if (!d) {
      ::Error("DrawCPU", "could not find directory '%s' ...", dirn);
      fout->Close();
      delete fout;
      return;
   }
   d->cd();

   TString hprofn, hmaxn;
   const char *lx = (isX) ? "_x" : "";
   const char *ln = (isNorm) ? "Norm" : "Prof";
   hprofn.Form("%s%s_CPU_QR_Evts", ln, lx);
   hmaxn.Form("%s%s_CPU_PS_MaxEvts", ln, lx);

   Double_t xmin = -1., xmax = -1.;
   Double_t ami = -1., amx = -1., mmi = -1., mmx = -1.;
   Int_t kamx = -1, kmmx = -1, nbins = -1;
   Double_t ymx = -1., ymi = -1.;

   TProfile *pf = 0;
   Int_t kmx = -1;

   TProfile *pfav = 0;
   TGraphErrors *grav = 0;
   if (doAvg) {
      if (!(grav = GetGraph(d, hprofn, nbins, xmin, xmax, ami, amx, kamx, pfav))) {
         ::Error("DrawCPU", "could not find '%s' ...", hprofn.Data());
         fout->Close();
         delete fout;
         return;
      }
      ymx = amx;
      ymi = ami;
      pf = pfav;
      kmx = kamx;
   }
   TProfile *pfmx = 0;
   TGraphErrors *grmx = 0;
   if (doMax) {
      if (!(grmx = GetGraph(d, hmaxn, nbins, xmin, xmax, mmi, mmx, kmmx, pfmx))) {
         ::Warning("DrawCPU", "could not find '%s': feature added in 5.34/11", hmaxn.Data());
         if (!grav) {
            // Nothing to do if not asked for the average
            fout->Close();
            delete fout;
            return;
         }
         doMax = kFALSE;
      }
      if (mmx > ymx) ymx = mmx;
      if ((ymi > 0 && mmi < ymi) || (ymi < 0.)) ymi = mmi;
      pf = pfmx;
      kmx = kmmx;
   }

   // Create the canvas
   TCanvas *cpu = new TCanvas("cpu", "Rate vs wrks",204,69,1050,502);
   cpu->Range(-3.106332,0.7490716,28.1362,1.249867);

   TH1F *hgr = new TH1F("Graph-CPU"," CPU speed-up", nbins*4, xmin, xmax);
   hgr->SetMaximum(ymx + (ymx-ymi)*0.2);
   hgr->SetMinimum(0);
   if (isNorm) hgr->SetMaximum(ymx*1.2);
   hgr->SetDirectory(0);
   hgr->SetStats(0);
   hgr->GetXaxis()->SetTitle(pf->GetXaxis()->GetTitle());
   hgr->GetXaxis()->CenterTitle(true);
   hgr->GetXaxis()->SetLabelSize(0.05);
   hgr->GetXaxis()->SetTitleSize(0.06);
   hgr->GetXaxis()->SetTitleOffset(0.62);
   hgr->GetYaxis()->SetTitleSize(0.08);
   hgr->GetYaxis()->SetTitleOffset(0.52);
   hgr->GetYaxis()->SetTitle("Rate (events/s)");

   TLegend *leg = 0;
   if (isNorm) {
      leg = new TLegend(0.7, 0.8, 0.9, 0.9);
   } else {
      leg = new TLegend(0.1, 0.8, 0.3, 0.9);
   }

   gStyle->SetOptTitle(0);
   TGraphErrors *gr = 0;
   if (doAvg) {
      grav->SetFillColor(1);
      grav->SetLineColor(13);
      grav->SetMarkerColor(4);
      grav->SetMarkerStyle(21);
      grav->SetMarkerSize(1.2);
      grav->SetHistogram(hgr);

      if (verbose) grav->Print();
      grav->Draw("alp");
      leg->AddEntry(grav, "Average", "P");
      gr = grav;
   }
   if (doMax) {
      grmx->SetFillColor(1);
      grmx->SetLineColor(13);
      grmx->SetMarkerColor(2);
      grmx->SetMarkerStyle(29);
      grmx->SetMarkerSize(1.8);
      grmx->SetHistogram(hgr);

      if (verbose) grmx->Print();
      if (doAvg) {
         grmx->Draw("lpSAME");
      } else {
         grmx->Draw("alp");
      }
      leg->AddEntry(grmx, "Maximum", "P");
      gr = grmx;
   }
   leg->Draw();
   gPad->Update();

   if (dofit > 0) {
      // Make sure the fitting functions are defined
      Double_t xmi = 0.9;
      if (nbins > 5) xmi = 1.5;
      AssertFittingFun(xmi, nbins + .1);

      // Starting point for the parameters and fit
      Double_t normrate = -1.;
      if (dofit == 1) {
         if (isNorm) {
            fgFp1n->SetParameter(0, pf->GetBinContent(1));
            fgFp1n->SetParameter(1, pf->GetBinContent(nbins-1));
            gr->Fit(fgFp1n);
            if (verbose) fgFp1n->Print();
            normrate = fgFp1n->GetParameter(1);
         } else {
            fgFp1->SetParameter(0, 0.);
            fgFp1->SetParameter(1, pf->GetBinContent(1));
            gr->Fit(fgFp1);
            if (verbose) fgFp1->Print();
            normrate = fgFp1->Derivative(1.);
         }
      } else if (dofit == 2) {
         if (isNorm) {
            fgFp2n->SetParameter(0, pf->GetBinContent(1));
            fgFp2n->SetParameter(1, pf->GetBinContent(nbins-1));
            fgFp2n->SetParameter(2, 0.);
            gr->Fit(fgFp2n);
            if (verbose) fgFp2n->Print();
            normrate = fgFp2n->GetParameter(1);
         } else {
            fgFp2->SetParameter(0, 0.);
            fgFp2->SetParameter(1, pf->GetBinContent(1));
            fgFp2->SetParameter(2, 0.);
            gr->Fit(fgFp2);
            if (verbose) fgFp2->Print();
            normrate = fgFp2->Derivative(1.);
         }
      } else {
         // Starting point for the parameters and fit
         gFioVn0 = (n0 > 0) ? n0 : (Int_t) (nbins + .1)/2.;
         gFioVn1 = (n1 > 0) ? n1 : (Int_t) (nbins + .1);
         if (isNorm) {
            fgFp3n->SetParameter(0, 0.);
            fgFp3n->SetParameter(1, pf->GetBinContent(1));
            fgFp3n->SetParameter(2, pf->GetBinContent(nbins-1));
            gr->Fit(fgFp3n);
            if (verbose) fgFp3n->Print();
            normrate = pf->GetBinContent(1);
         } else {
            fgFp3->SetParameter(0, 0.);
            fgFp3->SetParameter(1, 0.);
            fgFp3->SetParameter(2, pf->GetBinContent(1));
            gr->Fit(fgFp3);
            if (verbose) fgFp3->Print();
            normrate = fgFp3->Derivative(1.);
         }
      }

      // Notify the cluster performance parameters
      if (!isNorm) {
         printf("* ************************************************************ *\n");
         printf("*                                                              *\r");
         printf("* Cluster: %s\n", description.Data());
         printf("* Performance measurement from scalability plot:               *\n");
         printf("*                                                              *\r");
         printf("*    rate max:         %.3f\tmegaRNGPS (@ %d workers)\n", ymx/1000000, kmx);
         printf("*                                                              *\r");
         printf("*    per-worker rate:  %.3f\tmegaRNGPS \n", normrate/1000000);
         printf("* ************************************************************ *\n");
      } else {
         printf("* ************************************************************ *\n");
         printf("*                                                              *\r");
         printf("* Cluster: %s\n", description.Data());
         printf("*                                                              *\r");
         printf("* Per-worker rate from normalized plot:  %.3f\tmegaRNGPS\n", normrate/1000000);
         printf("* ************************************************************ *\n");
      }
   }
   // Close the file
   fout->Close();
   if (grav) fgGraphs->Add(grav);
   if (grmx) fgGraphs->Add(grmx);
}

////////////////////////////////////////////////////////////////////////////////
/// Get from TDirectory 'd' the TProfile named 'pfn' and create the graph.
/// Return also the max y in mx.

TGraphErrors *TProofBench::GetGraph(TDirectory *d, const char *pfn, Int_t &nb,
                                    Double_t &xmi, Double_t &xmx,
                                    Double_t &ymi, Double_t &ymx, Int_t &kmx, TProfile *&pf)
{
   // Sanity checks
   if (!d || !pfn || (pfn && strlen(pfn) <= 0)) {
      ::Error("TProofBench::GetGraph", "directory or name not defined!");
      return (TGraphErrors *)0;
   }

   TList *keylist = d->GetListOfKeys();
   TKey *key = 0;
   TIter nxk(keylist);
   while ((key = (TKey *) nxk())) {
      if (TString(key->GetName()).BeginsWith(pfn)) {
         pf = (TProfile *) d->Get(key->GetName());
         break;
      }
   }
   // Sanity checks
   if (!pf) {
      ::Error("TProofBench::GetGraph", "TProfile for '%s' not found in directory '%s'", pfn, d->GetName());
      return (TGraphErrors *)0;
   }

   nb = pf->GetNbinsX();
   TGraphErrors *gr = new TGraphErrors(nb);
   gr->SetName(TString::Format("Graph_%s", pfn));
   Double_t xx, ex, yy, ey;
   ymi = pf->GetBinContent(1);
   ymx = ymi;
   xmi = pf->GetBinCenter(1) - pf->GetBinWidth(1)/2. ;
   xmx = pf->GetBinCenter(nb) + pf->GetBinWidth(nb)/2. ;
   kmx = -1;
   for (Int_t k = 1;k <= nb; k++) {
      xx = pf->GetBinCenter(k);
      ex = pf->GetBinWidth(k) * .001;
      yy = pf->GetBinContent(k);
      ey = pf->GetBinError(k);
      if (k == 1) {
         ymi = yy;
         ymx = yy;
         kmx = k;
      } else {
         if (yy < ymi) ymi = yy;
         if (yy > ymx) { ymx = yy; kmx = k; }
      }
      gr->SetPoint(k-1, xx, yy);
      gr->SetPointError(k-1, ex, ey);
   }

   // Done
   return gr;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that the fitting functions are defined

void TProofBench::AssertFittingFun(Double_t mi, Double_t mx)
{
   if (!fgFp1) {
      fgFp1 = new TF1("funp1", funp1, mi, mx, 2);
      fgFp1->SetParNames("offset", "slope");
   }

   if (!fgFp1n) {
      fgFp1n = new TF1("funp1n", funp1n, mi, mx, 2);
      fgFp1n->SetParNames("decay", "norm rate");
   }

   if (!fgFp2) {
      fgFp2 = new TF1("funp2", funp2, mi, mx, 3);
      fgFp2->SetParNames("offset", "slope", "deviation");
   }

   if (!fgFp2n) {
      fgFp2n = new TF1("funp2n", funp2n, mi, mx, 3);
      fgFp2n->SetParNames("decay", "norm rate", "deviation");
   }

   if (!fgFp3) {
      fgFp3 = new TF1("funcpuv", funcpuv, mi, mx, 3);
      fgFp3->SetParNames("offset", "slope real", "slope hyper");
   }

   if (!fgFp3n) {
      fgFp3n = new TF1("funcpuvn", funcpuvn, mi, mx, 3);
      fgFp3n->SetParNames("offset", "slope real", "slope hyper");
   }

   if (!fgFio) {
      fgFio = new TF1("funio", funio, mi, mx, 3);
      fgFio->SetParNames("R1", "RIO", "TotIO");
   }
   if (!fgFioV) {
      fgFioV = new TF1("funiov", funiov, mi, mx, 5);
      fgFioV->SetParNames("rio", "b1", "b2", "nc", "ri");
   }

}

////////////////////////////////////////////////////////////////////////////////

class fileDesc : public TNamed {
public:
   Long_t  fMtime; // Modification time
   TString fDesc; // Test description string, if any
   fileDesc(const char *n, const char *o,
            Long_t t, const char *d) : TNamed(n, o), fMtime(t), fDesc(d) { }
   Int_t   Compare(const TObject *o) const {
      const fileDesc *fd = static_cast<const fileDesc *>(o);
      if (!fd || (fd && fd->fMtime == fMtime)) return 0;
      if (fMtime < fd->fMtime) return -1;
      return 1;
   }
};

////////////////////////////////////////////////////////////////////////////////
/// Get performance specs. Check file 'path', or files in directory 'path'
/// (default current directory).
/// The degree of the polynomial used for the fit is 'degfit' (default 1).

void TProofBench::GetPerfSpecs(const char *path, Int_t degfit)
{
   // Locate the file (ask if many)
   TString pp(path), fn, oo;
   if (pp.IsNull()) pp = gSystem->WorkingDirectory();
   FileStat_t st;
   if (gSystem->GetPathInfo(pp.Data(), st) != 0) {
      ::Error("TProofBench::GetPerfSpecs", "path '%s' could not be stat'ed - abort", pp.Data());
      return;
   }
   TSortedList filels;
   if (R_ISDIR(st.fMode)) {
      // Scan the directory
      void *dirp = gSystem->OpenDirectory(pp.Data());
      if (!dirp) {
         ::Error("TProofBench::GetPerfSpecs", "directory path '%s' could nto be open - abort", pp.Data());
         return;
      }
      const char *ent = 0;
      while ((ent = gSystem->GetDirEntry(dirp))) {
         if (!strcmp(ent, ".") || !strcmp(ent, "..")) continue;
         fn.Form("%s/%s", pp.Data(), ent);
         if (gSystem->GetPathInfo(fn.Data(), st) != 0) continue;
         if (!R_ISREG(st.fMode)) continue;
         fn += "?filetype=raw";
         TFile *f = TFile::Open(fn);
         if (!f) continue;
         char rr[5] = {0};
         if (!f->ReadBuffer(rr, 4)) {
            if (!strncmp(rr, "root", 4)) {
               SafeDelete(f);
               fn.ReplaceAll("?filetype=raw", "");
               if ((f = TFile::Open(fn))) {
                  TString desc("<no decription>");
                  TNamed *nmdesc = (TNamed *) f->Get("PB_description");
                  if (nmdesc) desc = nmdesc->GetTitle();
                  if (f->GetListOfKeys()->FindObject("RunCPU"))
                     filels.Add(new fileDesc(fn, "std:", st.fMtime, desc.Data()));
                  if (f->GetListOfKeys()->FindObject("RunCPUx"))
                     filels.Add(new fileDesc(fn, "stdx:", st.fMtime, desc.Data()));
               } else {
                  ::Warning("TProofBench::GetPerfSpecs", "problems opening '%s'", fn.Data());
               }
            }
         }
         SafeDelete(f);
      }
   } else if (!R_ISREG(st.fMode)) {
      ::Error("TProofBench::GetPerfSpecs",
              "path '%s' not a regular file nor a directory - abort", pp.Data());
      return;
   } else {
      // This is the file
      fn = pp;
      // Check it
      TString emsg;
      Bool_t isOk = kFALSE;
      if (gSystem->GetPathInfo(fn.Data(), st) == 0) {
         fn += "?filetype=raw";
         TFile *f = TFile::Open(fn);
         if (f) {
            char rr[5] = {0};
            if (!(f->ReadBuffer(rr, 4))) {
               if (!strncmp(rr, "root", 4)) {
                  fn.ReplaceAll("?filetype=raw", "");
                  if ((f = TFile::Open(fn))) {
                     if (f->GetListOfKeys()->FindObject("RunCPU")) oo = "std:";
                     if (f->GetListOfKeys()->FindObject("RunCPUx")) oo = "stdx:";
                     SafeDelete(f);
                     if (!oo.IsNull()) {
                        isOk = kTRUE;
                     } else {
                        emsg.Form("path '%s' does not contain the relevant dirs - abort", fn.Data());
                     }
                  } else {
                     emsg.Form("path '%s' cannot be open - abort", fn.Data());
                  }
               } else {
                  emsg.Form("'%s' is not a ROOT file - abort", fn.Data());
               }
            } else {
               emsg.Form("could not read first 4 bytes from '%s' - abort", fn.Data());
            }
            SafeDelete(f);
         } else {
            emsg.Form("path '%s' cannot be open in raw mode - abort", fn.Data());
         }
      } else {
         emsg.Form("path '%s' cannot be stated - abort", fn.Data());
      }
      if (!isOk) {
         ::Error("TProofBench::GetPerfSpecs", "%s", emsg.Data());
         return;
      }
   }

   fileDesc *nm = 0;
   // Ask the user, if more then 1
   if (filels.GetSize() == 1) {
      nm = (fileDesc *) filels.First();
      fn = nm->GetName();
      oo = nm->GetTitle();
   } else if (filels.GetSize() > 1) {
      TIter nxf(&filels);
      Int_t idx = 0;
      Printf("Several possible files found:");
      while ((nm = (fileDesc *) nxf())) {
         Printf("  %d\t%s\t%s\t%s (file: %s)", idx++, nm->GetTitle(),
                TTimeStamp(nm->fMtime).AsString("s"), nm->fDesc.Data(), nm->GetName());
      }
      TString a(Getline(TString::Format("Make your choice [%d] ", idx-1)));
      if (a.IsNull() || a[0] == '\n') a.Form("%d", idx-1);
      idx = a.Atoi();
      if ((nm = (fileDesc *) filels.At(idx))) {
         fn = nm->GetName();
         oo = nm->GetTitle();
      } else {
         ::Error("TProofBench::GetPerfSpecs", "chosen index '%d' does not exist - abort", idx);
         return;
      }
   } else {
      if (fn.IsNull()) {
         ::Error("TProofBench::GetPerfSpecs",
                 "path '%s' is a directory but no ROOT file found in it - abort", pp.Data());
         return;
      }
   }

   // Now get the specs
   TProofBench::DrawCPU(fn.Data(), oo.Data(), kFALSE, degfit);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a test using dataset 'dset'
/// Return 0 on success, -1 on error
/// Open the file for the results

Int_t TProofBench::RunDataSet(const char *dset,
                              Int_t start, Int_t stop, Int_t step)
{
   if (OpenOutFile(kTRUE) != 0) {
      Error("RunDataSet", "problems opening '%s' to save the result", fOutFileName.Data());
      return -1;
   }
   fUnlinkOutfile = kFALSE;

   if (fReleaseCache) ReleaseCache(dset);
   SafeDelete(fRunDS);
   TPBReadType *readType = fReadType;
   if (!readType) readType = new TPBReadType(TPBReadType::kReadOpt);
   fRunDS = new TProofBenchRunDataRead(fDS, readType, fOutFile);
   if (!fDataSel.IsNull()) fRunDS->SetSelName(fDataSel);
   if (!fSelOption.IsNull()) fRunDS->SetSelOption(fSelOption);
   if (!fDataPar.IsNull()) fRunDS->SetParList(fDataPar);
   fRunDS->SetReleaseCache(fReleaseCache);
   fRunDS->Run(dset, start, stop, step, fNtries, fDebug, -1);
   SafeDelete(readType);

   // Close the file
   if (SetOutFile(0) != 0)
      Warning("RunDataSet", "problems closing '%s'", fOutFileName.Data());

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a test using dataset 'dset' scanning over the number of workers
/// per node.
/// Return 0 on success, -1 on error
/// Open the file for the results

Int_t TProofBench::RunDataSetx(const char *dset, Int_t start, Int_t stop)
{
   if (OpenOutFile(kTRUE) != 0) {
      Error("RunDataSetx", "problems opening '%s' to save the result", fOutFileName.Data());
      return -1;
   }
   fUnlinkOutfile = kFALSE;

   ReleaseCache(dset);
   SafeDelete(fRunDS);
   TPBReadType *readType = fReadType;
   if (!readType) readType = new TPBReadType(TPBReadType::kReadOpt);
   fRunDS = new TProofBenchRunDataRead(fDS, readType, fOutFile);
   if (!fDataSel.IsNull()) fRunDS->SetSelName(fDataSel);
   if (!fSelOption.IsNull()) fRunDS->SetSelOption(fSelOption);
   if (!fDataPar.IsNull()) fRunDS->SetParList(fDataPar);
   fRunDS->Run(dset, start, stop, -2, fNtries, fDebug, -1);
   SafeDelete(readType);

   // Close the file
   if (SetOutFile(0) != 0)
      Warning("RunDataSetx", "problems closing '%s'", fOutFileName.Data());

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the CPU speedup plot.
///  opt = 'typewhat', e.g. 'std:max:'
///    type = 'std:'      draw standard plot
///           'stdx:'     draw standard plot, 1 worker per node
///           'norm:'     draw normalized plot
///           'normx:'    draw normalized plot, 1 worker per node
///    what = 'max:'      draw max rate
///           'avg:'      draw average rate
///           'all:'      draw max and average rate on same plot (default)
/// type =    'mbs'           MB/s scaling plots (default)
///           'evts'          Event/s scaling plots
///  dofit =  0           no fit
///           1           fit with default 3 parameter saturated I/O formula
///           2           fit with 4 parameter saturated I/O formula (varying Rcpu)
///     n0 = for dofit == 2, number of real cores
///     n1 = for dofit == 2, number of total cores (real + hyperthreaded)
///

void TProofBench::DrawDataSet(const char *outfile,
                              const char *opt, const char *type, Bool_t verbose,
                              Int_t dofit, Int_t n0, Int_t n1)
{
   // Get the TProfile an create the graphs
   TFile *fout = TFile::Open(outfile, "READ");
   if (!fout || (fout && fout->IsZombie())) {
      ::Error("DrawDataSet", "could not open file '%s' ...", outfile);
      return;
   }

   // Get description
   TString description("<not available>");
   TNamed *nmdesc = (TNamed *) fout->Get("PB_description");
   if (nmdesc) description = nmdesc->GetTitle();

   // Parse option
   TString oo(opt);
   Bool_t isNorm = (oo.Contains("norm")) ? kTRUE : kFALSE;
   Bool_t isX = (oo.Contains("stdx:") || oo.Contains("normx:")) ? kTRUE : kFALSE;
   Bool_t doAvg = (oo.Contains("all:") || oo.Contains("avg:")) ? kTRUE : kFALSE;
   Bool_t doMax = (oo.Contains("all:") || oo.Contains("max:")) ? kTRUE : kFALSE;

   const char *dirn = (isX) ? "RunDataReadx" : "RunDataRead";
   TDirectory *d = (TDirectory *) fout->Get(dirn);
   if (!d) {
      ::Error("DrawCPU", "could not find directory '%s' ...", dirn);
      fout->Close();
      delete fout;
      return;
   }
   d->cd();

   TString hprofn, hmaxn;
   const char *lx = (isX) ? "_x" : "";
   const char *ln = (isNorm) ? "Norm" : "Prof";
   Bool_t isIO = kTRUE;
   if (type && !strcmp(type, "evts")) {
      hprofn.Form("%s%s_DataRead_QR_Evts", ln, lx);
      hmaxn.Form("%s%s_DataRead_PS_MaxEvts", ln, lx);
      isIO = kFALSE;
   } else {
      hprofn.Form("%s%s_DataRead_QR_IO", ln, lx);
      hmaxn.Form("%s%s_DataRead_PS_MaxIO", ln, lx);
   }

   Double_t xmin = -1., xmax = -1.;
   Double_t ami = -1., amx = -1., mmi = -1., mmx = -1.;
   Int_t kamx = -1, kmmx = -1, nbins = -1;
   Double_t ymx = -1., ymi = -1.;

   TProfile *pf = 0;
   Int_t kmx = -1;

   TProfile *pfav = 0;
   TGraphErrors *grav = 0;
   if (doAvg) {
      if (!(grav = GetGraph(d, hprofn, nbins, xmin, xmax, ami, amx, kamx, pfav))) {
         ::Error("DrawCPU", "could not find '%s' ...", hprofn.Data());
         fout->Close();
         delete fout;
         return;
      }
      ymx = amx;
      ymi = ami;
      pf = pfav;
      kmx = kamx;
   }

   TProfile *pfmx = 0;
   TGraphErrors *grmx = 0;
   if (doMax) {
      if (!(grmx = GetGraph(d, hmaxn, nbins, xmin, xmax, mmi, mmx, kmmx, pfmx))) {
         ::Warning("DrawCPU", "could not find '%s': feature added in 5.34/11", hmaxn.Data());
         if (!grav) {
            // Nothing to do if not asked for the average
            fout->Close();
            delete fout;
            return;
         }
         doMax = kFALSE;
      }
      if (mmx > ymx) ymx = mmx;
      if ((ymi > 0 && mmi < ymi) || (ymi < 0.)) ymi = mmi;
      pf = pfmx;
      kmx = kmmx;
   }

   // Create the canvas
   TCanvas *cpu = new TCanvas("dataset", "Rate vs wrks",204,69,1050,502);
   cpu->Range(-3.106332,0.7490716,28.1362,1.249867);

   TH1F *hgr = new TH1F("Graph-DataSet"," Data Read speed-up", nbins*4, xmin, xmax);
   hgr->SetMaximum(ymx + (ymx-ymi)*0.2);
   hgr->SetMinimum(0);
   if (isNorm) hgr->SetMaximum(ymx*1.2);
   hgr->SetDirectory(0);
   hgr->SetStats(0);
   hgr->GetXaxis()->SetTitle(pf->GetXaxis()->GetTitle());
   hgr->GetXaxis()->CenterTitle(true);
   hgr->GetXaxis()->SetLabelSize(0.05);
   hgr->GetXaxis()->SetTitleSize(0.06);
   hgr->GetXaxis()->SetTitleOffset(0.62);
   hgr->GetYaxis()->SetLabelSize(0.06);
   hgr->GetYaxis()->SetTitleSize(0.08);
   hgr->GetYaxis()->SetTitleOffset(0.52);
   if (isIO) {
      hgr->GetYaxis()->SetTitle("Rate (MB/s)");
   } else {
      hgr->GetYaxis()->SetTitle("Rate (events/s)");
   }

   TLegend *leg = 0;
   if (isNorm) {
      leg = new TLegend(0.7, 0.8, 0.9, 0.9);
   } else {
      leg = new TLegend(0.1, 0.8, 0.3, 0.9);
   }

   TGraphErrors *gr = 0;
   if (doAvg) {
      grav->SetFillColor(1);
      grav->SetLineColor(13);
      grav->SetMarkerColor(4);
      grav->SetMarkerStyle(21);
      grav->SetMarkerSize(1.2);
      grav->SetHistogram(hgr);

      if (verbose) grav->Print();
      grav->Draw("alp");
      leg->AddEntry(grav, "Average", "P");
      gr = grav;
   }
   if (doMax) {
      grmx->SetFillColor(1);
      grmx->SetLineColor(13);
      grmx->SetMarkerColor(2);
      grmx->SetMarkerStyle(29);
      grmx->SetMarkerSize(1.8);
      grmx->SetHistogram(hgr);

      if (verbose) grmx->Print();
      if (doAvg) {
         grmx->Draw("lpSAME");
      } else {
         grmx->Draw("alp");
      }
      leg->AddEntry(grmx, "Maximum", "P");
      gr = grmx;
   }
   leg->Draw();
   gPad->Update();

   Double_t normrate = -1.;
   if (dofit > 0) {
      // Make sure the fitting functions are defined
      Double_t xmi = 0.9;
      if (nbins > 5) xmi = 1.5;
      AssertFittingFun(xmi, nbins + .1);

      if (dofit == 1) {
         // Starting point for the parameters and fit
         fgFio->SetParameter(0, pf->GetBinContent(1));
         fgFio->SetParameter(1, pf->GetBinContent(nbins-1));
         fgFio->SetParameter(2, pf->GetBinContent(nbins-1));
         gr->Fit(fgFio);
         if (verbose) fgFio->Print();
         normrate = fgFio->Derivative(1.);
      } else if (dofit > 1) {
         // Starting point for the parameters and fit
         gFioVn0 = (n0 > 0) ? n0 : (Int_t) (nbins + .1)/2.;
         gFioVn1 = (n1 > 0) ? n1 : (Int_t) (nbins + .1);
         fgFioV->SetParameter(0, 20.);
         fgFioV->SetParameter(1, pf->GetBinContent(1));
         fgFioV->SetParameter(2, pf->GetBinContent(1));
         fgFioV->SetParameter(3, 4.);
         fgFioV->SetParameter(4, 1000.);

         gr->Fit(fgFioV);
         if (verbose) fgFio->Print();
         normrate = fgFioV->Derivative(1.);
      }
   }

   // Notify the cluster performance parameters
   if (!isNorm) {
      printf("* ************************************************************ *\n");
      printf("*                                                              *\r");
      printf("* Cluster: %s\n", description.Data());
      printf("* Performance measurement from scalability plot:               *\n");
      printf("*                                                              *\r");
      if (isIO) {
         printf("*    rate max:         %.3f\tMB/s (@ %d workers)\n", ymx, kmx);
         printf("*                                                              *\r");
         printf("*    per-worker rate:  %.3f\tMB/s \n", normrate);
      } else {
         printf("*    rate max:         %.3f\tevts/s (@ %d workers)\n", ymx, kmx);
      }
      printf("* ************************************************************ *\n");
   }
   // Close the file
   fout->Close();
   if (grav) fgGraphs->Add(grav);
   if (grmx) fgGraphs->Add(grmx);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the efficiency plot.
///  opt = 'cpu' or 'data' (default the first found)
///

void TProofBench::DrawEfficiency(const char *outfile,
                                 const char *opt, Bool_t verbose)
{
   // Get the TProfile an create the graphs
   TFile *fout = TFile::Open(outfile, "READ");
   if (!fout || (fout && fout->IsZombie())) {
      ::Error("DrawEfficiency", "could not open file '%s' ...", outfile);
      return;
   }

   // Get description
   TString description("<not available>");
   TNamed *nmdesc = (TNamed *) fout->Get("PB_description");
   if (nmdesc) description = nmdesc->GetTitle();

   // Parse option
   TString oo(opt), ln("CPU");
   const char *dirs[4] = { "RunCPU", "RunCPUx", "RunDataRead", "RunDataReadx"};
   const char *labs[4] = { "CPU", "CPU", "DataRead", "DataRead"};
   Int_t fst = 0, lst = 3;
   if (oo == "cpu") {
      lst = 0;
   } else if (oo == "cpux") {
      fst = 1;
      lst = 1;
   } else if (oo.BeginsWith("data")) {
      if (oo.EndsWith("x")) {
         fst = 3;
         lst = 3;
      } else {
         fst = 2;
         lst = 2;
      }
   }
   const char *dirn = 0;
   TDirectory *d = 0;
   for (Int_t i = fst; i <= lst; i++) {
      if ((d = (TDirectory *) fout->Get(dirs[i]))) {
         dirn = dirs[i];
         ln = labs[i];
         break;
      }
   }
   if (!d && !dirn) {
      ::Error("DrawEfficiency", "could not find directory ...");
      fout->Close();
      delete fout;
      return;
   }
   d->cd();

   TString hprof;
   hprof.Form("Prof_%s_CPU_eff", ln.Data());

   Double_t xmin = -1., xmax = -1.;
   Int_t kmx = -1, nbins = -1;
   Double_t ymx = -1., ymi = -1.;

   TProfile *pf = 0;
   TGraphErrors *gr = 0;
   if (!(gr = GetGraph(d, hprof, nbins, xmin, xmax, ymi, ymx, kmx, pf))) {
      ::Error("DrawEfficiency", "could not find '%s' ...", hprof.Data());
      fout->Close();
      delete fout;
      return;
   }

   // Create the canvas
   TCanvas *cpu = new TCanvas("efficiency", "efficiency vs wrks",204,69,1050,502);
   cpu->Range(-3.106332,0.7490716,28.1362,1.249867);

   TH1F *hgr = new TH1F("Graph-Efficiency","CPU effectiveness", nbins*4, xmin, xmax);
   hgr->SetMaximum(1.2);
   hgr->SetMinimum(0);
   hgr->SetDirectory(0);
   hgr->SetStats(0);
   hgr->GetXaxis()->SetTitle(pf->GetXaxis()->GetTitle());
   hgr->GetXaxis()->CenterTitle(true);
   hgr->GetXaxis()->SetLabelSize(0.05);
   hgr->GetXaxis()->SetTitleSize(0.06);
   hgr->GetXaxis()->SetTitleOffset(0.62);
   hgr->GetYaxis()->SetLabelSize(0.06);
   hgr->GetYaxis()->SetTitleSize(0.08);
   hgr->GetYaxis()->SetTitleOffset(0.52);
   hgr->GetYaxis()->SetTitle("CPU effectiveness");

   gr->SetFillColor(1);
   gr->SetLineColor(13);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetMarkerSize(1.2);
   gr->SetHistogram(hgr);

   if (verbose) gr->Print();
   gr->Draw("alp");

   // Notify the cluster performance parameters
   printf("* ************************************************************ *\n");
   printf("*                                                              *\r");
   printf("* Cluster: %s\n", description.Data());
   printf("* CPU effectiveness measurement:                               *\n");
   printf("*                                                              *\r");
   printf("*    effectiveness max:     %.3f (@ %d workers)\n", ymx, kmx);
   printf("*                                                              *\r");
   printf("* ************************************************************ *\n");
   // Close the file
   fout->Close();
   if (gr) fgGraphs->Add(gr);
}

////////////////////////////////////////////////////////////////////////////////
/// Release memory cache for dataset 'dset'
/// Return 0 on success, -1 on error

Int_t TProofBench::ReleaseCache(const char *dset)
{
   // Do it via the dataset handler
   if (!fDS) fDS = new TProofBenchDataSet(fProofDS);
   return fDS ? fDS->ReleaseCache(dset) : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Physically remove the dataset 'dset', i.e. remove the dataset and the files
/// it describes
/// Return 0 on success, -1 on error

Int_t TProofBench::RemoveDataSet(const char *dset)
{
   // Do it via the dataset handler
   if (!fDS) fDS = new TProofBenchDataSet(fProofDS);
   return fDS ? fDS->RemoveFiles(dset) : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the largest dataset for the run.
/// Defaults for
///          dataset name, filename root
/// are
///          "BenchDataSet", "event"
/// respectively.
/// These can be changed via dset and fnroot, respectively.
/// The string 'fnroot' defines the location of the files, interpreted as an URL.
/// Examples:
///          fnroot             files
///          'event'            <datadir>/event_<ord>_<#>.root
///          '/mss/event'       /mss/event_<ord>_<#>.root
///          'root://srv//mss/event?remote=1'
///                             root://srv//mss/event_<ord>_<#>?remote=1.root
/// Default selector is TSelEventGen. Use SetDataGenSel and SetDataGenPar to change it
/// and to pass the list of PARs defining the alternative selector.
/// The argument 'nevt' controls the number of events per file (-1 for the default,
/// which is 30000).
/// Return 0 on success, -1 on error

Int_t TProofBench::MakeDataSet(const char *dset, Long64_t nevt, const char *fnroot,
                               Bool_t regenerate)
{
   if (dset && strlen(dset) > 0) fDataSet = dset;

   // Load the selector, if needed
   if (!TClass::GetClass(fDataGenSel)) {
      // Is it the default selector?
      if (fDataGenSel == kPROOF_BenchSelDataGenDef) {
         // Load the parfile
         TString par = TString::Format("%s/%s%s.par", TROOT::GetEtcDir().Data(), kPROOF_BenchParDir, kPROOF_BenchDataSelPar);
         Info("MakeDataSet", "uploading '%s' ...", par.Data());
         if (fProof->UploadPackage(par) != 0) {
            Error("MakeDataSet", "problems uploading '%s' - cannot continue", par.Data());
            return -1;
         }
         Info("MakeDataSet", "enabling '%s' ...", kPROOF_BenchDataSelPar);
         if (fProof->EnablePackage(kPROOF_BenchDataSelPar) != 0) {
            Error("MakeDataSet", "problems enabling '%s' - cannot continue", kPROOF_BenchDataSelPar);
            return -1;
         }
      } else {
         if (fDataGenPar.IsNull()) {
            Error("MakeDataSet", "you should load the class '%s' before running the benchmark", fDataGenSel.Data());
            return -1;
         }
      }
      // Load additional PAR files, if any or required by the alternative selector
      TString par;
      Int_t from = 0;
      while (fDataGenPar.Tokenize(par, from, ",")) {
         Info("MakeDataSet", "Uploading '%s' ...", par.Data());
         if (fProof->UploadPackage(par) != 0) {
            Error("MakeDataSet", "problems uploading '%s' - cannot continue", par.Data());
            return -1;
         }
         Info("MakeDataSet", "Enabling '%s' ...", par.Data());
         if (fProof->EnablePackage(par) != 0) {
            Error("MakeDataSet", "problems enabling '%s' - cannot continue", par.Data());
            return -1;
         }
      }
      // Check
      if (!TClass::GetClass(fDataGenSel)) {
         Error("MakeDataSet", "failed to load '%s'", fDataGenSel.Data());
         return -1;
      }
   }

   // For files, 30000 evst each (about 600 MB total) per worker
   TString fn, fnr("event");
   Bool_t remote = kFALSE;
   if (fnroot && strlen(fnroot) > 0) {
      TUrl ur(fnroot, kTRUE);
      if (!strcmp(ur.GetProtocol(), "file") &&
         !gSystem->IsAbsoluteFileName(ur.GetFile())) {
         fnr = fnroot;
      } else {
         fnr = gSystem->BaseName(ur.GetFile());
         // We need to set the basedir
         TString bdir = gSystem->GetDirName(fnroot);
         bdir += "/<fn>";
         fProof->SetParameter("PROOF_BenchmarkBaseDir", bdir.Data());
         // Flag as remote, if so
         if (strcmp(ur.GetProtocol(), "file")) remote = kTRUE;
      }
   }
   TProofNodes pn(fProof);
   TMap *filesmap = new TMap;
   TMap *nodesmap = pn.GetMapOfNodes();
   TIter nxnd(nodesmap);
   TList *wli = 0;
   TObject *obj = 0;
   Int_t kf = 1;
   while ((obj = nxnd()) != 0) {
      if ((wli = dynamic_cast<TList *>(nodesmap->GetValue(obj)))) {
         THashList *fli = new THashList;
         Int_t nf = wli->GetSize() * fNFilesWrk;
         TSlaveInfo *wi = (TSlaveInfo *) wli->First();
         while (nf--) {
            fn.Form("%s-%s-%d.root", fnr.Data(), wi->GetName(), kf++);
            // Add to the node list for generation
            fli->Add(new TObjString(fn));
         }
         filesmap->Add(new TObjString(obj->GetName()), fli);
      }
   }
   filesmap->Print();
   // Prepare for file generation ... add map in the input list
   filesmap->SetName("PROOF_FilesToProcess");
   fProof->AddInput(filesmap);

   // Set parameters for processing
   TString oldpack;
   if (TProof::GetParameter(fProof->GetInputList(), "PROOF_Packetizer", oldpack) != 0) oldpack = "";
   fProof->SetParameter("PROOF_Packetizer", "TPacketizerFile");
   Int_t oldnotass = -1;
   if (TProof::GetParameter(fProof->GetInputList(), "PROOF_ProcessNotAssigned", oldnotass) != 0) oldnotass = -1;
   fProof->SetParameter("PROOF_ProcessNotAssigned", (Int_t)0);

   // Process
   Long64_t ne = (nevt > 0) ? nevt : 30000;
   fProof->SetParameter("PROOF_BenchmarkNEvents", ne);
   fProof->SetParameter("PROOF_BenchmarkRegenerate", Int_t(regenerate));
   fProof->Process(fDataGenSel, (Long64_t) 1);
   fProof->DeleteParameters("PROOF_BenchmarkNEvents");
   fProof->DeleteParameters("PROOF_BenchmarkRegenerate");
   fProof->DeleteParameters("PROOF_BenchmarkBaseDir");

   // Restore parameters
   if (!oldpack.IsNull())
      fProof->SetParameter("PROOF_Packetizer", oldpack);
   else
      fProof->DeleteParameters("PROOF_Packetizer");
   if (oldnotass != -1)
      fProof->SetParameter("PROOF_ProcessNotAssigned", oldnotass);
   else
      fProof->DeleteParameters("PROOF_ProcessNotAssigned");

   // Cleanup
   if (fProof->GetInputList()) fProof->GetInputList()->Remove(filesmap);
   filesmap->SetOwner(kTRUE);
   delete filesmap;

   // The dataset to be registered in the end with proper port
   TFileCollection *fc = new TFileCollection("dum", "dum");

   if (fProof->GetOutputList()) {
      fProof->GetOutputList()->Print();
      TIter nxout(fProof->GetOutputList());
      while ((obj = nxout())) {
         TList *fli = dynamic_cast<TList *>(obj);
         if (fli && TString(fli->GetName()).BeginsWith("PROOF_FilesGenerated_")) {
            TIter nxfg(fli);
            TFileInfo *fi = 0;
            while ((fi = (TFileInfo *) nxfg()))
               fc->Add(fi);
            fli->SetOwner(kFALSE);
         }
      }
      // Register the new dataset, overwriting any existing dataset wth the same name
      // trusting the existing information
      fc->Update();
      if (fc->GetNFiles() > 0) {
         if (remote) fc->SetBit(TFileCollection::kRemoteCollection);
         if (!(fProof->RegisterDataSet(fDataSet, fc, "OT")))
            Warning("MakeDataSet", "problems registering '%s'", dset);
      } else {
         Warning("MakeDataSet", "dataset '%s' is empty!", dset);
      }
   } else {
      Warning("MakeDataSet", "PROOF output list is empty!");
   }

   SafeDelete(fc);

   // Get updated information
   fc = fProof->GetDataSet(fDataSet);
   if (fc) {
      fc->Print("F");
   } else {
      Warning("MakeDataSet", "dataset '%s' was not generated!", fDataSet.Data());
   }

   SafeDelete(fc);

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the files of dataset 'dset' to 'destdir' and create a new dataset named 'dsetdst'
/// decribing them.
/// Return 0 on success, -1 on error

Int_t TProofBench::CopyDataSet(const char *dset, const char *dsetdst, const char *destdir)
{
   // Make some checks
   if (!fProof) {
      Error("CopyDataSet", "no PROOF found - cannot continue");
      return -1;
   }
   if (!dset || (dset && !fProof->ExistsDataSet(dset))) {
      Error("CopyDataSet", "dataset '%s' does not exist", dset);
      return -1;
   }
   if (!dsetdst || (dsetdst && fProof->ExistsDataSet(dsetdst))) {
      if (isatty(0) != 0 && isatty(1) != 0) {
         Printf("Target dataset '%s' exists already:"
                                          " do you want to remove it first?", dsetdst);
         const char *a = Getline("[Y,n] ");
         Printf("a: %s", a);
         if (a[0] == 'Y' || a[0] == 'y' || a[0] == '\n') {
            Info("CopyDataSet", "removing dataset '%s' ...", dsetdst);
            RemoveDataSet(dsetdst);
         } else {
            return -1;
         }
      } else {
         Error("CopyDataSet", "destination dataset '%s' does already exist: remove it first", dsetdst);
         return -1;
      }
   }

   // The TFileCollection object for the new dataset
   TFileCollection *fc = fProof->GetDataSet(dset);
   if (!fc) {
      Error("CopyDataSet", "problems retrieving TFileCollection for dataset '%s'", dset);
      return -1;
   }
   TFileCollection *fcn = new TFileCollection(dsetdst, "");
   TString fn;
   TFileInfo *fi = 0;
   TIter nxfi(fc->GetList());
   while ((fi = (TFileInfo *) nxfi())) {
      fn.Form("%s/%s", destdir, gSystem->BaseName(fi->GetCurrentUrl()->GetFile()));
      Info("CopyDataSet", "adding info for file '%s'", fn.Data());
      fcn->Add(new TFileInfo(fn));
   }
   delete fc;

   // Do it via the dataset handler
   if (!fDS) fDS = new TProofBenchDataSet(fProofDS);
   if (fDS->CopyFiles(dset, destdir) != 0) {
      Error("CopyDataSet", "problems copying files of dataset '%s' to dest dir '%s'", dset, destdir);
      delete fcn;
      return -1;
   }

   // Register the new dataset, overwriting any existing dataset wth the same name
   // trusting the existing information
   Int_t rc = 0;
   if (!(fProof->RegisterDataSet(dsetdst, fcn, "OT"))) {
      Error("CopyDataSet", "problems registering and verifying '%s'", dsetdst);
      rc = -1;
   }
   delete fcn;

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the PROOF instance to be used for dataset operations, like releasing
/// cache ...
/// Use SetProofDS(0) to reset and using the default PROOF

void TProofBench::SetProofDS(TProof *pds)
{
   if (pds && !pds->IsValid()) {
      Error("SetProofDS", "trying to set an invalid PROOF instance");
      return;
   }
   fProofDS = pds ? pds : fProof;
   if (fProofDS) {
      SafeDelete(fDS);
      fDS = new TProofBenchDataSet(fProofDS);
   }
   // Done
   return;
}

