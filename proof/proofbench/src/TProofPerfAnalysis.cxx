// @(#)root/proofx:$Id$
// Author: G.Ganis Nov 2011

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TProofPerfAnalysis
\ingroup proofbench

Set of tools to analyse the performance tree

*/

#include <errno.h>

#include "TProofPerfAnalysis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THashList.h"
#include "TKey.h"
#include "TList.h"
#include "TSortedList.h"
#include "TObjString.h"
#include "TPerfStats.h"
#include "TRegexp.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TMath.h"

//
// Auxilliary internal classes
////////////////////////////////////////////////////////////////////////////////

class TProofPerfAnalysis::TWrkInfo : public TNamed {
public:
   TWrkInfo(const char *ord, const char *name) :
      TNamed(ord, name), fPackets(0), fRemotePackets(0), fEventsProcessed(0),
      fBytesRead(0), fLatency(0), fProcTime(0), fCpuTime(0), fStart(0), fStop(-1),
      fRateT(0), fRateRemoteT(0), fMBRateT(0), fMBRateRemoteT(0), fLatencyT(0) { }
   virtual ~TWrkInfo() { SafeDelete(fRateT); SafeDelete(fRateRemoteT);
                         SafeDelete(fMBRateT); SafeDelete(fMBRateRemoteT);
                         SafeDelete(fLatencyT); }

   Int_t     fPackets;          // Number of packets processed
   Int_t     fRemotePackets;    // Number of processed packet from non-local files
   Long64_t  fEventsProcessed;  // Tot events processed
   Long64_t  fBytesRead;        // Tot bytes read
   Double_t  fLatency;          // Tot latency
   Double_t  fProcTime;         // Tot proc time
   Double_t  fCpuTime;          // Tot CPU time

   Float_t   fStart;            // Start time
   Float_t   fStop;             // Stop time

   TGraph   *fRateT;             // Event processing rate vs time
   TGraph   *fRateRemoteT;       // Event processing rate of remote packets vs time
   TGraph   *fMBRateT;            // Byte processing rate vs time
   TGraph   *fMBRateRemoteT;      // Byte processing rate of remote packets vs time
   TGraph   *fLatencyT;          // Packet latency vs time

   Double_t  AvgRate() { if (fProcTime > 0) return (fEventsProcessed/fProcTime); return -1.; }
   Double_t  AvgIO() { if (fProcTime > 0) return (fBytesRead/fProcTime); return -1.; }

   void Print(Option_t * = "") const {
      Printf(" +++ TWrkInfo ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
      Printf(" +++ Worker:             %s, %s", GetName(), GetTitle());
      Printf(" +++ Activity interval:  %f -> %f", fStart, fStop);
      Printf(" +++ Amounts processed:  %d packets (%d remote), %lld evts, %lld bytes",
                                       fPackets, fRemotePackets, fEventsProcessed, fBytesRead);
      if (fProcTime) {
         Printf(" +++ Processing time:    %f s (CPU: %f s)", fProcTime, fCpuTime);
         Printf(" +++ Averages:           %f evts/s, %f MB/s", (Double_t)fEventsProcessed / fProcTime, (Double_t)fBytesRead /1024./1024./fProcTime);
      }
      Printf(" +++ Total latency:      %f", fLatency);
      Printf(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
   }

   Int_t Compare(const TObject *o) const { TWrkInfo *wi = (TWrkInfo *)o;
                                           if (wi) {
                                             if (fStop < wi->fStop) {
                                                return -1;
                                             } else if (fStop == wi->fStop) {
                                                return 0;
                                             }
                                          }
                                          return 1; }
};


class TProofPerfAnalysis::TPackInfo : public TNamed {
public:
   TPackInfo(const char *ord, const char *host) :  TNamed(ord, host), fStart(0), fStop(-1), fSize(0), fMBRate(0.) { }
   TPackInfo(const char *ord, const char *host, Float_t start, Float_t stop, Long64_t sz, Double_t mbr)
            :  TNamed(ord, host), fStart(start), fStop(stop), fSize(sz), fMBRate(mbr) { }
   Float_t   fStart;            // When the packet has been assigned
   Float_t   fStop;             // When the packet has been finished
   Long64_t  fSize;             // Packet size
   Double_t  fMBRate;           // Processing rate MB/s
   void Print(Option_t *opt= "") const {
      if (!strcmp(opt, "S")) {
         Printf("       \t%10lld evts, \t%12.2f MB/s, \t%12.3f -> %12.3f s", fSize, fMBRate, fStart, fStop);
      } else {
         Printf("   %s:\t%s  \t%10lld evts, \t%12.2f MB/s, \t%12.3f -> %12.3f s", GetTitle(), GetName(), fSize, fMBRate, fStart, fStop);
      }
   }
};

class TProofPerfAnalysis::TWrkInfoFile : public TNamed {
public:
   TWrkInfoFile(const char *ord, const char *name) :  TNamed(ord, name) { }
   ~TWrkInfoFile() {fPackets.SetOwner(kFALSE); fPackets.Clear("nodelete");}
   TList     fPackets;          // Packest from this file processed by this worker
   void Print(Option_t *opt= "") const {
      if (!strcmp(opt, "R")) {
         Printf(" Worker: %s,\tpacket(s): %d", GetName(), fPackets.GetSize());
      } else {
         Printf(" Worker: %s,\t%d packet(s) from file: %s", GetName(), fPackets.GetSize(), GetTitle());
      }
      TIter nxp(&fPackets);
      TObject *o = 0;
      while ((o = nxp())) { o->Print("S"); }
   }
};

class TProofPerfAnalysis::TWrkEntry : public TObject {
public:
   TWrkEntry(Double_t xx, Double_t er, Double_t mbr, Double_t pt) : fXx(xx), fEvtRate(er), fMBRate(mbr), fProcTime(pt) { }
   Double_t fXx;             // Bin center
   Double_t fEvtRate;        // Event processing rate from this worker for this packet
   Double_t fMBRate;         // I/O processing rate from this worker for this packet
   Double_t fProcTime;       // Processing time
   void Print(Option_t * = "") const { Printf("%.4f \t%.3f evt/s \t%.3f MB/s \t%.3f s ", fXx, fEvtRate, fMBRate, fProcTime); }
};

////////////////////////////////////////////////////////////////////////////////

class TProofPerfAnalysis::TFileInfo : public TNamed {
public:
   TFileInfo(const char *name, const char *srv) :
      TNamed(name, srv), fPackets(0), fRPackets(0), fStart(0), fStop(-1),
      fSizeAvg(0), fSizeMax(-1.), fSizeMin(-1.),
      fMBRateAvg(0), fMBRateMax(-1.), fMBRateMin(-1.), fSizeP(0),
      fRateP(0), fRatePRemote(0), fMBRateP(0), fMBRatePRemote(0) { }
   virtual ~TFileInfo() {SafeDelete(fSizeP);
                         SafeDelete(fRateP); SafeDelete(fRatePRemote);
                         SafeDelete(fMBRateP); SafeDelete(fMBRatePRemote);
                         fPackList.SetOwner(kTRUE); fPackList.Clear();
                         fWrkList.SetOwner(kTRUE); fWrkList.Clear();
                         fRWrkList.SetOwner(kTRUE); fRWrkList.Clear();}

   Int_t     fPackets;          // Number of packets from this file
   Int_t     fRPackets;         // Number of different remote workers processing this file

   TList     fPackList;          // List of packet info
   TList     fWrkList;          // List of worker names processing this packet
   TList     fRWrkList;         // List of remote worker names processing this packet

   Float_t   fStart;            // When the first packet has been assigned
   Float_t   fStop;             // When the last packet has been finished

   Long64_t  fSizeAvg;          // Average Packet size
   Long64_t  fSizeMax;          // Max packet size
   Long64_t  fSizeMin;          // Min packet size

   Double_t  fMBRateAvg;        // Average MB rate
   Double_t  fMBRateMax;        // Max MB rate
   Double_t  fMBRateMin;        // Min MB rate

   TGraph   *fSizeP;             // Packet size vs packet (all)
   TGraph   *fRateP;             // Event processing rate vs packet (all)
   TGraph   *fRatePRemote;       // Event processing rate vs packet (remote workers)
   TGraph   *fMBRateP;           // Byte processing rate vs packet (all)
   TGraph   *fMBRatePRemote;     // Byte processing rate vs packet (remote workers)

   void Print(Option_t *opt = "") const {
      Printf(" +++ TFileInfo ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
      Printf(" +++ Server:         %s", GetTitle());
      Printf(" +++ File:           %s", GetName());
      Printf(" +++ Processing interval:  %f -> %f", fStart, fStop);
      Printf(" +++ Packets:         %d (%d remote)", fPackets, fRPackets);
      Printf(" +++ Processing wrks: %d (%d remote)", fWrkList.GetSize(), fRWrkList.GetSize());
      if (!strcmp(opt, "P")) fPackList.Print();
      if (!strcmp(opt, "WP")) fWrkList.Print("R");
      if (fPackets > 0) {
         Printf(" +++ MB rates:       %f MB/s (avg), %f MB/s (min), %f MB/s (max)",
                fMBRateAvg / fPackets, fMBRateMin, fMBRateMax);
         Printf(" +++ Sizes:          %lld  (avg), %lld (min), %lld (max)",
                fSizeAvg / fPackets, fSizeMin, fSizeMax);
      }
      Printf(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
   }

   Int_t Compare(const TObject *o) const { TFileInfo *wi = (TFileInfo *)o;
                                           if (wi) {
                                             if (fStop < wi->fStop) {
                                                return -1;
                                             } else if (fStop == wi->fStop) {
                                                return 0;
                                             }
                                          }
                                          return 1; }
};

Bool_t TProofPerfAnalysis::fgDebug = kTRUE;
////////////////////////////////////////////////////////////////////////////////
/// Constructor: open the file and attach to the tree

TProofPerfAnalysis::TProofPerfAnalysis(const char *perffile,
                               const char *title, const char *treename)
               : TNamed(perffile, title), fFile(0), fTreeName(treename), fTree(0),
                 fInitTime(-1.), fMergeTime(-1.), fMaxTime(-1.),
                 fEvents(0), fPackets(0),
                 fEvtRateMax(-1.), fMBRateMax(-1.), fLatencyMax(-1.),
                 fEvtRate(0), fEvtRateRun(0), fMBRate(0), fMBRateRun(0),
                 fEvtRateAvgMax(-1.), fMBRateAvgMax(-1.),
                 fEvtRateAvg(-1.), fMBRateAvg(0),
                 fFileResult(""), fSaveResult(kFALSE),
                 fDebug(0)
{
   // Use default title, if not specified
   if (!title) SetTitle("PROOF Performance Analysis");

   fTree = 0;
   fFile = TFile::Open(perffile);
   if (!fFile || (fFile && fFile->IsZombie())) {
      SafeDelete(fFile);
      Error("TProofPerfAnalysis", "problems opening file '%s'",
                              perffile ? perffile : "<undef>");
      SetBit(TObject::kInvalidObject);
      return;
   }

   // Set the subdirectory name, if any
   if (fTreeName.Contains("/")) {
      fDirName = gSystem->GetDirName(fTreeName);
      fTreeName = gSystem->BaseName(fTreeName);
   }

   // Adjust the name, if requested
   if (fTreeName.BeginsWith("+"))
      fTreeName.Replace(0, 1, "PROOF_PerfStats");

   // Point to the right TDirectory
   TDirectory *dir = fFile;
   if (!fDirName.IsNull()) {
      if (!(dir = dynamic_cast<TDirectory *>(fFile->Get(fDirName)))) {
         Error("TProofPerfAnalysis", "directory '%s' not found or not loadable", fDirName.Data());
         fFile->Close();
         SafeDelete(fFile);
         SetBit(TObject::kInvalidObject);
         return;
      }
   }

   // Load the performance tree
   LoadTree(dir);
   if (!fTree) {
      Error("TProofPerfAnalysis", "tree '%s' not found or not loadable", fTreeName.Data());
      fFile->Close();
      SafeDelete(fFile);
      SetBit(TObject::kInvalidObject);
      return;
   }
   if (fgDebug)
      Printf(" +++ TTree '%s' has %lld entries", fTreeName.Data(), fTree->GetEntries());

   // Init worker information
   FillWrkInfo();

   // Init file information
   FillFileInfo();

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor: open the file and attach to the tree

TProofPerfAnalysis::TProofPerfAnalysis(TTree *tree, const char *title)
               : TNamed("", title), fFile(0), fTree(0),
                 fInitTime(-1.), fMergeTime(-1.), fMaxTime(-1.),
                 fEvents(0), fPackets(0),
                 fEvtRateMax(-1.), fMBRateMax(-1.), fLatencyMax(-1.),
                 fEvtRate(0), fEvtRateRun(0), fMBRate(0), fMBRateRun(0),
                 fEvtRateAvgMax(-1.), fMBRateAvgMax(-1.),
                 fEvtRateAvg(-1.), fMBRateAvg(0),
                 fDebug(0)
{
   // The tree must be defined
   if (!tree) {
      SetBit(TObject::kInvalidObject);
      return;
   }

   // Use default title, if not specified
   if (!title) SetTitle("PROOF Performance Analysis");

   fTree = tree;
   fTreeName = fTree->GetName();
   SetName(TString::Format("heap_%s", fTreeName.Data()));

   // Adjust the name, if requested
   if (fTreeName.BeginsWith("+"))
      fTreeName.Replace(0, 1, "PROOF_PerfStats");

   if (fgDebug)
      Printf(" +++ TTree '%s' has %lld entries", fTreeName.Data(), fTree->GetEntries());

   // Init worker information
   FillWrkInfo();

   // Init file information
   FillFileInfo();

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor: detach the tree and close the file

TProofPerfAnalysis::~TProofPerfAnalysis()
{
   SafeDelete(fEvents);
   SafeDelete(fPackets);
   if (fFile) fFile->Close();
   SafeDelete(fFile);
}

////////////////////////////////////////////////////////////////////////////////
/// If defined, add '- <this title>' to the canvas title 't'

TString TProofPerfAnalysis::GetCanvasTitle(const char *t)
{
   if (fTitle.IsNull()) return TString(t);

   TString newt;
   if (t && strlen(t) > 0) {
      newt.Form("%s - %s", t, GetTitle());
   } else {
      newt = GetTitle();
   }
   // Done
   return newt;
}

////////////////////////////////////////////////////////////////////////////////
/// Load tree fTreeName from directory 'dir'. If not found, look for the
/// first TTree in the directory (and sub-directories) with the name containing
/// fTreeName.
/// The tree pointer is saved in fTree.

void TProofPerfAnalysis::LoadTree(TDirectory *dir)
{
   fTree = 0;
   if (!dir) return;

   // Try first the full name in the top directory
   if ((fTree = dynamic_cast<TTree *>(dir->Get(fTreeName)))) return;

   TRegexp re(fTreeName);
   // Now look inside: iter on the list of keys first
   TIter nxk(dir->GetListOfKeys());
   TKey *k = 0;
   while ((k = (TKey *) nxk())) {
      if (!strcmp(k->GetClassName(), "TDirectoryFile")) {
         TDirectory *kdir = (TDirectory *) dir->Get(k->GetName());
         LoadTree(kdir);
         if (fTree) return;
      } else if (!strcmp(k->GetClassName(), "TTree")) {
         TString tn(k->GetName());
         if (tn.Index(re) != kNPOS) {
            if ((fTree = dynamic_cast<TTree *>(dir->Get(tn)))) {
               fTreeName = tn;
               if (fgDebug) Printf(" +++ Found and loaded TTree '%s'", tn.Data());
               return;
            }
         }
      }
   }

   // Nothing found
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Analyse the file distribution. If writedet, underling details are
/// written out to a text file.

void TProofPerfAnalysis::FileDist(Bool_t writedet)
{
   if (!IsValid()) {
      Error("FileDist","not a valid instance - do nothing");
      return;
   }

   // Fill file info
   TList *wrkList = new TList;
   TList *srvList = new TList;
   GetWrkFileList(wrkList, srvList);
   Info("FileDist", "%d workers were active during this query", wrkList->GetSize());
   Info("FileDist", "%d servers were active during this query", srvList->GetSize());

   // Fill the worker-data server mapping
   TIter nxs(srvList);
   TIter nxw(wrkList);
   TNamed *sn = 0, *wn = 0;
   while ((sn = (TNamed *)nxs())) {
      nxw.Reset();
      while ((wn = (TNamed *) nxw())) {
         if (!strcmp(TUrl(sn->GetName()).GetHostFQDN(), wn->GetTitle())) {
            sn->SetTitle(wn->GetName());
         }
      }
   }

   // Reorder the lists following the title
   TList *nwl = new TList;
   TList *nsl = new TList;
   nxw.Reset();
   while ((wn = (TNamed *) nxw())) {
      TIter nnxw(nwl);
      TNamed *nwn = 0;
      while ((nwn = (TNamed *) nnxw())) {
         if (CompareOrd(wn->GetName(), nwn->GetName()) < 0) {
            nwl->AddBefore(nwn, wn);
            break;
         }
      }
      if (!nwn) nwl->Add(wn);
      // Find the server name, if any
      nxs.Reset();
      while ((sn = (TNamed *)nxs())) {
         if (!strcmp(sn->GetTitle(), wn->GetName())) {
            TIter nnxs(nsl);
            TNamed *nsn = 0;
            while ((nsn = (TNamed *) nnxs())) {
               if (CompareOrd(sn->GetTitle(), nsn->GetTitle()) < 0) {
                  nsl->AddBefore(nsn, sn);
                  break;
               }
            }
            if (!nsn) nsl->Add(sn);
            break;
         }
      }
      if (sn) srvList->Remove(sn);
   }
   // Add remaining servers at the end
   nxs.Reset();
   while ((sn = (TNamed *)nxs())) {
      nsl->Add(sn);
   }
   // Clean the orginal lists
   wrkList->SetOwner(kFALSE);
   srvList->SetOwner(kFALSE);
   delete wrkList;
   delete srvList;
   wrkList = nwl;
   srvList = nsl;

   // Notify
   wrkList->ls();
   srvList->ls();

   // Separate out the case with only one file server
   if (srvList->GetSize() == 1) {

      Printf("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ");
      Printf(" + Only one data server found: full analysis meaningful  + ");
      Printf(" + only when there are more file servers                 + ");
      Printf(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");


      // Create a 1D histo for cross packets
      TH1F *hxpak = new TH1F("hxpak", "MBytes / Worker",
                                    wrkList->GetSize(), 0., (Double_t)wrkList->GetSize());
      hxpak->SetDirectory(0);
      hxpak->SetMinimum(0.);
      hxpak->GetXaxis()->SetTitle("Worker");

      // Set the labels
      Int_t j = 1;
      TIter nxsw(wrkList);
      while ((wn = (TNamed *)nxsw())) {
         hxpak->GetXaxis()->SetBinLabel(j++, wn->GetName());
      }

      // Fill the histograms
      FillFileDistOneSrv(hxpak, writedet);

      // Display histos
      gStyle->SetOptStat(0);

      TCanvas *c2 = new TCanvas("cv-hxpak",  GetCanvasTitle(hxpak->GetTitle()), 800,350,700,700);
      c2->cd();
      DoDraw(hxpak);
      c2->Update();

   } else {
      // Create a 1D histo for file distribution
      TH1F *hfdis = new TH1F("hfdis", "Packet dist",
                           srvList->GetSize(), 0., (Double_t)srvList->GetSize());
      hfdis->SetDirectory(0);
      hfdis->SetMinimum(0);
      hfdis->GetXaxis()->SetTitle("Server");
      TH1F *hbdis = new TH1F("hbdis", "MBytes dist",
                           srvList->GetSize(), 0., (Double_t)srvList->GetSize());
      hbdis->SetDirectory(0);
      hbdis->SetMinimum(0);
      hbdis->GetXaxis()->SetTitle("Server");
      // Create a 2D histo for cross packets
      TH2F *hxpak = new TH2F("hxpak", "MBytes / {Worker,Server}",
                                    srvList->GetSize(), 0., (Double_t)srvList->GetSize(),
                                    wrkList->GetSize(), 0., (Double_t)wrkList->GetSize());
      hxpak->SetDirectory(0);
      hxpak->GetYaxis()->SetTitle("Worker");
      hxpak->GetXaxis()->SetTitle("Server");
      hxpak->GetXaxis()->SetTitleOffset(1.4);
      hxpak->GetYaxis()->SetTitleOffset(1.7);

      // Set the labels
      Int_t j = 1;
      TIter nxsw(wrkList);
      while ((wn = (TNamed *)nxsw())) {
         hxpak->GetYaxis()->SetBinLabel(j++, wn->GetName());
      }
      j = 1;
      TIter nxss(srvList);
      while ((sn = (TNamed *)nxss())) {
         hfdis->GetXaxis()->SetBinLabel(j, sn->GetName());
         hbdis->GetXaxis()->SetBinLabel(j, sn->GetName());
         hxpak->GetXaxis()->SetBinLabel(j++, sn->GetName());
      }

      // Fill the histograms
      FillFileDist(hfdis, hbdis, hxpak, writedet);

      j = 1;
      nxss.Reset();
      while ((sn = (TNamed *)nxss())) {
         TString lab(sn->GetName());
         lab = TUrl(sn->GetName()).GetHost();
         if (strcmp(sn->GetTitle(), "remote") && lab.Index(".") != kNPOS) lab.Remove(lab.Index("."));
         hfdis->GetXaxis()->SetBinLabel(j, lab);
         hbdis->GetXaxis()->SetBinLabel(j, lab);
         hxpak->GetXaxis()->SetBinLabel(j++, lab);
      }

      // Display histos
      gStyle->SetOptStat(0);

      TCanvas *c1 = new TCanvas("cv-hfdis", GetCanvasTitle(hfdis->GetTitle()), 800,50,700,700);
      c1->Divide(1,2);
      TPad *pad1 = (TPad *) c1->GetPad(1);
      TPad *pad2 = (TPad *) c1->GetPad(2);
      pad1->cd();
      DoDraw(hfdis);
      pad2->cd();
      DoDraw(hbdis);
      c1->Update();

      TCanvas *c2 = new TCanvas("cv-hxpak",  GetCanvasTitle(hxpak->GetTitle()), 500,350,700,700);
      c2->cd();
      DoDraw(hxpak, "lego");
      c2->Update();
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill file info

void TProofPerfAnalysis::GetWrkFileList(TList *wl, TList *sl)
{
   if (!wl || !sl) return;

   // Extract information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   fTree->SetBranchAddress("PerfEvents", &pep);
   Long64_t entries = fTree->GetEntries();
   TNamed *wn = 0, *sn = 0;
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType != TVirtualPerfStats::kPacket) continue;
      // Find out the worker instance
      TString wrk(TUrl(pe.fSlaveName.Data()).GetHostFQDN());
      wn = (TNamed *) wl->FindObject(pe.fSlave.Data());
      if (!wn) {
         wn = new TNamed(pe.fSlave.Data(), wrk.Data());
         wl->Add(wn);
      }
      // Find out the file server instance
      TUrl uf(pe.fFileName);
      TString srv(uf.GetUrl());
      Int_t ifn = srv.Index(uf.GetFile());
      if (ifn != kNPOS) srv.Remove(ifn);
      sn = (TNamed *) sl->FindObject(srv.Data());
      if (!sn) {
         sn = new TNamed(srv.Data(), "remote");
         sl->Add(sn);
      }
   }

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return -1 if ord1 comes before ord2, 0 i they are equal,
/// 1 if ord1 comes after ord2

Int_t TProofPerfAnalysis::CompareOrd(const char *ord1, const char *ord2)
{
   TString o1(ord1), o2(ord2), p1, p2;
   Int_t o1d = 0, o2d = 0;
   if ((o1d = o1.CountChar('.')) > (o2d = o2.CountChar('.'))) {
      return 1;
   } else if (o1d < o2d) {
      return -1;
   } else {
      o1.ReplaceAll(".", " ");
      o2.ReplaceAll(".", " ");
      Bool_t b1 = o1.Tokenize(p1, o1d, " ");
      Bool_t b2 = o2.Tokenize(p2, o2d, " ");
      while (b1 && b2) {
         if (p1.Atoi() > p2.Atoi()) {
            return 1;
         } else if (p1.Atoi() < p2.Atoi()) {
            return -1;
         } else {
            b1 = o1.Tokenize(p1, o1d, " ");
            b2 = o2.Tokenize(p2, o2d, " ");
         }
      }
      if (b1 && !b2) {
         return 1;
      } else if (b2 && !b1) {
         return -1;
      } else {
         return 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill file info

void TProofPerfAnalysis::FillFileDist(TH1F *hf, TH1F *hb, TH2F *hx, Bool_t wdet)
{
   if (!hf || !hb || !hx) return;

   TString fnout;
   FILE *fout = 0;
   if (wdet) {
      fnout.Form("%s-FileDist-Details.txt", GetName());
      if (!(fout = fopen(fnout.Data(), "w"))) {
         Warning("FillFileDist", "asked to save details in '%s' but file could"
                                 " not be open (errno: %d)", fnout.Data(), (int)errno);
      } else {
         Info("FillFileDist", "saving details to '%s'", fnout.Data());
      }
   }
   // Extract information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   fTree->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = fTree->GetEntries();
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType != TVirtualPerfStats::kPacket) continue;
      // Find out the labels ...
      TString wrk(pe.fSlave.Data());
      TUrl uf(pe.fFileName);
      TString srv(uf.GetUrl());
      Int_t ifn = srv.Index(uf.GetFile());
      if (ifn != kNPOS) srv.Remove(ifn);
      // ... and the bins
      Double_t xhf = hf->GetXaxis()->GetBinCenter(hf->GetXaxis()->FindBin(srv.Data()));
      Double_t xhx = hx->GetXaxis()->GetBinCenter(hx->GetXaxis()->FindBin(srv.Data()));
      Double_t yhx = hx->GetYaxis()->GetBinCenter(hx->GetYaxis()->FindBin(wrk.Data()));
      // Save details, if asked
      if (fout)
         fprintf(fout, "%s,%s -> %f,%f (%f)\n",
                       srv.Data(), wrk.Data(), xhx, yhx, pe.fBytesRead / 1024.);
      // Fill now
      hf->Fill(xhf);
      hb->Fill(xhf, pe.fBytesRead / 1024. / 1024.);
      hx->Fill(xhx, yhx, pe.fBytesRead / 1024. / 1024.);
   }
   if (fout) fclose(fout);
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill file info when there is only one file server

void TProofPerfAnalysis::FillFileDistOneSrv(TH1F *hx, Bool_t wdet)
{
   if (!hx) return;

   TString fnout;
   FILE *fout = 0;
   if (wdet) {
      fnout.Form("%s-FileDist-Details.txt", GetName());
      if (!(fout = fopen(fnout.Data(), "w"))) {
         Warning("FillFileDistOneSrv", "asked to save details in '%s' but file could"
                                       " not be open (errno: %d)", fnout.Data(), (int)errno);
      } else {
         Info("FillFileDistOneSrv", "saving details to '%s'", fnout.Data());
      }
   }
   // Extract information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   fTree->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = fTree->GetEntries();
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType != TVirtualPerfStats::kPacket) continue;
      // Find out the labels ...
      TString wrk(pe.fSlave.Data());
      TUrl uf(pe.fFileName);
      TString srv(uf.GetUrl());
      Int_t ifn = srv.Index(uf.GetFile());
      if (ifn != kNPOS) srv.Remove(ifn);
      // ... and the bins
      Double_t xhx = hx->GetXaxis()->GetBinCenter(hx->GetXaxis()->FindBin(wrk.Data()));
      // Save details, if asked
      if (fout)
         fprintf(fout, "%s,%s -> %f (%f)\n",
                       srv.Data(), wrk.Data(), xhx, pe.fBytesRead / 1024.);
      // Fill now
      hx->Fill(xhx, pe.fBytesRead / 1024. / 1024.);
   }
   if (fout) fclose(fout);
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Measure the worker activity

void TProofPerfAnalysis::WorkerActivity()
{
   if (!IsValid()) {
      Error("WorkerActivity","not a valid instance - do nothing");
      return;
   }

   // Fill basic worker info
   if (!WrkInfoOK()) FillWrkInfo();
   if (!WrkInfoOK()) {
      Error("WorkerActivity", "workers information not available - do nothing");
      return;
   }

   TObject *o = 0;
   // Create the histograms with activity vs time
   if ((o = gDirectory->FindObject("act10"))) delete o;
   Float_t t0 = fMergeTime - 2.* (fMaxTime - fMergeTime);
   Float_t t1 = 2.*fInitTime;
   if (t1 > t0) t1 = t0;
   TH1F *hact10 = new TH1F("act10", "Worker activity start (seconds)", 50, 0., t1);
   hact10->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("act11"))) delete o;
   TH1F *hact11 = new TH1F("act11", "Worker activity stop (seconds)", 50, t0, fMaxTime);
   hact11->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("act2"))) delete o;
   TH1F *hact2 = new TH1F("act2", "End of activity (seconds)", 50, t0, fMaxTime);
   hact2->GetXaxis()->SetTitle("Query Processing Time (s)");

   // Fine-tune stat printing
   Int_t curoptstat = gStyle->GetOptStat();
   gStyle->SetOptStat(1100);

   // Create the sorted list
   TIter nxw(&fWrksInfo);
   TWrkInfo *wi = 0;
   while ((wi = (TWrkInfo *)nxw())) {
      Int_t j = 0;
      for (j = 1; j < hact10->GetNbinsX()+1 ; j++) {
         if (wi->fStart < hact10->GetBinLowEdge(j))
            hact10->Fill(hact10->GetBinCenter(j));
      }
      for (j = 1; j < hact11->GetNbinsX()+1 ; j++) {
         if (wi->fStop > hact11->GetBinLowEdge(j))
            hact11->Fill(hact11->GetBinCenter(j));
      }
      hact2->Fill(wi->fStop);
   }

   // Display histos
   TCanvas *c1 = new TCanvas("perf", GetCanvasTitle("Activity histos"), 800,10,700,780);
   c1->Divide(1,2);
   TPad *pad1 = (TPad *) c1->GetPad(1);
   pad1->Divide(2,1);
   TPad *pad10 = (TPad *) pad1->GetPad(1);
   TPad *pad11 = (TPad *) pad1->GetPad(2);
   pad10->cd();
   DoDraw(hact10);
   pad11->cd();
   DoDraw(hact11);
   TPad *pad2 = (TPad *) c1->GetPad(2);
   pad2->cd();
   DoDraw(hact2);
   c1->cd();
   c1->Update();

   // Restore stat options
   gStyle->SetOptStat(curoptstat);

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Print information for all or the slowest showlast workers.
/// Use showlast < 0 to print all

void TProofPerfAnalysis::PrintWrkInfo(Int_t showlast)
{
   // Create the sorted list
   Int_t k = fWrksInfo.GetSize();
   TIter nxw(&fWrksInfo);
   TWrkInfo *wi = 0;
   while ((wi = (TWrkInfo *)nxw())) {
      // Print info about slowest workers
      k--;
      if (showlast < 0 || k < showlast) wi->Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print information for worker 'wn' (ordinal) or on the machine whose
/// ordinal or fqdn matches 'wn'. Multiple specifications separated by ','
/// or ' ' are supported, as well as wildcards '*', e.g. '0.2*,lxb10* lxf2323.doma.in"

void TProofPerfAnalysis::PrintWrkInfo(const char *wn)
{
   if (!wn || (wn && strlen(wn) <= 0)) {
      Error("PrintWrkInfo", "worker name or host must be defined!");
      return;
   }

   // Check exact name
   TWrkInfo *wi = (TWrkInfo *) fWrksInfo.FindObject(wn);
   if (wi) {
      wi->Print();
   } else {
      // Check matching
      TString ww(wn), w;
      TIter nxw(&fWrksInfo);
      while ((wi = (TWrkInfo *)nxw())) {
         TString n(wi->GetName()), t(wi->GetTitle());
         Ssiz_t from = 0;
         while (ww.Tokenize(w, from, "[, ]")) {
            TRegexp re(w, kTRUE);
            if (n.Index(re) != kNPOS || t.Index(re) != kNPOS) wi->Print();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print information for all or the slowest showlast workers.
/// Use showlast < 0 to print all

void TProofPerfAnalysis::PrintFileInfo(Int_t showlast, const char *opt, const char *out)
{
   RedirectHandle_t rh;
   if (out && strlen(out) > 0) gSystem->RedirectOutput(out, "w", &rh);

   // Create the sorted list
   Int_t k = fFilesInfo.GetSize();
   TIter nxf(&fFilesInfo);
   TFileInfo *fi = 0;
   while ((fi = (TFileInfo *)nxf())) {
      // Print info about files processed last
      k--;
      if (showlast < 0 || k < showlast) fi->Print(opt);
   }

   if (out && strlen(out) > 0) gSystem->RedirectOutput(0, 0, &rh);
}

////////////////////////////////////////////////////////////////////////////////
/// Print information for file 'fn' (path including directory) or server 'fn'.
/// Multiple specifications separated by ','
/// or ' ' are supported, as well as wildcards '*', e.g. 'pippo.root, h4mu*,lxb10*"

void TProofPerfAnalysis::PrintFileInfo(const char *fn, const char *opt, const char *out)
{
   if (!fn || (fn && strlen(fn) <= 0)) {
      Error("PrintFileInfo", "file path must be defined!");
      return;
   }

   RedirectHandle_t rh;
   if (out && strlen(out) > 0) gSystem->RedirectOutput(out, "w", &rh);

   // Check exact name
   TFileInfo *fi = (TFileInfo *) fFilesInfo.FindObject(fn);
   if (fi) {
      fi->Print(opt);
   } else {
      // Check matching
      TString fw(fn), f;
      TIter nxf(&fFilesInfo);
      while ((fi = (TFileInfo *)nxf())) {
         TString n(fi->GetName()), s(fi->GetTitle());
         Ssiz_t from = 0;
         while (fw.Tokenize(f, from, "[, ]")) {
            TRegexp re(f, kTRUE);
            if (n.Index(re) != kNPOS || s.Index(re) != kNPOS) fi->Print(opt);
         }
      }
   }

   if (out && strlen(out) > 0) gSystem->RedirectOutput(0, 0, &rh);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill basic worker info; if 'force' rescan the TTree even already done

void TProofPerfAnalysis::FillWrkInfo(Bool_t force)
{
   // Nothing to do if already called
   if (fWrksInfo.GetSize() > 0 && !force) return;

   // Cleanup existing information
   fWrksInfo.SetOwner(kTRUE);
   fWrksInfo.Clear();
   fInitTime = -1.;
   fMergeTime = -1.;
   fMaxTime = -1.;
   fEvtRateMax = -1.;
   fMBRateMax = -1.;
   fLatencyMax = -1.;

   TList *wl = new TList;
   // Extract worker information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   fTree->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = fTree->GetEntries();

   // First determine binning for global rates
   Int_t nraw = entries * 2, jj = 0, kk = 0;
   Double_t *xraw = new Double_t[nraw];
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType == TVirtualPerfStats::kPacket) {
         Float_t stop = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         Float_t start = stop - pe.fProcTime;
         // Bins
         xraw[jj++] = start;
         xraw[jj++] = stop;
      }
   }
   Int_t nbins = jj;
   Int_t *jidx = new Int_t[nbins];
   memset(jidx, 0, nbins * sizeof(Int_t));
   TMath::Sort(nbins, xraw, jidx, kFALSE);
   Double_t *xbins = new Double_t[nbins];
   jj = 0;
   for (kk = 0; kk < nbins; kk++) {
      Double_t xtmp = xraw[jidx[kk]];
      if (jj == 0 || xtmp > xbins[jj - 1] + .5) {
         xbins[jj] = xtmp;
         jj++;
      }
   }
   nbins = jj;
   delete [] xraw;
   delete [] jidx;

   // Create the global histograms
   Int_t nbin = nbins - 1;
   TObject *o = 0;
   if ((o = gDirectory->FindObject("gEvtRate"))) delete o;
   fEvtRate = new TH1F("gEvtRate", "Total event processing rate (evt/s)", nbin, xbins);
   fEvtRate->SetMinimum(0.);
   fEvtRate->SetStats(kFALSE);
   fEvtRate->SetFillColor(kCyan-8);
   fEvtRate->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("gEvtRateAvg"))) delete o;
   fEvtRateRun = new TH1F("gEvtRateAvg", "Event processing rate running average (evt/s)", nbin, xbins);
   fEvtRateRun->SetMinimum(0.);
   fEvtRateRun->SetStats(kFALSE);
   fEvtRateRun->SetLineColor(kBlue);
   fEvtRateRun->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("gMBRate"))) delete o;
   fMBRate = new TH1F("gMBRate", "Total processing rate (MB/s)", nbin, xbins);
   fMBRate->SetMinimum(0.);
   fMBRate->SetStats(kFALSE);
   fMBRate->SetFillColor(kCyan-8);
   fMBRate->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("gMBRateAvg"))) delete o;
   fMBRateRun = new TH1F("gMBRateAvg", "Processing rate running average (MB/s)", nbin, xbins);
   fMBRateRun->SetMinimum(0.);
   fMBRateRun->SetStats(kFALSE);
   fMBRateRun->SetLineColor(kBlue);
   fMBRateRun->GetXaxis()->SetTitle("Query Processing Time (s)");
   // Not needed any longer
   delete [] xbins;

   THashList gBins;
   TList *gwl = 0, *gbl = 0;

   // Extract the worker info now
   TWrkInfo *wi = 0;
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType == TVirtualPerfStats::kPacket) {
         // Find out the worker instance
         if (!(wi = (TWrkInfo *) wl->FindObject(pe.fSlave.Data()))) {
            wi = new TWrkInfo(pe.fSlave.Data(), pe.fSlaveName.Data());
            wl->Add(wi);
            wi->fRateT = new TGraph(100);
            wi->fRateRemoteT = new TGraph(100);
            wi->fMBRateT = new TGraph(100);
            wi->fMBRateRemoteT = new TGraph(100);
            wi->fLatencyT = new TGraph(100);
         }
         // Add Info now
         Float_t stop = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         Float_t start = stop - pe.fProcTime;
         if (wi->fPackets <= 0) {
            wi->fStart = start;
         } else {
            wi->fStop = stop;
         }
         TUrl uf(pe.fFileName), uw(pe.fSlaveName);
         fMaxTime = stop;
         wi->fEventsProcessed += pe.fEventsProcessed;
         wi->fBytesRead += pe.fBytesRead;
         wi->fLatency += pe.fLatency;
         wi->fProcTime += pe.fProcTime;
         wi->fCpuTime += pe.fCpuTime;
         // Fill graphs
         Double_t tt = stop;
         Double_t ert = pe.fEventsProcessed / pe.fProcTime ;
         Double_t brt = pe.fBytesRead / pe.fProcTime / 1024. / 1024. ;
         wi->fRateT->SetPoint(wi->fPackets, tt, ert);
         if (brt > 0.) wi->fMBRateT->SetPoint(wi->fPackets, tt, brt);
         wi->fLatencyT->SetPoint(wi->fPackets, tt, pe.fLatency);
         if (!pe.fFileName.IsNull() && strcmp(uf.GetHostFQDN(), uw.GetHostFQDN())) {
            wi->fRateRemoteT->SetPoint(wi->fRemotePackets, tt, ert);
            wi->fMBRateRemoteT->SetPoint(wi->fRemotePackets, tt, brt);
            wi->fRemotePackets++;
         }
         wi->fPackets++;
         if (ert > fEvtRateMax) fEvtRateMax = ert;
         if (brt > fMBRateMax) fMBRateMax = brt;
         if (pe.fLatency > fLatencyMax) fLatencyMax = pe.fLatency;

         // Fill global rate histos
         for (kk = 1; kk <= nbins; kk++) {
            Double_t mi = fEvtRate->GetBinLowEdge(kk);
            if (mi > stop) break;
            Double_t wd = fEvtRate->GetBinWidth(kk);
            Double_t mx = mi + wd;
            Double_t xx = fEvtRate->GetBinCenter(kk);
            // Overlap length
            Double_t olap = stop - mi;
            if (start > mi) olap = mx - start;
            if (olap >= 0) {
               TString sb = TString::Format("%d", kk);
               if (!(gbl = (TList *) gBins.FindObject(sb))) {
                  gbl = new TList;
                  gbl->SetName(sb);
                  gBins.Add(gbl);
               }
               if (!(gwl = (TList *) gbl->FindObject(pe.fSlave))) {
                  gwl = new TList;
                  gwl->SetName(pe.fSlave);
                  gbl->Add(gwl);
               }
               gwl->Add(new TWrkEntry(xx, ert, brt, pe.fProcTime));
            }
         }

         // Notify
         if (fDebug > 1) {
            if (pe.fProcTime > 0.) {
               Printf(" +++ %s #:%d at:%fs lat:%fs proc:%fs evts:%lld bytes:%lld (rates:%f evt/s, %f MB/s)",
                     wi->GetName(), wi->fPackets, fMaxTime - pe.fProcTime,
                     pe.fLatency, pe.fProcTime, pe.fEventsProcessed, pe.fBytesRead,
                     ert, brt);
            } else {
               Printf(" +++ %s #:%d at:%fs lat:%fs proc:%fs rate:-- evt/s (-- bytes/s)",
                     wi->GetName(), wi->fPackets, fMaxTime, pe.fLatency, pe.fProcTime);
            }
         }
      } else if (pe.fType == TVirtualPerfStats::kStart) {
         Float_t start = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         if (fDebug > 1) Printf(" +++ %s Start: %f s", pe.fEvtNode.Data(), start);
      } else if (pe.fType == TVirtualPerfStats::kStop) {
         Float_t stop = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         if (fDebug > 1) Printf(" +++ %s Stop: %f s", pe.fEvtNode.Data(), stop);
      } else {
         if (fDebug > 2) Printf(" +++ Event type: %d", pe.fType);
      }
   }

   TIter nxb(&gBins);
   gbl = 0;
   while ((gbl = (TList *) nxb())) {
      gwl = 0;
      TIter nxw(gbl);
      while ((gwl = (TList *) nxw())) {
         Double_t er = 0, br = 0, pt = 0, xx = 0;
         TIter nxp(gwl);
         TWrkEntry *we = 0;
         while ((we = (TWrkEntry *) nxp())) {
            if (we->fProcTime > 0) {
               er += we->fEvtRate * we->fProcTime;
               br += we->fMBRate * we->fProcTime;
               pt += we->fProcTime;
            }
            xx = we->fXx;
         }
         if (pt > 0.) {
            er /= pt;
            br /= pt;
            fEvtRate->Fill(xx, er);
            if (br > 0.) fMBRate->Fill(xx, br);
         }
      }
   }

   // Running averages
   Double_t er = 0, br = 0, pt = 0;
   for (kk = 1; kk < nbins; kk++) {
      Double_t wd = fEvtRate->GetBinWidth(kk);
      Double_t wx = fEvtRate->GetBinCenter(kk);
      Double_t wer = fEvtRate->GetBinContent(kk);
      Double_t wbr = fMBRate->GetBinContent(kk);

      if (kk == 1) {
         er = wer;
         br = wbr;
         pt = wd;
      } else {
         er *= pt;
         br *= pt;
         pt += wd;
         er += wer * wd;
         br += wbr * wd;
         er /= pt;
         br /= pt;
      }
      if (er > fEvtRateAvgMax) fEvtRateAvgMax = er;
      if (br > fMBRateAvgMax) fMBRateAvgMax = br;
      fEvtRateAvg = er;
      fMBRateAvg = br;
      // Fill
      fEvtRateRun->Fill(wx, er);
      fMBRateRun->Fill(wx, br);
   }


   // Final analysis to find relevant times
   TIter nxw(wl);
   while ((wi = (TWrkInfo *) nxw())) {
      fWrksInfo.Add(wi);
      if (wi->fStart > fInitTime) fInitTime = wi->fStart;
      // Resize the graphs
      wi->fRateT->Set(wi->fPackets);
      wi->fRateRemoteT->Set(wi->fRemotePackets);
      wi->fLatencyT->Set(wi->fPackets);
      wi->fMBRateT->Set(wi->fPackets);
      wi->fMBRateRemoteT->Set(wi->fRemotePackets);
   }
   wl->SetOwner(kFALSE);
   delete wl;

   // Final analysis to find relevant times
   fMergeTime = fMaxTime;
   Int_t rsw = (fWrksInfo.GetSize() > 1) ? 2 : 1, ksw = 0;
   TIter nxsw(&fWrksInfo);
   while ((wi = (TWrkInfo *) nxsw())) {
      if (wi->fStop > 0.) ksw++;
      if (ksw == rsw) break;
   }
   if (wi) fMergeTime = wi->fStop;

   // (Re-)create the event and packet distribution histograms
   SafeDelete(fEvents);
   SafeDelete(fPackets);
   fEvents = new TH1F("hevents", "Events per worker", fWrksInfo.GetSize(), -.5, fWrksInfo.GetSize()-.5);
   fEvents->SetDirectory(0);
   fPackets = new TH1F("hpackets", "Packets per worker", fWrksInfo.GetSize(), -.5, fWrksInfo.GetSize()-.5);
   fPackets->SetDirectory(0);
   Int_t j = 0;
   TIter nxwi(&fWrksInfo);
   while ((wi = (TWrkInfo *)nxwi())) {
      fEvents->GetXaxis()->SetBinLabel(j+1, wi->GetName());
      fEvents->Fill(j, wi->fEventsProcessed);
      fPackets->GetXaxis()->SetBinLabel(j+1, wi->GetName());
      fPackets->Fill(j++, wi->fPackets);
   }
   fEvents->SetMinimum(0.);
   fPackets->SetMinimum(0.);
   fEvents->SetFillColor(38);
   fPackets->SetFillColor(38);
   fEvents->GetYaxis()->SetTitle("Events");
   fEvents->GetXaxis()->SetTitle("Worker");
   fPackets->GetYaxis()->SetTitle("Packets");
   fPackets->GetXaxis()->SetTitle("Worker");

   // Print summary
   if (fgDebug) Summary();
}

////////////////////////////////////////////////////////////////////////////////
/// Print summary of query. Use opt = 'S' for compact version.
/// Output to 'out' or to screen.

void TProofPerfAnalysis::Summary(Option_t *opt, const char *out)
{
   TString o(out);
   RedirectHandle_t rh;
   if (!o.IsNull()) {
      const char *m = (o.BeginsWith("+")) ? "a" : "w";
      o.Remove(TString::kLeading, '+');
      gSystem->RedirectOutput(o, m, &rh);
   }

   // Print summary
   if (!strcmp(opt, "S")) {
      // Short version
      Printf("%d %f %f %f %f %f %f %f",
              fWrksInfo.GetSize(), fMaxTime, fInitTime, fMaxTime - fMergeTime,
              fEvtRateAvg, fEvtRateAvgMax, fMBRateAvg, fMBRateAvgMax);
   } else {
      // Long version
      Printf(" +++ %d workers were active during this query", fWrksInfo.GetSize());
      Printf(" +++ Total query time: %f secs (init: %f secs, merge: %f secs)",
             fMaxTime, fInitTime, fMaxTime - fMergeTime);
      Printf(" +++ Avg processing rates: %.4f evts/s, %.4f MB/s", fEvtRateAvg, fMBRateAvg);
      Printf(" +++ Max processing rates: %.4f evts/s, %.4f MB/s", fEvtRateAvgMax, fMBRateAvgMax);
   }

   if (!o.IsNull()) gSystem->RedirectOutput(0, 0, &rh);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill basic worker info; if 'force' rescan the TTree even already done

void TProofPerfAnalysis::FillFileInfo(Bool_t force)
{
   // Nothing to do if already called
   if (fFilesInfo.GetSize() > 0 && !force) return;

   // Cleanup existing information
   fFilesInfo.SetOwner(kTRUE);
   fFilesInfo.Clear();

   TList *fl = new TList;
   // Extract worker information
   TPerfEvent pe;
   TPerfEvent* pep = &pe;
   fTree->SetBranchAddress("PerfEvents",&pep);
   Long64_t entries = fTree->GetEntries();
   TFileInfo *fi = 0;
   for (Long64_t k=0; k<entries; k++) {
      fTree->GetEntry(k);
      // Analyse only packets
      if (pe.fType == TVirtualPerfStats::kPacket) {
         TUrl uf(pe.fFileName);
         TString srv(uf.GetUrl());
         Int_t ifn = srv.Index(uf.GetFile());
         if (ifn != kNPOS) srv.Remove(ifn);
         // Find out the file instance
         fi = (TFileInfo *) fl->FindObject(uf.GetFile());
         if (!fi) {
            fi = new TFileInfo(uf.GetFile(), srv.Data());
            fl->Add(fi);
            fi->fSizeP = new TGraph(10);
            fi->fRateP = new TGraph(10);
            fi->fRatePRemote = new TGraph(10);
            fi->fMBRateP = new TGraph(10);
            fi->fMBRatePRemote = new TGraph(10);
         }
         // Add Info now
         Float_t stop = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         Float_t start = stop - pe.fProcTime;
         if (fi->fPackets <= 0) {
            fi->fStart = start;
         } else {
            fi->fStop = stop;
         }
         TUrl uw(pe.fSlaveName);

         // Fill size graphs
         fi->fSizeP->SetPoint(fi->fPackets, (Double_t) fi->fPackets, (Double_t) pe.fEventsProcessed);
         fi->fSizeAvg += pe.fEventsProcessed;
         if (pe.fEventsProcessed > fi->fSizeMax || fi->fSizeMax < 0.) fi->fSizeMax = pe.fEventsProcessed;
         if (pe.fEventsProcessed < fi->fSizeMin || fi->fSizeMin < 0.) fi->fSizeMin = pe.fEventsProcessed;

         // Fill rate graphs
         Double_t tt = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         Double_t ert = pe.fEventsProcessed / pe.fProcTime ;
         Double_t brt = pe.fBytesRead / pe.fProcTime / 1024. / 1024. ;
         fi->fRateP->SetPoint(fi->fPackets, tt, ert);
         if (brt > 0.) fi->fMBRateP->SetPoint(fi->fPackets, tt, brt);
         if (!pe.fFileName.IsNull() && strcmp(uf.GetHostFQDN(), uw.GetHostFQDN())) {
            if (!(fi->fRWrkList.FindObject(pe.fSlave))) fi->fRWrkList.Add(new TNamed(pe.fSlave, pe.fSlaveName));
            fi->fRatePRemote->SetPoint(fi->fRPackets, tt, ert);
            fi->fMBRatePRemote->SetPoint(fi->fRPackets, tt, brt);
            fi->fRPackets++;
         }
         fi->fPackets++;
         if (brt > 0) {
            fi->fMBRateAvg += brt;
            if (brt > fi->fMBRateMax || fi->fMBRateMax < 0.) fi->fMBRateMax = brt;
            if (brt < fi->fMBRateMin || fi->fMBRateMin < 0.) fi->fMBRateMin = brt;
         }

         // Packet info
         TPackInfo *pi = new TPackInfo(pe.fSlave, pe.fSlaveName, start, stop, pe.fEventsProcessed, brt);
         fi->fPackList.Add(pi);
         TWrkInfoFile *wif = 0;
         if (!(wif = (TWrkInfoFile *) fi->fWrkList.FindObject(pe.fSlave))) {
            wif = new TWrkInfoFile(pe.fSlave, uf.GetFile());
            fi->fWrkList.Add(wif);
         }
         wif->fPackets.Add(pi);

         // Notify
         if (fDebug > 1) {
            if (pe.fProcTime > 0.) {
               Printf(" +++ %s #:%d at:%fs lat:%fs proc:%fs evts:%lld bytes:%lld (rates:%f evt/s, %f MB/s)",
                     fi->GetName(), fi->fPackets, fMaxTime - pe.fProcTime,
                     pe.fLatency, pe.fProcTime, pe.fEventsProcessed, pe.fBytesRead,
                     ert, brt);
            } else {
               Printf(" +++ %s #:%d at:%fs lat:%fs proc:%fs rate:-- evt/s (-- bytes/s)",
                     fi->GetName(), fi->fPackets, fMaxTime, pe.fLatency, pe.fProcTime);
            }
         }
      } else if (pe.fType == TVirtualPerfStats::kStart) {
         Float_t start = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         if (fDebug > 1) Printf(" +++ %s Start: %f s", pe.fEvtNode.Data(), start);
      } else if (pe.fType == TVirtualPerfStats::kStop) {
         Float_t stop = pe.fTimeStamp.GetSec() + 1e-9*pe.fTimeStamp.GetNanoSec();
         if (fDebug > 1) Printf(" +++ %s Stop: %f s", pe.fEvtNode.Data(), stop);
      } else {
         if (fDebug > 2) Printf(" +++ Event type: %d", pe.fType);
      }
   }
   // Final analysis to find relevant times
   TIter nxf(fl);
   while ((fi = (TFileInfo *) nxf())) {
      fFilesInfo.Add(fi);
      // Resize the graphs
      fi->fRateP->Set(fi->fPackets);
      fi->fRatePRemote->Set(fi->fRPackets);
      fi->fMBRateP->Set(fi->fPackets);
      fi->fMBRatePRemote->Set(fi->fRPackets);
   }
   fl->SetOwner(kFALSE);
   delete fl;

   // Print summary
   if (fgDebug)
      Printf(" +++ %d files were processed during this query", fFilesInfo.GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Static setter for the verbosity level

void TProofPerfAnalysis::SetDebug(Int_t d)
{
   fDebug = d;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw object 'o' with options 'opt'
/// Save it with 'name' if in saving mode (see SetSaveResult)

void TProofPerfAnalysis::DoDraw(TObject *o, Option_t *opt, const char *name)
{
   // Draw
   o->Draw(opt);

   // Save the result
   if (fSaveResult) {
      // Preparation is done in SetSaveResult, here we just update
      TDirectory *curdir = gDirectory;
      TFile *f = TFile::Open(fFileResult, "UPDATE");
      if (f && !f->IsZombie()) {
         const char *n = (name && strlen(name) > 0) ? name : 0;
         o->Write(n);
         f->Close();
      }
      if (f) delete f;
      gDirectory = curdir;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set save result mode and validate 'file' according to 'mode'.
/// Return 0 on success, -1 if any problem with the file is encountered
/// (save result mode is not enabled in such a case).
/// If 'file' is null saving is disabled.

Int_t TProofPerfAnalysis::SetSaveResult(const char *file, Option_t *mode)
{
   // A null 'file' indicates the will to disable
   if (!file) {
      fFileResult = "";
      fSaveResult = kFALSE;
      // Notify
      Printf("Drawn objects saving disabled");
      return 0;
   }

   // Check if there is a change
   if (!fFileResult.IsNull() && fFileResult == file) {
      // No change
      fSaveResult = kTRUE;
      return 0;
   }
   // New or changed file: validate
   fFileResult = "";
   fSaveResult = kFALSE;
   TDirectory *curdir = gDirectory;
   TFile *f = TFile::Open(file, mode);
   if (!f || f->IsZombie()) {
      if (f) delete f;
      fFileResult = "";
      Error("SetSaveResult", "could not open file '%s' in mode '%s'", file, mode);
      gDirectory = curdir;
      return -1;
   }
   f->Close();
   delete f;
   gDirectory = curdir;
   // Ok
   fFileResult = file;
   fSaveResult = kTRUE;
   // Notify
   Printf("Drawn objects will be saved in file '%s'", file);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static setter for the verbosity level

void TProofPerfAnalysis::SetgDebug(Bool_t on)
{
   fgDebug = on;
}
////////////////////////////////////////////////////////////////////////////////
/// Display event and packet distribution

void TProofPerfAnalysis::EventDist()
{
   if (!fEvents || !fPackets) {
      Error("EventDist", "distributions not initialized - do nothing");
   }

   // Display histos
   TCanvas *c1 = new TCanvas("evtdist", GetCanvasTitle("Event distributions"),800,10,700,780);
   c1->Divide(1,2);
   TPad *pad1 = (TPad *) c1->GetPad(1);
   pad1->cd();
   fEvents->SetStats(kFALSE);
   DoDraw(fEvents);
   TPad *pad2 = (TPad *) c1->GetPad(2);
   pad2->cd();
   fPackets->SetStats(kFALSE);
   DoDraw(fPackets);
   c1->cd();
   c1->Update();

}

////////////////////////////////////////////////////////////////////////////////
/// Show event processing or MB processing rate plot vs time

void TProofPerfAnalysis::RatePlot(const char *wrks)
{
   Bool_t global = (wrks && !strcmp(wrks, "global")) ? kTRUE : kFALSE;

   TH1F *hrt1 = 0, *hrt2 = 0;
   if (global) {
      hrt1 = fEvtRate;
      hrt2 = fMBRate;
   } else {
      // Create the histograms
      TObject *o = 0;
      if ((o = gDirectory->FindObject("rt1"))) delete o;
      hrt1 = new TH1F("rt1", "Evt processing rate (evt/s)", 100, 0., fMaxTime);
      hrt1->SetMinimum(0.);
      hrt1->SetMaximum(1.05*fEvtRateMax);
      hrt1->SetStats(kFALSE);
      hrt1->GetXaxis()->SetTitle("Query Processing Time (s)");
      if ((o = gDirectory->FindObject("rt2"))) delete o;
      hrt2 = new TH1F("rt2", "MB processing rate (MB/s)", 100, 0., fMaxTime);
      hrt2->SetMinimum(0.);
      hrt2->SetMaximum(1.05*fMBRateMax);
      hrt2->SetStats(kFALSE);
      hrt2->GetXaxis()->SetTitle("Query Processing Time (s)");
   }

   // Display histo frames
   TCanvas *c1 = new TCanvas("rates", GetCanvasTitle("Processing rates"), 800,10,700,780);
   c1->Divide(1,2);
   TPad *pad1 = (TPad *) c1->GetPad(1);
   pad1->cd();
   hrt1->Draw();
   if (global) DoDraw(fEvtRateRun, "SAME", "EvtRateRun");
   TPad *pad2 = (TPad *) c1->GetPad(2);
   pad2->cd();
   hrt2->Draw();
   if (global) DoDraw(fMBRateRun, "SAME", "MBRateRun");
   c1->cd();
   c1->Update();

   // Done if global
   if (global) return;

   // Which workers?
   THashList *wl = 0;
   TString ww(wrks);
   if (!ww.IsNull() && ww != "*" && ww != "all") {
      TString w;
      Ssiz_t from = 0;
      while ((ww.Tokenize(w, from, ","))) {
         if (!wl) wl = new THashList();
         wl->Add(new TObjString(w.Data()));
      }
   }

   // Now plot the graphs per worker
   Int_t ci = 40, cir = 30, ic = 0;
   TIter nxw(&fWrksInfo);
   TWrkInfo *wi = 0;
   while ((wi = (TWrkInfo *) nxw())) {
      if (wl && !wl->FindObject(wi->GetName())) continue;
      if (wi->fRateT && wi->fRateT->GetN() > 0) {
         wi->fRateT->SetNameTitle(wi->GetName(), wi->GetTitle());
         pad1->cd();
         wi->fRateT->SetLineColor(ci);
         DoDraw(wi->fRateT, "L", TString::Format("RateT-%s", wi->fRateT->GetName()));
      }
      if (wi->fRateRemoteT && wi->fRateRemoteT->GetN() > 0) {
         wi->fRateRemoteT->SetNameTitle(wi->GetName(), wi->GetTitle());
         pad1->cd();
         wi->fRateRemoteT->SetLineColor(cir);
         DoDraw(wi->fRateRemoteT, "L", TString::Format("RateRemoteT-%s", wi->fRateRemoteT->GetName()));
      }
      if (wi->fMBRateT && wi->fMBRateT->GetN() > 0) {
         wi->fMBRateT->SetNameTitle(wi->GetName(), wi->GetTitle());
         pad2->cd();
         wi->fMBRateT->SetLineColor(ci);
         DoDraw(wi->fMBRateT, "L", TString::Format("MBRateT-%s", wi->fMBRateT->GetName()));
      }
      if (wi->fMBRateRemoteT && wi->fMBRateRemoteT->GetN() > 0) {
         wi->fMBRateRemoteT->SetNameTitle(wi->GetName(), wi->GetTitle());
         pad2->cd();
         wi->fMBRateRemoteT->SetLineColor(cir);
         DoDraw(wi->fMBRateRemoteT, "L", TString::Format("MBRateRemoteT-%s", wi->fMBRateRemoteT->GetName()));
      }
      ic++;
      ci = ic%10 + 40;
      cir = ic%10 + 30;
      c1->cd();
      c1->Update();
   }

   // Cleanup
   if (wl) {
      wl->SetOwner(kTRUE);
      delete wl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show event processing or MB processing rate plot vs time
/// Create the histograms

void TProofPerfAnalysis::LatencyPlot(const char *wrks)
{
   TObject *o = 0;
   if ((o = gDirectory->FindObject("lt1"))) delete o;
   TH1F *hlt1 = new TH1F("lt1", "Packet retrieval latency", 100, 0., fMaxTime);
   hlt1->SetMinimum(0.);
   hlt1->SetMaximum(1.05*fLatencyMax);
   hlt1->SetStats(kFALSE);
   hlt1->GetXaxis()->SetTitle("Query Processing Time (s)");
   hlt1->GetYaxis()->SetTitle("Latency (s)");

   // Display histo frames
   TCanvas *c1 = new TCanvas("latency", GetCanvasTitle("Packet Retrieval Latency"), 800,10,700,780);
   hlt1->Draw();
   c1->cd();
   c1->Update();

   // Which workers?
   THashList *wl = 0;
   TString ww(wrks);
   if (!ww.IsNull() && ww != "*" && ww != "all") {
      TString w;
      Ssiz_t from = 0;
      while ((ww.Tokenize(w, from, ","))) {
         if (!wl) wl = new THashList();
         wl->Add(new TObjString(w.Data()));
      }
   }

   // Now plot the graphs per worker
   Int_t ci = 40, ic = 0;
   TIter nxw(&fWrksInfo);
   TWrkInfo *wi = 0;
   while ((wi = (TWrkInfo *) nxw())) {
      if (wl && !wl->FindObject(wi->GetName())) continue;
      if (wi->fLatencyT) {
         wi->fLatencyT->SetNameTitle(wi->GetName(), wi->GetTitle());
         wi->fLatencyT->SetLineColor(ci);
         DoDraw(wi->fLatencyT, "L", TString::Format("LatencyT-%s", wi->fLatencyT->GetName()));
      }
      ic++;
      ci = ic%10 + 40;
      c1->cd();
      c1->Update();
   }

   // Cleanup
   if (wl) {
      wl->SetOwner(kTRUE);
      delete wl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show event processing or MB processing rate plot vs time

void TProofPerfAnalysis::FileProcPlot(const char *fn, const char *out)
{
   if (!fn || strlen(fn) <= 0) {
      Error("FileRatePlot", "file name is mandatory!");
      return;
   }
   // Get the file info object
   TFileInfo *fi = (TFileInfo *) fFilesInfo.FindObject(fn);
   if (!fi) {
      Error("FileRatePlot", "TFileInfo object for '%s' not found!", fn);
      return;
   }

   // Output text file, if required
   FILE *fo = stdout;
   if (out && strlen(out) > 0) {
      if (!(fo = fopen(out, "w"))) {
         Warning("FileRatePlot", "problems creating '%s': logging to stdout", out);
         fo = stdout;
      } else {
         Printf(" Details logged to %s", out);
      }
   }

   // Get bins
   Int_t nbins = fi->fPackList.GetSize() * 2;
   Double_t *xraw = new Double_t[nbins];
   Int_t jj = 0;
   TPackInfo *pi = 0;
   TIter nxp(&(fi->fPackList));
   while ((pi = (TPackInfo *) nxp())) {
      // Bins
      xraw[jj++] = pi->fStart;
      xraw[jj++] = pi->fStop;
   }
   Int_t *jidx = new Int_t[nbins];
   memset(jidx, 0, nbins * sizeof(Int_t));
   TMath::Sort(nbins, xraw, jidx, kFALSE);
   Double_t *xbins = new Double_t[nbins];
   Int_t kk =0;
   for (kk = 0; kk < nbins; kk++) {
      xbins[kk] = xraw[jidx[kk]];
   }
   delete [] xraw;
   delete [] jidx;

   // Create the histograms
   Int_t nbin = nbins - 1;
   TObject *o = 0;
   if ((o = gDirectory->FindObject("rt1"))) delete o;
   TH1F *hrt1 = new TH1F("rt1", "Total processing rate (MB/s)", nbins - 1, xbins);
   hrt1->SetMinimum(0.);
   hrt1->SetStats(kFALSE);
   hrt1->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("rt2"))) delete o;
   TH1F *hrt2 = new TH1F("rt2", "Number of processing workers", nbins - 1, xbins);
   hrt2->SetMinimum(0.);
   hrt2->SetMaximum(1.2*fWrksInfo.GetSize());
   hrt2->SetStats(kFALSE);
   hrt2->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("rt3"))) delete o;
   TH1F *hrt3 = new TH1F("rt3", "Total processing events", nbins - 1, xbins);
   hrt3->SetMinimum(0.);
   hrt3->SetStats(kFALSE);
   hrt3->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("rt4"))) delete o;
   TH1F *hrt4 = new TH1F("rt4", "Weighted processing rate (MB/s)", nbins - 1, xbins);
   hrt4->SetMinimum(0.);
   hrt4->SetStats(kFALSE);
   hrt4->GetXaxis()->SetTitle("Query Processing Time (s)");
   // Not needed any longer
   delete [] xbins;

   // Fill histos now
   Int_t ii = 0;
   for (ii = 1; ii <= nbin; ii++) {
      Double_t mi = hrt1->GetBinLowEdge(ii);
      Double_t wd = hrt1->GetBinWidth(ii);
      Double_t mx = mi + wd;
      Double_t xx = hrt1->GetBinCenter(ii);
      fprintf(fo, " Bin: %d/%d [%f, %f]\n", ii, nbin, mi, mx);
      pi = 0;
      kk = 0;
      nxp.Reset();
      while ((pi = (TPackInfo *) nxp())) {
         // Overlap length
         Double_t olap = pi->fStop - mi;
         if (pi->fStart > mi) olap = mx - pi->fStart;
         if (olap >= 0) {
            hrt1->Fill(xx, pi->fMBRate);
            hrt2->Fill(xx, 1.);
            hrt3->Fill(xx, pi->fSize);
            hrt4->Fill(xx, pi->fMBRate * pi->fSize);
            fprintf(fo, "    %d: %s \t%lld \tevts \t%f \tMB/s\n", kk++, pi->GetName(), pi->fSize, pi->fMBRate);
         }
      }
   }
   if (fo != stdout) fclose(fo);

   // Display histo frames
   TCanvas *c1 = new TCanvas("rates", GetCanvasTitle("File processing info"), 800,10,700,780);
   c1->Divide(1,3);
   TPad *pad1 = (TPad *) c1->GetPad(1);
   pad1->cd();
   DoDraw(hrt1);
   TPad *pad2 = (TPad *) c1->GetPad(2);
   pad2->cd();
   DoDraw(hrt2);
   TPad *pad4 = (TPad *) c1->GetPad(3);
   pad4->cd();
   hrt4->Divide(hrt3);
   DoDraw(hrt4);
   c1->cd();
   c1->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Show MB processing rate plot per file vs time

void TProofPerfAnalysis::FileRatePlot(const char *fns)
{
   // Create the histograms
   TObject *o = 0;
   if ((o = gDirectory->FindObject("rt1"))) delete o;
   TH1F *hrt1 = new TH1F("rt1", "Event processing rate per packet (evt/s)", 100, 0., fMaxTime);
   hrt1->SetMinimum(0.);
   hrt1->SetMaximum(1.05*fEvtRateMax);
   hrt1->SetStats(kFALSE);
   hrt1->GetXaxis()->SetTitle("Query Processing Time (s)");
   if ((o = gDirectory->FindObject("rt2"))) delete o;
   TH1F *hrt2 = new TH1F("rt2", "I/O processing rate per packet (MB/s)", 100, 0., fMaxTime);
   hrt2->SetMinimum(0.);
   hrt2->SetMaximum(1.05*fMBRateMax);
   hrt2->SetStats(kFALSE);
   hrt2->GetXaxis()->SetTitle("Query Processing Time (s)");

   // Display histo frames
   TCanvas *c1 = new TCanvas("rates", GetCanvasTitle("Processing rates"), 800,10,700,780);
   c1->Divide(1,2);
   TPad *pad1 = (TPad *) c1->GetPad(1);
   pad1->cd();
   hrt1->Draw();
   TPad *pad2 = (TPad *) c1->GetPad(2);
   pad2->cd();
   hrt2->Draw();
   c1->cd();
   c1->Update();

   // Which workers?
   THashList *fl = 0;
   TString fw(fns);
   if (!fw.IsNull() && fw != "*" && fw != "all") {
      TString w;
      Ssiz_t from = 0;
      while ((fw.Tokenize(w, from, ","))) {
         if (!fl) fl = new THashList();
         fl->Add(new TObjString(w.Data()));
      }
   }

   // Now plot the graphs per worker
   Int_t ci = 40, cir = 30, ic = 0;
   TIter nxf(&fFilesInfo);
   TFileInfo *fi = 0;
   while ((fi = (TFileInfo *) nxf())) {
      if (fl && !fl->FindObject(fi->GetName())) continue;
      if (fi->fRateP && fi->fRateP->GetN() > 0) {
         fi->fRateP->SetNameTitle(fi->GetName(), fi->GetTitle());
         pad1->cd();
         fi->fRateP->SetLineColor(ci);
         DoDraw(fi->fRateP, "L", TString::Format("RateP-%d", ic));
      }
      if (fi->fRatePRemote && fi->fRatePRemote->GetN() > 0) {
         fi->fRatePRemote->SetNameTitle(fi->GetName(), fi->GetTitle());
         pad1->cd();
         fi->fRatePRemote->SetLineColor(cir);
         DoDraw(fi->fRatePRemote, "L", TString::Format("RatePRemote-%d", ic));
      }
      if (fi->fMBRateP && fi->fMBRateP->GetN() > 0) {
         fi->fMBRateP->SetNameTitle(fi->GetName(), fi->GetTitle());
         pad2->cd();
         fi->fMBRateP->SetLineColor(ci);
         DoDraw(fi->fMBRateP, "L", TString::Format("MBRateP-%d", ic));
      }
      if (fi->fMBRatePRemote && fi->fMBRatePRemote->GetN() > 0) {
         fi->fMBRatePRemote->SetNameTitle(fi->GetName(), fi->GetTitle());
         pad2->cd();
         fi->fMBRatePRemote->SetLineColor(cir);
         DoDraw(fi->fMBRatePRemote, "L", TString::Format("MBRatePRemote-%d", ic));
      }
      ic++;
      ci = ic%10 + 40;
      cir = ic%10 + 30;
      c1->cd();
      c1->Update();
   }

   // Cleanup
   if (fl) {
      fl->SetOwner(kTRUE);
      delete fl;
   }
}
