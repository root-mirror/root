/// \file
/// \ingroup tutorial_http
///  This program creates :
///    - a one dimensional histogram
///    - a two dimensional histogram
///    - a profile histogram
///    - a memory-resident ntuple
///
///  These objects are filled with some random numbers and saved on a in-memory file.
///  All objects can be seen in web browser is open url:
/// ~~~
///      http://localhost:8080
/// ~~~
///
/// \macro_code
///
/// \author Sergey Linev

#include <TFile.h>
#include <TMemFile.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <THttpServer.h>

void httpserver(const char* jobname = "job1", Long64_t maxcnt = 0)
{
   TString filename = Form("%s.root", jobname);
   TFile *hfile = new TMemFile(filename,"RECREATE","Demo ROOT file with histograms");

   // Create some histograms, a profile histogram and an ntuple
   TH1F *hpx = new TH1F("hpx","This is the px distribution",100,-4,4);
   hpx->SetFillColor(48);
   TH2F *hpxpy = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
   TProfile *hprof = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
   TNtuple *ntuple = new TNtuple("ntuple","Demo ntuple","px:py:pz:random:i");
   hfile->Write();


   // http server with port 8080, use jobname as top-folder name
   THttpServer* serv = new THttpServer(Form("http:8080?top=%s", jobname));

   // fastcgi server with port 9000, use jobname as top-folder name
   // THttpServer* serv = new THttpServer(Form("fastcgi:9000?top=%s_fastcgi", jobname));

   // dabc agent, connects to DABC master_host:1237, works only when DABC configured
   // THttpServer* serv = new THttpServer(Form("dabc:master_host:1237?top=%s_dabc", jobname));

   // when read-only mode disabled one could execute object methods like TTree::Draw()
   serv->SetReadOnly(kFALSE);

   // One could specify location of newer version of JSROOT
   // serv->SetJSROOT("https://root.cern/js/latest/");
   // serv->SetJSROOT("http://jsroot.gsi.de/latest/");

   gBenchmark->Start(jobname);

   // Create a new canvas.
   TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",200,10,700,500);
   c1->SetFillColor(42);
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(6);
   c1->GetFrame()->SetBorderMode(-1);


   // Fill histograms randomly
   TRandom3 random;
   Float_t px, py, pz;
   const Int_t kUPDATE = 1000;
   Long64_t i = 0;

   while (true) {
      random.Rannor(px,py);
      pz = px*px + py*py;
      Float_t rnd = random.Rndm(1);
      hpx->Fill(px);
      hpxpy->Fill(px,py);
      hprof->Fill(px,pz);
      // fill only first 25000 events in NTuple
      if (i<25000) ntuple->Fill(px,py,pz,rnd,i);
      if (i && (i%kUPDATE) == 0) {
         if (i == kUPDATE) hpx->Draw();
         c1->Modified();
         c1->Update();
         if (i == kUPDATE) hfile->Write();

         if (gSystem->ProcessEvents()) break;
      }
      i++;
      if ((maxcnt>0) && (i>=maxcnt)) break;
   }

   gBenchmark->Show(jobname);
}
