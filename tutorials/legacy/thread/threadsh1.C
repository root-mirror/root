/// \file
/// \ingroup tutorial_thread
/// Example of a simple script creating 3 threads.
/// This script can only be executed via ACliC .x threadsh1.C++.
///
/// \macro_code
///
/// \author Victor Perevovchikov

#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TThread.h"


TCanvas *c[4];
TH1F    *hpx[4];
TThread *t[5];
TRandom * rng[4];
Bool_t finished;

void *handle(void *ptr)
{
   long nr = (long) ptr;
   Long64_t nfills = 25000000;
   int upd = 50000;

   char name[32];
   sprintf(name,"hpx%ld",nr);
   TThread::Lock();
   hpx[nr] = new TH1F(name,"This is the px distribution",100,-4,4);
   hpx[nr]->SetFillColor(48);
   TThread::UnLock();
   Float_t px, py, pz;
   rng[nr]->SetSeed(0);
   for (Int_t i = 0; i < nfills; i++) {
      rng[nr]->Rannor(px,py);
      pz = px*px + py*py;
      hpx[nr]->Fill(px);
      if (i && (i%upd) == 0) {
         if (i == upd) {
            TThread::Lock();
            c[nr]->cd();
            hpx[nr]->Draw();
            TThread::UnLock();
         }
         if (c[nr]) c[nr]->Modified();
         gSystem->Sleep(10);
      }
   }
   return 0;
}

void *joiner(void *)
{
   t[0]->Join();
   t[1]->Join();
   t[2]->Join();
   t[3]->Join();

   finished = kTRUE;

   return 0;
}

void closed(Int_t id)
{
   // kill the thread matching the canvas being closed
   t[id]->Kill();
   // and set the canvas pointer to 0
   c[id] = 0;
}

void threadsh1()
{

   finished = kFALSE;
   //gDebug = 1;

   c[0] = new TCanvas("c0","Dynamic Filling Example",100,20,400,300);
   c[0]->SetFillColor(42);
   c[0]->GetFrame()->SetFillColor(21);
   c[0]->GetFrame()->SetBorderSize(6);
   c[0]->GetFrame()->SetBorderMode(-1);
   c[1] = new TCanvas("c1","Dynamic Filling Example",510,20,400,300);
   c[1]->SetFillColor(42);
   c[1]->GetFrame()->SetFillColor(21);
   c[1]->GetFrame()->SetBorderSize(6);
   c[1]->GetFrame()->SetBorderMode(-1);
   c[2] = new TCanvas("c2","Dynamic Filling Example",100,350,400,300);
   c[2]->SetFillColor(42);
   c[2]->GetFrame()->SetFillColor(21);
   c[2]->GetFrame()->SetBorderSize(6);
   c[2]->GetFrame()->SetBorderMode(-1);
   c[3] = new TCanvas("c3","Dynamic Filling Example",510,350,400,300);
   c[3]->SetFillColor(42);
   c[3]->GetFrame()->SetFillColor(21);
   c[3]->GetFrame()->SetBorderSize(6);
   c[3]->GetFrame()->SetBorderMode(-1);

   // connect to the Closed() signal to kill the thread when a canvas is closed
   c[0]->Connect("Closed()", 0, 0, "closed(Int_t=0)");
   c[1]->Connect("Closed()", 0, 0, "closed(Int_t=1)");
   c[2]->Connect("Closed()", 0, 0, "closed(Int_t=2)");
   c[3]->Connect("Closed()", 0, 0, "closed(Int_t=3)");

   rng[0] = new TRandom3(1);
   rng[1] = new TRandom3(2);
   rng[2] = new TRandom3(3);
   rng[3] = new TRandom3(4);

   printf("Starting Thread 0\n");
   t[0] = new TThread("t0", handle, (void*) 0);
   t[0]->Run();
   printf("Starting Thread 1\n");
   t[1] = new TThread("t1", handle, (void*) 1);
   t[1]->Run();
   printf("Starting Thread 2\n");
   t[2] = new TThread("t2", handle, (void*) 2);
   t[2]->Run();
   printf("Starting Thread 3\n");
   t[3] = new TThread("t3", handle, (void*) 3);
   t[3]->Run();
   printf("Starting Thread 4\n");
   t[4] = new TThread("t4", joiner, (void*) 3);
   t[4]->Run();

   TThread::Ps();

   while (!finished) {
      for (int i = 0; i < 4; i++) {
         if (c[i] && c[i]->IsModified()) {
            //printf("Update canvas %d\n", i);
            c[i]->Update();
         }
      }
      gSystem->Sleep(100);
      gSystem->ProcessEvents();
   }

   t[4]->Join();
   TThread::Ps();

   delete t[0];
   delete t[1];
   delete t[2];
   delete t[3];
   delete t[4];

   delete rng[0];
   delete rng[1];
   delete rng[2];
   delete rng[3];
}
