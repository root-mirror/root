#include "TApplication.h"
#include "TObject.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TThread.h"

/*
** Module name : threads
**
** Description :
**  This file implements four functions, which shall perform some
**  actions on objects like histograms and a canvas.
*/


// a global canvas for the drawing of histograms
TCanvas *c1;
// a pad for the histogram, will be used by thread mhs
TPad *pad1;
// a pad for the histogram, will be used by thread mhs2
TPad *pad2;

// a global histogram object, which will be accessed both by mhs and mhs1
TH1F *total0;

// thread function (filling histogram main and drawing total0 on pad1)
void *mhs(void *)
{
   TThread::Printf("Start of mhs");

   // change to the canvas and draw its correspondig pad 1
   c1->cd();
   c1->SetGrid();
   pad1->Draw();

   // creating a histogram is explicitly locked by a global mutex
   TThread::Lock();
   TH1F *main   = new TH1F("main","Main contributor",100,-4,4);
   TThread::UnLock();

   // other operations needn't be locked
   total0->Sumw2();
   total0->SetMarkerStyle(21);
   total0->SetMarkerSize(0.7);
   main->SetFillColor(16);

   // Fill histogram main randomly
   gRandom->SetSeed();
   const Int_t kUPDATE = 100;
   Float_t xmain;
   // run for a certain amount of events or
   //for ( Int_t i=0; i<=100000; i++) {
   // run forever
   for ( Int_t i=0; ; i++) {
      xmain = gRandom->Gaus(-1,1.5);
      // filling of a histogram is locked to avoid data loss
      TThread::Lock();
      main->Fill(xmain);
      total0->Fill(xmain);
      TThread::UnLock();
      // drawing of the histogram is locked automatically by ROOT
      if (!(i%kUPDATE)) {
         if (i == kUPDATE) {
            pad1->cd();
            total0->Draw("e1p");
         }
         pad1->Modified();
         pad1->Update();
         c1->Modified();
         c1->Update();
      }
      // sleep for 1 ms: sleep not necessary, slows things only a bit down
      // because the threads are actually doing nothing, which is eventually
      // very fast ;-)
      TThread::Sleep(0, 1000000);
   }
   TThread::Printf("End of mhs\n");

   c1->Modified();
   c1->Update();
   return 0;
}

// thread function (filling total0)
void *mhs1(void *)
{
   TThread::Printf("Start of mhs1");

   // instantiation of objects has to be locked explicitly
   TThread::Lock();
   TH1F *s1     = new TH1F("s1","This is the first signal",100,-4,4);
   TH1F *s2     = new TH1F("s2","This is the second signal",100,-4,4);
   TThread::UnLock();

   s1->SetFillColor(42);
   s2->SetFillColor(46);

   // Fill histograms randomly
   gRandom->SetSeed();
   Float_t xs1, xs2;
   //for ( Int_t i=0; i<=100000; i++) {
   for ( Int_t i=0; ; i++) {
      xs1   = gRandom->Gaus(-0.5,0.5);
      xs2   = gRandom->Gaus(1,0.3);
      // locking avoids data loss
      TThread::Lock();
      s1->Fill(xs1,0.3);
      s2->Fill(xs2,0.2);
      total0->Fill(xs1,0.3);
      total0->Fill(xs2,0.2);
      TThread::UnLock();
      // sleep for 6 ms: not necessary (see above)
      TThread::Sleep(0, 6000000);
   }
   TThread::Printf("End of mhs1\n");

   // wo is last draws the histogram a last time
   pad1->cd();
   total0->Draw("e1p");
   c1->Modified();
   c1->Update();
   return 0;
}

// thread function: plays with its own histograms, draw on canvas in pad2
void *mhs2(void *)
{
   TThread::Printf("Start of mhs2");

   c1->cd();
   pad2->Draw();

   // Create some histograms: explicitly locked.
   TThread::Lock();
   TH1F *total  = new TH1F("total2","This is the total distribution",100,-4,4);
   TH1F *main   = new TH1F("main2","Main contributor",100,-4,4);
   TH1F *s1     = new TH1F("s12","This is the first signal",100,-4,4);
   TH1F *s2     = new TH1F("s22","This is the second signal",100,-4,4);
   TThread::UnLock();

   total->Sumw2();
   total->SetMarkerStyle(21);
   total->SetMarkerSize(0.7);
   main->SetFillColor(16);
   s1->SetFillColor(42);
   s2->SetFillColor(46);

   // Fill histograms randomly
   gRandom->SetSeed();
   const Int_t kUPDATE = 100;
   Float_t xs1, xs2, xmain;
   // for ( Int_t i=0; i<=100000; i++) {
   for ( Int_t i=0; ; i++) {
      xmain = gRandom->Gaus(-1,1.5);
      xs1   = gRandom->Gaus(-0.5,0.5);
      xs2   = gRandom->Gaus(1,0.3);
      // no locking here because nobody else uses these objects
      main->Fill(xmain);
      s1->Fill(xs1,0.3);
      s2->Fill(xs2,0.2);
      total->Fill(xmain);
      total->Fill(xs1,0.3);
      total->Fill(xs2,0.2);
      if (!(i%kUPDATE)) {
         if (i == kUPDATE) {
            pad2->cd();
            total->Draw("e1p");
         }
         pad2->Modified();
         pad2->Update();
         c1->Modified();
         c1->Update();
      }
      // sleep for 1 ms: not necessary
      TThread::Sleep(0, 1000000);
   }
   TThread::Printf("End of mhs2\n");
   c1->Modified();
   c1->Update();
   return 0;
}

// thread to run Ps(): perform every 5 seconds a TThread::Ps()
void *top(void *)
{
   TThread::Printf("Start of top");

   for (Int_t i = 0; ;i++) {
      TThread::Ps();
      TThread::Sleep(5);
   }
   TThread::Printf("End of top");
   return 0;
}

int main(int argc, char **argv)
{
   TApplication theApp("h1Thread", &argc, argv);

   // a global canvas for the drawing of histograms
   c1 = new TCanvas("c1","The HSUM example",800,400);
   // a pad for the histogram, will be used by thread mhs
   pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.48,0.98,33);
   // a pad for the histogram, will be used by thread mhs2
   pad2 = new TPad("pad2","This is pad2",0.52,0.02,0.98,0.98,33);

   // a global histogram object, which will be accessed both by mhs and mhs1
   total0  = new TH1F("total","This is the total distribution",100,-4,4);

   // create the TThread instances
   TThread *th1 = new TThread("th1",mhs);
   TThread *th2 = new TThread("th2",mhs2);
   TThread *th3 = new TThread("th3",mhs1);

   // top thread
   TThread *thp = new TThread("top",top);

   th1->Run();
   th2->Run();
   th3->Run();
   thp->Run();

   printf("*** Exit the program by selecting Quit from the File menu ***\n");
   theApp.Run(kTRUE);

   th1->SetCancelAsynchronous();
   th2->SetCancelAsynchronous();
   th3->SetCancelAsynchronous();
   thp->SetCancelAsynchronous();

   th1->Kill();
   th2->Kill();
   th3->Kill();
   thp->Kill();

   printf("The END...\n");

   return 0;
}

