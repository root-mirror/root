// @(#)root/test:$name:  $:$id: stressSpectrum.cxx,v 1.15 2002/10/25 10:47:51 rdm exp $
// Author: Rene Brun 17/01/2006

/////////////////////////////////////////////////////////////////
//
//    TSPectrum test suite
//    ====================
//
// This stress program tests many elements of the TSpectrum, TSpectrum2 classes.
//
// To run in batch, do
//   stressSpectrum        : run 100 experiments with graphics (default)
//   stressSpectrum 1000   : run 1000 experiments with graphics
//   stressSpectrum -b 200 : run 200 experiments in batch mode
//   stressSpectrum -b     : run 100 experiments in batch mode
//
// To run interactively, do
// root -b
//  Root > .x stressSpectrum.cxx      : run 100 experiments with graphics (default)
//  Root > .x stressSpectrum.cxx(20)  : run 20 experiments
//  Root > .x stressSpectrum.cxx+(30) : run 30 experiments via ACLIC
//
// Several tests are run sequentially. Each test will produce one line (Test OK or Test FAILED) .
// At the end of the test a table is printed showing the global results
// Real Time and Cpu Time.
// One single number (ROOTMARKS) is also calculated showing the relative
// performance of your machine compared to a reference machine
// a Pentium IV 3.0 Ghz) with 512 MBytes of memory
// and 120 GBytes IDE disk.
//
// An example of output when all the tests run OK is shown below:
//
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//****************************************************************************
//*  Starting  stress S P E C T R U M                                        *
//****************************************************************************
//Peak1 : found = 70.21/ 73.75, good = 65.03/ 68.60, ghost = 8.54/ 8.39,--- OK
//Peak2 : found =163/300, good =163, ghost =8,----------------------------  OK
//****************************************************************************
//stressSpectrum: Real Time =  19.86 seconds Cpu Time =  19.04 seconds
//****************************************************************************
//*  ROOTMARKS = 810.9   *  Root5.09/01   20051216/1229
//****************************************************************************

#include <cstdlib>
#include <iostream>
#include "snprintf.h"
#include "TApplication.h"
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"

Int_t npeaks;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}
Double_t fpeaks2(Double_t *x, Double_t *par) {
   Double_t result = 0.1;
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm   = par[5*p+0];
      Double_t mean1  = par[5*p+1];
      Double_t sigma1 = par[5*p+2];
      Double_t mean2  = par[5*p+3];
      Double_t sigma2 = par[5*p+4];
      result += norm*TMath::Gaus(x[0],mean1,sigma1)*TMath::Gaus(x[1],mean2,sigma2);
   }
   return result;
}
void findPeaks(Int_t pmin, Int_t pmax, Int_t &nfound, Int_t &ngood, Int_t &nghost) {
   npeaks = (Int_t)gRandom->Uniform(pmin,pmax);
   Int_t nbins = 500;
   Double_t dxbins = 2;
   TH1F *h = new TH1F("h","test",nbins,0,nbins*dxbins);
   //generate n peaks at random
   Double_t par[3000];
   par[0] = 0.8;
   par[1] = -0.6/1000;
   Int_t p,pf;
   for (p=0;p<npeaks;p++) {
      par[3*p+2] = 1;
      par[3*p+3] = 10+gRandom->Rndm()*(nbins-20)*dxbins;
      par[3*p+4] = 3+2*gRandom->Rndm();
   }
   TF1 *f = new TF1("f",fpeaks,0,nbins*dxbins,2+3*npeaks);
   f->SetNpx(1000);
   f->SetParameters(par);
   h->FillRandom("f",200000);
   TSpectrum *s = new TSpectrum(4*npeaks);
   nfound = s->Search(h,2,"goff");
   //Search found peaks
   ngood = 0;
   Double_t *xpeaks = s->GetPositionX();
   for (p=0;p<npeaks;p++) {
      for (pf=0;pf<nfound;pf++) {
         Double_t dx = TMath::Abs(xpeaks[pf] - par[3*p+3]);
         if (dx <dxbins) ngood++;
      }
   }
   //Search ghost peaks
   nghost = 0;
   for (pf=0;pf<nfound;pf++) {
      Int_t nf=0;
      for (p=0;p<npeaks;p++) {
         Double_t dx = TMath::Abs(xpeaks[pf] - par[3*p+3]);
         if (dx <dxbins) nf++;
      }
      if (nf == 0) nghost++;
   }
   delete f;
   delete h;
   delete s;
}

void stress1(Int_t ntimes) {
   Int_t pmin = 5;
   Int_t pmax = 55;
   TCanvas *c1 = new TCanvas("c1","Spectrum results",10,10,800,800);
   c1->Divide(2,2);
   gStyle->SetOptFit();
   TH1F *hpeaks = new TH1F("hpeaks","Number of peaks",pmax-pmin,pmin,pmax);
   TH1F *hfound = new TH1F("hfound","% peak founds",100,0,100);
   TH1F *hgood  = new TH1F("hgood", "% good peaks",100,0,100);
   TH1F *hghost = new TH1F("hghost","% ghost peaks",100,0,100);
   Int_t nfound,ngood,nghost;
   for (Int_t i=0;i<ntimes;i++) {
      findPeaks(pmin,pmax,nfound,ngood,nghost);
      hpeaks->Fill(npeaks);
      hfound->Fill(100*Double_t(nfound)/Double_t(npeaks));
      hgood->Fill(100*Double_t(ngood)/Double_t(npeaks));
      hghost->Fill(100*Double_t(nghost)/Double_t(npeaks));
      //printf("npeaks = %d, nfound = %d, ngood = %d, nghost = %d\n",npeaks,nfound,ngood,nghost);
   }
   c1->cd(1);
   hpeaks->Fit("pol1","lq");
   c1->cd(2);
   hfound->Fit("gaus","lq");
   c1->cd(3);
   hgood->Fit("gaus","lq");
   c1->cd(4);
   hghost->Fit("gaus","lq","",0,30);
   c1->cd();
   Double_t p1  = hfound->GetFunction("gaus")->GetParameter(1);
   Double_t ep1 = hfound->GetFunction("gaus")->GetParError(1);
   Double_t p2  = hgood->GetFunction("gaus")->GetParameter(1);
   Double_t ep2 = hgood->GetFunction("gaus")->GetParError(1);
   Double_t p3  = hghost->GetFunction("gaus")->GetParameter(1);
   Double_t ep3 = hghost->GetFunction("gaus")->GetParError(1);
   Double_t p1ref = 70.21; //ref numbers obtained with ntimes=1000
   Double_t p2ref = 65.03;
   Double_t p3ref =  8.54;

   //printf("p1=%g+-%g, p2=%g+-%g, p3=%g+-%g\n",p1,ep1,p2,ep2,p3,ep3);

   char sok[20];
   if (TMath::Abs(p1ref-p1) < 2*ep1 && TMath::Abs(p2ref-p2) < 2*ep2  && TMath::Abs(p3ref-p3) < 2*ep3 ) {
      snprintf(sok,20,"OK");
   } else {
      snprintf(sok,20,"failed");
   }
   printf("Peak1 : found =%6.2f/%6.2f, good =%6.2f/%6.2f, ghost =%5.2f/%5.2f,--- %s\n",
          p1,p1ref,p2,p2ref,p3,p3ref,sok);
}
void stress2(Int_t np2) {
   npeaks = np2;
   TRandom r;
   Int_t nbinsx = 200;
   Int_t nbinsy = 200;
   Double_t xmin   = 0;
   Double_t xmax   = (Double_t)nbinsx;
   Double_t ymin   = 0;
   Double_t ymax   = (Double_t)nbinsy;
   Double_t dx = (xmax-xmin)/nbinsx;
   Double_t dy = (ymax-ymin)/nbinsy;
   TH2F *h2 = new TH2F("h2","test",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   h2->SetStats(0);
   //generate n peaks at random
   Double_t par[3000];
   Int_t p;
   for (p=0;p<npeaks;p++) {
      par[5*p+0] = r.Uniform(0.2,1);
      par[5*p+1] = r.Uniform(xmin,xmax);
      par[5*p+2] = r.Uniform(dx,5*dx);
      par[5*p+3] = r.Uniform(ymin,ymax);
      par[5*p+4] = r.Uniform(dy,5*dy);
   }
   TF2 *f2 = new TF2("f2",fpeaks2,xmin,xmax,ymin,ymax,5*npeaks);
   f2->SetNpx(100);
   f2->SetNpy(100);
   f2->SetParameters(par);
   h2->FillRandom("f2",500000);
   //now the real stuff
   TSpectrum2 *s = new TSpectrum2(2*npeaks);
   Int_t nfound = s->Search(h2,2,"goff noMarkov");

   //searching good and ghost peaks (approximation)
   Int_t pf,ngood = 0;
   Double_t *xpeaks = s->GetPositionX();
   Double_t *ypeaks = s->GetPositionY();
   for (p=0;p<npeaks;p++) {
      for (pf=0;pf<nfound;pf++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 2*dx && diffy < 2*dy) ngood++;
      }
   }
   if (ngood > nfound) ngood = nfound;
   //Search ghost peaks (approximation)
   Int_t nghost = 0;
   for (pf=0;pf<nfound;pf++) {
      Int_t nf=0;
      for (p=0;p<npeaks;p++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 2*dx && diffy < 2*dy) nf++;
      }
      if (nf == 0) nghost++;
   }

   delete s;
   delete f2;
   delete h2;
   Int_t nfoundRef = 163;
   Int_t ngoodRef  = 163;
   Int_t nghostRef = 8;
   char sok[20];
   if (  TMath::Abs(nfound - nfoundRef) < 5
      && TMath::Abs(ngood - ngoodRef) < 5
      && TMath::Abs(nghost - nghostRef) < 5)  {
      snprintf(sok,20,"OK");
   } else {
      snprintf(sok,20,"failed");
   }
   printf("Peak2 : found =%d/%d, good =%d, ghost =%2d,---------------------------- %s\n",
          nfound,npeaks,ngood,nghost,sok);
}

void stressSpectrum(Int_t ntimes=100) {
   std::cout << "****************************************************************************" <<std::endl;
   std::cout << "*  Starting  stress S P E C T R U M                                        *" <<std::endl;
   std::cout << "****************************************************************************" <<std::endl;
   gBenchmark->Start("stressSpectrum");
   stress1(ntimes);
   stress2(300);
   gBenchmark->Stop ("stressSpectrum");
   Double_t reftime100 = 19.04; //pcbrun compiled
   Double_t ct = gBenchmark->GetCpuTime("stressSpectrum");
   const Double_t rootmarks = 800*reftime100*ntimes/(100*ct);
   printf("****************************************************************************\n");

   gBenchmark->Print("stressSpectrum");
   printf("****************************************************************************\n");
   printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
   printf("****************************************************************************\n");
}

#ifndef __CINT__

int main(int argc, char **argv)
{
   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   gBenchmark = new TBenchmark();
   Int_t ntimes = 100;
   if (argc > 1)  ntimes = atoi(argv[1]);
   stressSpectrum(ntimes);
   return 0;
}

#endif
