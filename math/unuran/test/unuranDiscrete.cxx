// test using discrete distribution.
// Generate numbers from a given probability vector or from a discrete distribution like
// the Poisson distribution.
// Compare also the Unuran method for generating Poisson numbers with TRandom::Poisson
//
// run within ROOT (.x unuranDiscrete.cxx+) or pass any extra parameter in the command line to get
// a graphics output (./unuranDiscrete 1 )


#include "TUnuran.h"
#include "TUnuranDiscrDist.h"

#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TError.h"


#ifdef WRITE_OUTPUT
#include "TFile.h"
#include "TGraph.h"
#endif

#include "Math/Functor.h"
#include "Math/DistFunc.h"
#include "Math/Util.h"

#include <iostream>
#include <iomanip>

//#define DEBUG
#ifdef DEBUG
int n = 100;
#else
int n = 1000000;
#endif
TCanvas * c1;
int icanv = 1;

bool useRandomSeed = false;   // to use a random seed different every time

double poisson_pmf(double * x, double * p) {

   double y = ROOT::Math::poisson_pdf(int(x[0]),p[0]);
//   std::cout << x[0] << " f(x) = " << y << std::endl;
   return y;
}
double binomial_pmf(double * x, double * p) {

   double y = ROOT::Math::binomial_pdf(static_cast<unsigned int>(x[0]),p[1], static_cast<unsigned int>(p[0]));
//   std::cout << x[0] << " f(x) = " << y << std::endl;
   return y;
}

int testUnuran(TUnuran & unr, double & time, TH1 * h1, const TH1 * href,  bool weightHist = false ) {


   // test first the time
   TStopwatch w;

   w.Start();
   for (int i = 0; i < n; ++i)
      unr.SampleDiscr();

   w.Stop();
   time = w.CpuTime()*1.E9/n;


   // test quality (use cdf to avoid zero bins)
   h1->Reset();
   int n2 = n/10;
   int x;
   for (int i = 0; i<n2; ++i) {
      x =  unr.SampleDiscr();
      h1->Fill( double( x) );
   }
   double prob;
   if (weightHist)
      prob = h1->Chi2Test(href,"UW");
   else
      prob = h1->Chi2Test(href,"UU");
   std::string s = "Time using Unuran  " +  unr.MethodName();
   std::cout << std::left << std::setw(40) << s << "\t=\t " << time << "\tns/call \t\tChi2 Prob = "<< prob << std::endl;
   if (prob < 1E-06) {
      std::cout << "Chi2 Test failed for method " << unr.MethodName() << std::endl;
      if (weightHist)
         h1->Chi2Test(href,"UWP"); // print all chi2 test info
      else
         h1->Chi2Test(href,"UUP"); // print all chi2 test info

      return 1;
   }
   return 0;
}

int testRootPoisson(double mu, double &time, TH1 * h2) {
   TStopwatch w;

   w.Start();
   for (int i = 0; i < n; ++i)
      gRandom->Poisson(mu);

   w.Stop();
   time = w.CpuTime()*1.E9/n;

   // make ref histo
   int n2 = n/10;
   for (int i = 0; i< n2; ++i) {
      h2->Fill ( gRandom->Poisson(mu) );
   }

   std::string s = "Time using TRandom::Poisson";
   std::cout << std::left << std::setw(40) << s << "\t=\t " << time << std::endl;
   return 0;
}

int testRootBinomial(int m, double p, double &time, TH1 * h2) {

   TStopwatch w;

   w.Start();
   for (int i = 0; i < n; ++i)
      gRandom->Binomial(m,p);

   w.Stop();
   time = w.CpuTime()*1.E9/n;

   // make ref histo
   int n2 = n/10;
   for (int i = 0; i< n2; ++i) {
      h2->Fill ( gRandom->Binomial(m,p) );
   }

   std::string s = "Time using TRandom::Binomial";
   std::cout << std::left << std::setw(40) << s << "\t=\t " << time << std::endl;
   return 0;
}

int testProbVector() {

   int iret = 0;

   TH1D * h0 = new TH1D("h0","reference prob",10,-0.5,9.5);
   TH1D * h1 = new TH1D("h1","UNURAN gen prob",10,-0.5,9.5);


   double p[10] = {1.,2.,3.,5.,3.,2.,1.,0.5,0.3,0.5 };
   for (int i = 0; i< 10; ++i) {
      h0->SetBinContent(i+1,p[i]);
      h0->SetBinError(i+1,0.);
   }
   double sum = h0->GetSumOfWeights();
   std::cout << " prob sum = " << sum << std::endl;

   TUnuran unr(gRandom,2);

   TUnuranDiscrDist dist(p,p+10);

//    TUnuranDiscrDist fDDist = dist;
//    std::cout << fDDist.ProbVec().size() << std::endl;

   std::cout << "Test  generation with a PV :\n\n";

   double time;

   bool ret = unr.Init(dist,"method=dau");
   if (!ret) iret = -1;
   iret |= testUnuran(unr,time,h1,h0,true);

   ret = unr.Init(dist,"method=dgt");
   if (!ret) iret = -1;
   iret |= testUnuran(unr,time,h1,h0,true);

   // dss require the sum
   dist.SetProbSum(sum);
   ret = unr.Init(dist,"method=dss");
   if (!ret) iret = -1;
   iret |= testUnuran(unr,time,h1,h0,true);





   c1->cd(icanv++);
   h1->Sumw2();
   h1->Scale( h0->GetSumOfWeights()/(h1->GetSumOfWeights() ) );
   h0->Draw();
   h1->SetLineColor(kBlue);
   h1->Draw("Esame");
   c1->Update();

   return iret;
}

int testPoisson() {
   // test with a function
   // Poisson distribution

   int iret = 0;
   std::cout << "\nTest  generation with a Probability function (Poisson) :\n\n";

   TF1 * f = new TF1("f",poisson_pmf,1,0,1);

   // loop on mu values for Nmu times

   const int Nmu = 5;
   double muVal[Nmu] = {5,10,20,50,100};
   double tR[Nmu];
   double tU[Nmu];
   double tUdari[Nmu];
   double tUdsrou[Nmu];

  TUnuran unr(gRandom,2);

  for (int imu = 0; imu < Nmu; ++imu) {

     const double mu = muVal[imu];

     int nmax = static_cast<int>(mu*3);
     TH1D * h2 = new TH1D("h2","reference Poisson prob",nmax,0,nmax);
     TH1D * h3 = new TH1D("h3","UNURAN gen prob",nmax,0,nmax);
     TH1D * h4 = new TH1D("h4","UNURAN gen prob 2",nmax,0,nmax);


     std::cout << "\nPoisson mu = " << mu << std::endl << std::endl;

     testRootPoisson(mu,tR[imu],h2);

     // test UNURAN with standard method
     bool ret = unr.InitPoisson(mu);
     if (!ret) iret = -1;
     testUnuran(unr,tU[imu],h3,h2);


     // test changing all the time the mu
     // use re-init for a fast re-initialization
     TStopwatch w;
     unr.InitPoisson(mu,"dstd");
     double p[1]; p[0] = mu;
     w.Start();
     for (int i = 0; i < n; ++i) {
        unr.ReInitDiscrDist(1,p);
        int k = unr.SampleDiscr();
        if (n % 10 == 0) h4->Fill(k);
     }
     w.Stop();
     double time = w.CpuTime()*1.E9/n;
     double prob = h2->Chi2Test(h4,"UU");
     std::string s = "Time using Unuran w/ re-init  method=" + unr.MethodName();
     std::cout << std::left << std::setw(40) << s << "\t=\t " << time
               << "\tns/call \t\tChi2 Prob = "<< prob << std::endl;

     if (prob < 1E-06) {
        std::cout << "Chi2 Test failed for re-init " << std::endl;
        iret = -2;
     }

     f->SetParameter(0,mu);

#ifdef USE_FUNCTOR
     ROOT::Math::Functor1D f2(f);
     TUnuranDiscrDist dist2 = TUnuranDiscrDist(f2);
#else
     TUnuranDiscrDist dist2 = TUnuranDiscrDist(f);
#endif

     // dari method (needs mode and pdf sum)
     dist2.SetMode(int(mu) );
     dist2.SetProbSum(1. );
     ret = unr.Init(dist2,"method=dari");
     if (!ret) iret = -1;

     iret |= testUnuran(unr,tUdari[imu],h3,h2);

     // dsrou method (needs mode and sum)

     dist2.SetProbSum(1);
     ret = unr.Init(dist2,"method=dsrou");
     if (!ret) iret = -1;

     iret |= testUnuran(unr,tUdsrou[imu],h3,h2);


     c1->cd(icanv++);
     h2->DrawCopy();
     h3->SetLineColor(kBlue);
     h3->DrawCopy("Esame");
     c1->Update();

     delete h2;
     delete h3;
     delete h4;
  }

#ifdef WRITE_OUTPUT

   TFile * file = new TFile("unuranPoisson.root","RECREATE");
   // create graphs with results
   TGraph * gR = new TGraph(Nmu,muVal,tR);
   TGraph * gU = new TGraph(Nmu,muVal,tU);
   TGraph * gU2 = new TGraph(Nmu,muVal,tUdari);
   TGraph * gU3 = new TGraph(Nmu,muVal,tUdsrou);
   gR->Write();
   gU->Write();
   gU2->Write();
   gU3->Write();
   file->Close();

#endif

   return iret;
}

int testBinomial() {
   // test using binomial distribution
   int iret = 0;

   std::cout << "\nTest  generation with a Probability function (Binomimal) :\n\n";

   TF1 * f = new TF1("f",binomial_pmf,1,0,2);

   // loop on binomual values

   const int NBin = 3;
   double pVal[NBin] = {0.5,0.1,0.01};
   double NVal[NBin] = {20,100,1000};
   double tR[NBin];
   double tU[NBin];
   double tUdari[NBin];
   double tUdsrou[NBin];


  TUnuran unr(gRandom,2);


  for (int ib = 0; ib < NBin; ++ib) {


     double par[2];
     par[0] = NVal[ib];
     par[1] = pVal[ib];

     int nmax = static_cast<int>(par[0]*par[1]*3);
     TH1D * h2 = new TH1D("h2","reference Binomial prob",nmax,0,nmax);
     TH1D * h3 = new TH1D("h3","UNURAN gen prob",nmax,0,nmax);
     TH1D * h4 = new TH1D("h4","UNURAN gen prob 2",nmax,0,nmax);


     std::cout << "\nBinomial n = " << par[0] << "  " << par[1] << std::endl;

     testRootBinomial(static_cast<int>(par[0]),par[1],tR[ib],h2);


     // test UNURAN with standard method
     bool ret = unr.InitBinomial(static_cast<int>(par[0]),par[1]);
     if (!ret) iret = -1;
     testUnuran(unr,tU[ib],h3,h2);


     // test changing all the time the mu
     // use re-init for a fast re-initialization

     TStopwatch w;
     unr.InitBinomial(static_cast<int>(par[0]), par[1],"dstd");
     w.Start();
     for (int i = 0; i < n; ++i) {
        unr.ReInitDiscrDist(2,par);
        int k = unr.SampleDiscr();
        if (n % 10 == 0) h4->Fill(k);
     }
     w.Stop();
     double time = w.CpuTime()*1.E9/n;
     double prob = h2->Chi2Test(h4,"UU");
     std::string s = "Time using Unuran w/ re-init  method=" + unr.MethodName();
     std::cout << std::left << std::setw(40) << s << "\t=\t " << time
               << "\tns/call \t\tChi2 Prob = "<< prob << std::endl;

     if (prob < 1E-06) {
        std::cout << "Chi2 Test failed for re-init " << std::endl;
        iret = -2;
     }

     // test the universal methods

     f->SetParameters(par);
     TUnuranDiscrDist dist = TUnuranDiscrDist(f);

     // dari method (needs mode and pdf sum)
     dist.SetMode(int(par[0]*par[1]) );
     dist.SetProbSum(1. );
     ret = unr.Init(dist,"method=dari");
     if (!ret) iret = -1;

     iret |= testUnuran(unr,tUdari[ib],h3,h2);

     // dsrou method (needs mode and sum)

     ret = unr.Init(dist,"method=dsrou");
     if (!ret) iret = -1;

     iret |= testUnuran(unr,tUdsrou[ib],h3,h2);


     c1->cd(icanv++);
     h2->DrawCopy();
     h3->SetLineColor(kBlue);
     h3->DrawCopy("Esame");
     c1->Update();

     delete h2;
     delete h3;
     delete h4;
  }


   return iret;
}

int unuranDiscrete() {

   int iret = 0;

   c1 = new TCanvas("c1_unuranDiscr_PV","Discrete distribution from PV",10,10,800,800);
   c1->Divide(3,3);

   // switch off printing of  info messages from chi2 test
   gErrorIgnoreLevel = 1001;

   // check if using a random seed
   if (useRandomSeed) gRandom->SetSeed(0);

   iret |= testProbVector();
   iret |= testPoisson();
   iret |= testBinomial();

   if (iret != 0)
      std::cerr <<"\n\nUnuRan Discrete Distribution Test:\t  Failed !!!!!!!" << std::endl;
   else
      std::cerr << "\n\nUnuRan Discrete Distribution Test:\t OK" << std::endl;

   return iret;

}


#ifndef __CINT__
int main(int argc, char **argv)
{
   int iret = 0;
   if (argc > 1) {
      TApplication theApp("App",&argc,argv);
      iret =  unuranDiscrete();
      theApp.Run();
   }
   else
      iret =  unuranDiscrete();

   return iret;
}
#endif
