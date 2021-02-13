//test  Unuran using as  empirical distribution.
// Use as input an histogram (using its buffer) or a TGraph and TGraph2D for multi-dimensional data
//
// run the test within ROOT (.x unuranHist.cxx+) or pass any extra parameter in the command line to get
// a graphics output  (./unuranHist 1)
//

#include "TUnuran.h"
#include "TUnuranEmpDist.h"
#include "TUnuranMultiContDist.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TF2.h"
#include "TF3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "TStopwatch.h"
#include "TError.h"

#include <iterator>

#include <iostream>


int unuranHist() {

   // switch off printing of  info messages from chi2 test
   gErrorIgnoreLevel = 1001;


   int nbin = 100;
   double xmin = -5;
   double xmax = 20;

   int n = 100000; // number of events generated and of reference histo
   int ns = 10000;   // number of events from starting histo

   // h0 is reference histo
   TH1D * h0 = new TH1D("h0","Landau ref data",nbin,xmin,xmax);
   h0->FillRandom("landau",n);

   // h1 is low statistic histo from where we generate numbers
   TH1D * h1 = new TH1D("h1","Landau source data",nbin,xmin,xmax);
   h1->SetBuffer(ns);


   h1->FillRandom("landau",ns);

//    const double * bf = h1->GetBuffer();
//    std::ostream_iterator<double> oi(std::cout," ,  ");
//    std::copy(bf, bf+ 2*h1->GetBufferLength() + 2, oi);
//    std::cout << std::endl << std::endl;

   int iret = 0;

   std::cout << "Test Using 1D UnBinned data\n " << std::endl;


   TUnuran unr;
   TUnuranEmpDist dist(h1);

   if (!unr.Init(dist)) return -1;

   TStopwatch w;
   w.Start();
   TH1D * h2 = new TH1D("h2","Landau from Unuran unbin generation",nbin,xmin,xmax);
   for (int i = 0; i < n; ++i) {
      h2->Fill( unr.Sample() );
   }
   w.Stop();
   double time = w.CpuTime()*1.E9/n;
   std::cout << "Time using Unuran  " << unr.MethodName() << "   \t=\t " << time << "\tns/call" << std::endl;

   w.Start();
   TH1D * h3 = new TH1D("h3","h3",nbin,xmin,xmax);
   for (int i = 0; i < n; ++i) {
      h3->Fill( h1->GetRandom() );
   }
   w.Stop();
   time = w.CpuTime()*1.E9/n;
   std::cout << "Time using TH1::GetRandom()  \t=\t " <<  time << "\tns/call" << std::endl;


   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << " (with h0)\t:\t";
   double prob = h2->Chi2Test(h1,"UUP");
   if (prob < 1.E-6 ) {
      std::cerr << "Chi2 Test failed for UNURAN method " <<  unr.MethodName() << std::endl;
      iret = -1;
   }
   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << " (with ref)\t:\t";
   h2->Chi2Test(h0,"UUP");

   std::cout << "Test quality TH1::GetRandom (with h1) \t:\t";
   h3->Chi2Test(h1,"UUP");
   std::cout << "Test quality TH1::GetRandom (with ref) \t:\t";
   h3->Chi2Test(h0,"UUP");
   std::cout << "Comparison UnuRan-TH1::GetRandom \t:\t";
   h2->Chi2Test(h3,"UUP");



   TCanvas * c1 = new TCanvas("c1_unuranEmp","Empirical 1D distribution",10,10,800,800);
   c1->Divide(1,2);
   c1->cd(1);
   h2->SetLineColor(kBlue);
   h2->Draw();
   h0->Draw("same");
   h3->SetLineColor(kRed);
   h3->Draw("same");



   // generate using binned data from h1
   std::cout << "\nTest Using Binned data\n " << std::endl;

   TUnuranEmpDist dist2(h1,false);


   if (!unr.Init(dist2) ) return -1;



   w.Start();
   TH1D * h4 = new TH1D("h4","Landau from Unuran binned generation",nbin,xmin,xmax);
   for (int i = 0; i < n; ++i) {
      h4->Fill( unr.Sample() );
   }
   w.Stop();
   time = w.CpuTime()*1.E9/n;
   std::cout << "Time using Unuran  " << unr.MethodName() << "   \t=\t " << time << "\tns/call" << std::endl;

   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << " (with h1)\t:\t";
   double prob2 = h4->Chi2Test(h1,"UUP");
   if (prob2 < 1.E-6) {
      std::cerr << "Chi2 Test failed for UNURAN method " <<  unr.MethodName() << std::endl;
      iret = -2;
   }
   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << " (with ref)\t:\t";
   h4->Chi2Test(h0,"UUP");

   c1->cd(2);
   h4->SetLineColor(kBlue);
   h4->Draw();
   h0->Draw("same");

   return iret;
}

double gaus2d(double *x, double *p) {

   double sigma_x = p[0];
   double sigma_y = p[1];
   double rho = p[2];
   double u = x[0] / sigma_x ;
   double v = x[1] / sigma_y ;
   double c = 1 - rho*rho ;
   double result = (1 / (2 * TMath::Pi() * sigma_x * sigma_y * sqrt(c)))
      * exp (-(u * u - 2 * rho * u * v + v * v ) / (2 * c));
   return result;
}

double gaus3d(double *x, double *p) {

   double sigma_x = p[0];
   double sigma_y = p[1];
   double sigma_z = p[2];
   double rho = p[3];
   double u = x[0] / sigma_x ;
   double v = x[1] / sigma_y ;
   double w = x[2] / sigma_z ;
   double c = 1 - rho*rho ;
   double result = (1 / (2 * TMath::Pi() * sigma_x * sigma_y * sigma_z * sqrt(c)))
      * exp (-(u * u - 2 * rho * u * v + v * v + w*w) / (2 * c));
   return result;
}


int unuranGraf() {

   int iret = 0;

   std::cout << "\nTest Using 2D UnBinned data\n " << std::endl;

   // generate graf with x-y data
   const int Ndata = 100000;
   const int n = 100000;

   TF2 * f2 = new TF2("g2d",gaus2d,-10,10,-10,10,3);
   double p2[3] = {1,1,0.5};
   f2->SetParameters(p2);

   TUnuran unr(gRandom,2);
   if (!unr.Init(TUnuranMultiContDist(f2),"vnrou") ) return -1;
   // generate 2d data
   double xx[2];
   // generate ref data
   TH2D * href = new TH2D("href2d","UNURAN reference 2D gauss data",100,-5,5,100,-5,5);
   for (int i = 0; i < n; ++i) {
      unr.SampleMulti(xx);
      href->Fill(xx[0],xx[1]);
   }

   // create a graph with source xy data
   TGraph * gr = new TGraph();

   for (int i = 0; i< Ndata; ++i) {
      unr.SampleMulti( xx);
      gr->SetPoint(i,xx[0],xx[1]);
   }

   TH2D * h2 = new TH2D("h2d","UNURAN generated 2D gauss data",100,-5,5,100,-5,5);

   TH1D * hx = new TH1D("hx","x gen variable",100,-5,5);
   TH1D * hy = new TH1D("hy","y gen variable",100,-5,5);

   // now generate random points from the graph
   TUnuranEmpDist dist(Ndata,gr->GetX(), gr->GetY() );
   if (!unr.Init(dist) ) {
      std::cerr << "Initialization failed for method " << unr.MethodName() << std::endl;
      return -1;
   }


   TStopwatch w;
   w.Start();
   for (int i = 0; i < n; ++i) {
      unr.SampleMulti(xx);
      h2->Fill(xx[0],xx[1]);
      hx->Fill(xx[0]);
      hy->Fill(xx[1]);
   }
   w.Stop();
   double time = w.CpuTime()*1.E9/n;
   std::cout << "Time using Unuran  " << unr.MethodName() << "   \t=\t " << time << "\tns/call" << std::endl;

   TCanvas * c1 = new TCanvas("c1_unuranEmp2D","Empirical Multidim distribution",300,10,800,400);

   c1->Divide(2,1);
   c1->cd(1);
   gr->Draw("AP");
   c1->cd(2);
   h2->Draw("surf");
   //f2->Draw("same");

   TCanvas * c2 = new TCanvas("c2d","Empirical Multidim distribution",300,100,800,400);
   c2->Divide(2,1);
   c2->cd(1);
   hx->Draw();
   c2->cd(2);
   hy->Draw();

   // apply chi2 test to href
   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << "\t\t:\t";
   double prob = href->Chi2Test(h2,"UUP");
   if (prob < 1.E-6) {
      std::cerr << "Chi2 Test failed for UNURAN method " <<  unr.MethodName() << std::endl;
      iret = -2;
   }

   return iret;
}

int unuranGraf2D() {
   int iret = 0;

   std::cout << "\nTest Using 3D UnBinned data\n " << std::endl;

   // generate graf with x-y data
   const int Ndata = 100000;
   const int n = 100000;

   TF3 * f3 = new TF3("g3d",gaus3d,-10,10,-10,10,-10,10,4);
   double p[4] = {1,1,1,0.5};
   f3->SetParameters(p);

   TUnuran unr(gRandom,2);
   if (!unr.Init(TUnuranMultiContDist(f3),"vnrou") ) return -1;
   // generate 3d data
   double xx[3];
   // generate ref data
   TH3D * href = new TH3D("href3d","UNURAN reference 3D gauss data",50,-5,5,50,-5,5,50,-5,5);

   for (int i = 0; i < n; ++i) {
      unr.SampleMulti(xx);
      href->Fill(xx[0],xx[1],xx[2]);
   }

   // create a graph with xy data
   TGraph2D * gr = new TGraph2D();
   for (int i = 0; i< Ndata; ++i) {
      unr.SampleMulti( xx);
      gr->SetPoint(i,xx[0],xx[1],xx[2]);
   }

   TH3D * h3 = new TH3D("h3d","UNURAN generated 3D gauss data",50,-5,5,50,-5,5,50,-5,5);

   TH1D * hx = new TH1D("hx3","x gen variable",100,-5,5);
   TH1D * hy = new TH1D("hy3","y gen variable",100,-5,5);
   TH1D * hz = new TH1D("hz3","z gen variable",100,-5,5);

   // now generate random points from the graph
   TUnuranEmpDist dist(Ndata,gr->GetX(), gr->GetY(), gr->GetZ() );
   if (!unr.Init(dist) ) {
      std::cerr << "Initialization failed for method " << unr.MethodName() << std::endl;
      return -1;
   }


   TStopwatch w;
   w.Start();
   for (int i = 0; i < n; ++i) {
      unr.SampleMulti(xx);
      h3->Fill(xx[0],xx[1],xx[2]);
      hx->Fill(xx[0]);
      hy->Fill(xx[1]);
      hz->Fill(xx[2]);
   }
   w.Stop();
   double time = w.CpuTime()*1.E9/n;
   std::cout << "Time using Unuran  " << unr.MethodName() << "   \t=\t " << time << "\tns/call" << std::endl;

   TCanvas * c3 = new TCanvas("c3d","Empirical Multidim distribution",300,600,800,400);
   c3->Divide(3,1);
   c3->cd(1);
   hx->Draw();
   c3->cd(2);
   hy->Draw();
   c3->cd(3);
   hz->Draw();

   TCanvas * c1 = new TCanvas("c1_unuranEmp3D","Empirical Multidim distribution",400,500,900,400);
   c1->Divide(2,1);
   c1->cd(1);
   gr->Draw("AP");
   c1->cd(2);
   h3->Draw("surf1");
   f3->Draw("same");

   // apply chi2 test to href
   std::cout << "\nTest quality UNURAN  " << unr.MethodName() << "\t\t:\t";
   double prob = href->Chi2Test(h3,"UUP");
   if (prob < 1.E-6) {
      std::cerr << "Chi2 Test failed for UNURAN method " <<  unr.MethodName() << std::endl;
      iret = -2;
   }

   return iret;
}

#ifndef __CINT__
int main(int argc, char **argv)
{
   int iret = 0;
   if (argc > 1) {
      TApplication theApp("App",&argc,argv);
      iret |=  unuranHist();
      iret |=  unuranGraf();
      iret |=  unuranGraf2D();
      theApp.Run();
   }
   else {
      iret |=  unuranHist();
      iret |=  unuranGraf();
      iret |=  unuranGraf2D();
   }

   if (iret != 0)
      std::cerr <<"\n\nUnuRan Empirical Distribution Test:\t  Failed !!!!!!!\n" << std::endl;
   else
      std::cout << "\n\nUnuRan  Empirical Distribution Test:\t OK\n" << std::endl;
   return iret;
}
#endif
