/// \file
/// \ingroup tutorial_unuran
/// This program must be compiled and executed with Aclic as follows
///
/// ~~~{.cpp}
/// .x unuranFoamTest.C+
/// ~~~
///
/// it is an extension of tutorials foam_kanwa.C to compare
/// generation of a 2D distribution with unuran and Foam
///
/// \macro_code
///
/// \author Lorenzo Moneta

#include "TH2.h"
#include "TF2.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TStopwatch.h"
#include "TROOT.h"


#include "TUnuran.h"
#include "TUnuranMultiContDist.h"

#include <iostream>

//_____________________________________________________________________________
Double_t sqr(Double_t x){return x*x;};
//_____________________________________________________________________________
//_____________________________________________________________________________


Double_t Camel2(Int_t nDim, Double_t *Xarg){
// 2-dimensional distribution for Foam, normalized to one (within 1e-5)
  Double_t x=Xarg[0];
  Double_t y=Xarg[1];
  Double_t GamSq= sqr(0.100e0);
  Double_t Dist= 0;
  Dist +=exp(-(sqr(x-1./3) +sqr(y-1./3))/GamSq)/GamSq/TMath::Pi();
  Dist +=exp(-(sqr(x-2./3) +sqr(y-2./3))/GamSq)/GamSq/TMath::Pi();
  return 0.5*Dist;
}// Camel2

class FoamFunction : public TFoamIntegrand {
   public:
   virtual ~FoamFunction() {}
   double Density(int nDim, double * x) {
      return Camel2(nDim,x);
   }
   ClassDef(FoamFunction,1);

};

TH2 * hFoam;
TH2 * hUnr;


Int_t run_foam(int nev){
  cout<<"--- kanwa started ---"<<endl;
  gSystem->Load("libFoam.so");
  TH2D  *hst_xy = new TH2D("foam_hst_xy" ,  "FOAM x-y plot", 50,0,1.0, 50,0,1.0);
  hFoam = hst_xy;

  Double_t MCvect[2]; // 2-dim vector generated in the MC run
  //
  TRandom     *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(2);         // No. of dimensions, obligatory!
  FoamX->SetnCells(500);     // Optionally No. of cells, default=2000
  FoamX->SetRho(new FoamFunction() );  // Set 2-dim distribution, included below
  FoamX->SetPseRan(PseRan);  // Set random number generator
  //
  // From now on FoamX is ready to generate events

   // test first the time
   TStopwatch w;

  w.Start();
  FoamX->Initialize();       // Initialize simulator, may take time...

  //int nshow=5000;
  int nshow=nev;

  for(long loop=0; loop<nev; loop++){
    FoamX->MakeEvent();            // generate MC event
    FoamX->GetMCvect( MCvect);     // get generated vector (x,y)
    Double_t x=MCvect[0];
    Double_t y=MCvect[1];
    //if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    hst_xy->Fill(x,y);
    // live plot
    if(loop == nshow){
      nshow += 5000;
      hst_xy->Draw("lego2");
      //cKanwa->Update();
    }
  }// loop
  w.Stop();

  double time = w.CpuTime()*1.E9/nev;
  cout << "Time using FOAM \t\t " << "   \t=\t " << time << "\tns/call" << endl;

  //
  hst_xy->Draw("lego2");  // final plot
  //
  Double_t MCresult, MCerror;
  FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral, should be one
  cout << " MCresult= " << MCresult << " +- " << MCerror <<endl;
  cout<<"--- kanwa ended ---"<<endl;

  return 0;
}//kanwa



double UCamel2(double * x, double *) {
   return Camel2(2,x);
}

int run_unuran(int nev, std::string method = "hitro") {
   // use unuran

   std::cout << "run unuran " << std::endl;

   gSystem->Load("libUnuran.so");

   TH2D  *h1 = new TH2D("unr_hst_xy" ,  "UNURAN x-y plot", 50,0,1.0, 50,0,1.0);
   hUnr= h1;

   TF2 * f = new TF2("f",UCamel2,0,1,0,1,0);

   TUnuranMultiContDist dist(f);

   TRandom3 r;

   TUnuran unr(&r,2);  // 2 is debug level


   // test first the time
   TStopwatch w;

   w.Start();

   // init unuran
   bool ret =   unr.Init(dist,method);
   if (!ret) {
      std::cerr << "Error initializing unuran with method " << unr.MethodName() << endl;
      return -1;
   }

   double x[2];
   for (int i = 0; i < nev; ++i) {
      unr.SampleMulti(x);
      h1->Fill(x[0],x[1]);
//       if (method == "gibbs" && i < 100)
//          std::cout << x[0] << " , " << x[1] << std::endl;
   }

   w.Stop();
   double time = w.CpuTime()*1.E9/nev;
   cout << "Time using Unuran  " << unr.MethodName() << "   \t=\t " << time << "\tns/call" << endl;
   h1->Draw("lego2");
   return 0;
}

Int_t unuranFoamTest(){

  // visualising generated distribution
  TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,1000);
  cKanwa->Divide(1,2);
  cKanwa->cd(1);
  int n = 100000;


  run_unuran(n,"hitro");
  cKanwa->Update();

  cKanwa->cd(2);

  run_foam(n);
  cKanwa->Update();


  std::cout <<"\nChi2 Test Results (UNURAN-FOAM):\t";
  // test chi2
  hFoam->Chi2Test(hUnr,"UUP");

  return 0;
}
