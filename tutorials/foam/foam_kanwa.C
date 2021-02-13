/// \file
/// \ingroup tutorial_FOAM
/// \notebook -js
/// This program can be execute from the command line as folows:
///
/// ~~~{.cpp}
///      root -l foam_kanwa.C
/// ~~~
///
/// \macro_code
///
/// \author Stascek Jadach

#include "Riostream.h"
#include "TFoam.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TMath.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"

//_____________________________________________________________________________
Double_t sqr(Double_t x){
   return x*x;
}
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
}
//_____________________________________________________________________________

Int_t foam_kanwa(){
   cout<<"--- kanwa started ---"<<endl;
   TH2D  *hst_xy = new TH2D("hst_xy" ,  "x-y plot", 50,0,1.0, 50,0,1.0);
   Double_t MCvect[2]; // 2-dim vector generated in the MC run
   //
   TRandom     *PseRan   = new TRandom3();  // Create random number generator
   PseRan->SetSeed(4357);
   TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
   FoamX->SetkDim(2);         // No. of dimensions, obligatory!
   FoamX->SetnCells(500);     // Optionally No. of cells, default=2000
   FoamX->SetRhoInt(Camel2);  // Set 2-dim distribution, included below
   FoamX->SetPseRan(PseRan);  // Set random number generator
   FoamX->Initialize();       // Initialize simulator, may take time...
   //
   // visualising generated distribution
   TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,600);
   cKanwa->cd();
   // From now on FoamX is ready to generate events
   int nshow=5000;
   for(long loop=0; loop<100000; loop++){
      FoamX->MakeEvent();            // generate MC event
      FoamX->GetMCvect( MCvect);     // get generated vector (x,y)
      Double_t x=MCvect[0];
      Double_t y=MCvect[1];
      if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
      hst_xy->Fill(x,y);
      // live plot
      if(loop == nshow){
         nshow += 5000;
         hst_xy->Draw("lego2");
         cKanwa->Update();
      }
   }// loop
   //
   hst_xy->Draw("lego2");  // final plot
   cKanwa->Update();
   //
   Double_t MCresult, MCerror;
   FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral, should be one
   cout << " MCresult= " << MCresult << " +- " << MCerror <<endl;
   cout<<"--- kanwa ended ---"<<endl;

   return 0;
}
