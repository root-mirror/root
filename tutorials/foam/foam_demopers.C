/// \file
/// \ingroup tutorial_FOAM
/// \notebook -nodraw
/// This simple macro demonstrates persistency of FOAM object.
/// First run macro foam_demo.C to create file foam_demo.root with FOAM object.
///
/// Next type `root -l foam_demopers.C` from shell command line
///
/// \macro_code
///
/// \author Stascek Jadach

#include "Riostream.h"
#include "TFile.h"
#include "TFoam.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFoamIntegrand.h"


Int_t foam_demopers()
{

   // need to load the foam_demo tutorial for the definition of the function
   TString macroName = gROOT->GetTutorialDir();
   macroName.Append("/foam/foam_demo.C");
   gROOT->ProcessLine(TString::Format(".L %s+",macroName.Data()));

   //******************************************
   cout<<"====================== TestVector ================================"<<endl;
   TFile fileA("foam_demo.root");
   fileA.cd();
   cout<<"------------------------------------------------------------------"<<endl;
   fileA.ls();
   cout<<"------------------------------------------------------------------"<<endl;
   fileA.Map();
   cout<<"------------------------------------------------------------------"<<endl;
   fileA.ShowStreamerInfo();
   cout<<"------------------------------------------------------------------"<<endl;
   fileA.GetListOfKeys()->Print();
   cout<<"------------------------------------------------------------------"<<endl;
   //*******************************************
   TFoam  *FoamX = (TFoam*)fileA.Get("FoamX");
   //*******************************************
   //  FoamX->PrintCells();
   FoamX->CheckAll(1);

   //N.B. the integrand functions need to be reset
   // because cannot be made persistent
   TFoamIntegrand * rho = (TFoamIntegrand*) gROOT->ProcessLine("return new TFDISTR();");
   FoamX->SetRho(rho);

   Double_t MCvect[2]; // 2-dim vector generated in the MC run

   for(long loop=0; loop<50000; loop++){
      FoamX->MakeEvent();            // generate MC event
      FoamX->GetMCvect( MCvect);     // get generated vector (x,y)
      Double_t x=MCvect[0];
      Double_t y=MCvect[1];
      if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
      }// loop
   //
   Double_t IntNorm, Errel;
   FoamX->Finalize(   IntNorm, Errel);     // final printout
   Double_t MCresult, MCerror;
   FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral, should be one
   cout << " MCresult= " << MCresult << " +- " << MCerror <<endl;
   cout<<"===================== TestPers FINISHED ======================="<<endl;
   return 0;
}
//_____________________________________________________________________________
//

