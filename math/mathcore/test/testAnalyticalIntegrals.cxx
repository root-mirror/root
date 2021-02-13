//
//  TestAnalyticalIntegrals
//  
//
//  Created by Aurélie Flandi on 10.09.14.
//
//

#include <TCanvas.h>
#include <iostream>//for cout
#include <TROOT.h>
#include "TMath.h"
#include <TF1.h>
#include <TApplication.h>
#include <Math/PdfFuncMathCore.h>//for pdf
#include <Math/ProbFuncMathCore.h>//for cdf
#include <Math/WrappedTF1.h>
#include <Math/Integrator.h>


using namespace std;

void testAnalyticalIntegrals()
{

//compare analytical integration with numerical one

   ROOT::Math::Integrator ig; 

   TCanvas * c1 = new TCanvas("pol3","pol3",800,1000);
   c1->Divide(3,3);
   int ipad = 0;
   {
      TF1 *f = new TF1("TRY","pol3",-5.,5.);
      f -> SetParameters(1.,1.,1.,1.);
      c1->cd(++ipad);
      f->Draw();

      ROOT::Math::WrappedTF1 wf(*f);
      ig.SetFunction(wf);

      double numInt = ig.Integral(-5.,5.);
      double anaInt = f->Integral(-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);

   }

   {
      TF1 *f  = new TF1("MyExp","expo",-5.,5.);
      f -> SetParameters(0.2,-0.3);
      c1->cd(++ipad);
      f ->Draw();

      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);

   }  
   {
      TF1 *f  = new TF1("MyCrystalBall","crystalball",-5.,5.);
      f -> SetParameters(2,1,0.5,2.,0.9);
      c1->cd(++ipad);
      f ->Draw();

      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);
      
   }
   {
      // CB with alpha < 0
      TF1 *f  = new TF1("MyCrystalBall","crystalball",-5.,5.);
      f -> SetParameters(2,-1,0.5,-2.,0.9);
      c1->cd(++ipad);
      f ->Draw();

      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);
   }
   {
      TF1 *f  = new TF1("MyGauss","gaus",-5.,5.);
      f -> SetParameters(2.,0.,0.3);
      c1->cd(++ipad);
      f ->Draw();
      ROOT::Math::WrappedTF1 wf(*f);
      
      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);
}
   {
      TF1 *f  = new TF1("MyGauss","gausn",-5.,5.);
      f -> SetParameters(2.,0.,0.3);
      c1->cd(++ipad);
      f ->Draw();

      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-6))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %20.10f ana = %20.10f diff = %20.10f",f->GetTitle(),numInt,anaInt,numInt-anaInt);

   }

   {
      TF1 *f  = new TF1("MyExp","landau",-5.,5.);
      f -> SetParameters(2.,2,0.3);
      c1->cd(++ipad);
      f ->Draw();

      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;

      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);

   }  
   {
      TF1 *f  = new TF1("MyExp","landaun",-5.,5.);
      f -> SetParameters(2.,-2.,0.3);
      c1->cd(++ipad);
      f ->Draw();
      ROOT::Math::WrappedTF1 wf(*f);

      double anaInt = f->Integral(-5.,5.);
      double numInt = ig.Integral(wf,-5.,5.);

      std::cout<<"analytical integral for " << f->GetTitle()  << " = " << anaInt << std::endl;
      std::cout<<"numerical  integral for " << f->GetTitle()  << " = " << numInt << std::endl;
      
      if (!TMath::AreEqualAbs(numInt, anaInt, 1.E-8))
         Error("TestAnalyticalIntegral","Different integral value for %s num = %f ana = %f diff = %f",f->GetTitle(),numInt,anaInt,numInt-anaInt);

   }  
}

int main(int argc, char **argv)
{
   bool showGraphics = false;
   // Parse command line arguments
   for (Int_t i=1 ;  i<argc ; i++) {
      std::string arg = argv[i] ;
      if (arg == "-g") {
         showGraphics = true;
      }
      if (arg == "-h") {
         cerr << "Usage: " << argv[0] << " [-g] [-v]\n";
         cerr << "  where:\n";
         cerr << "     -g : graphics mode\n";
         cerr << "     -v : verbose  mode";
         cerr << endl;
         return -1; 
      }
   }
   
   TApplication* theApp = 0;
   if ( showGraphics )
      theApp = new TApplication("App",&argc,argv);
  
   testAnalyticalIntegrals();
  
   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return 0;
}
