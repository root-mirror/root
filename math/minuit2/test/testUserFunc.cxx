// @(#)root/minuit2:$Id$
// Author: L. Moneta    10/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#include "TApplication.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "Math/MinimizerOptions.h"
#include "TMath.h"

#include <iostream>

double myfunc(double *x, double *p)
{

   return p[0] * TMath::Gaus(x[0], p[1], p[2]);
}

void testUserFunc(std::string type = "Minuit2", int n = 1000)
{

   gRandom = new TRandom3();

   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(type.c_str());

   TH1D *h1 = new TH1D("h1", "fit histo 1", 100, -5, 5.);

   //   gStyle->SetOptStat(1111111);
   //   gStyle->SetOptFit(1111111);

   for (int i = 0; i < n; ++i) {
      h1->Fill(gRandom->Gaus(0, 1));
   }

   TF1 *f = new TF1("f", myfunc, -10, 10, 3);
   double p[3] = {100.0, 0.0, 1.0};
   f->SetParameters(p);

   h1->Fit(f);
   // try fix a parameter
   // TVirtualFitter * fitter = TVirtualFitter::GetFitter();
   // std::cout << typeid(*fitter).name() << std::endl;
   // fitter->FixParameter(2);
   f->FixParameter(2, 1.0);

   h1->Fit(f, "V");

   h1->Draw();
}

#ifndef __CINT__
int main(int argc, char **argv)
{
   if (argc > 1) {
      TApplication theApp("App", &argc, argv);
      testUserFunc();
      theApp.Run();
   } else
      testUserFunc();
   return 0;
}
#endif

//#ifndef __CINT__
// int main() {
//  testUserFunc( );
//}
//#endif
