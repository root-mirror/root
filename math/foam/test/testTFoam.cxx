// Tests for TFoam
// Author: Stephan Hageboeck, CERN  04/2020

#include "TFoam.h"
#include "TFile.h"
#include "TMath.h"

#include "gtest/gtest.h"

Double_t sqr(Double_t x){
   return x*x;
}

Double_t Camel2(Int_t/*nDim*/, Double_t *Xarg){
   // 2-dimensional distribution for Foam, normalized to one (within 1e-5)
   Double_t x=Xarg[0];
   Double_t y=Xarg[1];
   Double_t GamSq= sqr(0.100e0);
   Double_t Dist= 0;
   Dist +=exp(-(sqr(x-1./3) +sqr(y-1./3))/GamSq)/GamSq/TMath::Pi();
   Dist +=exp(-(sqr(x-2./3) +sqr(y-2./3))/GamSq)/GamSq/TMath::Pi();
   return 0.5*Dist;
}

// Read a simple v6.20 (and before) TFoam
// This is the example from the TFoam documentation.
TEST(TFoam, Readv1) {
  TFile file("testTFoam_1.root", "READ");
  ASSERT_TRUE(file.IsOpen());

  TFoam* foam = nullptr;
  file.GetObject("foam", foam);
  ASSERT_NE(foam, nullptr);

  foam->SetRhoInt(Camel2);
  EXPECT_EQ(foam->GetTotDim(), 2);

  double results[5][2] = {
      {0.7010131486, 0.6501394535},
      {0.3375651929, 0.334055375},
      {0.2137994088, 0.2531630374},
      {0.7420976609, 0.7784924952},
      {0.7099058364, 0.8037966862}};

  for (int i=0; i<5; ++i) {
    double x[2];
    foam->MakeEvent();
    foam->GetMCvect(x);
    EXPECT_NEAR(x[0], results[i][0], 1.E-9);
    EXPECT_NEAR(x[1], results[i][1], 1.E-9);
  }
}
