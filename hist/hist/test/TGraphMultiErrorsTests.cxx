#include "TFile.h"
#include "TSystem.h"
#include "TGraphMultiErrors.h"

#include "gtest/gtest.h"

TEST(TGraphMultiErrorsTests, tgraphmultierrorstest)
{
   const auto ofileName = "tgraphmultierrorstest.root";
   TFile *f = TFile::Open(ofileName, "RECREATE");

   Double_t x[]       = {0, 1, 2, 3, 4};
   Double_t y[]       = {0, 2, 4, 1, 3};
   Double_t exl[]     = {0.3, 0.3, 0.3, 0.3, 0.3};
   Double_t exh[]     = {0.3, 0.3, 0.3, 0.3, 0.3};
   Double_t eylstat[] = {1, 0.5, 1, 0.5, 1};
   Double_t eyhstat[] = {0.5, 1, 0.5, 1, 2.};
   Double_t eylsys[]  = {0.5, 0.4, 0.8, 0.3, 1.2};
   Double_t eyhsys[]  = {0.6, 0.7, 0.6, 0.4, 0.8};

   auto gme = new TGraphMultiErrors(5, x, y, exl, exh, eylstat, eyhstat, TGraphMultiErrors::kOnlyFirst);
   gme->AddYError(5, eylsys, eyhsys);
   gme->SetMarkerStyle(20);
   gme->SetLineColor(kRed);
   gme->GetAttLine(0)->SetLineColor(kRed);
   gme->GetAttLine(1)->SetLineColor(kBlue);
   gme->GetAttFill(1)->SetFillStyle(0);

   gme->Write("gme");
   delete gme;

   gme = (TGraphMultiErrors*) f->Get("gme");

   EXPECT_DOUBLE_EQ(gme->GetErrorY(0), 0.79056941504209488);
   EXPECT_DOUBLE_EQ(gme->GetErrorY(4), 1.5811388300841898);

   gme->SetSumErrorsMode(TGraphMultiErrors::kSquareSum);

   EXPECT_DOUBLE_EQ(gme->GetErrorY(1), 0.97467943448089633);
   EXPECT_DOUBLE_EQ(gme->GetErrorY(3), 0.86602540378443871);

   gme->SetSumErrorsMode(TGraphMultiErrors::kAbsSum);

   EXPECT_DOUBLE_EQ(gme->GetErrorY(2), 1.4916433890176297);

   gme->DeleteYError(1);

   EXPECT_EQ(gme->GetNYErrors(), 1);

   gme->DeleteYError(0); // Must not work!

   EXPECT_EQ(gme->GetNYErrors(), 1);

   gme->AddYError(5, eylsys, eyhsys);

   EXPECT_EQ(gme->GetNYErrors(), 2);
   EXPECT_DOUBLE_EQ(gme->GetErrorY(0), 1.3152946437965904);

   gme->Set(4);
   gme->Set(5);
   gme->SetPoint(4, x[4], y[4]);
   gme->SetPointEX(4, exl[4], exh[4]);
   gme->SetPointEY(4, 0, eylstat[4], eyhstat[4]);
   gme->SetPointEY(4, 1, eylsys[4], eyhsys[4]);

   EXPECT_DOUBLE_EQ(gme->GetErrorY(4), 2.5179356624028344);

   delete gme;
   f->Close();
   delete f;
   gSystem->Unlink(ofileName);
}
